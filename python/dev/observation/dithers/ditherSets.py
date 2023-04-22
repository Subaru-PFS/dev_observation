from importlib import reload
import os.path
import copy

import numpy as np

from pfs.drp.stella.utils import raster
reload(raster)

class DitherSets:
    def __init__(self, config, loadFrom=None, spectrographs=[1,3]):
        self.config = config
        self.spectrographs = spectrographs

        self.ditherSets = dict()
        if loadFrom is not None:
            for ds in loadFrom.ditherSets.values():
                self.newDitherSet(ds.name, ds.visits, ds.comment)
                self[ds.name].setRawDithers(ds.rawDithers)

    def __str__(self):
        return f'DitherSet("{list(self.ditherSets.keys())}"'

    def __getitem__(self, key):
        return self.ditherSets[key]

    def newDitherSet(self, name, visits, comment=None, force=False):
        if name in self.ditherSets:
            if force:
                print(f'overwriting ditherSet {name}')
            else:
                print(f'ditherSet {name} exists, not overwriting')
                return self[name]

        self.ditherSets[name] = DitherSet(self.config, name, visits, 
                                          spectrographs=self.spectrographs,
                                          comment=comment)
        return self[name]

class DitherConfig:
    def __init__(self, butler, dbConn, dataId, fiberTraces,
                 camera, lmin, lmax):
        self.butler = butler
        self.dbConn = dbConn
        self.dataId = dataId
        self.fiberTraces = fiberTraces
        self.camera = camera
        self.lmin = lmin
        self.lmax = lmax

class DitherSet:
    def __init__(self, config, name, visits, comment=None, 
                 useGuideOffsets=False, spectrographs=[1,3]):
        self.config = config

        self.name = name
        self.visits = list(visits)
        self.comment = comment
        self.useGuideOffsets = useGuideOffsets
        self.rawDithers = []
        self.dithers = []
        self.images = []
        self.imageCfg = dict()
        self.imageDithers = []
        self.spectrographs = spectrographs

    @property
    def pfsConfig(self):
        if len(self.dithers) == 0:
            return None

        return self.dithers[0].pfsConfig

    def setRawDithers(self, rawDithers):
        self.rawDithers = rawDithers
        self.dithers = []
        self.images = []
        self.extent_CI = []

    def updateVisits(self, newVisitList):
        if len(self.visits) > newVisitList:
            raise ValueError('fewer visits than already exist')
        if self.visits != newVisitList[:len(self.visit)]:
            raise ValueError('new visits do not start with existing visits')
        self.visits = newVisitList

    def flush(self):
        self.setRawDithers([])

    def fullName(self, suffix='', dir=None):
        if suffix:
            suffix = '_' + suffix

        name = f'{self.name}_{self.visits[0]}{suffix}'
        if dir is not None:
            name = os.path.join(dir, name)

        return name

    def loadAvailable(self, dataId, forceReload=False,
                      useDitherOffset=True):
        visit0 = self.visits[0]

        if forceReload:
            self.setRawDithers([])

        cfg = self.config
        b = cfg.butler
        ra0 = None
        for v in self.visits:
            if v in [d.visit for d in self.rawDithers]:
                print(f"{v} is cached")
                continue

            print(f"Processing {v}")

            ds = []
            for spectrograph in self.spectrographs:
                dataId.update(visit=v, arm='r', spectrograph=spectrograph)
                try:
                    ft = cfg.fiberTraces[f'r{spectrograph}']
                    d = raster.makeDither(b, dataId, fiberTraces=ft, lmin=cfg.lmin, lmax=cfg.lmax,
                                          camera=cfg.camera, useButler=False, visit0=visit0, usePfsArm=True)
                except Exception as e:
                    print(f'Failed to makeDither: {e}')
                    continue
                ds.append(d)

            d = None
            if len(ds) > 1:
                d = raster.concatenateDithers(b, ds)
            print(f"madeDithers {v}")

            try:
                dra, ddec = raster.getGuideOffset(self.config.dbConn, d.visit)[:2]
                d.dra = dra
                d.ddec = ddec

                if useDitherOffset:
                    if ra0 is None:
                        ra0, dec0 = d.ra, d.dec

                    d.dither_ra, d.dither_dec = raster.getDitherRaDec(self.config.dbConn, d.visit)[0:2]
                    raoff, decoff = d.dither_ra/np.cos(np.deg2rad(d.dec)), d.dither_dec

                    if np.isnan(raoff):
                        raoff = 0
                    if np.isnan(decoff):
                        decoff = 0

                    if self.useGuideOffsets:
                        raoff += d.dra
                        decoff += d.ddec

                    d.ra =  ra0  + raoff/3600
                    d.dec = dec0 + decoff/3600

            except Exception as e:
                print(f'failed to calculate guide offsets: {e}')

            if d is not None:
                self.rawDithers.append(d)
            self.freeze()
        print(f"Done: {len(self.rawDithers)}/{len(self.visits)}")

    def updateImageCfg(self, side=None, pixelScale=None, R=None,
                       usePFImm=True, icrosstalk=False, correctExtinction=True):

        imageCfgUpdates = dict(side=side, pixelScale=pixelScale, R=R, usePFImm=usePFImm,
                               icrosstalk=icrosstalk, correctExtinction=correctExtinction)
        newImageCfg = self.imageCfg.copy()
        newImageCfg.update(imageCfgUpdates)

        changed = newImageCfg != self.imageCfg
        if changed:
            print(f'forcing image rebuild: {self.imageCfg} != {newImageCfg}')
            self.imageCfg = newImageCfg
            self.imageDithers = []
            self.images = []

        return changed

    def freeze(self):
        """ Heinous. """
        self.dithers = copy.deepcopy(self.rawDithers)[:len(self.visits)]
        self.images = []
        return self.dithers

    def getCobraImages(self, fiberId=None, forceRebuild=False):
        """ Get current dither images. If anything has changed or the images are not available, (re-)construct them.


        This still needs to learn how to update images just for new dithers.
        """
        changed = len(self.dithers) != len(self.imageDithers)
        changed |= forceRebuild

        if changed:
            print(f'rebuilding images: {changed} force={forceRebuild} {len(self.dithers)} != {len(self.imageDithers)}')
            imCfg = self.imageCfg
            if imCfg['correctExtinction']:
                extinction = {}
                for d in self.dithers:
                    extinction[d.visit] = raster.estimateExtinction(self.config.dbConn, d.visit)
            else:
                extinction = None
            self.extinction = extinction

            images, extent_CI, visitImage = raster.makeCobraImages(self.dithers,
                                                                   extinction=extinction,
                                                                   side=imCfg['side'],
                                                                   pixelScale=imCfg['pixelScale'],
                                                                   R=imCfg['R'], usePFImm=imCfg['usePFImm'],
                                                                   fiberIds=None if fiberId is None else [fiberId],
                                                                   icrosstalk=imCfg['icrosstalk'],
                                                                   setUnimagedPixelsToNaN=True)
            self.images = images
            self.extent_CI = extent_CI
            self.visitImage = visitImage
            self.imageDithers = [d.visit for d in self.dithers]

        return self.images, self.extent_CI, self.visitImage
