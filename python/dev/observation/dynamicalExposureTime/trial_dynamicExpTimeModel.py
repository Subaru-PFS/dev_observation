#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
#from opdb import utils as op
import psycopg2
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import datetime as dt
from pfs.utils.coordinates.DistortionCoefficients import radec_to_subaru

# new database access 
#from pfs.utils.database import qadb

# Ignore warning abour space motion
import warnings
# An Astropy module that raises ErfaWarning
import erfa
warnings.filterwarnings('ignore', module='erfa')

# For option
from argparse import ArgumentParser

def get_option():

    bottom_text = 'You need to ssh with -Y option to send a window. And you may want configure .pgpass under your home directory'

    argparser = ArgumentParser(epilog=bottom_text)
    argparser.add_argument('expt', type=float, help='nominal exposure time [s]')
    #argparser.add_argument('visit',help='visit id for spectrograph')
    argparser.add_argument('-m', '--mcs', type=int, default=None,
                           help='visit id for convergence to adjust ecposure time')
    argparser.add_argument('-vi', '--visit', action="store_true",
                           help='Use previous visit instead of time')
    #argparser.add_argument('-t', '--date', type=str, default=None,
    #                       help='Date to show seeing/transparency. (e.g., 2025-01-01)')
    #argparser.add_argument('-r', '--rotRange', type=float, default=3.,
    #                       help='threshold of rotator error')
    #argparser.add_argument('-a', '--alazRange', type=float, default=0.5,
    #                       help='threshold of alt/az error')
    return argparser.parse_args()


def get_skycondition_for_visit(visit):

    """ Query seeing and transparency for a given visit
    ----------
    visit: visit id to get sky condition
    n: int , number of get
    """

    ## qadb
    conn = psycopg2.connect("dbname='qadb' host='pfsa-db' port=5436 user='pfs'") 
    # new method
    #qa_db=qadb.QaDB()
    #conn=qa_db.connect()


    items = f'*'

    #tables1 = f'seeing_agc_exposure'
    #tables2 = f'transparency_agc_exposure'
    tables1 = f'seeing'
    tables2 = f'transparency'

    condition = f"pfs_visit_id={visit}"
    
    que1 = f'select {items} from {tables1} WHERE {condition}'
    que2 = f'select {items} from {tables2} WHERE {condition}'

    with conn.cursor() as cur:
        cur.execute(que1)
        df1 = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])
        cur.execute(que2)
        df2 = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])

    # new method
    #df1.query_dataframe(que1, conn=conn)
    #df2.query_dataframe(que2, conn=conn)


    return df1, df2

def get_skycondition_time(ts, te):

    """ Query seeing and transparency for a given visit
    ----------
    ts: timestanp, the starting time to search window
    te: timestanp, the ending time to search window
    """

    ## qadb
    conn = psycopg2.connect("dbname='qadb' host='pfsa-db' port=5436 user='pfs'") 
    # new method
    #qa_db=qadb.QaDB()
    #conn=qa_db.connect()

    items = f'*'

    tables1 = f'seeing_agc_exposure'
    tables2 = f'transparency_agc_exposure'

    condition = f"taken_at BETWEEN '{ts}' AND '{te}'"
    
    que1 = f'select {items} from {tables1} WHERE {condition}'
    que2 = f'select {items} from {tables2} WHERE {condition}'

    with conn.cursor() as cur:
        cur.execute(que1)
        df1 = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])
        cur.execute(que2)
        df2 = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])

    # new method
    # df1.query_dataframe(que1, conn=conn)
    # df2.query_dataframe(que2, conn=conn)

    return df1, df2


def estimate_eet(seeing_past, transp_past, airmass):

    """
    https://github.com/KevinMacAstro/KevinMacAstro.github.io/blob/main/forecaster.html
    """
    
    coef_dicts = {
        'b': {'a1': 0.25205554669238134, 'a2': -1.428688585473366, 'a3': -0.5883696020777965,
              'C': 0.8360016758605285, 'R': 0.056749237191502454, 'sig': 0.420986524418459},
        'r': {'a1': 0.8509599727172921, 'a2': -1.1651223051928536,'a3': -1.0313878703838228, 
              'C': 1.0815585176161657, 'R': 0.060295422415749964, 'sig': 0.39646837711891375},
        'm': {'a1': 0.45446489233839404, 'a2': -1.554070290330937, 'a3': -1.470590358828559,
              'C': 1.272745155125946, 'R': 0.05550936668365635, 'sig': 0.5700657130611495},
        'n': {'a1': -0.10485118203041388, 'a2': -0.241044955115914, 'a3': -0.33116599395355245, 
              'C': 0.7485761248254111, 'R': 0.047163676545725314, 'sig': 0.24156412218308423}
    }

    factor=[]
    for arm in ['b', 'r', 'n', 'm']:
        coef = coef_dicts[arm]
        log_ratio = coef['a1'] * np.log(transp_past) + coef['a2'] * np.log(seeing_past) + coef['a3'] * np.log(airmass) + np.log(coef['C'])
        ratP = np.exp(log_ratio) + coef['R']
        ratP = 1/ratP
        #print(ratP)
        factor.append(ratP)

    return np.array(factor)


def show_dynamical_eet(expt, mcs=None, useVisit=False):

    # opdb
    conn2 = psycopg2.connect("dbname='opdb' host='pfsa-db' port=5432 user='pfs'") 
    # new method
    #op_db=DB(host='pfsa-db', user='pfs', port='5432', dbname='opdb')
    #conn2=op_db.connect()


    # Select convergence of interest
    items = '*'
    tables = 'tel_status'
    if mcs is None:
        condition= "caller='mcs' order BY pfs_visit_id DESC LIMIT 1"
    else:
        condition= f"caller='mcs' and pfs_visit_id= {mcs} order BY pfs_visit_id DESC LIMIT 1"

    que = f'select {items} from {tables} where {condition}'

    with conn2.cursor() as cur:
        cur.execute(que)
        df = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])
    # new method
    # df=op_db.query_dataframe(que, conn=conn2)

    df = df.drop_duplicates()

    visit_mcs = df.pfs_visit_id.values[0]
    ra=df.tel_ra.values[0]
    dec=df.tel_dec.values[0]
    pa=df.inst_pa.values[0]
    time=df.created_at.values[0]
    # 10 min (and HST -> UTC)
    time = time + np.timedelta64(dt.timedelta(minutes=10, hours=10))

    if useVisit:    # get sky condition from the last sps exposure
        items = '*'
        tables = 'tel_status JOIN sps_visit ON tel_status.pfs_visit_id=sps_visit.pfs_visit_id'
        condition= f"tel_status.pfs_visit_id <{visit_mcs} AND caller='iic' ORDER BY tel_status.pfs_visit_id DESC LIMIT 1"

        que = f'select {items} from {tables} where {condition}'

        with conn2.cursor() as cur:
            cur.execute(que)
            df = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])

        # new method
        #df=op_db.query_dataframe(que, conn=conn2)
        #df = df.drop_duplicates()

        visit_sps = df.pfs_visit_id.values[0][0]
        dfs, dft = get_skycondition_for_visit(visit_sps)
        ts = np.nan
        te = np.nan
        try:
            seeing_last=dfs.seeing_median.values[0]
            trans_last=dft.transparency_median.values[0]
        except IndexError:
            print(f"SpS data before {visit_mcs} doesn't have seeing/trnsparency information. I set it nan.")
            seeing_last=np.nan
            trans_last=np.nan
    else:     # get sky condition a certain time window
        items = "taken_at"
        tables = "mcs_exposure" 
        condition= f"mcs_frame_id = {visit_mcs}00"

        que = f'select {items} from {tables} where {condition}'

        with conn2.cursor() as cur:
            cur.execute(que)
            df = pd.DataFrame(cur.fetchall(), columns=[col.name for col in cur.description])

        # new method
        #df=op_db.query_dataframe(que, conn=conn2)
        #df = df.drop_duplicates()

        t0 = df.taken_at.values[0]
        ts = t0 - np.timedelta64(dt.timedelta(minutes=15, hours=0))
        te = t0 - np.timedelta64(dt.timedelta(minutes=5, hours=0))
        #print(pd.to_datetime(ts).strftime('%Y-%m-%d %H:%M:%S.%f'), te.astype(dt.datetime))
        dfs, dft = get_skycondition_time(ts, te)
        #dfs, dft = get_skycondition_time(pd.to_datetime(ts).strftime('%Y-%m-%d %H:%M:%S.%f'), 
        #                                 pd.to_datetime(te).strftime('%Y-%m-%d %H:%M:%S.%f'))
        try:
            seeing_last=dfs.median().seeing_median
            trans_last=dft.median().transparency_median
        except IndexError:
            print(f"SpS data before {visit_mcs} doesn't have seeing/trnsparency information. I set it nan.")
            seeing_last=np.nan
            trans_last=np.nan



    if trans_last > 1:
        trans_last = 1

    time_str=np.datetime_as_string(time)
    if not useVisit:
        ts=np.datetime_as_string(ts)
        te=np.datetime_as_string(te)
    else:
        ts='9999-99-99 99:99:99.000000000'
        te='9999-99-99 99:99:99.000000000'
    az, el, inr = radec_to_subaru(ra, dec, pa, time_str, 2016., 0., 0., 1e-5)
    airmass = 1/np.cos(np.deg2rad(90.-el))

    factor=estimate_eet(seeing_last, trans_last, airmass)
    eet=[f*expt for f in factor]
    #print(f"{visit_mcs},{ra:.2f}, {dec:.2f}, {pa:.2f}, {time}, {az:.2f}, el, inr, seeing_last, trans_last, airmass, eet)

    print(f"visit_mcs,     ra,    dec,     pa,    aimed_time_UTC,       sky_time1_HST,      sky_time2_HST,seeing,transp,airmass,   el,expt_b,expt_r,expt_n,expt_m")
    print(f"   {visit_mcs},{ra:7.2f},{dec:7.2f},{pa:7.2f},{time_str[:-10]},{ts[:-10]},{te[:-10]}, {seeing_last:.2f},  {trans_last:.2f},    {airmass:.2f},{el:.2f},  {eet[0]:4.0f},  {eet[1]:4.0f},  {eet[2]:4.0f},  {eet[3]:4.0f}")

    with open('log_dynamical_exposre.txt', 'a') as fout:
        print(f"{visit_mcs},{ra:f},{dec:f},{pa:f},{time_str[:-10]},{ts},{te},{seeing_last:f},{trans_last:f},{airmass:f},{el:f},{eet[0]:f},{eet[1]:f},{eet[2]:f},{eet[3]:f}", file=fout)
        


if __name__ == '__main__':
    
    args = get_option()
    expt=args.expt
    
    show_dynamical_eet(expt, mcs=args.mcs, useVisit=args.visit)