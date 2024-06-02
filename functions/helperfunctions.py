# a bunch of other functions that are useful and used in many scripts, 
# but not fundamental to the project

import pandas as pd
import time
import csv
import gzip
import math
import numpy as np
import warnings
from PairingVDJdb_GA import get_mode, check_correct_pairs_ID

def load_df(mypath, N):
    chunksize = 1000000

    # the list that contains all the dataframes
    list_of_dataframes = []

    s = time.time()
    for df in pd.read_csv(mypath, chunksize=chunksize):
        df['correct'] = df['correct'].astype(float)
        results = df.loc[(df['iteration'] == df['iteration'].max()) & df['paired'] == 1]
        list_of_dataframes.append(results)
    R = pd.concat(list_of_dataframes)
    R = R.loc[(R['iteration'] == R['iteration'].max()) & R['paired'] == 1]
    print(R['iteration'].unique())
    assert len(R['iteration'].unique()) == 1
    try:
        assert R['iteration'].unique()[0] == math.floor(N/6) + 1 # this is only true when no training set used (so not for best for example)
    except:
        assert R['iteration'].unique()[0] < math.floor(N/6) + 1
        warnings.warn('fewer iteations than max', UserWarning)
    print(R.shape)
    print(R['iteration'].unique())
    print('elapsed:', time.time()-s)
    for r in R['n_repeat'].unique():
        assert R.loc[R['n_repeat'] == r]['paired'].sum() == N
    
    return(R)

def load_df_fast(mypath, N):
    # find out which iteration number you want so that you only load useful rows
    its = math.floor(N/6)+1
    print(its)
    list_of_rows= []

    s = time.time()
    # only load the rows I need
    with gzip.open(mypath, mode='rt') as f:
        csvreader = csv.reader(f)
        header = next(csvreader)
        for row in csvreader:
            # print(row)
            if (row[7] == '1') & (row[6]==str(its)):
                list_of_rows.append(row)

    R = pd.DataFrame(list_of_rows, columns=header)
    # print(R.head())
    print(R.shape)
    print(R['iteration'].unique())
    R = R.loc[(R['iteration'] == R['iteration'].max()) & R['paired'] == 1]
    print(R.shape)
    assert R['paired'].astype(int).sum() == N*10 # 10 repeats
    print('elapsed:', time.time()-s)
    
    return(R)

def load_df_TPR(mypath):
    chunksize = 1000000

    # the list that contains all the dataframes
    list_of_dataframes = []
    N = 0
    # i = 0

    for df in pd.read_csv(mypath, chunksize=chunksize, usecols=['iteration', 'correct', 'n_repeat', 'paired']):
        df['correct'] = df['correct'].astype(float)
        df = df.loc[df['paired'] == 1]
        # results = df.groupby(['iteration', 'n_repeat']).sum().reset_index()
        list_of_dataframes.append(df)

    # if you want all the dataframes together, here it is
    result = pd.concat(list_of_dataframes)
    # print(result)
    N = result[['iteration','n_repeat', 'paired']].groupby(by=['iteration','n_repeat']).sum()['paired'].tolist()
    assert len(set(N)) == 1
    N = N[0]
    result = result.drop('paired', axis = 1)
    result = result[['iteration', 'correct', 'n_repeat']].groupby(['iteration', 'n_repeat']).sum().reset_index()
    R = result[['correct', 'n_repeat']].groupby('n_repeat').agg(list).reset_index()
    O = np.array(R['correct'])
    return(O, N)

def load_df_recall(mypath):
    # same function as above, but only keeps unique a/b/ID pairs to calculate number of correct
    chunksize = 1000000

    # the list that contains all the dataframes
    list_of_dataframes = []
    N = 0
    i = 0

    for df in pd.read_csv(mypath, chunksize=chunksize, usecols=['alpha','beta', 'iteration', 'correct', 'n_repeat', 'paired']):
        df['correct'] = df['correct'].astype(float)
        df = df.loc[df['paired'] == 1]
        list_of_dataframes.append(df)

    # if you want all the dataframes together, here it is
    result = pd.concat(list_of_dataframes)
    # print(result)
    N = result[['iteration','n_repeat', 'paired']].groupby(by=['iteration','n_repeat']).sum()['paired'].tolist()
    assert len(set(N)) == 1
    N = N[0]
    # use the line below to find out the unique number of correct pairs you find
    result = result.drop_duplicates()
    result = result.drop('paired', axis = 1)
    result = result[['iteration', 'correct', 'n_repeat']].groupby(['iteration', 'n_repeat']).sum().reset_index()
    R = result[['correct', 'n_repeat']].groupby('n_repeat').agg(list).reset_index()
    O = np.array(R['correct'])
    return(O, N)

def check_mode_correct_IPA(IPApairs, epdf):
    IPApairs['repeat'] = ['repeat' + str(x) for x in IPApairs['n_repeat']]
    IPAres = IPApairs[['alpha', 'alpha_ID', 'beta','repeat']].copy()
    IPAres['key']=IPAres.groupby(['alpha', 'alpha_ID', 'repeat']).cumcount() # I need this key because alpha can occur multiple times in the same ID
    IPAres = IPAres.pivot_table(index=['alpha', 'key', 'alpha_ID'], columns = 'repeat', values='beta',aggfunc='sum').reset_index()
    # print(IPAres.head())
    IPAres = IPAres.rename(columns={'alpha_ID':'subject-PMID'})
    info = epdf[['alpha', 'beta', 'subject-PMID']].copy()
    info['key']=info.groupby(['alpha', 'subject-PMID']).cumcount()
    IPAres = pd.merge(IPAres, info)
    IPAres = IPAres.drop('key', axis = 1)
    IPAres1 = get_mode(IPAres)
    IPAres1['alpha'] = IPAres['alpha']
    IPAres1['beta'] = IPAres['beta']
    IPAres1['subject-PMID'] = IPAres['subject-PMID']
    # print(IPAres1.head())
    check_correct = check_correct_pairs_ID(IPAres1, (epdf['alpha'].str.replace(', ','').str.replace('-','') + '::' \
                    + epdf['beta'].str.replace(', ','').str.replace('-','') + '::' \
                        + epdf['subject-PMID']).tolist(), IPAres1['subject-PMID'].tolist())
    return(check_correct)