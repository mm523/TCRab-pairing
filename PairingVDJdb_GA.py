import numpy as np
import pandas as pd
import time
from pathlib import Path
import argparse
import functions.pwdistances as pwdist
import warnings
from datetime import timedelta
from typing import Sequence
import PairingVDJdb_MI as MI
currentPath = Path('.').resolve()

# I define the function to get the input args
def parse_cli_args(myargs: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--input',
        help='Path to the location of the source data file.'
    )
    parser.add_argument(
        '--epitope',
    )
    parser.add_argument(
        '--Vgene', default='no'
    )
    parser.add_argument(
        '--small_ind', default='no'
    )
    parser.add_argument(
        '--n_repeats', default='100'
    )
    parser.add_argument(
        '--kNN', default='10'
    )
    parser.add_argument(
        '--output', default=''
    )
    parser.add_argument(
        '--distance_type', default='lev'
    )
    parser.add_argument(
        '--distance_dfa', default=''
    )
    parser.add_argument(
        '--distance_dfb', default=''
    )

    return parser.parse_args(myargs)

# I define the preprocessing steps
def get_inputs(myargs: Sequence[str] | None = None):
    argv = vars(parse_cli_args(myargs))

    data = Path(argv['input'])
    assert Path.exists(data.resolve()), 'file ' + argv['input'] + ' does not exist'
    folder = Path('data/output/pairing_GA/' + argv['output'])
    assert Path.exists(folder.resolve()), 'folder ' + 'data/output/pairing_GA/' + argv['output'] + ' does not exist'
    print('Saving results to ', folder)

    distance_type = argv['distance_type'].lower()
    assert distance_type in ['tcrdist', 'lev', 'weightedlev', 'triplet'], 'distance-type must be one of [triplet, lev, weightedlev, triplet]'
    
    if distance_type == 'tcrdist':
        assert (argv['distance_dfa'] != '') & (argv['distance_dfb'] != ''), 'distance matrix file must be provided when tcrdist is used'
        distanceDF_a = Path(argv['distance_dfa']).resolve()
        assert Path.exists(distanceDF_a.resolve()), 'folder ' + argv['distance_dfa'] + ' does not exist'
        distanceDF_b = Path(argv['distance_dfb']).resolve()
        assert Path.exists(distanceDF_b.resolve()), 'folder ' + argv['distance_dfb'] + ' does not exist'
    else:
        if argv['distance_dfa'] == '': 
            distanceDF_a = None
        else:
            distanceDF_a = None
            warnings.warn('path to precomputed distances are ignored when [distance-type] is not set to tcrdist')
        if argv['distance_dfb'] == '': 
            distanceDF_b = None
        else:
            distanceDF_b = None
            warnings.warn('path to precomputed distances are ignored when [distance-type] is not set to tcrdist')

    epitope = argv['epitope']
    ind = argv['small_ind'].lower()
    assert ind in ['yes', 'no']
    vgene = argv['Vgene'].lower()
    assert vgene in ['yes', 'no', 'only'], 'Vgene must be one of [yes, no, only]'
    n_repeats = int(argv['n_repeats'])
    assert n_repeats > 0
    kNN = int(argv['kNN'])
    assert kNN > 0

    return_dict = {
        'for_preprocessing': {
            'data':data, 
            'epitope':epitope, 
            'ind':ind, 
            'vgene':vgene, 
            },
        'for_GA':{
            'kNN':int(kNN),
            'n_repeats':int(n_repeats),
            'distance_type':distance_type,
            'distanceDF': [distanceDF_a, distanceDF_b]
            },
        'folder':folder,
    }

    return(return_dict)

def preprocessing(data, epitope, ind, vgene):
    print()
    print('-----------------INITIALISING-----------------')
    print('Running graph alignment algorithm on epitope ' + epitope +' with the following settings:')
    print('Using V gene? ', vgene)

    print()
    print('----------------PREPARING DATA----------------')
    vdj = pd.read_csv(data, index_col = 0)
    epdf = MI._define_individuals(vdj, epitope, ind)

    print('Number of entries to be paired = ', epdf.shape[0])
    print('Number of individuals = ', len(set(epdf['subject-PMID'])))
    print('Min number of sequences/id = ', min(epdf['subject-PMID'].value_counts()))
    print('Max number of sequences/id = ', max(epdf['subject-PMID'].value_counts()))
    print('Avg number of sequences/id = ', sum(epdf['subject-PMID'].value_counts())/len(set(epdf['subject-PMID'])))

    epdf = MI._generate_seqs_of_interest(epdf, vgene)
    epdf['alpha'] = epdf['alpha1'].str.replace(', ', '').str.replace('-','')
    epdf['beta'] = epdf['beta1'].str.replace(', ', '').str.replace('-','')
    # should I have a test vs golden split here? 
    # It's meaningless for GA, but might be important for integration between the two?
    test = epdf.sort_values(by = 'subject-PMID').reset_index(drop = True)
    test = test[['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID']]

    print('Min number of sequences/id in test after clean-up = ', min(test['subject-PMID'].value_counts()))
    print('Max number of sequences/id in test after clean-up = ', max(test['subject-PMID'].value_counts()))
    print('Avg number of sequences/id in test after clean-up = ', sum(test['subject-PMID'].value_counts())/len(set(test['subject-PMID'])))

    print()
    print('--------------PREPARE EXTRA INFO--------------')
    print('(1) subject list...')
    individuals = test['subject-PMID'].tolist()

    print('-------------PREPROCESSING DONE---------------')

    return_dict = {
        'individuals':individuals,
        'test':test.sort_values(by = 'subject-PMID').reset_index(drop = True)
    }

    return(return_dict)

def get_id_info(test):
    N = len(test['subject-PMID'].unique())
    M = test.shape[0]
    FirstSeqSpec = []
    LastSeqSpec = []
    MSeqSpec = []

    for s in test['subject-PMID'].unique():
        sid = test.loc[test['subject-PMID'] == s].index
        FirstSeqSpec.append(sid.min())
        LastSeqSpec.append(sid.max())
        MSeqSpec.append(len(sid))

    FirstSeqSpec = np.array(FirstSeqSpec)+1 # add 1 because Julia is not 0-indexed
    LastSeqSpec = np.array(LastSeqSpec)+1
    MSeqSpec = np.array(MSeqSpec)

    index_map = dict(zip(test['subject-PMID'].unique(), range(1, len(test['subject-PMID'].tolist())+1)))
    IndexSeqSpec = [index_map[x] for x in test['subject-PMID'].tolist()]
    assert IndexSeqSpec == sorted(IndexSeqSpec) # because subject-PMID should be sorted BEFORE you do this
    IndexSeqSpec = np.array(IndexSeqSpec)
    
    return(M,N,FirstSeqSpec,LastSeqSpec,MSeqSpec,IndexSeqSpec)

def _reorder_distance_frame(df, seqs):
    # remove duplicates
    # print(df)
    assert sorted(seqs) == sorted(df.index.tolist())
    df1 = df.reset_index(names='seqs').drop_duplicates()
    assert sorted(list(set(seqs))) == sorted(df1['seqs'].tolist())
    # print(df)
    # print('number of NANs:', df1.isnull().sum().sum())
    
    # assign order to sequences
    order = pd.DataFrame(seqs, columns=['seqs'])
    order = order.reset_index(names='order')

    # reorder x-axis
    X = pd.merge(df1, order, how='outer').drop_duplicates()
    # print(X[X.duplicated()])
    # print(X['seqs'].tolist())
    # print(X)
    # print('number of NANs:', X.isnull().sum().sum())
    X.index = X['order'].tolist()
    X = X.drop(['order','seqs'], axis=1).sort_index()
    # print(X)
    # print('number of NANs:', X.isnull().sum().sum())
    assert X.isnull().sum().sum() == 0
    # print('x-axis reordered')
    # transpose and reorder y-axis
    X1 = X.T
    X1 = X1.reset_index(names='seqs')
    X1 = pd.merge(X1, order, how='right').drop_duplicates()
    assert X1.isnull().sum().sum() == 0
    X1.index = X1['order'].tolist()
    # print(X1)
    X1 = X1.drop(['order','seqs'], axis=1).sort_index()
    # print(X1)
    # transpose and make array
    X2 = np.array(X1.T.fillna(-10).values)
    # print(X2)
    # print(np.where(X2==-10))
    # print(X2[208:211,208:211])
    return(X2)

def get_distances(seqs, distance_type, distanceDF):
    print()
    print('(3) get alpha and beta...')
    print('distance type: ', distance_type)
    As, Bs = seqs

    if distance_type == 'tcrdist':
        print('load distances...')
        dij_A = pd.read_csv(distanceDF[0], index_col=0)
        dij_B = pd.read_csv(distanceDF[1], index_col=0)
        print(len(As), dij_A.shape)
        assert dij_A.shape == (len(As), len(As))
        assert dij_B.shape == (len(Bs), len(Bs))
        # now reorder the df

        dij_A = _reorder_distance_frame(dij_A, As)
        print(dij_A.shape)
        dij_B = _reorder_distance_frame(dij_B, Bs)
        assert dij_A.shape == (len(As), len(As))
        assert dij_B.shape == (len(Bs), len(Bs))
        
    elif distance_type == 'lev':
        dij_A = pwdist.levenshtein_dist(As)
        dij_B = pwdist.levenshtein_dist(Bs)
    elif distance_type == 'weightedlev':
        dij_A = pwdist.weighted_levenshtein_dist(As)
        dij_B = pwdist.weighted_levenshtein_dist(Bs)
    elif distance_type == 'triplet':
        dij_A = pwdist.triplet_diversity(As).round(1)*10
        dij_B = pwdist.triplet_diversity(Bs).round(1)*10
    
    # I return the sequences back because if I have calculated distances outside of the script,
    # I want to make sure that the seqs I am replacing later are correct
    return(dij_A, dij_B)

def one_kNN_run(k, n_replicates, M, N, 
                FirstSeqSpec, LastSeqSpec, MSeqSpec, IndexSeqSpec, 
                dij_A, dij_B):
    # not testing this as not sure how to write a test for the GA...
    # it turns out you need to load julia on all processes for this to work: 
    # https://stackoverflow.com/questions/74438358/error-while-targeting-a-julia-function-into-multiprocessing-process-of-python
    
    # this first two lines are a workaround to make julia work in env
    print('loading packages...')
    from julia.api import Julia
    jl = Julia(compiled_modules=False)
    from julia import Main
    Main.include("./GraphAlignment/RunGraphAlignment.jl")
    print('all packages loaded, starting...')
    s = time.time()
    output = Main.robust_GA_SA(k, n_replicates, M, N,
                                FirstSeqSpec, LastSeqSpec, MSeqSpec, IndexSeqSpec,
                                dij_A,dij_B)

    print('k = ', k, ' completed, elapsed: ', time.time()-s)
    return(output)

def parse_output(OutputArray):
    num_seqs = OutputArray[0,:]
    assert len(set(num_seqs)) == 1
    seqs_per_species = OutputArray[1,:]
    assert len(set(seqs_per_species)) == 1
    true_pos_repeat = OutputArray[2,:]
    energy_final_match = OutputArray[3,:]
    final_match_seq = OutputArray[4:,:]

    return(num_seqs[0],seqs_per_species[0],true_pos_repeat,energy_final_match,final_match_seq)

def replace_with_seq(ParsedOutputArray, alphas, betas, individuals):
    df = pd.DataFrame(ParsedOutputArray, columns=['repeat' + str(x) for x in range(1, ParsedOutputArray.shape[1] + 1)])
    # print(df)
    mapping = dict(zip(range(1,len(betas)+1), betas))
    df = df.replace(mapping)
    df['alpha'] = alphas
    df['correct_beta'] = betas
    df['subject-PMID'] = individuals    
    print(df)
    return(df)

def get_mode(GAoutput):
    mode = []
    mode_freq = []
    for i in range(GAoutput.shape[0]):
        m = GAoutput.iloc[i][[c for c in GAoutput.columns if 'repeat' in c]].value_counts()
        mmode = m.loc[m == m.max()]
        if mmode.shape[0] == 1:
            cdr3, v = mmode.index[0], mmode.values[0]
        else:
            # pick a random one if there is a tie
            mmode1 = pd.DataFrame(mmode.sample(1))
            # print(mmode1)
            # print(mmode1.index, mmode1.values)
            cdr3, v = mmode1.index[0], mmode1.values[0][0]
            # print(cdr3, v)
        mode_freq.append(v)
        mode.append(cdr3)
    
    GAoutput['mode'] = mode
    GAoutput['mode_freq'] = mode_freq

    return(GAoutput)

def check_correct_pairs_ID(GAoutput, correctPairs, ids = None):

    if ids:
        pairs = (GAoutput['alpha'].str.replace(', ','').str.replace('-','') + '::' \
                    + GAoutput['mode'].str.replace(', ','').str.replace('-','') + '::' \
                        + GAoutput['subject-PMID']).tolist()
    else:
        pairs = (GAoutput['alpha'].str.replace(', ','').str.replace('-','') + '::' \
                    + GAoutput['mode'].str.replace(', ','').str.replace('-','')).tolist()
    # print(pairs)
    GAoutput['correct'] = [p in correctPairs for p in pairs]
    return(GAoutput)

def save_results_to_df(GAoutput, correctpairs, individuals):
    
    GAoutput = get_mode(GAoutput)
    GAoutput = check_correct_pairs_ID(GAoutput, correctpairs, individuals)
    num_correct_pairs = sum([int(x) for x in GAoutput['correct']])
    print('Number of correct pairs: ', num_correct_pairs, ' of ', GAoutput.shape[0])

    return(GAoutput)

def GraphAlignment(test, individuals, kNN, n_repeats, distance_type, distanceDF):
    # I define the function that will allow me to parallelise n_repeats
    print('(2) get all information about IDs to feed into Julia...')

    M,N,FirstSeqSpec,LastSeqSpec,MSeqSpec,IndexSeqSpec = get_id_info(test)    
    
    print()
    print('(4) getting distances...')
    As = test['alpha'].tolist()
    Bs = test['beta'].tolist()
    # print(As)
    # print(Bs)
    correctPairs = (test['alpha'].str.replace(', ','').str.replace('-','') + '::' \
                    + test['beta'].str.replace(', ','').str.replace('-','') + '::' \
                        + test['subject-PMID']).tolist()
    # print(correctPairs)
    seqs = [As, Bs]
    dij_A, dij_B = get_distances(seqs, distance_type, distanceDF)
    # A1s, B1s = seqs 
    # assert As == A1s # to make sure the right sequences are replaced later
    # assert Bs == B1s
    # print(As)
    dij_A = dij_A.astype('int64')
    dij_B = dij_B.astype('int64')
    print(dij_A.min(), dij_A.max())
    assert dij_A.max() > 1
    assert dij_B.max() > 1
    assert dij_A.min() >= 0
    assert dij_B.min() >= 0

    print()
    print('(5) run graph alignment...')
    print('K = ', kNN)
    results = one_kNN_run(kNN, n_repeats, M, N, 
                          FirstSeqSpec, LastSeqSpec, MSeqSpec, IndexSeqSpec, 
                          dij_A, dij_B)
    
    num_seqs,seqs_per_species,true_pos_repeat,energy_final_match,final_match_seq = parse_output(results)
    print(final_match_seq)
    print(As)
    res = replace_with_seq(final_match_seq, As, Bs, individuals)
    print(res.head())
    # print(res)

    res = save_results_to_df(res, correctPairs, individuals)
    return(res)

if __name__ == "__main__":
    start = time.time()
    input_args = get_inputs()
    preprocessed = preprocessing(**input_args['for_preprocessing'])

    res = GraphAlignment(**preprocessed, **input_args['for_GA'])
        
    print()
    print('Pairing finished, elapsed time = ', str(timedelta(seconds=time.time() - start)))
    print('Saving...')
    mylist = ['epitope', 'ind', 'vgene']
    part1 = [k + '-' + str(input_args['for_preprocessing'][k]) for k in mylist]
    mylist1 = ['distance_type', 'kNN', 'n_repeats']
    part2 = [k + '-' + str(input_args['for_GA'][k]) for k in mylist1]
    fileinfo = part1+part2
    filename = '_'.join(fileinfo) + '.csv.gz'
    print(filename)
    res.to_csv(str(input_args['folder']) + '/GA-' + filename)
    