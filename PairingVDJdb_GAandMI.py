'''
This script combines pairing with GA first, 
then getting the more stable results and using them as a golden set for the MI-IPA
'''

import PairingVDJdb_GA as GA
import PairingVDJdb_MI as MI
import argparse
from pathlib import Path
import warnings
from typing import Sequence
import pandas as pd
from time import time
from concurrent.futures import ProcessPoolExecutor
from functions.groupAA import translate_aa_into_groups as translate
from warnings import warn
import numpy as np

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
        '--prop_test', default='all'
    )
    parser.add_argument(
        '--ones', default='keep'
    )
    parser.add_argument(
        '--output', default=''
    )

    # GA-related arguments
    parser.add_argument(
        '--GA_repeats', default='100'
    )
    parser.add_argument(
        '--GA_kNN', default='10'
    )
    parser.add_argument(
        '--GA_thresh', default='0.95'
    )
    parser.add_argument(
        '--GA_precomputed', default=None
    )
    parser.add_argument(
        '--GA_distance_type', default='lev'
    )
    parser.add_argument(
        '--GA_distance_dfa', default=''
    )
    parser.add_argument(
        '--GA_distance_dfb', default=''
    )
    # MI-related arguments
    parser.add_argument(
        '--MI_repair', default='yes'
    )
    parser.add_argument(
        '--MI_weights', default='0.6'
    )
    parser.add_argument(
        '--MI_method', default='None'
    )
    parser.add_argument(
        '--MI_confidence', default='none'
    )
    parser.add_argument(
        '--MI_correlation', default='no'
    )
    parser.add_argument(
        '--MI_translation', default='no'
    )
    # note sequences are only translated for MI
    parser.add_argument(
        '--MI_L', default='0.6'
    )
    parser.add_argument(
        '--MI_step', default='6'
    )
    parser.add_argument(
        '--MI_repeats', default='1'
    )
    parser.add_argument(
        '--MI_output', default=''
    )

    return parser.parse_args(myargs)

def get_inputs(myargs: Sequence[str] | None = None):
    argv = vars(parse_cli_args(myargs))

    data = Path(argv['input'])
    assert Path.exists(data.resolve()), 'file ' + argv['input'] + ' does not exist'
    folder = Path('data/output/pairing_GAMI/' + argv['output'])
    assert Path.exists(folder.resolve()), 'folder ' + 'data/output/pairing_GAMI/' + argv['output'] + ' does not exist'
    print('Saving results to ', folder)

    epitope = argv['epitope']

    # GA-related
    GAprecomputed = argv['GA_precomputed']
    try:
        assert GAprecomputed == None, 'Pre-computed GA is provided'
        # this is when GA is not provided so we need to make a new one
        myGA = None
        
        distance_type = argv['GA_distance_type'].lower()
        assert distance_type in ['tcrdist', 'lev', 'weightedlev', 'triplet'], 'distance-type must be one of [tcrdist, lev, weightedlev, triplet]'
        
        kNN = int(argv['GA_kNN'])
        assert kNN > 0
        
        if distance_type == 'tcrdist':
            assert (argv['GA_distance_dfa'] != '') & (argv['GA_distance_dfb'] != ''), 'distance matrix file must be provided when tcrdist is used'
            distanceDF_a = Path(argv['GA_distance_dfa']).resolve()
            assert Path.exists(distanceDF_a.resolve()), 'folder ' + argv['GA_distance_dfa'] + ' does not exist'
            distanceDF_b = Path(argv['GA_distance_dfb']).resolve()
            assert Path.exists(distanceDF_b.resolve()), 'folder ' + argv['GA_distance_dfb'] + ' does not exist'
        else:
            if argv['GA_distance_dfa'] == '': 
                distanceDF_a = None
            else:
                distanceDF_a = None
                warnings.warn('path to precomputed distances are ignored when [distance-type] is not set to tcrdist')
            if argv['GA_distance_dfb'] == '': 
                distanceDF_b = None
            else:
                distanceDF_b = None
                warnings.warn('path to precomputed distances are ignored when [distance-type] is not set to tcrdist')

    except AssertionError as m:
        print(m)
        if str(m) == 'Pre-computed GA is provided':
            # load a precomputed GA to use results of
            print('Path: ', GAprecomputed)
            GAprecomputed = Path(GAprecomputed).resolve()
            assert Path.exists(GAprecomputed.resolve()), 'Path to precomputed GA does not exist'
            distance_type = kNN = ''
            distanceDF_a = distanceDF_b = None
        else:
            raise

    if GAprecomputed:
        GArepeats=None
    else:   
        GArepeats = int(argv['GA_repeats'])
        assert GArepeats > 0
    
    GAthresh = float(argv['GA_thresh'])
    assert (GAthresh <= 1) & (GAthresh >= 0)

    # MI-related
    repair = argv['MI_repair'].lower()
    assert repair in ['yes', 'no']
    if repair == 'yes':
        repair = True
    else:
        repair = False
    weights = argv['MI_weights'].lower()
    try:
        assert (weights == 'no')
    except:
        weights = float(weights)
        assert (0<weights) & (weights<1), 'weights needs to be either no or a float between 0 and 1'

    method = argv['MI_method']
    assert method in ['None']
    if method == 'None':
        method = None
    
    prop_test = argv['prop_test'].lower()
    try:
        assert (prop_test == 'all')
    except:
        assert (0<float(prop_test)) & (float(prop_test)<1), 'prop_test must be either all or a float between 0 and 1'
 
    confidence = argv['MI_confidence'].lower()
    assert confidence in ['hungarian', 'greedy', 'none'], 'confidence must be one of [hungarian, greedy, none]'
    
    correlation = argv['MI_correlation'].lower()
    assert correlation in ['no'], 'only correlation = no implemented'
    
    translation = argv['MI_translation'].lower()
    assert translation in ['yes', 'no'], 'translation must be one of [yes, no]'
    
    ind = argv['small_ind'].lower()
    assert ind in ['yes', 'no']

    vgene = argv['Vgene'].lower()
    assert vgene in ['yes', 'no', 'only'], 'Vgene must be one of [yes, no, only]'
    
    step = int(argv['MI_step'])
    
    L = float(argv['MI_L'])
    assert (L>=0) & (L<=1), 'L must be 0<=L<=1'
    
    ones = argv['ones'].lower()
    assert ones in ['throw', 'keep', 'use'], 'ones must be one of [throw, keep, use]'
    
    MIrepeats = int(argv['MI_repeats'])
    assert MIrepeats > 0

    AAs = list("ARNDCEQGHILKMFPSTWYV-")

    return_dict = {
        'for_preprocessing': {
            'data':data, 
            'epitope':epitope, 
            'weights':weights, 
            'L':L,
            'method':method, 
            'prop_test':prop_test, 
            'confidence':confidence, 
            'correlation':correlation, 
            'translation':translation, 
            'ind':ind, 
            'vgene':vgene, 
            'step':step, 
            'ones':ones, 
            'AAs':AAs},
        'for_MI':{
            'n_repeats':MIrepeats
        },
        'for_GA':{
            'kNN':kNN,
            'n_repeats':GArepeats,
            'distance_type':distance_type,
            'distanceDF': [distanceDF_a, distanceDF_b]
            },
        'repair':repair,
        'GA_precomputed':GAprecomputed,
        'GA_thresh':GAthresh,
        'folder':folder
    }

    return(return_dict)

def _derive_GA_set_from_MI_preprocessing(df):
    test = df.copy()
    test['alpha'] = test['alpha1'].str.replace('-','').str.replace(', ','')
    test['beta'] = test['beta1'].str.replace('-','').str.replace(', ','')
    test=test.sort_values(by = 'subject-PMID').reset_index(drop=True)

    return(test)

def _get_stable_GA(GAoutput, GAthresh, GA_repeats):
    # now we want to get the most stable repeats
    # first, try to get what is stable > thresh
    GAoutput['mode_freq'] = GAoutput['mode_freq'].astype(int)
    GA_stable = GAoutput.loc[GAoutput['mode_freq'] >= GAthresh*GA_repeats].reset_index(drop=True)
    # print(GA_stable)
    # if no stable ones found, sort by stability and take top 5
    if GA_stable.shape[0] < 5:
        warn('<3 pairs were stable enough, picking the top 5')
        GA_stable = GAoutput.sort_values(by = 'mode_freq', ascending = False).head().reset_index(drop=True)
    print('correct in GA above threshold')
    print(GA_stable['correct'].value_counts()) # I can do this because I know the answers

    return(GA_stable[['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct']])

def _make_new_golden_set(GAstable, MItest):

    # sequence map to change the seqs back to have spaces
    MItest1 = MItest.copy()
    MItest1['alpha1'] = MItest1['alpha1'].str.replace('-','').str.replace(', ','')
    MItest1['beta1'] = MItest1['beta1'].str.replace('-','').str.replace(', ','')
    seq_maps_a = dict(zip(MItest1['alpha1'], MItest1['alpha']))
    seq_maps_b = dict(zip(MItest1['beta1'], MItest1['beta']))

    golden2 = GAstable.copy().rename(columns={'alpha':'alpha1', 'mode':'beta1'})
    golden2['alpha'] = golden2['alpha1'].apply(lambda x: seq_maps_a[x])
    golden2['beta'] = golden2['beta1'].apply(lambda x: seq_maps_b[x])
    golden2 = golden2[['alpha1', 'alpha', 'beta1', 'beta', 'subject-PMID', 'mode_freq', 'correct', 'correct_beta']]
    # print(golden2)

    return(golden2)

def _prepare_new_test_set(dict_MI_run, golden2, translation, MI_repair):
    startingA = np.array(golden2['alpha'].str.split(', ', expand=True))
    startingB = np.array(golden2['beta'].str.split(', ', expand=True))
    if translation == 'yes':
        startingA = translate(startingA)
        startingB = translate(startingB)

    current_alphas = dict_MI_run['alphas']
    current_betas = dict_MI_run['betas']
    to_remove_a = []
    to_remove_b = []
    if MI_repair == False:
        for j in range(len(startingA)):
            matches = [i for i in range(len(current_alphas)) if (current_alphas[i]==startingA[j]).all()]
            matches1 = [i for i in range(len(dict_MI_run['individuals'][0])) if dict_MI_run['individuals'][0][i] == golden2.iloc[j]['subject-PMID']]
            intersect = [x for x in set(matches1) & set(matches)]
            # when a sequence is found multiple times, I add a new iteration of it each time
            intersect1 = [x for x in intersect if x not in to_remove_a]
            newidx = intersect1[0]
            to_remove_a.append(newidx)  
        for j in range(len(startingB)):
            matches = [i for i in range(len(current_betas)) if (current_betas[i]==startingB[j]).all()]
            matches1 = [i for i in range(len(dict_MI_run['individuals'][1])) if dict_MI_run['individuals'][1][i] == golden2.iloc[j]['subject-PMID']]
            intersect = [x for x in set(matches1) & set(matches)]
            # when a sequence is found multiple times, I add a new iteration of it each time
            intersect1 = [x for x in intersect if x not in to_remove_b]
            newidx = intersect1[0]
            to_remove_b.append(newidx)

        to_remove_a = list(set(to_remove_a))
        to_remove_b = list(set(to_remove_b))
        assert len(to_remove_a) == startingA.shape[0]
        assert len(to_remove_b) == startingB.shape[0]

    dict_MI_run['alphas'] = np.delete(current_alphas, to_remove_a, axis=0)
    dict_MI_run['betas'] = np.delete(current_betas, to_remove_b, axis=0)
    # print(dict_MI_run['alphas'].shape)
    # print(dict_MI_run['betas'].shape)
    dict_MI_run['individuals'] = [np.delete(np.array(dict_MI_run['individuals'][0]), to_remove_a).tolist(), 
                                    np.delete(np.array(dict_MI_run['individuals'][1]), to_remove_b).tolist()]

    # finally concatenate the new golden set to the old golden set
    if dict_MI_run['is_golden'] == 'yes':
        # there is a previous training set, so we concatenate
        dict_MI_run['startingA'] = np.concatenate((dict_MI_run['startingA'], startingA), axis=0)
        dict_MI_run['startingB'] = np.concatenate((dict_MI_run['startingB'], startingB), axis=0)
    else: # the training set is just a random shuffle, so we overwrite
        dict_MI_run['startingA'] = startingA
        dict_MI_run['startingB'] = startingB
        dict_MI_run['is_golden'] = 'yes'

    return(dict_MI_run)


if __name__ == "__main__":
    start = time()
    input_args = get_inputs()
    MIpreprocessed = MI.preprocessing(**input_args['for_preprocessing'])
    # print(test)
    if input_args['GA_precomputed'] is None:
        GAtest = _derive_GA_set_from_MI_preprocessing(MIpreprocessed['test'])
        GApreprocessed = {
            'individuals':GAtest['subject-PMID'].tolist(),
            'test':GAtest
        }
        GAres = GA.GraphAlignment(**GApreprocessed, **input_args['for_GA'])
        mylist = ['epitope', 'prop_test', 'ind', 'vgene', 'ones']
        part1 = [k + '-' + str(input_args['for_preprocessing'][k]) for k in mylist]
        mylist1 = ['distance_type', 'kNN', 'n_repeats']
        part2 = [k + '-' + str(input_args['for_GA'][k]) for k in mylist1]
        fileinfo = part1+part2
        filename = '_'.join(fileinfo) + '.csv.gz'
        print(filename)
        GAres.to_csv(str(input_args['folder']) + '/GA-' + filename)
    else:
        GAres = pd.read_csv(input_args['GA_precomputed'], index_col=0).copy()
        input_args['for_GA']['n_repeats'] = len([x for x in GAres.columns if "repeat" in x])

    GA_stable = _get_stable_GA(GAres, input_args['GA_thresh'],input_args['for_GA']['n_repeats'])

    golden2 = _make_new_golden_set(GA_stable, MIpreprocessed['test'])

    dict_MI_run = MI.prepare_info_for_OnePairingRepeat(MIpreprocessed)
    # print(dict_MI_run)

    # I now need to replace startingA, startingB, alphas and betas.
    # I also need to remove the corresponding individuals.
    # the list of correct pairs will remain the same (I don't think there is any point in updating)
    new_dict_MI_run = _prepare_new_test_set(dict_MI_run, golden2, translation = MIpreprocessed['translation'], MI_repair=input_args['repair'])

    if input_args['for_MI']['n_repeats']>1:
        executor = ProcessPoolExecutor(3) # to control total number of processes when parallelised (big computer has 36, so using 10 in IPAfunctions and 3 here)
        future = [executor.submit(MI.OnePairingRepeat, **new_dict_MI_run, n_repeat = repeat) for repeat in range(input_args['for_MI']['n_repeats'])] 
        r = [f.result() for f in future]
        executor.shutdown()
        results = pd.concat(r)
    else:
        results = MI.OnePairingRepeat(**new_dict_MI_run, n_repeat = 1)

    folder = input_args['folder']
    print(results)
    filename = '_'.join([k + '-' + str(input_args['for_preprocessing'][k]) for k in input_args['for_preprocessing'].keys() if (k != 'data') and (k != 'AAs')] + 
                        ['GAthresh-' + str(input_args['GA_thresh'])] + ['repairing-' + str(input_args['repair'])] + 
                        [k + '-' + str(input_args['for_GA'][k]) for k in input_args['for_GA'].keys() if (k != 'distanceDF') and (k != 'n_repeats')] + 
                        ['GA_repeats-' + str(input_args['for_GA']['n_repeats'])] + 
                        ['MI_repeats-' + str(input_args['for_MI']['n_repeats'])]) + '.csv.gz'
    print(filename)
    results.to_csv(str(folder) + '/' + filename, chunksize = 1e6)