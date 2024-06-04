import functions.IPAfunctions as ipa
import functions.IPAscoring as IPAscoring
from typing import Sequence
import pandas as pd
import numpy as np
import functions.myfunctions as mf
from pathlib import Path
from functions.groupAA import translate_aa_into_groups as translate
from concurrent.futures import ProcessPoolExecutor
from sklearn.model_selection import train_test_split
import time
import argparse
from datetime import timedelta
from itertools import chain
import openturns as ot
from scipy.stats import gausshyper, genextreme, weibull_max, johnsonsu
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
        '--weights', default='0.4'
    )
    parser.add_argument(
        '--method', default='None'
    )
    parser.add_argument(
        '--prop_test', default='all'
    )
    parser.add_argument(
        '--confidence', default='none'
    )
    parser.add_argument(
        '--correlation', default='no'
    )
    parser.add_argument(
        '--translation', default='no'
    )
    parser.add_argument(
        '--small_ind', default='no'
    )
    parser.add_argument(
        '--L', default='0.6'
    )
    parser.add_argument(
        '--step', default='3'
    )
    parser.add_argument(
        '--ones', default='keep'
    )
    parser.add_argument(
        '--n_repeats', default='10'
    )
    parser.add_argument(
        '--output', default=''
    )
    parser.add_argument(
        '--save_train', default='no'
    )

    return parser.parse_args(myargs)

# I define the preprocessing steps
def get_inputs(myargs: Sequence[str] | None = None):
    argv = vars(parse_cli_args(myargs))
    print(argv)

    data = Path(argv['input'])
    assert Path.exists(data.resolve()), 'file ' + argv['input'] + ' does not exist'
    folder = Path('data/output/pairing_MI-IPA/' + argv['output'])
    assert Path.exists(folder.resolve()), 'folder ' + 'data/output/pairing_MI-IPA/' + argv['output'] + ' does not exist'
    print('Saving results to ', folder)

    save_train = argv['save_train'].lower()
    assert save_train in ['no'], 'save_train not currently implemented'
    
    epitope = argv['epitope']
    
    weights = argv['weights'].lower()
    try:
        assert (weights == 'no')
    except:
        weights = float(weights)
        assert (0<weights) & (weights<1), 'weights needs to be either no or a float between 0 and 1'
    
    method = argv['method']
    assert method in ['None']
    if method == 'None':
        method = None
    
    prop_test = argv['prop_test'].lower()
    try:
        assert (prop_test == 'all')
    except:
        assert (0<=float(prop_test)) & (float(prop_test)<1), 'prop_test must be either all or a float between 0 and 1'
    
    confidence = argv['confidence'].lower()
    assert confidence in ['hungarian', 'greedy', 'none'], 'confidence must be one of [hungarian, greedy, none]'
    
    correlation = argv['correlation'].lower()
    assert correlation in ['no'], 'only correlation = no implemented'
    
    translation = argv['translation'].lower()
    assert translation in ['yes', 'no'], 'translation must be one of [yes, no]'
    
    ind = argv['small_ind'].lower()
    assert ind in ['yes', 'no']
    
    vgene = argv['Vgene'].lower()
    assert vgene in ['yes', 'no', 'only'], 'Vgene must be one of [yes, no, only]'
    
    step = int(argv['step'])
    
    L = float(argv['L'])
    assert (L>=0) & (L<=1), 'L must be 0<=L<=1'
    
    ones = argv['ones'].lower()
    assert ones in ['throw', 'keep', 'use'], 'ones must be one of [throw, keep, use]'
    
    n_repeats = int(argv['n_repeats'])
    assert n_repeats > 0

    AAs = list('ARNDCEQGHILKMFPSTWYV-')

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
        'n_repeats':n_repeats,
        'folder':folder,
        'save_train':save_train
    }

    return(return_dict)

def _define_individuals(df, epitope, ind):
    if 'Unnamed: 0' in df.columns:
        df = df.drop('Unnamed: 0', axis = 1)
    df = df.replace('https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#', '10xGenomics')
    df = df.replace('https://github.com/antigenomics/vdjdb-db/issues/*', 'GitHubIssue', regex=True)
    epdf = df.loc[(df.Epitope == epitope)].reset_index(drop=True)
    if ind == 'no':
        print('Using individuals as provided')
        subject_id = mf.read_json_column(epdf, 'Meta', 'subject.id').tolist()
        epdf['subject_id'] = subject_id
        epdf['subject-PMID'] = epdf['Reference'].replace(pd.NA, 'unknown') + '_' + epdf['subject_id'].replace(pd.NA, 'unknown')
    else:
        print('Providing fake (simple) ID division')
        individuals = [[x] * 10 for x in range(1000)]
        individuals = list(chain.from_iterable(individuals))[0:epdf.shape[0]]
        epdf['subject-PMID'] = individuals

    epdf = epdf.sort_values(by = 'subject-PMID').reset_index(drop=True)

    return(epdf)

def _generate_seqs_of_interest(epdf, vgene):

    if vgene == 'no':
        epdf1 = mf.prepare_data(epdf, col1='cdr3a_IMGTgaps', col2 = 'cdr3b_IMGTgaps')
        epdf1 = epdf1.rename(columns={'cdr3a_IMGTgaps_padded':'alpha1', 'cdr3b_IMGTgaps_padded':'beta1'})
    elif vgene == 'yes':
        epdf1 = mf.prepare_data(epdf, col1='cdr1a_IMGTgaps', col2 = 'cdr1b_IMGTgaps')
        epdf1 = mf.prepare_data(epdf, col1='cdr2a_IMGTgaps', col2 = 'cdr2b_IMGTgaps')
        epdf1 = mf.prepare_data(epdf, col1='cdr3a_IMGTgaps', col2 = 'cdr3b_IMGTgaps')
        epdf1['alpha1'] = epdf1['cdr1a_IMGTgaps_padded'] + ', ' + epdf1['cdr2a_IMGTgaps_padded'] + ', ' + epdf1['cdr3a_IMGTgaps_padded']
        epdf1['beta1'] = epdf1['cdr1b_IMGTgaps_padded'] + ', ' + epdf1['cdr2b_IMGTgaps_padded'] + ', ' + epdf1['cdr3b_IMGTgaps_padded']
    elif vgene == 'only':
        epdf1 = mf.prepare_data(epdf, col1='cdr1a_IMGTgaps', col2 = 'cdr1b_IMGTgaps')
        epdf1 = mf.prepare_data(epdf, col1='cdr2a_IMGTgaps', col2 = 'cdr2b_IMGTgaps')
        epdf1['alpha1'] = epdf1['cdr1a_IMGTgaps_padded'] + ', ' + epdf1['cdr2a_IMGTgaps_padded']
        epdf1['beta1'] = epdf1['cdr1b_IMGTgaps_padded'] + ', ' + epdf1['cdr2b_IMGTgaps_padded']
    
    epdf2 = epdf1[['alpha1', 'beta1', 'subject-PMID']]

    return(epdf2)

def _remove_constant_res_positions(epdf1):
        # from alpha and beta columns, we want to remove the positions that are constant
    # this is so that we can skip a few positions and optimise performance
    # I need to do this here, as the train/test split might end up different otherwise
    a = np.array(epdf1['alpha1'].str.split(', ', expand=True))
    a = a[:, ~np.all(a[1:] == a[:-1], axis=0)] # remove columns where all values are equal to first value
    epdf1['alpha'] = [', '.join(x) for x in a] # then substitute original column
    b = np.array(epdf1['beta1'].str.split(', ', expand=True))
    b = b[:, ~np.all(b[1:] == b[:-1], axis=0)] # remove columns where all values are equal to first value
    epdf1['beta'] = [', '.join(x) for x in b] # then substitute original column
    return(epdf1)

def _deal_with_small_individuals(test, golden, ones):
    if ones == 'throw':
        print('Throwing away IDs that have only 1 sequence...')
        more_than_one = [x for x in set(test['subject-PMID']) if test['subject-PMID'].value_counts().loc[x] > 1]
        keep = test['subject-PMID'].isin(more_than_one)
        test = test.loc[keep].sort_values(by = 'subject-PMID').reset_index(drop = True)
    elif ones == 'use':
        print('Adding sequences from IDs that have only 1 sequence to golden set...')
        onlyone = [x for x in set(test['subject-PMID']) if test['subject-PMID'].value_counts().loc[x] == 1]
        keep = test['subject-PMID'].isin(onlyone)
        golden = pd.concat([golden, test.loc[keep].sort_values(by = 'subject-PMID').reset_index(drop = True)])
        print('...and removing them from the other set')
        more_than_one = [x for x in set(test['subject-PMID']) if test['subject-PMID'].value_counts().loc[x] > 1]
        keep = test['subject-PMID'].isin(more_than_one)
        test = test.loc[keep].sort_values(by = 'subject-PMID').reset_index(drop = True)
        golden = golden[['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID']]
    elif ones == 'keep':
        test = test.sort_values(by = 'subject-PMID').reset_index(drop = True)
    
    test = test[['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID']]
    return(test, golden)

def _calculate_correlation(test, correlation):
    # currently not tested because not really using
    # legacy code which was going to test the effect of simulating correlation between chains
    # never fully implemented
    if correlation == 'no':
        print('...no correlation calculated as per instructions')
        corr = None
    
    return(corr)


def preprocessing(data, epitope, weights, method, prop_test, confidence, correlation, translation, ind, vgene, step, L, ones, AAs):
    # this takes all outputs from get_input() - it doesn't need all but makes code easier!
    print()
    print('-----------------INITIALISING-----------------')
    print('Running pairing algorithm on epitope ' + epitope +' with the following settings:')
    print('Using V gene? ', vgene)
    print('theta = ', weights)
    print('lambda = ', L)
    print('method = ', method)
    print('confidence scoring = ', confidence)
    print('correlation = ', correlation)
    print('translate = ', translation)
    if prop_test == 'all':
        print('proportion to use as training = 0')
    elif float(prop_test) == 0:
        print('use all as training')
    else:
        print('proportion to use as training = ', 1-float(prop_test))
    print('training step = ', step)

    print()
    print('----------------PREPARING DATA----------------')
    
    df = pd.read_csv(data, index_col = 0)
    epdf = _define_individuals(df, epitope, ind)

    print('Number of entries to be paired = ', epdf.shape[0])
    print('Number of individuals = ', len(set(epdf['subject-PMID'])))
    print('Min number of sequences/id = ', min(epdf['subject-PMID'].value_counts()))
    print('Max number of sequences/id = ', max(epdf['subject-PMID'].value_counts()))
    print('Avg number of sequences/id = ', sum(epdf['subject-PMID'].value_counts())/len(set(epdf['subject-PMID'])))

    epdf1 = _generate_seqs_of_interest(epdf, vgene)
    epdf1 = _remove_constant_res_positions(epdf1)

    # prepare train / test
    if prop_test == 'all':
        golden = pd.DataFrame()
        test = epdf1.copy()
    elif float(prop_test) == 0:
        golden = epdf1.copy()
        test = epdf1.copy()
        step = None
    else:
        golden, test = train_test_split(epdf1, test_size=float(prop_test))

    test, golden = _deal_with_small_individuals(test, golden, ones)

    print('Golden set size = ', golden.shape[0])
    print('Test set size = ', test.shape[0])
    print('Min number of sequences/id in test after clean-up = ', min(test['subject-PMID'].value_counts()))
    print('Max number of sequences/id in test after clean-up = ', max(test['subject-PMID'].value_counts()))
    print('Avg number of sequences/id in test after clean-up = ', sum(test['subject-PMID'].value_counts())/len(set(test['subject-PMID'])))

    print()
    print('--------------PREPARE EXTRA INFO--------------')
    print('(1) subject list...')
    individuals = test['subject-PMID'].tolist()
    individuals = [individuals, individuals]

    print('(2) background correlation...')
    corr = _calculate_correlation(test, correlation)
    
    print('-------------PREPROCESSING DONE---------------')

    return_dict = {
        'weights':weights, 
        'method':method, 
        'confidence':confidence, 
        'correlation':correlation, 
        'corr':corr,
        'translation':translation, 
        'step':step, 
        'L':L, 
        'AAs':AAs,
        'individuals':individuals,
        'golden':golden,
        'test':test
    }

    return(return_dict)

def _create_starting_cache_set(test, golden):
    
    alphas = np.array(test['alpha'].str.split(', ', expand=True))
    betas = np.array(test['beta'].str.split(', ', expand=True))

    if golden.shape[0] > 0:
        startingA = np.array(golden['alpha'].str.split(', ', expand=True))
        startingB = np.array(golden['beta'].str.split(', ', expand=True))
    else:
        startingA = alphas.copy()
        startingB = betas.copy()
        np.random.shuffle(startingB)
    
    return(alphas, betas, startingA, startingB)

def _translation(alphas, betas, startingA, startingB, AAs, translation):

    if translation == 'yes':
        AAs = list('abcdefghijkl-')
        startingA = translate(startingA)
        startingB = translate(startingB)
        alphas = translate(alphas)
        betas = translate(betas)
    
    return(alphas, betas, startingA, startingB, AAs)

def _apply_corr_and_confidence(matrix_0, correlation, corr, confidence, individuals):
    if correlation == 'yes':
        # correlation == 'yes' is not actually implemented
        matrix_0 = ipa.multiply_corr(matrix_0, corr)

    if confidence in ['hungarian', 'greedy']:
        confScores = IPAscoring.IPAscoring(matrix_0, confidence, IDs = individuals)
    else:
        confScores = matrix_0.copy()
    
    return(confScores)

def prepare_info_for_OnePairingRepeat(preprocessed):

    alphas, betas, startingA, startingB = _create_starting_cache_set(preprocessed['test'], preprocessed['golden'])
    alphas, betas, startingA, startingB, AAs = _translation(alphas, betas, startingA, startingB, preprocessed['AAs'], preprocessed['translation'])

    # I should translate these here when translation is on
    correctPairs = (preprocessed['test']['alpha'].str.replace(', ','').str.replace('-','') + '::' \
                    + preprocessed['test']['beta'].str.replace(', ','').str.replace('-','') + '::' \
                        + preprocessed['test']['subject-PMID']).tolist()
    
    if preprocessed['translation'] == 'yes':
        print('translate correct')
        correctPairs = [
            '::'.join([''.join(translate(np.array(list(x.split('::')[0])))),
                       ''.join(translate(np.array(list(x.split('::')[1])))), 
                       x.split('::')[2]]) for x in correctPairs
        ]

    assert len(correctPairs) == alphas.shape[0]
    print(len(correctPairs))

    OnePairingRepeat_vars = ['alphas', 'betas', 'startingA', 'startingB', 'is_golden', 'correctPairs', 
                            'weights', 'method', 'confidence', 'correlation', 'corr', 'step', 'L', 
                            'AAs', 'individuals', 'n_repeat']
    OnePairingRepeat_dict = {k:preprocessed[k] for k in preprocessed.keys() if k in OnePairingRepeat_vars}

    OnePairingRepeat_dict['AAs'] = AAs
    OnePairingRepeat_dict['is_golden'] = 'yes' if preprocessed['golden'].shape[0] >0  else 'no'
    OnePairingRepeat_dict['alphas'] = alphas
    OnePairingRepeat_dict['betas'] = betas
    OnePairingRepeat_dict['startingA'] = startingA
    OnePairingRepeat_dict['startingB'] = startingB
    OnePairingRepeat_dict['correctPairs'] = correctPairs

    return(OnePairingRepeat_dict)

def OnePairingRepeat(alphas, betas, startingA, startingB, is_golden, correctPairs, weights, method, confidence, correlation, corr, step, L, AAs, individuals, n_repeat):
    # I define the function that will allow me to parallelise n_repeats
    print('(3) alpha/beta aligned arrays...')
    start = time.time()
    r = []

    print()
    print('----------------INITIALISE PMI----------------')
    pmi0 = ipa.calculate_PMI_AllPositions(startingA, startingB, method = method, AAs = AAs, L = L, weights=weights)

    print()
    print('--------------STARTING ITERATIONS-------------')
    print('(0) Iteration 0...')
    iteration = 0
    matrix_0 = ipa.create_score_matrix(alphas, betas, individuals, pmi0, AAs)
    confScores = _apply_corr_and_confidence(matrix_0, correlation, corr, confidence, individuals)

    res = ipa.save_results_to_df(confScores, alphas, betas, correctPairs, iteration, n_repeat, individuals)
    r.append(res.dropna(subset = 'confScore'))

    if step: # step is none when I want to iterate only once by training on all
        trainSize = step
        while trainSize < alphas.shape[0]:
            iteration += 1  
            print('(' + str(iteration) + ') Iteration ' + str(iteration) + '...' )
            print('training set size: ', trainSize)
            if is_golden == 'yes':
                matrix_1 = ipa.do_a_training_loop(trainSize, confScores, alphas, betas, IDs = individuals, method=method,
                                                    background_alpha=startingA, background_beta=startingB, 
                                                    AAs = AAs, L = L, weights=weights)
            else:
                matrix_1 = ipa.do_a_training_loop(trainSize, confScores, alphas, betas, IDs = individuals, method=method,
                                                    AAs = AAs, L = L, weights=weights)
            
            confScores = _apply_corr_and_confidence(matrix_1, correlation, corr, confidence, individuals)
            
            res = ipa.save_results_to_df(confScores, alphas, betas, correctPairs, iteration, n_repeat, individuals)
            r.append(res.dropna(subset = 'confScore'))
            trainSize += step

        print()
        print('----------------FINAL ITERATION---------------')
        trainSize = alphas.shape[0]
        iteration += 1  
        if is_golden == 'yes':
            matrix_1 = ipa.do_a_training_loop(trainSize, confScores, alphas, betas, IDs = individuals, method=method,
                                                background_alpha=startingA, background_beta=startingB, 
                                                AAs = AAs, L = L, weights=weights)
        else:
            matrix_1 = ipa.do_a_training_loop(trainSize, confScores, alphas, betas, IDs = individuals, method=method,
                                                AAs = AAs, L = L, weights=weights)
        
        confScores = _apply_corr_and_confidence(matrix_1, correlation, corr, confidence, individuals)

        res = ipa.save_results_to_df(confScores, alphas, betas, correctPairs, iteration, n_repeat, individuals)
        r.append(res.dropna(subset = 'confScore'))
    print('Repeat #', n_repeat, 'finished, total elapsed time = ', str(timedelta(seconds=time.time() - start)))
    R = pd.concat(r)
    print(R)
    return(R)

def save(input_args, r):
    folder = input_args['folder']
    results = pd.concat(r)
    print(results)
    filename = '_'.join([k + '-' + str(input_args['for_preprocessing'][k]) for k in input_args['for_preprocessing'].keys() if (k != 'data') and (k != 'AAs')]) + '_test.csv.gz'
    print(filename)
    results.to_csv(str(folder) + '/' + filename, chunksize = 1e6)

if __name__ == '__main__':
    start = time.time()
    input_args = get_inputs()
    preprocessed = preprocessing(**input_args['for_preprocessing'])

    # Note that I need to call prepare_info_for_OnePairingRepeat within each iteration or they all start with the same initial shuffle
    executor = ProcessPoolExecutor(3) # to control total number of processes when parallelised (big computer has 36, so using 10 in IPAfunctions and 3 here)
    future = [executor.submit(OnePairingRepeat, **prepare_info_for_OnePairingRepeat(preprocessed), n_repeat = repeat) for repeat in range(input_args['n_repeats'])] 
    r = [f.result() for f in future]
    executor.shutdown()

    # r = []
    # for repeat in range(input_args['n_repeats']):
    #     res = OnePairingRepeat(**preprocessed, n_repeat = repeat)
    #     r.append(res)
        
    print()
    print('Pairing finished, elapsed time = ', str(timedelta(seconds=time.time() - start)))
    print('Saving...')
    save(input_args, r)
    