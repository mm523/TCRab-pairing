import PairingVDJdb_MI as MI
from pathlib import Path
import pytest
import copy
import pandas as pd
import numpy as np


class TestInputs:
    def test_get_inputs(self):
        output = {
            'for_preprocessing': {
                'data':Path('data/vdj_cleaned_subset_for_MI.csv'), 
                'epitope':'myep', 
                'weights':0.4, 
                'L':0.6,
                'method':None, 
                'prop_test':'all', 
                'confidence':'none', 
                'correlation':'no', 
                'translation':'no', 
                'ind':'no', 
                'vgene':'no', 
                'step':3, 
                'ones':'keep', 
                'AAs':list('ARNDCEQGHILKMFPSTWYV-')},
            'n_repeats':10,
            'folder':Path('data/output/pairing_MI-IPA/'),
            'save_train':'no'
            }
        
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep']) == output

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['weights'] = 0.3
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--weights', '0.3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['prop_test'] = '0.75'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--prop_test', '0.75']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['confidence'] = 'hungarian'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--confidence', 'HUNGARIAN']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['correlation'] = 'no'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--correlation', 'NO']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['translation'] = 'yes'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--translation', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['ind'] = 'yes'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--small_ind', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['vgene'] = 'yes'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--Vgene', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['step'] = 6
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--step', '6']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['L'] = 0.5
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--L', '0.5']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['ones'] = 'throw'
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--ones', 'THROW']) == output1

        output1 = copy.deepcopy(output)
        output1['n_repeats'] = 3
        assert MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--n_repeats', '3']) == output1



        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/aaa.csv', '--epitope', 'myep'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','aaa'])
        
        with pytest.raises(ValueError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--weights','yes'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--weights','2'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--method','sthg'])
        
        with pytest.raises(ValueError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--prop_test','no'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--prop_test','3'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--confidence','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--correlation','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--translation','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--small_ind','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--Vgene','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--L','3'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--ones','STHG'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--n_repeats','0'])
        
        with pytest.raises(AssertionError):
            MI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','./crazypath'])
        
class TestPreprocessing:
    def test_define_individuals(self):
        inputdf = pd.read_csv('tests/pairing_vdjdb/inputdf_test.csv', index_col=0)
        outputdf1 = pd.read_csv('tests/pairing_vdjdb/outputdf1_test.csv', index_col=0)
        pd.testing.assert_frame_equal(MI._define_individuals(inputdf, 'ELAGIGILTV', 'no'), outputdf1)

        # test small ind
        inputdf1 = pd.concat([inputdf, inputdf, inputdf, inputdf]).sort_values(by = 'complex.id').reset_index(drop=True)
        outputdf2 = pd.concat([outputdf1, outputdf1, outputdf1, outputdf1]).sort_values(by = 'complex.id').reset_index(drop=True)
        outputdf2['subject-PMID'] = [0,0,0,0,0,0,0,0,0,0,1,1]
        outputdf2 = outputdf2.drop('subject_id', axis=1)
        pd.testing.assert_frame_equal(MI._define_individuals(inputdf1, 'ELAGIGILTV', 'yes'), outputdf2)

    def test_generate_seqs_of_interest(self):
        # mf.prepare_data is tested within myfunctions.py
        # using the same df to test
        df = pd.DataFrame.from_dict(
             {'cdr1a_IMGTgaps':'CaS---VTF',
              'cdr2a_IMGTgaps':'CAs---VTF',
              'cdr3a_IMGTgaps':'CAS---VTF',
              'cdr1b_IMGTgaps':'CaST---VTF',
              'cdr2b_IMGTgaps':'CAsT---VTF',
              'cdr3b_IMGTgaps':'CAST---VTF',
              'subject-PMID':1},
              orient = 'index').T

        # if not set, max_len is the max for each column
        output1 = pd.DataFrame.from_dict(
                   {'alpha1':'C, A, S, -, -, -, V, T, F', 
                    'beta1':'C, A, S, T, -, -, -, V, T, F',
                    'subject-PMID':1},
                    orient = 'index').T
        
        pd.testing.assert_frame_equal(MI._generate_seqs_of_interest(df, 'no'), output1)

        output2 = pd.DataFrame.from_dict(
                   {'alpha1':'C, a, S, -, -, -, V, T, F, C, A, s, -, -, -, V, T, F, C, A, S, -, -, -, V, T, F', 
                    'beta1':'C, a, S, T, -, -, -, V, T, F, C, A, s, T, -, -, -, V, T, F, C, A, S, T, -, -, -, V, T, F',
                    'subject-PMID':1},
                    orient = 'index').T
        
        pd.testing.assert_frame_equal(MI._generate_seqs_of_interest(df, 'yes'), output2)

        output3 = pd.DataFrame.from_dict(
                   {'alpha1':'C, a, S, -, -, -, V, T, F, C, A, s, -, -, -, V, T, F', 
                    'beta1':'C, a, S, T, -, -, -, V, T, F, C, A, s, T, -, -, -, V, T, F',
                    'subject-PMID':1},
                    orient = 'index').T
        
        pd.testing.assert_frame_equal(MI._generate_seqs_of_interest(df, 'only'), output3)
        
    def test_remove_constant_res_positions(self):
        df = pd.DataFrame.from_dict(
             {'alpha1':['A, A, A, X, X, X, A, A',
                        'B, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'beta1': ['A, A, A, X, X, X, A, A',
                        'A, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'subject-PMID':[1,1,1]},
              orient = 'index').T

        output = pd.DataFrame.from_dict(
             {'alpha1':['A, A, A, X, X, X, A, A',
                        'B, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'beta1': ['A, A, A, X, X, X, A, A',
                        'A, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'subject-PMID':[1,1,1],
              'alpha': ['A, A, A, A, A',
                        'B, B, B, B, A', 
                        'A, A, A, A, B'],
              'beta': ['A, A, A, A',
                        'B, B, B, A', 
                        'A, A, A, B']
                        },
              orient = 'index').T
        print(output)
        pd.testing.assert_frame_equal(MI._remove_constant_res_positions(df), output)
    
    def test_deal_with_small_individuals(self):
        df = pd.DataFrame.from_dict(
             {'alpha': ['A, A, A, A, A',
                        'B, B, B, B, A', 
                        'A, A, A, A, B'],
              'alpha1':['A, A, A, X, X, X, A, A',
                        'B, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'beta': ['A, A, A, A',
                        'B, B, B, A', 
                        'A, A, A, B'],
              'beta1': ['A, A, A, X, X, X, A, A',
                        'A, B, B, X, X, X, B, A', 
                        'A, A, A, X, X, X, A, B'],
              'subject-PMID':[1,2,1]},
              orient = 'index').T
        golden = pd.DataFrame()
        
        output1 = pd.DataFrame.from_dict(
             {'alpha': ['A, A, A, A, A',
                        'A, A, A, A, B'],
              'alpha1':['A, A, A, X, X, X, A, A',
                        'A, A, A, X, X, X, A, B'],
              'beta': ['A, A, A, A',
                        'A, A, A, B'],
              'beta1': ['A, A, A, X, X, X, A, A',
                        'A, A, A, X, X, X, A, B'],
              'subject-PMID':[1,1],
                        },
              orient = 'index').T
        t, g = MI._deal_with_small_individuals(df, golden, 'throw')
        pd.testing.assert_frame_equal(t, output1)
        pd.testing.assert_frame_equal(g, golden)

        t, g = MI._deal_with_small_individuals(df, golden, 'keep')
        pd.testing.assert_frame_equal(t, df.sort_values(by='subject-PMID').reset_index(drop=True))
        pd.testing.assert_frame_equal(g, golden)

        golden1 = pd.DataFrame.from_dict(
             {'alpha': ['B, B, B, B, A'],
              'alpha1':['B, B, B, X, X, X, B, A'],
              'beta': ['B, B, B, A'],
              'beta1': ['A, B, B, X, X, X, B, A'],
              'subject-PMID':[2],},
              orient = 'index').T
        t, g = MI._deal_with_small_individuals(df, golden, 'use')
        pd.testing.assert_frame_equal(t, output1)
        pd.testing.assert_frame_equal(g, golden1)
    
    def test_entire_preprocessing(self):
        dict_out = MI.preprocessing(data = './tests/pairing_vdjdb/inputdf_test.csv', 
                                    epitope='ELAGIGILTV', weights=0.4, 
                                    method = None, prop_test='all', confidence = None, 
                                    correlation='no', translation='no', ind='no', 
                                    vgene='no', step=3, L=0.6, ones='keep', AAs=list('ARNDCEQGHILKMFPSTWYV-'))
        
        expected_df = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        assert dict_out['weights'] == 0.4
        assert dict_out['method'] == None
        assert dict_out['confidence'] == None
        assert dict_out['correlation'] == 'no'
        assert dict_out['corr'] == None
        assert dict_out['translation'] == 'no'
        assert dict_out['step'] == 3
        assert dict_out['L'] == 0.6
        assert dict_out['AAs'] == list('ARNDCEQGHILKMFPSTWYV-')
        assert dict_out['individuals'] == [['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138'], ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']]
        pd.testing.assert_frame_equal(dict_out['golden'], pd.DataFrame())
        pd.testing.assert_frame_equal(dict_out['test'], expected_df)

class TestRunning:
    def test_initial_cache(self):
        golden=pd.DataFrame()

        test = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        alphas, betas, startingA, startingB = MI._create_starting_cache_set(test, golden)

        alphas_exp = np.array([
            ['Y','T','V','L','G','-','-','N','E','K','T'],
            ['V','A','G','Y','G','G','S','Q','G','N','I'],
            ['V','S','F','G','-','-','-','N','E','K','T']
        ])

        betas_exp = np.array([
            ['S','S','F','T','P','-','Y','N','E','Q','F'],
            ['S','S','P','Q','G','L','G','T','E','A','F'],
            ['E','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        np.testing.assert_equal(alphas_exp, alphas)
        np.testing.assert_equal(betas_exp, betas)
        np.testing.assert_equal(alphas_exp, startingA) # alphas used as are
        assert startingB is not None
        with np.testing.assert_raises(AssertionError):
            np.testing.assert_equal(betas_exp, startingB) # I want to make sure betas are actually shuffled
    
    def test_initial_cache_with_golden(self):
        test = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        golden = pd.DataFrame(np.array(
            [['y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 's, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['v, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 's, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['v, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'e, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        alphas, betas, startingA, startingB = MI._create_starting_cache_set(test, golden)

        alphas_exp = np.array([
            ['Y','T','V','L','G','-','-','N','E','K','T'],
            ['V','A','G','Y','G','G','S','Q','G','N','I'],
            ['V','S','F','G','-','-','-','N','E','K','T']
        ])

        betas_exp = np.array([
            ['S','S','F','T','P','-','Y','N','E','Q','F'],
            ['S','S','P','Q','G','L','G','T','E','A','F'],
            ['E','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        startingA_exp = np.array([
            ['y','T','V','L','G','-','-','N','E','K','T'],
            ['v','A','G','Y','G','G','S','Q','G','N','I'],
            ['v','S','F','G','-','-','-','N','E','K','T']
        ])

        startingB_exp = np.array([
            ['s','S','F','T','P','-','Y','N','E','Q','F'],
            ['s','S','P','Q','G','L','G','T','E','A','F'],
            ['e','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        np.testing.assert_equal(alphas_exp, alphas)
        np.testing.assert_equal(betas_exp, betas)
        np.testing.assert_equal(startingA_exp, startingA)
        np.testing.assert_equal(startingB_exp, startingB)
    
    def test_translate(self):
        AAs=list('ARNDCEQGHILKMFPSTWYV-')
        alphas = np.array([
            ['Y','T','V','L','G','-','-','N','E','K','T'],
            ['V','A','G','Y','G','G','S','Q','G','N','I'],
            ['V','S','F','G','-','-','-','N','E','K','T']
        ])

        betas = np.array([
            ['S','S','F','T','P','-','Y','N','E','Q','F'],
            ['S','S','P','Q','G','L','G','T','E','A','F'],
            ['E','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        sA = np.array([
            ['D','T','V','L','G','-','-','N','E','K','T'],
            ['D','A','G','Y','G','G','S','Q','G','N','I'],
            ['D','S','F','G','-','-','-','N','E','K','T']
        ])

        sB = np.array([
            ['D','S','F','T','P','-','Y','N','E','Q','F'],
            ['D','S','P','Q','G','L','G','T','E','A','F'],
            ['D','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        a,b,s1,s2,aas = MI._translation(alphas, betas, sA, sB, AAs, 'no')

        np.testing.assert_equal(alphas, a)
        np.testing.assert_equal(betas, b)
        np.testing.assert_equal(sA, s1)
        np.testing.assert_equal(sB, s2)
        np.testing.assert_equal(aas, AAs)

        a,b,s1,s2,aas = MI._translation(alphas, betas, sA, sB, AAs, 'yes')
        
        ta = np.array([
            ['k','e','d','d','a','-','-','h','f','j','e'],
            ['d','d','a','k','a','a','e','h','a','h','d'],
            ['d','e','c','a','-','-','-','h','f','j','e']
        ])

        tb = np.array([
            ['e','e','c','e','b','-','k','h','f','h','c'],
            ['e','e','b','h','a','d','a','e','f','d','c'],
            ['f','a','h','a','c','d','a','h','b','h','i']
        ])

        tsa = np.array([
            ['f','e','d','d','a','-','-','h','f','j','e'],
            ['f','d','a','k','a','a','e','h','a','h','d'],
            ['f','e','c','a','-','-','-','h','f','j','e']
        ])

        tsb = np.array([
            ['f','e','c','e','b','-','k','h','f','h','c'],
            ['f','e','b','h','a','d','a','e','f','d','c'],
            ['f','a','h','a','c','d','a','h','b','h','i']
        ])

        np.testing.assert_equal(a, ta)
        np.testing.assert_equal(b, tb)
        np.testing.assert_equal(tsa, s1)
        np.testing.assert_equal(tsb, s2)
        np.testing.assert_equal(aas, list('abcdefghijkl-'))
    
    def test_apply_corr_and_confidence(self):
        matrix_0 = np.array([[1,2],
                             [2,1]])

        corr = np.array([[1,.5],
                         [.5,1]])
        
        exp = np.array([[1,1],
                        [1,1]])
        
        np.testing.assert_equal(MI._apply_corr_and_confidence(matrix_0, 'yes', corr, 'none',None),exp)
        np.testing.assert_equal(MI._apply_corr_and_confidence(matrix_0, 'no', corr, 'none',None),matrix_0)

    def test_dict_for_OnePairingRepeat(self):
        test = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        golden = pd.DataFrame(np.array(
            [['D, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'D, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['D, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'D, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['D, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'D, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        alphas_exp = np.array([
            ['Y','T','V','L','G','-','-','N','E','K','T'],
            ['V','A','G','Y','G','G','S','Q','G','N','I'],
            ['V','S','F','G','-','-','-','N','E','K','T']
        ])

        betas_exp = np.array([
            ['S','S','F','T','P','-','Y','N','E','Q','F'],
            ['S','S','P','Q','G','L','G','T','E','A','F'],
            ['E','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        startingA_exp = np.array([
            ['D','T','V','L','G','-','-','N','E','K','T'],
            ['D','A','G','Y','G','G','S','Q','G','N','I'],
            ['D','S','F','G','-','-','-','N','E','K','T']
        ])

        startingB_exp = np.array([
            ['D','S','F','T','P','-','Y','N','E','Q','F'],
            ['D','S','P','Q','G','L','G','T','E','A','F'],
            ['D','G','Q','G','F','V','G','Q','P','Q','H']
        ])

        correctpairs = ['YTVLGNEKT::SSFTPYNEQF::10xGenomics_M180',
                        'VAGYGGSQGNI::SSPQGLGTEAF::GitHubIssuesthg_M138',
                        'VSFGNEKT::EGQGFVGQPQH::PMID:12555663_M138']

        preprocess_output = {'weights' : 0.4,
                    'method': None,
                    'confidence': None,
                    'correlation' : 'no',
                    'corr' : None,
                    'translation' : 'no',
                    'step' : 3,
                    'L' : 0.6,
                    'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                    'individuals' : [['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138'], ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']],
                    'golden':golden,
                    'test':test}

        expected_output = {'weights' : 0.4,
                    'method': None,
                    'confidence': None,
                    'correlation' : 'no',
                    'corr' : None,
                    'step' : 3,
                    'L' : 0.6,
                    'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                    'individuals' : [['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138'], ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']],
                    'alphas' : alphas_exp,
                    'betas' : betas_exp,
                    'is_golden' : 'yes',
                    'startingA' : startingA_exp,
                    'startingB' : startingB_exp,
                    'correctPairs' : correctpairs}
        
        myoutput = MI.prepare_info_for_OnePairingRepeat(preprocess_output)

        assert sorted(list(myoutput.keys())) == sorted(list(expected_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == expected_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == expected_output[k]).all()
        

        # test with translation
        preprocess_output = {'weights' : 0.4,
                    'method': None,
                    'confidence': None,
                    'correlation' : 'no',
                    'corr' : None,
                    'translation' : 'yes',
                    'step' : 3,
                    'L' : 0.6,
                    'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                    'individuals' : [['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138'], ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']],
                    'golden':golden,
                    'test':test}
        
        correctpairs = ['keddahfje::eecebkhfhc::10xGenomics_M180',
                        'ddakaaehahd::eebhadaefdc::GitHubIssuesthg_M138',
                        'decahfje::fahacdahbhi::PMID:12555663_M138']
        
        ta = np.array([
            ['k','e','d','d','a','-','-','h','f','j','e'],
            ['d','d','a','k','a','a','e','h','a','h','d'],
            ['d','e','c','a','-','-','-','h','f','j','e']
        ])

        tb = np.array([
            ['e','e','c','e','b','-','k','h','f','h','c'],
            ['e','e','b','h','a','d','a','e','f','d','c'],
            ['f','a','h','a','c','d','a','h','b','h','i']
        ])

        tsa = np.array([
            ['f','e','d','d','a','-','-','h','f','j','e'],
            ['f','d','a','k','a','a','e','h','a','h','d'],
            ['f','e','c','a','-','-','-','h','f','j','e']
        ])

        tsb = np.array([
            ['f','e','c','e','b','-','k','h','f','h','c'],
            ['f','e','b','h','a','d','a','e','f','d','c'],
            ['f','a','h','a','c','d','a','h','b','h','i']
        ])
        
        expected_output = {'weights' : 0.4,
                    'method': None,
                    'confidence': None,
                    'correlation' : 'no',
                    'corr' : None,
                    'step' : 3,
                    'L' : 0.6,
                    'AAs' : list('abcdefghijkl-'),
                    'individuals' : [['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138'], ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']],
                    'alphas' : ta,
                    'betas' : tb,
                    'is_golden' : 'yes',
                    'startingA' : tsa,
                    'startingB' : tsb,
                    'correctPairs' : correctpairs}

        myoutput = MI.prepare_info_for_OnePairingRepeat(preprocess_output)

        assert sorted(list(myoutput.keys())) == sorted(list(expected_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == expected_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == expected_output[k]).all()
