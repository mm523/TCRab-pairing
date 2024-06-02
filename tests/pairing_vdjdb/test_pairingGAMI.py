import PairingVDJdb_GAandMI as GAMI
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
                'weights':0.6, 
                'L':0.6,
                'method':None, 
                'prop_test':'all', 
                'confidence':'none', 
                'correlation':'no', 
                'translation':'no', 
                'ind':'no', 
                'vgene':'no', 
                'step':6, 
                'ones':'keep', 
                'AAs':list('ARNDCEQGHILKMFPSTWYV-')},
            'for_MI':{
                'n_repeats':1 # GA is normally deterministic, so default to 1
            },
            'for_GA':{
                'kNN':10,
                'n_repeats':100,
                'distance_type':'lev',
                'distanceDF':[None, None]},
            'folder':Path('data/output/pairing_GAMI/'),
            'GA_precomputed':None,
            'GA_thresh':0.95,
            'repair':True
            }
        
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep']) == output
        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['weights'] = 0.3
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_weights', '0.3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['method'] = 'RP'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_method', 'RP']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['prop_test'] = '0.75'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--prop_test', '0.75']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['confidence'] = 'hungarian'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_confidence', 'HUNGARIAN']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['correlation'] = 'no'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_correlation', 'NO']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['translation'] = 'yes'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_translation', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['ind'] = 'yes'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--small_ind', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['vgene'] = 'yes'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--Vgene', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['step'] = 6
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_step', '6']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['L'] = 0.5
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_L', '0.5']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['ones'] = 'throw'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--ones', 'THROW']) == output1

        output1 = copy.deepcopy(output)
        output1['for_MI']['n_repeats'] = 3
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_repeats', '3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['kNN'] = 3
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--GA_kNN', '3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['n_repeats'] = 3
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--GA_repeats', '3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['distance_type'] = 'weightedlev'
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--GA_distance_type', 'weightedLev']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['distance_type'] = 'tcrdist'
        output1['for_GA']['distanceDF'] = [Path('./tests/pairing_vdjdb/inputdf_test.csv').resolve(), 
                                           Path('./tests/pairing_vdjdb/inputdf_test.csv').resolve()]
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', 
            '--epitope', 'myep', '--GA_distance_type', 'TCRDIST', 
            '--GA_distance_dfa', './tests/pairing_vdjdb/inputdf_test.csv', '--GA_distance_dfb', './tests/pairing_vdjdb/inputdf_test.csv']) == output1

        output1 = copy.deepcopy(output)
        output1['GA_thresh'] = 0.1
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--GA_thresh', '0.1']) == output1

        output1 = copy.deepcopy(output)
        output1['repair'] = False
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--MI_repair', 'NO']) == output1


        output1 = copy.deepcopy(output)
        output1['GA_precomputed'] = Path('./tests/pairing_vdjdb/inputdf_test.csv').resolve()
        output1['for_GA']['kNN'] = ''
        output1['for_GA']['n_repeats'] = None
        output1['for_GA']['distance_type'] = ''
        output1['for_GA']['distanceDF'] = [None, None]
        assert GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--GA_precomputed', './tests/pairing_vdjdb/inputdf_test.csv']) == output1


        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/aaa.csv', '--epitope', 'myep'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','aaa'])
        
        with pytest.raises(ValueError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_weights','yes'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_weights','2'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_method','sthg'])
        
        with pytest.raises(ValueError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--prop_test','no'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--prop_test','3'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_confidence','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_correlation','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_translation','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--small_ind','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--Vgene','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_L','3'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--ones','STHG'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_repeats','0'])

        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','aaa', '--GA_distance_type', 'sthg'])
        
        with pytest.warns():
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_distance_dfa','sthg'])

        with pytest.warns():
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_distance_dfb','sthg'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_distance_type','tcrdist'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_repeats','0'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_thresh','3'])

        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--GA_precomputed','./crazypath'])
        
        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','./crazypath'])

        with pytest.raises(AssertionError):
            GAMI.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--MI_repair','Nah'])

class TestGA:
    def test_derive_GA_set_from_MI_preprocessing(self):
        test_MI = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        test_GA = pd.DataFrame(np.array(
            [['CAYTVLGNEKLTF','C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'CASSFTPYNEQFF', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['CAVAGYGGSQGNLIF','C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'CASSPQGLGTEAFF', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['CAVSFGNEKLTF','C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'CAEGQGFVGQPQHF', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        pd.testing.assert_frame_equal(GAMI._derive_GA_set_from_MI_preprocessing(test_MI), test_GA)

class TestGAtoMI:
    def test_get_stable_GA(self):
        GA_output = pd.DataFrame(np.array([
            # I have messed about with the mode numbers, so this is not realistic
            # should not matter because I am only interested in mode and mode_freq
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1, 'b0', 1, True],
            ['b4','b4','b1','b4','b4', 'a1', 'b1', 1, 'b4', 4, False],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2, 'b2', 3, True],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2, 'b3', 2, True],
            ['b1','b0','b4','b1','b1', 'a4', 'b4', 1, 'b1', 3, False],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
    
        GAstable = pd.DataFrame(np.array([
            # I have messed about with the mode numbers, so this is not realistic
            # should not matter because I am only interested in mode and mode_freq
            ['a1', 'b1', 1, 'b4', 4, False],
            ['a2', 'b2', 2, 'b2', 3, True],
            ['a4', 'b4', 1, 'b1', 3, False],
        ]), columns=['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
        GAstable['mode_freq'] = GAstable['mode_freq'].astype('int')

        pd.testing.assert_frame_equal(GAMI._get_stable_GA(GA_output, 0.5, 5), GAstable)

    def test_make_new_golden_set(self):
        GAstable = pd.DataFrame(np.array([
            # I have messed about with the mode numbers, so this is not realistic
            # should not matter because I am only interested in mode and mode_freq
            ['CAYTVLGNEKLTF', 'CASSFTPYNEQFF', 1, 'CAEGQGFVGQPQHF', 4, False],
            ['CAVAGYGGSQGNLIF', 'CASSPQGLGTEAFF', 2, 'CASSPQGLGTEAFF', 3, True],
            ['CAVSFGNEKLTF', 'CAEGQGFVGQPQHF', 1, 'CASSFTPYNEQFF', 3, False],
        ]), columns=['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
        GAstable['mode_freq'] = GAstable['mode_freq'].astype('int')

        MItest = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', 1],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 2],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 1],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 1], # these two lines repeat but w different IDs
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 2], # these two lines repeat but w different IDs

            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])
    
        expected_outcome = pd.DataFrame(np.array(
            [['CAYTVLGNEKLTF', 'Y, T, V, L, G, -, -, N, E, K, T', 'CAEGQGFVGQPQHF', 'E, G, Q, G, F, V, G, Q, P, Q, H',  1, 4, False, 'CASSFTPYNEQFF'],
             ['CAVAGYGGSQGNLIF', 'V, A, G, Y, G, G, S, Q, G, N, I','CASSPQGLGTEAFF', 'S, S, P, Q, G, L, G, T, E, A, F', 2, 3, True, 'CASSPQGLGTEAFF'],
             ['CAVSFGNEKLTF', 'V, S, F, G, -, -, -, N, E, K, T', 'CASSFTPYNEQFF', 'S, S, F, T, P, -, Y, N, E, Q, F', 1, 3, False, 'CAEGQGFVGQPQHF'],
            ]
        ), columns=['alpha1', 'alpha', 'beta1', 'beta', 'subject-PMID', 'mode_freq', 'correct', 'correct_beta'])
        expected_outcome['mode_freq'] = expected_outcome['mode_freq'].astype('int')

        golden2 = GAMI._make_new_golden_set(GAstable, MItest)
        
        pd.testing.assert_frame_equal(golden2, expected_outcome)
 
    def test_prepare_new_test_set(self):
        MItest = pd.DataFrame(np.array(
            [['Y, T, V, L, G, -, -, N, E, K, T', 'C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'S, S, F, T, P, -, Y, N, E, Q, F', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', 1],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 2],
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 1],
             ['V, A, G, Y, G, G, S, Q, G, N, I', 'C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'S, S, P, Q, G, L, G, T, E, A, F', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 1], # these two lines repeat but w different IDs
             ['V, S, F, G, -, -, -, N, E, K, T', 'C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'E, G, Q, G, F, V, G, Q, P, Q, H', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 2], # these two lines repeat but w different IDs

            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        correctpairs = ['YTVLGNEKT::SSFTPYNEQF::1',
                        'VAGYGGSQGNI::SSPQGLGTEAF::2',
                        'VSFGNEKT::EGQGFVGQPQH::1',
                        'VAGYGGSQGNI::SSPQGLGTEAF::1',
                        'VSFGNEKT::EGQGFVGQPQH::2']

        input_alphas = np.array(
            [['Y', 'T', 'V', 'L', 'G', '-', '-', 'N', 'E', 'K', 'T'],
             ['V', 'A', 'G', 'Y', 'G', 'G', 'S', 'Q', 'G', 'N', 'I'],
             ['V', 'S', 'F', 'G', '-', '-', '-', 'N', 'E', 'K', 'T'],
             ['V', 'A', 'G', 'Y', 'G', 'G', 'S', 'Q', 'G', 'N', 'I'], # these two lines repeat but w different IDs
             ['V', 'S', 'F', 'G', '-', '-', '-', 'N', 'E', 'K', 'T'], # these two lines repeat but w different IDs
            ])

        input_betas = np.array(
            [['S', 'S', 'F', 'T', 'P', '-', 'Y', 'N', 'E', 'Q', 'F'],
             ['S', 'S', 'P', 'Q', 'G', 'L', 'G', 'T', 'E', 'A', 'F'],
             ['E', 'G', 'Q', 'G', 'F', 'V', 'G', 'Q', 'P', 'Q', 'H'],
             ['S', 'S', 'P', 'Q', 'G', 'L', 'G', 'T', 'E', 'A', 'F'], # these two lines repeat but w different IDs
             ['E', 'G', 'Q', 'G', 'F', 'V', 'G', 'Q', 'P', 'Q', 'H'], # these two lines repeat but w different IDs
            ])

        IDs = [['1','2','1','1','2'],['1','2','1','1','2']]

        golden_new = pd.DataFrame(np.array(
            [['CAYTVLGNEKLTF', 'Y, T, V, L, G, -, -, N, E, K, T', 'CAEGQGFVGQPQHF', 'E, G, Q, G, F, V, G, Q, P, Q, H',  '1', 4, False, 'CASSFTPYNEQFF'],
             ['CAVAGYGGSQGNLIF', 'V, A, G, Y, G, G, S, Q, G, N, I','CASSPQGLGTEAFF', 'S, S, P, Q, G, L, G, T, E, A, F', '2', 3, True, 'CASSPQGLGTEAFF'],
             ['CAVSFGNEKLTF', 'V, S, F, G, -, -, -, N, E, K, T', 'CASSFTPYNEQFF', 'S, S, F, T, P, -, Y, N, E, Q, F', '1', 3, False, 'CAEGQGFVGQPQHF'],
            ]
        ), columns=['alpha1', 'alpha', 'beta1', 'beta', 'subject-PMID', 'mode_freq', 'is_correct', 'correct_beta'])
        
        dict_input = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'no',
                # starting here would be alphas + random shuffle
                'startingA' : input_alphas, 
                'startingB' : np.random.shuffle(input_betas.copy()),
                'correctPairs' : correctpairs}

        expected_alpha_output = np.array([
            ['V', 'A', 'G', 'Y', 'G', 'G', 'S', 'Q', 'G', 'N', 'I'],
            ['V', 'S', 'F', 'G', '-', '-', '-', 'N', 'E', 'K', 'T'],
        ])

        expected_startingA_output = np.array([
            ['Y', 'T', 'V', 'L', 'G', '-', '-', 'N', 'E', 'K', 'T'],
            ['V', 'A', 'G', 'Y', 'G', 'G', 'S', 'Q', 'G', 'N', 'I'],
            ['V', 'S', 'F', 'G', '-', '-', '-', 'N', 'E', 'K', 'T']
        ])

        expected_beta_output = np.array([
            ['S', 'S', 'P', 'Q', 'G', 'L', 'G', 'T', 'E', 'A', 'F'],
            ['E', 'G', 'Q', 'G', 'F', 'V', 'G', 'Q', 'P', 'Q', 'H']
        ])

        expected_startingB_output = np.array([
            ['E', 'G', 'Q', 'G', 'F', 'V', 'G', 'Q', 'P', 'Q', 'H'],
            ['S', 'S', 'P', 'Q', 'G', 'L', 'G', 'T', 'E', 'A', 'F'],
            ['S', 'S', 'F', 'T', 'P', '-', 'Y', 'N', 'E', 'Q', 'F'],
        ])

        expected_IDs_output = [['1','2'],['1','2']]

        dict_output = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                'individuals' : expected_IDs_output,
                'alphas' : expected_alpha_output,
                'betas' : expected_beta_output,
                'is_golden' : 'yes',
                'startingA' : expected_startingA_output,
                'startingB' : expected_startingB_output,
                'correctPairs' : correctpairs}
        
        myoutput = GAMI._prepare_new_test_set(dict_input, golden_new, 'no',False)

        assert sorted(list(myoutput.keys())) == sorted(list(dict_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == dict_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == dict_output[k]).all()
        

        # now test correct behaviour when previous golden set provided
        previous_A_training = np.array([
            ['s', 'o', 'm', 'e', 's', 'e', 'q', 'u', 'e', 'n', 'c'],
        ])

        previous_B_training = np.array([
            ['a', 'n', 'o', 't', 'h', 'e', 'r', '-', 's', 'e', 'q'],
        ])

        expected_startingA_output = np.array([
            ['s', 'o', 'm', 'e', 's', 'e', 'q', 'u', 'e', 'n', 'c'],
            ['Y', 'T', 'V', 'L', 'G', '-', '-', 'N', 'E', 'K', 'T'],
            ['V', 'A', 'G', 'Y', 'G', 'G', 'S', 'Q', 'G', 'N', 'I'],
            ['V', 'S', 'F', 'G', '-', '-', '-', 'N', 'E', 'K', 'T']
        ])

        expected_startingB_output = np.array([
            ['a', 'n', 'o', 't', 'h', 'e', 'r', '-', 's', 'e', 'q'],
            ['E', 'G', 'Q', 'G', 'F', 'V', 'G', 'Q', 'P', 'Q', 'H'],
            ['S', 'S', 'P', 'Q', 'G', 'L', 'G', 'T', 'E', 'A', 'F'],
            ['S', 'S', 'F', 'T', 'P', '-', 'Y', 'N', 'E', 'Q', 'F'],
        ])

        dict_input = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'yes',
                # starting here would be alphas + random shuffle
                'startingA' : previous_A_training, 
                'startingB' : previous_B_training,
                'correctPairs' : correctpairs}
        
        dict_output = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('ARNDCEQGHILKMFPSTWYV-'),
                'individuals' : expected_IDs_output,
                'alphas' : expected_alpha_output,
                'betas' : expected_beta_output,
                'is_golden' : 'yes',
                'startingA' : expected_startingA_output,
                'startingB' : expected_startingB_output,
                'correctPairs' : correctpairs}
        
        myoutput = GAMI._prepare_new_test_set(dict_input, golden_new, 'no', False)

        assert sorted(list(myoutput.keys())) == sorted(list(dict_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == dict_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == dict_output[k]).all()

        # now test correct behaviour when translating

        input_alphas = np.array(
            [['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'], # these two lines repeat but w different IDs
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'], # these two lines repeat but w different IDs
            ])

        input_betas = np.array(
            [['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'], # these two lines repeat but w different IDs
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'], # these two lines repeat but w different IDs
            ])
            
        expected_alpha_output = np.array([
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e']
        ])

        expected_beta_output = np.array([
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'] 
        ])

        expected_startingA_output = np.array([
            ['s', 'o', 'm', 'e', 's', 'e', 'q', 'u', 'e', 'n', 'c'],
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
        ])

        expected_startingB_output = np.array([
            ['a', 'n', 'o', 't', 'h', 'e', 'r', '-', 's', 'e', 'q'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],            
        ])

        
        dict_input = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'yes',
                # starting here would be alphas + random shuffle
                'startingA' : previous_A_training, 
                'startingB' : previous_B_training,
                'correctPairs' : correctpairs}
        
        dict_output = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                'individuals' : expected_IDs_output,
                'alphas' : expected_alpha_output,
                'betas' : expected_beta_output,
                'is_golden' : 'yes',
                'startingA' : expected_startingA_output,
                'startingB' : expected_startingB_output,
                'correctPairs' : correctpairs}
        
        myoutput = GAMI._prepare_new_test_set(dict_input, golden_new, 'yes', False)

        assert sorted(list(myoutput.keys())) == sorted(list(dict_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == dict_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == dict_output[k]).all()
        
        # check that the result is still correct when you have duplicates
        # but want to put only one of them in golden set and keep the other in test

        input_alphas = np.array(
            [['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'], # these two lines repeat but w different IDs
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'], # these two lines repeat but w different IDs
            ])

        input_betas = np.array(
            [['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'], # these two lines repeat but w different IDs
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'], # these two lines repeat but w different IDs
            ])
            
        expected_alpha_output = np.array([
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e']
        ])

        expected_beta_output = np.array([
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'] 
        ])

        expected_startingA_output = np.array([
            ['s', 'o', 'm', 'e', 's', 'e', 'q', 'u', 'e', 'n', 'c'],
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
        ])

        expected_startingB_output = np.array([
            ['a', 'n', 'o', 't', 'h', 'e', 'r', '-', 's', 'e', 'q'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],            
        ])

        IDs = [['1','2', '1','2', '1','2', '1','1','2'],
               ['1','2', '1','2', '1','2', '1','1','2']]
        
        expected_IDs_output = [['1','2', '1','2', '1','2'],
                               ['1','2', '1','2', '1','2']]

        
        dict_input = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'yes',
                # starting here would be alphas + random shuffle
                'startingA' : previous_A_training, 
                'startingB' : previous_B_training,
                'correctPairs' : correctpairs}
        
        dict_output = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                'individuals' : expected_IDs_output,
                'alphas' : expected_alpha_output,
                'betas' : expected_beta_output,
                'is_golden' : 'yes',
                'startingA' : expected_startingA_output,
                'startingB' : expected_startingB_output,
                'correctPairs' : correctpairs}
        
        myoutput = GAMI._prepare_new_test_set(dict_input, golden_new, 'yes',False)

        assert sorted(list(myoutput.keys())) == sorted(list(dict_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == dict_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == dict_output[k]).all()
        
        #### If I DO want to pair again after GA
        # check that the result is still correct when you have duplicates
        # but want to put only one of them in golden set and keep the other in test

        input_alphas = np.array(
            [['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
             ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'], # these two lines repeat but w different IDs
             ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'], # these two lines repeat but w different IDs
            ])

        input_betas = np.array(
            [['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
             ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'], # these two lines repeat but w different IDs
             ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'], # these two lines repeat but w different IDs
            ])
            
        expected_alpha_output = np.array([
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e']
        ])

        expected_beta_output = np.array([
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'] 
        ])

        expected_startingA_output = np.array([
            ['s', 'o', 'm', 'e', 's', 'e', 'q', 'u', 'e', 'n', 'c'],
            ['k', 'e', 'd', 'd', 'a', '-', '-', 'h', 'f', 'j', 'e'],
            ['d', 'd', 'a', 'k', 'a', 'a', 'e', 'h', 'a', 'h', 'd'],
            ['d', 'e', 'c', 'a', '-', '-', '-', 'h', 'f', 'j', 'e'],
        ])

        expected_startingB_output = np.array([
            ['a', 'n', 'o', 't', 'h', 'e', 'r', '-', 's', 'e', 'q'],
            ['f', 'a', 'h', 'a', 'c', 'd', 'a', 'h', 'b', 'h', 'i'],
            ['e', 'e', 'b', 'h', 'a', 'd', 'a', 'e', 'f', 'd', 'c'],
            ['e', 'e', 'c', 'e', 'b', '-', 'k', 'h', 'f', 'h', 'c'],            
        ])

        IDs = [['1','2', '1','2', '1','2', '1','1','2'],
               ['1','2', '1','2', '1','2', '1','1','2']]
        
        dict_input = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'yes',
                # starting here would be alphas + random shuffle
                'startingA' : previous_A_training, 
                'startingB' : previous_B_training,
                'correctPairs' : correctpairs}
        
        dict_output = {'weights' : 0.4,
                'method': None,
                'confidence': None,
                'correlation' : 'no',
                'corr' : None,
                'step' : 3,
                'L' : 0.6,
                'AAs' : list('abcdefghijkl-'),
                # because I am allowing re-pairing, inputs to MI stay the same
                # but the starting set is changed as when MI_repair == False
                'individuals' : IDs,
                'alphas' : input_alphas,
                'betas' : input_betas,
                'is_golden' : 'yes',
                'startingA' : expected_startingA_output,
                'startingB' : expected_startingB_output,
                'correctPairs' : correctpairs}
        
        myoutput = GAMI._prepare_new_test_set(dict_input, golden_new, 'yes',True)

        assert sorted(list(myoutput.keys())) == sorted(list(dict_output.keys()))
        for k in myoutput.keys():
            print(k)
            try:
                assert myoutput[k] == dict_output[k]
            except Exception as e:
                print(e)
                assert (myoutput[k] == dict_output[k]).all()