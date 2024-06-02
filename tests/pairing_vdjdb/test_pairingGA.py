import PairingVDJdb_GA as GA
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
                'ind':'no', 
                'vgene':'no'},
            'for_GA':{
                'kNN':10,
                'n_repeats':100,
                'distance_type':'lev',
                'distanceDF':[None, None]},
            'folder':Path('data/output/pairing_GA/')
            }
        
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep']) == output

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['ind'] = 'yes'
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--small_ind', 'YES']) == output1

        output1 = copy.deepcopy(output)
        output1['for_preprocessing']['vgene'] = 'only'
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--Vgene', 'ONLY']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['kNN'] = 50
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--kNN', '50']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['n_repeats'] = 3
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--n_repeats', '3']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['distance_type'] = 'triplet'
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep', '--distance_type', 'TRIPLET']) == output1

        output1 = copy.deepcopy(output)
        output1['for_GA']['distance_type'] = 'tcrdist'
        output1['for_GA']['distanceDF'] = [Path('./tests/pairing_vdjdb/inputdf_test.csv').resolve(), 
                                           Path('./tests/pairing_vdjdb/inputdf_test.csv').resolve()]
        assert GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', 
            '--epitope', 'myep', '--distance_type', 'TCRDIST', 
            '--distance_dfa', './tests/pairing_vdjdb/inputdf_test.csv', '--distance_dfb', './tests/pairing_vdjdb/inputdf_test.csv']) == output1

        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/aaa.csv', '--epitope', 'myep'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','aaa'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','aaa', '--distance_type', 'sthg'])
        
        with pytest.warns():
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--distance_dfa','sthg'])

        with pytest.warns():
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--distance_dfb','sthg'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--distance_type','tcrdist'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--small_ind','STHG'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--Vgene','STHG'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--n_repeats','0'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--kNN','0'])
        
        with pytest.raises(AssertionError):
            GA.get_inputs(
            ['--input','data/vdj_cleaned_subset_for_MI.csv', '--epitope', 'myep',
            '--output','./crazypath'])

class TestPreprocessing:
    def test_preprocessing_function(self):
        dict_out = GA.preprocessing('./tests/pairing_vdjdb/inputdf_test.csv', 
                                    'ELAGIGILTV', 'no', 'no')

        expected_df = pd.DataFrame(np.array(
            [['CAYTVLGNEKLTF','C, A, Y, T, V, L, G, -, -, N, E, K, L, T, F', 'CASSFTPYNEQFF', 'C, A, S, S, F, T, P, -, -, Y, N, E, Q, F, F', '10xGenomics_M180'],
             ['CAVAGYGGSQGNLIF','C, A, V, A, G, Y, G, G, S, Q, G, N, L, I, F', 'CASSPQGLGTEAFF', 'C, A, S, S, P, Q, G, -, L, G, T, E, A, F, F', 'GitHubIssuesthg_M138'],
             ['CAVSFGNEKLTF','C, A, V, S, F, G, -, -, -, N, E, K, L, T, F', 'CAEGQGFVGQPQHF', 'C, A, E, G, Q, G, F, -, V, G, Q, P, Q, H, F', 'PMID:12555663_M138'],
            ]
        ), columns=['alpha', 'alpha1', 'beta', 'beta1', 'subject-PMID'])

        assert dict_out['individuals'] == ['10xGenomics_M180', 'GitHubIssuesthg_M138', 'PMID:12555663_M138']
        pd.testing.assert_frame_equal(dict_out['test'], expected_df)

    def test_reorder_distance_frame(self):

        seqs = ['a','b','b','c','a']
        inputdf = pd.DataFrame(
            data = [['aa', 'ac', 'aa', 'ab','ab'],
                    ['ca', 'cc', 'ca', 'cb','cb'],
                    ['aa', 'ac', 'aa', 'ab','ab'],
                    ['ba', 'bc', 'ba', 'bb','bb'],
                    ['ba', 'bc', 'ba', 'bb','bb']],
            index = ['a','c','a','b','b'],
            columns=['a','c','a','b','b']
        )
        outputdf = np.array([
                    ['aa', 'ab', 'ab', 'ac', 'aa'],
                    ['ba', 'bb', 'bb', 'bc', 'ba'],
                    ['ba', 'bb', 'bb', 'bc', 'ba'],
                    ['ca', 'cb', 'cb', 'cc', 'ca'],
                    ['aa', 'ab', 'ab', 'ac', 'aa']])
        
        np.testing.assert_array_equal(GA._reorder_distance_frame(inputdf, seqs), outputdf)

class TestOutputParsing:
    def test_parse_output(self):
        sample_input = np.array([
            [37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0],
            [37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0,37.0],
            [0.0,0.02702702702702703,0.0,0.02702702702702703,0.02702702702702703,0.05405405405405406,0.0,0.02702702702702703,0.02702702702702703,0.02702702702702703],
            [-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0],
            [18.0,5.0,15.0,22.0,15.0,30.0,27.0,26.0,5.0,6.0],
            [20.0,36.0,34.0,7.0,16.0,33.0,18.0,19.0,34.0,16.0],
            [11.0,35.0,31.0,32.0,22.0,4.0,28.0,7.0,14.0,24.0],
            [37.0,27.0,5.0,33.0,30.0,35.0,36.0,1.0,30.0,25.0],
            [22.0,10.0,9.0,18.0,10.0,29.0,1.0,2.0,8.0,36.0],
            [2.0,34.0,27.0,29.0,32.0,27.0,16.0,35.0,6.0,35.0]
        ])

        expected_output = np.array([
            [18.0,5.0,15.0,22.0,15.0,30.0,27.0,26.0,5.0,6.0],
            [20.0,36.0,34.0,7.0,16.0,33.0,18.0,19.0,34.0,16.0],
            [11.0,35.0,31.0,32.0,22.0,4.0,28.0,7.0,14.0,24.0],
            [37.0,27.0,5.0,33.0,30.0,35.0,36.0,1.0,30.0,25.0],
            [22.0,10.0,9.0,18.0,10.0,29.0,1.0,2.0,8.0,36.0],
            [2.0,34.0,27.0,29.0,32.0,27.0,16.0,35.0,6.0,35.0]
        ])

        outputs = GA.parse_output(sample_input)
        assert outputs[0] == 37
        outputs[1] == 37
        outputs[2] == np.array([0.0,0.02702702702702703,0.0,0.02702702702702703,0.02702702702702703,0.05405405405405406,0.0,0.02702702702702703,0.02702702702702703,0.02702702702702703])
        outputs[3] == np.array([-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0])
        np.testing.assert_array_equal(outputs[4], expected_output)
    
    def test_replace_with_seqs(self):
        sample_input = np.array([
            [1,2,1,1,1],
            [2,5,2,5,2],
            [4,3,3,3,4],
            [3,4,4,4,3],
            [5,1,5,2,5],
        ])

        alphas = ['a0', 'a1', 'a2', 'a3', 'a4']
        betas = ['b0', 'b1', 'b2', 'b3', 'b4']
        individuals = list([1,1,2,2,1])
        expected_output = pd.DataFrame(np.array([
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1],
            ['b1','b4','b1','b4','b1', 'a1', 'b1', 1],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2],
            ['b4','b0','b4','b1','b4', 'a4', 'b4', 1],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID'])
        expected_output['subject-PMID'] = expected_output['subject-PMID'].astype('int64') #hck to pass the test
        print(expected_output)
        pd.testing.assert_frame_equal(expected_output, 
                                      GA.replace_with_seq(sample_input, alphas, betas, individuals), 
                                      check_column_type=False, check_index_type=False)

    def test_get_mode(self):
        GA_output = pd.DataFrame(np.array([
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1],
            ['b1','b4','b1','b4','b1', 'a1', 'b1', 1],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2],
            ['b4','b0','b4','b1','b4', 'a4', 'b4', 1],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID'])

        expected_result = pd.DataFrame(np.array([
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1, 'b0', 4],
            ['b1','b4','b1','b4','b1', 'a1', 'b1', 1, 'b1', 3],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2, 'b2', 3],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2, 'b3', 3],
            ['b4','b0','b4','b1','b4', 'a4', 'b4', 1, 'b4', 3],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq'])
        expected_result['mode_freq'] = expected_result['mode_freq'].astype('int64') #hck to pass the test

        pd.testing.assert_frame_equal(GA.get_mode(GA_output), expected_result)

    def test_check_correct(self):
        myinput = pd.DataFrame(np.array([
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1, 'b0', 4],
            ['b4','b4','b1','b4','b4', 'a1', 'b1', 1, 'b4', 4],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2, 'b2', 3],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2, 'b3', 3],
            ['b1','b0','b4','b1','b1', 'a4', 'b4', 1, 'b1', 3],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq'])

        expected_output = pd.DataFrame(np.array([
            ['b0','b1','b0','b0','b0', 'a0', 'b0', 1, 'b0', 4, True],
            ['b4','b4','b1','b4','b4', 'a1', 'b1', 1, 'b4', 4, False],
            ['b3','b2','b2','b2','b3', 'a2', 'b2', 2, 'b2', 3, True],
            ['b2','b3','b3','b3','b2', 'a3', 'b3', 2, 'b3', 3, True],
            ['b1','b0','b4','b1','b1', 'a4', 'b4', 1, 'b1', 3, False],
        ]), columns=['repeat1', 'repeat2', 'repeat3', 'repeat4', 'repeat5', 'alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
        expected_output['correct'] = expected_output['correct'].map({'False':False, 'True':True})
        print(expected_output)

        correct_pairs = ['a0::b0::1', 'a1::b1::1', 'a2::b2::2', 'a3::b3::2', 'a4::b4::1']
        correct_pairs1 = ['a0::b0', 'a1::b1', 'a2::b2', 'a3::b3', 'a4::b4']

        pd.testing.assert_frame_equal(GA.check_correct_pairs_ID(myinput, correct_pairs1), expected_output)
        pd.testing.assert_frame_equal(GA.check_correct_pairs_ID(myinput, correct_pairs, myinput['subject-PMID'].tolist()), expected_output)

        # I also want a test where the same pair is correct in one ID but not in the other
        myinput = pd.DataFrame(np.array([
            ['aa', 'bb', 1, 'bb', 4],
            ['aaa', 'b1', 1, 'b1', 4],
            ['aa', 'b0', 2, 'bb', 3],
            ['aaa', 'bb', 2, 'b0', 3],
            ['aaaa', 'b4', 1, 'b4', 3],
        ]), columns=['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq'])

        correct_pairs = ['aa::bb::1', 'aaa::b1::1', 'aa::b0::2', 'aaa::bb::2', 'aaaa::b4::1']
        correct_pairs1 = ['aa::bb', 'aaa::b1', 'aa::b0', 'aaa::bb', 'aaaa::b4']

        output_wo_id = pd.DataFrame(np.array([
            ['aa', 'bb', 1, 'bb', 4, True],
            ['aaa', 'b1', 1, 'b1', 4, True],
            ['aa', 'b0', 2, 'bb', 3, True], # because not considering id assignment
            ['aaa', 'bb', 2, 'b0', 3, False],
            ['aaaa', 'b4', 1, 'b4', 3, True],
        ]), columns=['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
        output_wo_id['correct'] = output_wo_id['correct'].map({'False':False, 'True':True})

        output_w_id = pd.DataFrame(np.array([
            ['aa', 'bb', 1, 'bb', 4, True],
            ['aaa', 'b1', 1, 'b1', 4, True],
            ['aa', 'b0', 2, 'bb', 3, False], # because considering id assignment
            ['aaa', 'bb', 2, 'b0', 3, False],
            ['aaaa', 'b4', 1, 'b4', 3, True],
        ]), columns=['alpha', 'correct_beta', 'subject-PMID', 'mode', 'mode_freq', 'correct'])
        output_w_id['correct'] = output_w_id['correct'].map({'False':False, 'True':True})

        pd.testing.assert_frame_equal(GA.check_correct_pairs_ID(myinput, correct_pairs1, ids=None), output_wo_id)
        pd.testing.assert_frame_equal(GA.check_correct_pairs_ID(myinput, correct_pairs, ids=myinput['subject-PMID'].tolist()), output_w_id)

class TestGetGAInfo:
    def test_get_id_info(self):
        myinput = pd.DataFrame(np.array([
            ['a0', 'b0', 1],
            ['a1', 'b1', 1],
            ['a2', 'b2', 2],
            ['a3', 'b3', 2],
            ['a4', 'b4', 1],
        ]), columns=['alpha', 'beta', 'subject-PMID'])

        with pytest.raises(AssertionError):
             GA.get_id_info(myinput) # because not sorted

        M,N,FirstSeqSpec,LastSeqSpec,MSeqSpec,IndexSeqSpec = GA.get_id_info(myinput.sort_values(by='subject-PMID').reset_index(drop = True))

        assert M==5 # number of seqs
        assert N==2 # number of species
        assert FirstSeqSpec.tolist() == [1,4] # first seq of each species
        assert LastSeqSpec.tolist() == [3,5] # last seq of each species
        assert MSeqSpec.tolist() == [3,2] # number of seqs in each species
        assert IndexSeqSpec.tolist() == [1,1,1,2,2]