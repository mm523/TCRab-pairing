import functions.IPAfunctions as ipa
import functions.myfunctions as mf
import numpy as np
import pandas as pd
from itertools import product
import time
class TestCoreFunctions:
    def test_calculate_weights(self):
        A = np.array([
            ['m','a','r','t','i'],
            ['m','a','r','t','a'],
            ['m','a','m','m','a']
        ])
        B = np.array([
            ['t','o','m'],
            ['t','o','m'],
            ['m','a','x']
        ])

        assert ipa._calculate_weights(A,B,'no') == None
        np.testing.assert_equal(ipa._calculate_weights(A,B,0.15), np.array([.5, .5, 1.]))
        np.testing.assert_equal(ipa._calculate_weights(A,B,0.1), np.array([1., 1., 1.]))
        np.testing.assert_equal(ipa._calculate_weights(A,B,0.9), np.array([1/3, 1/3, 1/3]))


    def test_pseudocount_single(self):
        X = np.array(list("MARTINA"))
        pc = ipa.pseudocounter_single(X, L = 0.15, q = 21)
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
        AAdict = dict(zip(AAs, range(len(AAs))))
        weights = np.array([0.1, 0.2, 0.1, 0.3, 0.05, 0.1, 0.15])

        assert round(pc[AAdict["A"]], 3) == 0.25
        assert pc[AAdict["M"]] == pc[AAdict["R"]] == pc[AAdict["T"]] == pc[AAdict["I"]] == pc[AAdict["N"]]
        assert round(pc[AAdict["M"]], 3) == 0.129
        assert pc[AAdict["F"]] == 0.15/21 # one that is not found

        pc1 = ipa.pseudocounter_single(X, L = 0.15, q = 21, weights=weights)
        assert pc1[AAdict["F"]] == 0.15/21 # one that is not found
        pc1 = ipa.pseudocounter_single(X, L = 0, q = 21, weights=weights)
        assert pc1[AAdict["F"]] == 0. # one that is not found
        assert pc1[AAdict["M"]] == 0.1/1
        assert pc1[AAdict["A"]] == 0.35

        pc1 = ipa.pseudocounter_single(X, L = 1, q = 21, weights=weights)
        assert pc1[AAdict["F"]] == 1/21 # one that is not found
        assert pc1[AAdict["M"]] == 1/21
        assert pc1[AAdict["A"]] == 1/21

    def test_calculate_coincidence(self):
        X = np.array([list("MARTINAA"), list("ARIANNAA")]).T 
        mytuple = ['A', 'A']
        weights = np.array([0.1, 0.2, 0.1, 0.3, 0.05, 0.05, 0.15, 0.05])

        assert ipa._calculate_coincidence(X, mytuple) == 2/8
        assert ipa._calculate_coincidence(X, mytuple, weights) == 0.2/1


    def test_pseudocount_double(self):

        s = time.time()
        X = np.array([list("MARTINAA"), list("ARIANNAA")]).T 
        # here MARTINA acts as a column down alpha, 
        # ARIANNA is a column down beta
        pc = ipa.pseudocounter_double(X, L = 0.15, q = 21)
        print('pseudocounter_double elapsed: ', time.time() - s)
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
        AAdict = dict(zip(AAs, range(len(AAs))))

        assert round(pc[AAdict["A"], AAdict["A"]], 3) == 0.213
        assert pc[AAdict["M"], AAdict["A"]] == pc[AAdict["A"], AAdict["R"]] == pc[AAdict["R"], AAdict["I"]] == pc[AAdict["T"], AAdict["A"]] == pc[AAdict["I"], AAdict["N"]] == pc[AAdict["N"], AAdict["N"]]
        assert round(pc[AAdict["M"], AAdict["A"]], 3) == 0.107
        assert round(pc[AAdict["M"], AAdict["C"]], 5) == 0.00034 # one that is not found

        pc = ipa.pseudocounter_double(X, L = 1, q = 21)
        assert pc[AAdict["A"], AAdict["A"]] == 1/21**2
        assert pc[AAdict["M"], AAdict["A"]] == pc[AAdict["A"], AAdict["R"]] == pc[AAdict["R"], AAdict["I"]] == pc[AAdict["T"], AAdict["A"]] == pc[AAdict["I"], AAdict["N"]] == pc[AAdict["N"], AAdict["N"]] == 1/21**2
        assert pc[AAdict["M"], AAdict["A"]] == 1/21**2
        assert pc[AAdict["M"], AAdict["C"]] == 1/21**2 # one that is not found

    def test_PMIij(self):
        A = np.array(list("MARTINAA"))
        B = np.array(list("ARIANNAA"))
        AB = np.array([list("MARTINAA"), list("ARIANNAA")]).T
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
        AAdict = dict(zip(AAs, range(len(AAs)))) 

        pcA = ipa.pseudocounter_single(A, L = 0.15, q = 21)
        pcB = ipa.pseudocounter_single(B, L = 0.15, q = 21)
        pcAB = ipa.pseudocounter_double(AB)

        PMIij = ipa.calculate_PMIij(A, B)

        np.testing.assert_almost_equal(PMIij[AAdict["M"], AAdict["A"]], np.log(pcAB[AAdict["M"], AAdict["A"]]/(pcA[AAdict["M"]]*pcB[AAdict["A"]])))
        np.testing.assert_almost_equal(PMIij[AAdict["A"], AAdict["A"]], np.log(pcAB[AAdict["A"], AAdict["A"]]/(pcA[AAdict["A"]]*pcB[AAdict["A"]])))
        np.testing.assert_almost_equal(PMIij[AAdict["F"], AAdict["A"]], np.log(pcAB[AAdict["F"], AAdict["A"]]/(pcA[AAdict["F"]]*pcB[AAdict["A"]])))
        np.testing.assert_almost_equal(PMIij[AAdict["F"], AAdict["F"]], np.log(pcAB[AAdict["F"], AAdict["F"]]/(pcA[AAdict["F"]]*pcB[AAdict["F"]])))

    def test_PMIij_RP(self):
        A = np.array(list("MARTINAA"))
        B = np.array(list("ARIANNAA"))
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
        AAdict = dict(zip(AAs, range(len(AAs)))) 

        PMIij = ipa.calculate_PMIij(A, B, method='RP', L = 0)
        # print(PMIij)

        # please note this behaviour: all goes to 0 because happens as frequently as expected
        assert PMIij[AAdict["M"], AAdict["A"]].round(14) == 0.
        assert PMIij[AAdict["A"], AAdict["A"]].round(14) == 0.
        assert PMIij[AAdict["F"], AAdict["A"]].round(14) == 0.
        assert PMIij[AAdict["F"], AAdict["F"]].round(14) == 0.

        A = np.array(list("MARTINAMAAAAAAM"))
        B = np.array(list("ARIANNAAAAAAAAA"))
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
        AAdict = dict(zip(AAs, range(len(AAs)))) 

        PMIij = ipa.calculate_PMIij(A, B, method='RP', L = 0)
        # print(PMIij)

        assert PMIij[AAdict["M"], AAdict["A"]].round(14) == round(0.38659662380378124,14)
        assert PMIij[AAdict["A"], AAdict["A"]].round(14) == round(0.25306523117925883,14)
        assert PMIij[AAdict["F"], AAdict["A"]].round(14) == round(-2.6579258139196407,14)
        assert PMIij[AAdict["F"], AAdict["F"]].round(14) == round(1.0886619578149421,14)


    def test_PMI(self):
        s = time.time()
        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).head()
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps") # this is the way I prepare the data in the other script, to keep it consistent

        alphas = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))

        totPMI = ipa.calculate_PMI_AllPositions(alphas, betas)
        print('calculate_PMI_AllPositions elapsed: ', time.time() - s)

        assert (totPMI[(0,0)] == ipa.calculate_PMIij(alphas[:,0], betas[:,0])).all()
        assert (totPMI[(5,10)] == ipa.calculate_PMIij(alphas[:,5], betas[:,10])).all()
        assert (totPMI[(12,10)] == ipa.calculate_PMIij(alphas[:,12], betas[:,10])).all()

    def test_MI(self):
        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0)#.head()
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps") # this is the way I prepare the data in the other script, to keep it consistent

        alphas = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))

        s = time.time()
        MI = ipa.calculate_MI_Pairwise(alphas, betas)
        print('calculate_MI_Pairwise elapsed: ', time.time() - s)

        Ai = alphas[:,1] 
        Bj = betas[:, 3]
        X = np.stack([np.array(Ai), np.array(Bj)]).T

        fis = ipa.pseudocounter_single(Ai, L = 0.15, q = 21)
        fjs = ipa.pseudocounter_single(Bj, L = 0.15, q = 21)
        fijs = ipa.pseudocounter_double(X, L = 0.15, q = 21)

        MIab = (fijs*np.log(fijs/(fis[:,None]*fjs[None,:]))).sum()

        assert MI[1,3].round(15) == round(MIab,15)

    def test_Sab(self):
        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).head()
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 18, max_len2 = 18)

        test = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).iloc[0:10]
        test = mf.prepare_data(test, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 18, max_len2 = 18)

        alphas = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))

        alphas_t = np.array(test.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas_t = np.array(test.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))

        totPMI = ipa.calculate_PMI_AllPositions(alphas, betas)

        allpmis = {}

        for x in [1,2,5,8,9]:
            a = alphas_t[x, :]
            b = betas_t[x, :]

            allpmis[x] = []

            AAs = list("ARNDCEQGHILKMFPSTWYV-")
            AAdict = dict(zip(AAs, range(len(AAs)))) 

            for i in range(len(a)):
                for j in range(len(b)):
                    allpmis[x].append(totPMI[i, j][AAdict[a[i]], AAdict[b[j]]])

        s = time.time()
        assert ipa.calculate_Sab(alphas_t[1,:], betas_t[1,:], totPMI) == sum(allpmis[1])
        print('calculate_Sab elapsed: ', time.time() - s)
        assert ipa.calculate_Sab(alphas_t[2,:], betas_t[2,:], totPMI) == sum(allpmis[2])
        assert ipa.calculate_Sab(alphas_t[5,:], betas_t[5,:], totPMI) == sum(allpmis[5])
        assert ipa.calculate_Sab(alphas_t[8,:], betas_t[8,:], totPMI) == sum(allpmis[8])
        assert ipa.calculate_Sab(alphas_t[9,:], betas_t[9,:], totPMI) == sum(allpmis[9])

    def test_training_set(self):
        scoreMatrix = np.array([
            [0, 5, 3, 8],
            [7, 4, 9, 5],
            [1, 2, 7, 4],
            [1, 0, 5, 7]
        ])

        training_size = 2

        a, b = ipa.get_new_training_set(scoreMatrix, training_size)

        assert a == [1, 0]
        assert b == [2, 3]
class TestOneTrainingLoop:
    def test_create_score_matrix(self):
        alphas = np.array([
            ['Y','T'],
            ['V','A']
        ])

        betas = np.array([
            ['S','S'],
            ['T','S']
        ])
        AAs = ['Y', 'V', 'S', 'T', 'A']

        pmi = {(0,0):np.array([ #Y  V  S  T  A
                                [0, 0, 1, 0, 0], #Y
                                [0, 0, 0, 1, 0], #V
                                [0, 0, 0, 0, 0], #S
                                [0, 0, 0, 0, 0], #T
                                [0, 0, 0, 0, 0]  #A
                                ]), 
                (0,1):np.array([ #Y  V  S  T  A
                                [0, 0, 1, 0, 0], #Y
                                [0, 0, 1, 0, 0], #V
                                [0, 0, 0, 0, 0], #S
                                [0, 0, 0, 0, 0], #T
                                [0, 0, 0, 0, 0]  #A
                                ]), 
                (1,0):np.array([ #Y  V  S  T  A
                                [0, 0, 0, 0, 0], #Y
                                [0, 0, 0, 0, 0], #V
                                [0, 0, 0, 0, 0], #S
                                [0, 0, 1, 0, 0], #T
                                [0, 0, 0, 1, 0]  #A
                                ]), 
                (1,1):np.array([ #Y  V  S  T  A
                                [0, 0, 0, 0, 2], #Y
                                [0, 0, 0, 0, 2], #V
                                [0, 0, 0, 0, 0], #S
                                [0, 0, 1, 0, 0], #T
                                [0, 0, 1, 0, 0]  #A
                                ])}

        expected_output = np.array([
                                    [4,2],
                                    [2,4]
                                    ])
        
        X = ipa.create_score_matrix(alphas, betas, IDs=None, pmi = pmi, AAs = AAs)
        np.testing.assert_equal(X, expected_output)


    def test_training_loop(self):
        scoreMatrix = np.array([
            [0, 5, 3, 8],
            [7, 4, 9, 5],
            [1, 2, 7, 4],
            [1, 0, 5, 7]
        ])

        training_size = 2

        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).head(4)
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps") # this is the way I prepare the data in the other script, to keep it consistent

        alphas = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))

        expected_matrix = np.array([
            [-373.58, -438.08, -677.84, 106.61],
            [-373.58, -366.77, 106.61, -677.84],
            [-269.67, -304.96, -538.28, -146.06],
            [-347.69, -348.02, -17.27, -566.38]
        ])

        s = time.time()
        m1 = ipa.do_a_training_loop(training_size, scoreMatrix, alphas, betas)
        print('do_a_training_loop elapsed: ', time.time() - s)

        np.testing.assert_array_equal(m1.round(2), expected_matrix)

    def test_training_loop_with_background(self):
        scoreMatrix = np.array([
            [0, 5, 3, 8],
            [7, 4, 9, 5],
            [1, 2, 7, 4],
            [1, 0, 5, 7]
        ])

        training_size = 2

        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).head(4)
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps") # this is the way I prepare the data in the other script, to keep it consistent

        alphas = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))

        expected_matrix = np.array([
            [-373.58, -438.08, -677.84, 106.61],
            [-373.58, -366.77, 106.61, -677.84],
            [-269.67, -304.96, -538.28, -146.06],
            [-347.69, -348.02, -17.27, -566.38]
        ])

        # setting background to be exactly what the chosen training is so that my final results don't change
        backgroundA = alphas[[0,1]]
        backgroundB = betas[[3,2]]

        s = time.time()
        m1 = ipa.do_a_training_loop(training_size, scoreMatrix, alphas, betas, background_alpha=backgroundA, background_beta=backgroundB)
        print('do_a_training_loop elapsed: ', time.time() - s)

        np.testing.assert_array_equal(m1.round(2), expected_matrix)


    def test_training_loop_ID(self):
        scoreMatrix = np.array([
            [0, np.nan, 3, 8, np.nan, np.nan, np.nan, np.nan],
            [7, np.nan, 9, 5, np.nan, np.nan, np.nan, np.nan],
            [np.nan, 2, np.nan, np.nan, 0,0,0,0],
            [1, np.nan, 5, 7, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, 0, np.nan, 0, 5, 3, 8],
            [np.nan, np.nan, 0, np.nan, 7, 4, 9, 5],
            [np.nan, np.nan, 0, np.nan, 1, 2, 7, 4],
            [np.nan, np.nan, 0, np.nan, 1, 0, 5, 7]
        ])

        ids = [1,2,1,1,2,2,2,2]
        IDs = [ids, ids]

        training_size = 4 # making this 4 so I take top 2 from ID1 and top 2 from ID2 and it works out same numbers as test before

        vdj = pd.read_csv("data/vdj_cleaned_subset_for_MI.csv", index_col=0).head(4)
        vdj = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps") # this is the way I prepare the data in the other script, to keep it consistent

        a1 = np.array(vdj.cdr3a_IMGTgaps_padded.str.split(', ', expand=True))
        alphas = np.concatenate((a1, a1), axis = 0)
        # print(alphas)
        b1 = np.array(vdj.cdr3b_IMGTgaps_padded.str.split(', ', expand=True))
        betas = np.concatenate((b1, b1), axis = 0)

        expected_matrix = np.array([
            [-373.58, np.nan, -677.84, 106.61, np.nan, np.nan, np.nan, np.nan],
            [np.nan, -366.77, np.nan, np.nan, -373.58, -366.77, 106.61, -677.84],
            [-269.67, np.nan, -538.28, -146.06, np.nan, np.nan, np.nan, np.nan],
            [-347.69, np.nan, -17.27, -566.38, np.nan, np.nan, np.nan, np.nan],
            [np.nan, -438.08, np.nan, np.nan, -373.58, -438.08, -677.84, 106.61],
            [np.nan, -366.77, np.nan, np.nan, -373.58, -366.77, 106.61, -677.84],
            [np.nan, -304.96, np.nan, np.nan, -269.67, -304.96, -538.28, -146.06],
            [np.nan, -348.02, np.nan, np.nan, -347.69, -348.02, -17.27, -566.38]
        ])

        s = time.time()
        m1 = ipa.do_a_training_loop(training_size, scoreMatrix, alphas, betas, IDs=IDs)
        print('do_a_training_loop_ID elapsed: ', time.time() - s)
        print(m1.round(2))
        print(expected_matrix)

        np.testing.assert_array_equal(m1.round(2), expected_matrix)

class TestCheckCorrectPairs:
    def test_check_correct_pairs_IDisNone(self): # no ID provided
        alphas = ['martina', 'arianna', 'deborah', 'deborah', 'lucia']
        betas = ['tom', 'joe', 'marco', 'max', 'marco']

        paired = [(0, 1), (3,4), (1,0), (4,2), (2,3)]
        correctPairs = ['matrina::tom', 'arianna::joe', 'deborah::max', 'deborah::marco', 'lucia::marco']

        is_correct = np.array([False, True, False, True, True])

        np.testing.assert_array_equal(ipa.check_correct_pairs_ID(paired, alphas, betas, correctPairs), is_correct)

    def test_check_correct_pairs_ID(self):
        alphas = ['marti', 'ari', 'deb',   'deb', 'lucia', 'deb',   'marti', 'ari']
        betas =  ['tom',   'joe', 'marco', 'max', 'marco', 'marco', 'max',   'joe']
        ids =    ['1',      '2',  '1',     '1',    '2',     '2',    '2',     '2']
        IDs = [ids, ids]

        paired = [(0,0,'1'),(4,1,'2'),(2,3,'1'),(3,2,'1'),(5,4,'2'),(6,5,'2'),(1,6,'1'),(7,7,'2'),(5,6,'2')] 
        # we have an extra pair compared to a real thing because I want to check deb, max in second ID would be found wrong
        correctPairs = ['marti::tom::1', 'ari::joe::2', 'deb::marco::1', 'deb::max::1', 'lucia::marco::2', 'deb::marco::2', 'marti::max::2', 'ari::joe::2']

        is_correct = np.array([True, False, True, True, True, False, False, True, False])

        np.testing.assert_array_equal(ipa.check_correct_pairs_ID(paired, alphas, betas, correctPairs, IDs), is_correct)

class TestUseCorrelation:
    def test_multiply_scoring_corr(self):

        scorematrix = np.array([[0,1,-1],[0,0,0], [0,-1,1]])
        correlation = np.array([[0,1,-1],[0,0,0], [0,1,-1]])
        expected = np.array([[0, 1, -1], [0,0,0], [0,-1,-1]])

        np.testing.assert_array_equal(ipa.multiply_corr(scorematrix, correlation), expected)
class TestSaveResults:
    def test_save_results_to_df(self):
        alphas = [['a','-'], ['a', '1'], ['a','2']]
        betas = [['b','0'], ['b', '1'], ['b','-']]
        correctpairs = ['a::b0', 'a1::b1', 'a2::b']
        confScores = np.array([
            [1,0,0],
            [0,0,2],
            [0,3,0]
        ], dtype=float)

        expected = pd.DataFrame(data = [['a-', 'b0', 1., 5, 1, 1, 3, True],
                                        ['a1', 'b0', 0., 5, 0, 0, 3, np.nan],
                                        ['a2', 'b0', 0., 5, 0, 0, 3, np.nan],
                                        ['a-', 'b1', 0., 5, 0, 0, 3, np.nan],
                                        ['a1', 'b1', 0., 5, 0, 1, 3, np.nan],
                                        ['a2', 'b1', 3., 5, 1, 0, 3, False],
                                        ['a-', 'b-', 0., 5, 0, 0, 3, np.nan],
                                        ['a1', 'b-', 2., 5, 1, 0, 3, False],
                                        ['a2', 'b-', 0., 5, 0, 1, 3, np.nan]],
                                columns = ['alpha', 'beta', 'confScore', 'iteration', 'paired', 'real_pairs', 'n_repeat', 'correct'])
        
        pd.testing.assert_frame_equal(ipa.save_results_to_df(confScores, alphas, betas, correctpairs, 5, 3), expected)
    
    def test_save_results_to_df_withIDs(self):
        ids = [1,2,1]
        IDs = [ids, ids]
        alphas = [['a','-'], ['a', '1'], ['a','2']]
        betas = [['b','0'], ['b', '1'], ['b','-']]
        correctpairs = ['a::b0::1', 'a1::b1::2', 'a2::b::1']
        confScores = np.array([
            [0,0,1],
            [0,2,0],
            [3,0,0]
        ], dtype=float)

        expected = pd.DataFrame(data = [['a-', 1, 'b0', 1, 0., 5, 0, 1, 3, np.nan],
                                        ['a1', 2, 'b0', 1, 0., 5, 0, 0, 3, np.nan],
                                        ['a2', 1, 'b0', 1, 3., 5, 1, 0, 3, False],
                                        ['a-', 1, 'b1', 2, 0., 5, 0, 0, 3, np.nan],
                                        ['a1', 2, 'b1', 2, 2., 5, 1, 1, 3, True],
                                        ['a2', 1, 'b1', 2, 0., 5, 0, 0, 3, np.nan],
                                        ['a-', 1, 'b-', 1, 1., 5, 1, 0, 3, False],
                                        ['a1', 2, 'b-', 1, 0., 5, 0, 0, 3, np.nan],
                                        ['a2', 1, 'b-', 1, 0., 5, 0, 1, 3, np.nan]],
                                columns = ['alpha', 'alpha_ID', 'beta', 'beta_ID', 'confScore', 'iteration', 'paired', 'real_pairs', 'n_repeat', 'correct'])
        print(expected)
        print()
        pd.testing.assert_frame_equal(ipa.save_results_to_df(confScores, alphas, betas, correctpairs, 5, 3, IDs), expected)