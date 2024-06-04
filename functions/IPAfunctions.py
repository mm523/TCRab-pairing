from itertools import product
import numpy as np
import functions.IPAscoring as IPAscoring
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import warnings
from scipy.spatial.distance import hamming
import pandas as pd

def pseudocounter_single(X, L = 0.15, q = 21, AAs = None, weights = None):
    if AAs == None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
    k = L/q
    k1 = 1-L
    X = np.array(X)
    if weights is None:
        c = np.array([X[X == aa].shape[0]/X.shape[0] for aa in AAs])
    else:
        assert len(weights) == len(X)
        N = sum(weights)
        c = np.array([sum(weights[X == aa])/N for aa in AAs])

    pc = k + k1*c

    return pc

def _calculate_coincidence(X, mytuple, weights = None):
        aa1, aa2 = mytuple
        if weights is None:
            return(X[(X[:,0] == aa1) & (X[:,1] == aa2), :].shape[0]/X.shape[0])
        else:
            assert len(weights) == X.shape[0]
            N = sum(weights)
            return(sum(weights[(X[:,0] == aa1) & (X[:,1] == aa2)])/N)

def pseudocounter_double(X, L = 0.15, q = 21, AAs = None, weights = None):
    if AAs == None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
    k = L/(q**2)
    k1 = 1-L
    X = np.array(X)
    c = np.zeros(shape=(len(AAs), len(AAs)))

    for aa1, aa2 in product(AAs, AAs):
        x = AAs.index(aa1)
        y = AAs.index(aa2)
        c[x, y] = _calculate_coincidence(X, (aa1, aa2), weights=weights)
    
    pc = k + k1*c
    assert pc.shape == (len(AAs), len(AAs))

    return pc


def calculate_PMIij(Ai, Bj, AAs = None, method = None, L=0.15, weights = None):
    '''
    Pointwise mutual information, as defined in Bitbol, 2018.
    Produces a df, which we will use to store the PMIs from training set.
    Ai is a column in the alignment of alphas
    Bj is a column in the alignment of betas
    '''

    if AAs == None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
    assert method in [None], 'Method must be one of [None]'

    assert len(Ai) == len(Bj) # each will be paired, so this should be true
    X = np.stack([np.array(Ai), np.array(Bj)]).T

    if method == None:
        fis = pseudocounter_single(Ai, L = L, q = len(AAs), AAs = AAs, weights=weights)
        fjs = pseudocounter_single(Bj, L = L, q = len(AAs), AAs = AAs, weights=weights)
        fijs = pseudocounter_double(X, L = L, q = len(AAs), AAs = AAs, weights=weights)

    exp = np.dot(fis[:,None],fjs[None,:])
    PMIij = np.log(fijs) - np.log(exp)

    return PMIij

def _calculate_weights(A, B, weights):
    assert (weights in ['yes', None, 'no']) | (type(weights) == float), 'weights should be one of [yes, None, no] or a float'
    if weights == 'no':
        weights = None

    if isinstance(weights, float):
        X = np.concatenate([A, B], axis = 1)
        assert X.shape == (A.shape[0], A.shape[1] + B.shape[1])
        seqs = [''.join(row) for row in X]
        
        weights = [1/sum([int(hamming(list(x), list(y)) <= weights) for y in seqs]) for x in seqs]
        print('effective set size: ', sum(weights))
        weights = np.array(weights)
    
    return(weights)

def calculate_PMI_AllPositions(A, B, AAs = None, method = None, L=0.15, weights = None):
    '''
    Calculates the pairwise MI matrix starting from pairs (random or correct).
    I am assuming A and B to be arrays where each row is a sequence and each column is a position.
    All sequences within A and B need to be aligned and have the same length.
    For each position, we will calculate PMI
    '''

    weights = _calculate_weights(A, B, weights)

    A = np.array(A)
    B = np.array(B)
    pairlist = list(product(range(A.shape[1]), range(B.shape[1])))

    executor = ProcessPoolExecutor(10) # to control total number of processes when parallelised
    future = [executor.submit(calculate_PMIij, A[:,i], B[:,j], AAs, method = method, L=L, weights = weights) for i, j in pairlist] 
    PMI = {pairlist[x]: f.result() for x,f in enumerate(future)}
    executor.shutdown()

    return PMI 

def _PMI_calculation(A, B, i, j, AAs=None, L=0.15, method = None, weights = None):
    assert method in [None], 'Method must be one of [None]'
    if AAs == None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
    Ai = A[:,i]
    Bj = B[:,j]
    assert len(Ai) == len(Bj) # each will be paired, so this should be true
    X = np.stack([np.array(Ai), np.array(Bj)]).T

    if method == None:
        fijs = pseudocounter_double(X, L = L, q = len(AAs), AAs = AAs, weights=weights)
        fis = pseudocounter_single(Ai, L = L, q = len(AAs), AAs = AAs, weights=weights)
        fjs = pseudocounter_single(Bj, L = L, q = len(AAs), AAs = AAs, weights=weights)

    exp = np.dot(fis[:,None],fjs[None,:])
    assert fijs.shape == exp.shape
    where = np.where(fijs > 0)
    MI = (fijs*(np.log(fijs) - np.log(exp)))[where[0], where[1]].sum()
    if round(MI, 15) < 0: # np floating is accurate up to 15 digits, so I will only be checking those
        warnings.warn('MI was less than 0: ' + str(MI))

    return(MI)

def calculate_MI_Pairwise(A, B, AAs = None, L = 0.15, method = None, weights=None):
    assert (weights in ['yes', None, 'no']) | (type(weights) == float), 'weights should be one of [yes, None, no] or a float'
    if weights == 'no':
        weights = None

    if isinstance(weights, float):
        X = np.concatenate([A, B], axis = 1)
        assert X.shape == (A.shape[0], A.shape[1] + B.shape[1])
        seqs = [''.join(row) for row in X]
        
        weights = [1/sum([int(hamming(list(x), list(y)) <= weights) for y in seqs]) for x in seqs]
        print('effective set size: ', sum(weights))
        weights = np.array(weights)

    MIAllPositions = {}
    pairlist = list(product(range(A.shape[1]), range(B.shape[1])))

    executor = ProcessPoolExecutor(10) # to control total number of processes when parallelised
    future = [executor.submit(_PMI_calculation, A, B, i, j, L = L, AAs=AAs, method=method, weights=weights) for i, j in pairlist] 
    MI = [f.result() for f in future]
    executor.shutdown()

    MIAllPositions = np.array(MI).reshape(A.shape[1], B.shape[1])
    
    return(MIAllPositions)

def _calculate_partial_score(PMI, AAs, a, b, i, j):
    aa1 = AAs.index(a[i])
    aa2 = AAs.index(b[j])
    return(PMI[(i, j)][aa1, aa2])

def calculate_Sab(a, b, PMI, AAs = None):
    '''
    This is the pairing score, which is generated by looking at the PMI table.
    Produced by calculate_PMI_table
    '''

    if AAs is None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")

    Sab = 0
    for (i, j) in product(range(len(a)), range(len(b))):
        Sab += _calculate_partial_score(PMI, AAs, a, b, i, j)
    
    return Sab

def _ID_Sab(alphas, betas, new_pmi, IDs, id, AAs):
    '''
    Calculates Sab scores only within an individual
    '''

    a_indices = [i for i, x in enumerate(IDs[0]) if x == id]
    b_indices = [i for i, x in enumerate(IDs[1]) if x == id]
    mx = np.zeros(shape = (len(a_indices), len(b_indices)))
    for i, j in product(a_indices, b_indices):
        mx[a_indices.index(i), b_indices.index(j)] = calculate_Sab(alphas[i], betas[j], new_pmi, AAs = AAs)
    return(mx, np.array(a_indices), np.array(b_indices))

def get_new_training_set(scoreMatrix, trainSize):
    '''
    Get idx of the pairs that should be in training set.
    The score matrix can be both Sab scores, as well as confidence scores 
    '''
    
    pairs = IPAscoring.get_ranked_pairs(scoreMatrix)
    new_training = pairs[0:trainSize]
    assert len(new_training) == trainSize

    a = [x[0] for x in new_training] # index of alphas to be used
    b = [x[1] for x in new_training] # index of beta to be used

    return(a, b) 

def create_score_matrix(alphas, betas, IDs, pmi, AAs):
    matrix_1 = np.zeros(shape = (alphas.shape[0], betas.shape[0]))
    matrix_1[matrix_1 == 0] = np.nan

    if IDs == None:
        for i, j in product(range(alphas.shape[0]), range(betas.shape[0])):
            s = calculate_Sab(alphas[i], betas[j], pmi, AAs)
            matrix_1[i,j] = s
    else:
        executor = ProcessPoolExecutor(10)
        future = [executor.submit(_ID_Sab, alphas, betas, pmi, IDs, id, AAs) for id in set(IDs[0] + IDs[1])] 
        for f in future:
            matrix_1[np.ix_(f.result()[1], f.result()[2])] = f.result()[0]
        executor.shutdown()
    
    return(matrix_1)

def do_a_training_loop(training_size, scoreMatrix, alphas, betas, background_alpha = None, background_beta = None, AAs = None, IDs = None, method = None, L=0.15, weights = 'no'):
    '''
    Takes the model into 1 spin of learning. 
    It applies all the functions in the correct sequence and returns the new matrix with the scores. 
    It takes the pre-existing score matrix as input, as well as the list of alpha and beta sequences.
    The score matrix can be both Sab scores, as well as confidence scores, depending on what you want to do.
    '''

    # assert weights in ['yes', None, 'no'], 'weights should be one of [yes, None, no]'
    if IDs is not None:
        assert len(IDs) == 2, 'IDs should be a list of lists [alphaIDs, betaIDs]'

    if AAs is None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")

    a, b = get_new_training_set(scoreMatrix, training_size)

    alphas_1 = alphas[a]
    betas_1 = betas[b]

    if background_alpha is not None:
        assert background_beta is not None
        alphas2 = np.concatenate((background_alpha, alphas_1), axis = 0)
    else:
        alphas2 = alphas_1
    if background_beta is not None:
        assert background_alpha is not None
        betas2 = np.concatenate((background_beta, betas_1), axis = 0)
    else:
        betas2 = betas_1

    new_pmi = calculate_PMI_AllPositions(alphas2, betas2, AAs = AAs, method = method, L=L, weights = weights)

    matrix_1 = create_score_matrix(alphas, betas, IDs, new_pmi, AAs)
    
    return(matrix_1)  

def check_correct_pairs_ID(pairs, alphas, betas, correctPairs, IDs=None):
    if IDs:
        pairs = ['::'.join((''.join(alphas[pairs[i][0]]).replace('-',''), ''.join(betas[pairs[i][1]]).replace('-',''), str(pairs[i][2]))) for i in range(len(pairs))]
    else:
        pairs = ['::'.join((''.join(alphas[pairs[i][0]]).replace('-',''), ''.join(betas[pairs[i][1]]).replace('-',''))) for i in range(len(pairs))]
    correct = [p in correctPairs for p in pairs]
    return(correct)

def multiply_corr(scorematrix, correlation):
    '''
    I use this function to multiply the scoring matrix with the correlation matrix because I want to 
    have a negative results whenever either matrix is negative and I don't want two negatives to give me a positive
    '''

    negative_score = scorematrix < 0
    negative_corr = correlation < 0

    should_be_negative = negative_score | negative_corr

    scores = abs(scorematrix*correlation)

    scores[should_be_negative] = scores[should_be_negative]*-1
    
    return(scores)

def save_results_to_df(confScores, alphas, betas, correctPairs, iteration, n_repeat, IDs=None):
    pairs = IPAscoring.get_ranked_pairs(confScores, IDs)
    check_correct = check_correct_pairs_ID(pairs, alphas, betas, correctPairs, IDs)
    num_correct_pairs = sum([int(x) for x in check_correct])
    print('Number of correct pairs: ', num_correct_pairs)

    a = [''.join(alpha) for alpha in alphas]
    b = [''.join(beta) for beta in betas]

    if IDs is None:
        res = pd.melt(pd.DataFrame(confScores, index = a, columns = b), ignore_index=False).reset_index()
        res.columns = ['alpha', 'beta', 'confScore']
        res['pair_ids'] = res.apply(lambda row: 
                        row['alpha'].replace('-','') + '::' + row['beta'].replace('-',''), axis=1).tolist()
    else:
        res = pd.melt(pd.DataFrame(confScores, index = [a, IDs[0]], columns = [b, IDs[1]]), ignore_index=False).reset_index()
        res.columns = ['alpha', 'alpha_ID', 'beta', 'beta_ID', 'confScore']
        res['pair_ids'] = res.apply(lambda row: 
                        row['alpha'].replace('-','') + '::' + row['beta'].replace('-','') + '::' + str(row['beta_ID']) \
                            if row['beta_ID'] == row['alpha_ID'] else np.nan, axis=1).tolist()
    res['iteration'] = iteration

    paired = np.zeros(shape=(len(alphas), len(betas)), dtype=int)
    correct = np.empty(shape=(len(alphas), len(betas)), dtype=bool)
    for i, p in enumerate(pairs):
        a_pos = p[0]
        b_pos = p[1]
        paired[a_pos, b_pos] = 1
        correct[a_pos, b_pos] = check_correct[i]
        
    res['paired'] = pd.melt(pd.DataFrame(paired))['value'].tolist()
    res['real_pairs'] = res['pair_ids'].isin(correctPairs).astype('int')
    assert res['real_pairs'].sum() >= len(correctPairs)
    print(res['paired'].sum(), len(correctPairs))
    assert res['paired'].sum() <= len(correctPairs)
    if res['paired'].sum() != len(correctPairs):
        warnings.warn('There are fewer paired examples than correct pairs. This is expected with GA+MI-IPA when GA is not re-paired.')
    res['n_repeat'] = n_repeat
    res['correct'] = pd.melt(pd.DataFrame(correct))['value'].tolist()
    res.loc[res['paired'] == 0, 'correct'] = None
    assert res['correct'].sum() == num_correct_pairs
    assert res['correct'].sum() <= len(correctPairs)
    res = res.drop('pair_ids', axis = 1)
    
    if IDs is not None:
        # check that all pairs labelled correct are within the same ID
        X = res.loc[res['real_pairs'] == 1]
        if X['alpha_ID'].equals(X['beta_ID']):
            pass
        else:
            print(X[['alpha_ID', 'beta_ID']])
            raise ValueError('alpha and beta IDs do not correspond for correct pairs.')

    return(res)