import numpy as np
from concurrent.futures import ProcessPoolExecutor
from warnings import warn
import random
from munkres import Munkres

# Greedy algorithm implementation, as per Bitbol, 2016

def score(E1, E2, n) :
    score = (E1 - E2)/(n + 1)
    return(score)

def greedy_scoring(mx):
    if mx.shape[1]>1:
        confidenceScores = np.zeros(shape = (mx.shape[0], mx.shape[1]))
        confidenceScores[confidenceScores == 0] = np.nan

        alist = list(range(mx.shape[0])) # used to keep track of the idx in the real starting matrix
        blist = list(range(mx.shape[1]))
        mx0 = mx.copy()

        for it, m in enumerate(range(mx.shape[0])):
            if m < (mx.shape[0] - 1):
                indices = np.where(mx0 == mx0.max()) # find max score of table
                # print(indices)
                # use same procedure as in get_ranked_pairs to resolve ties to not inflate result
                if len(indices[0]) == len(indices[1]) == 1:
                    x, y = indices[0][0], indices[1][0]   
                else:
                    warn('Tie detected! Picking a random pair')
                    idx = random.sample(range(len(indices[0])),1)[0]
                    x, y = indices[0][idx], indices[1][idx]
                
                
                # find second best beta for alpha
                Xscores = mx0[x, :].tolist()
                Xscores.sort(reverse = True)

                assert Xscores[0] == mx0.max()
                assert Xscores[0] - Xscores[1] >= 0

                new_x = alist[x] # get alpha number  
                new_y = blist[y] # get beta number

                if it == 0:
                    PairScore = score(Xscores[0], Xscores[1], 0)
                else:
                    # figure out how many with higher scores I have already removed
                    originalScores = mx[new_x, :]
                    n = len(originalScores[originalScores > Xscores[0]]) # number of partners that would have been better, but have already been removed
                    PairScore = score(Xscores[0], Xscores[1], n)

                confidenceScores[new_x, new_y] = PairScore
                
                mx0 = np.delete(mx0, x, axis= 0)
                mx0 = np.delete(mx0, y, axis= 1)
                alist.remove(new_x) # remove used index so position is the same as in matrix I am sorting
                blist.remove(new_y)
            else:
                assert min(len(alist), len(blist)) == 1

                confidenceScores[alist[0], blist[0]] = np.nanmin(confidenceScores)
    else:
        confidenceScores = mx.copy()
    
    return(confidenceScores)

### Hungarian algorithm implementation, as per Bitbol, 2018

def hungarian_assignment(profit_matrix):
    # assigns pairs to maximise scores
    # the munkres function for this minimises the cost function
    # (not using scipy function because slower)
    # because I want to maximise my matrix, I call it a profit matrix
    # the cost matrix is the max value of this matrix - the matrix

    #Using the maximum value of the profit_matrix to get the corresponding cost_matrix
    max_value = np.nanmax(profit_matrix)
    #Using the cost matrix to find which positions are the answer
    cost_matrix = max_value - profit_matrix
    # add very small random numbers to the cost matrix to solve ties
    white_noise = np.random.uniform(0,10**(-7),len(cost_matrix.ravel())).reshape(cost_matrix.shape)
    cost_matrix = cost_matrix + white_noise
    # print(cost_matrix)

    # optimal_assignment = linear_sum_assignment(cost_matrix)
    m = Munkres()
    a = m.compute(cost_matrix)
    # print(a)
    optimal_assignment = [np.array([x[0] for x in a]), np.array([x[1] for x in a])]
    # print(optimal_assignment)

    return(optimal_assignment)

def _hungarian_score_suboptimal(mx, x, y):
    if mx.shape[1] > 1:
        mx1 = mx.copy()
        mx1[x, y] = -abs(np.nanmin(mx1)+0.1)*10**7 # this gives the pair the lowest possible score, so it is not chosen

        suboptimal_assignment = hungarian_assignment(mx1)
        assert (x, y) not in zip(suboptimal_assignment[0], suboptimal_assignment[1])
        score_suboptimal = mx1[suboptimal_assignment[0], suboptimal_assignment[1]].sum()
    else:
        score_suboptimal = 0
    
    return(score_suboptimal)

def hungarian_confidence(mx, optimal_assignment):

    confidenceScores = np.zeros(shape = (mx.shape[0], mx.shape[1]))
    confidenceScores[confidenceScores == 0] = np.nan

    score = mx[optimal_assignment[0], optimal_assignment[1]].sum()
    pairlist = list(zip(optimal_assignment[0], optimal_assignment[1]))

    executor = ProcessPoolExecutor(5)
    future = [executor.submit(_hungarian_score_suboptimal, mx, x, y) for x,y in pairlist] 
    scores = [score - f.result() for f in future]
    executor.shutdown()
    for n, pair in enumerate(pairlist):
        confidenceScores[pair[0], pair[1]] = scores[n]
    return(confidenceScores)

def hungarian_scoring(mx):

    optimal_assignment = hungarian_assignment(mx)
    confidenceScore = hungarian_confidence(mx, optimal_assignment)
    
    return(confidenceScore)

def get_ranked_pairs(mx, IDs=None):
    # Get ranked pairs from confidence tables
    # here mx will contain the confidence scores
    # this is the function that before was called "get_top_scoring"
    # this is not an efficient way to do this if the matrix is mostly NaN (as it is for greedy alg)
    # to be edited later if a format which takes the non-NaN from each row and column is enough

    pairs = []

    alist = list(range(mx.shape[0])) # used to keep track of the idx in the real starting matrix
    blist = list(range(mx.shape[1]))

    if IDs:
        individuals_a, individuals_b = IDs
        # print(individuals_a)

    for m in range(mx.shape[0]): # the loop does it once for each a
        indices = np.where(mx == np.nanmax(mx))

        # sometimes indices has a tie.
        # initially, I was taking the first each time. However, this artificially inflates benchmarking results
        # when I use sorted tables for alpha and beta
        if len(indices[0]) == len(indices[1]) == 1:
            idx = 0    
        else:
            warn('Tie detected! Picking a random pair')
            idx = random.sample(range(len(indices[0])),1)[0]
        new_a = alist[indices[0][idx]]
        new_b = blist[indices[1][idx]] 

        if IDs:
            a_id = individuals_a[new_a]
            b_id = individuals_b[new_b]
            assert a_id == b_id, 'Individual ID does not correspond!'
            pair_id = a_id
            pairs.append((new_a, new_b, pair_id))
        else:
            pairs.append((new_a, new_b))

        mx = np.delete(mx, indices[0][idx], axis= 0) # to avoid using the same sequence in 2 pairs
        mx = np.delete(mx, indices[1][idx], axis= 1)
        alist.remove(new_a) # remove used index so position is the same as in matrix I am sorting
        blist.remove(new_b)

    return(pairs)

def _ID_confidence(mx, IDs, id, scorefx):
    a_indices = [i for i, x in enumerate(IDs[0]) if x == id]
    b_indices = [i for i, x in enumerate(IDs[1]) if x == id]
    submatrix = mx[np.ix_(a_indices, b_indices)]
    scores = scorefx(submatrix)
    # print(id, scores)

    return(scores, a_indices, b_indices)

def IPAscoring(mx, scoringtype, IDs = None):
    if IDs is not None:
        assert len(IDs) == 2, 'IDs should be a list of lists [alphaIDs, betaIDs]'
    if scoringtype not in ['hungarian', 'greedy']:
        raise ValueError('scoringtype must be one of [hungarian, greedy]')

    if scoringtype == 'greedy':
        scorefx = greedy_scoring
    elif scoringtype == 'hungarian':
        scorefx = hungarian_scoring

    if IDs == None:
        confscores = scorefx(mx)        
        return(confscores)
    
    else:
        assert len(IDs[0]) == mx.shape[0]
        assert len(IDs[1]) == mx.shape[1]
        confScores = np.zeros(shape=(len(IDs[0]), len(IDs[1])))
        confScores[confScores == 0] = np.nan

        executor = ProcessPoolExecutor(5)
        future = [executor.submit(_ID_confidence, mx, IDs, sID, scorefx) for sID in set(IDs[0] + IDs[1])] 
        executor.shutdown()
        for f in future:
            r = f.result()
            confScores[np.ix_(r[1], r[2])] = r[0]
        # print(confScores)
        return(confScores)
