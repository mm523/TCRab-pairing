# Implementation of Seok and Kang, 2015

from scipy.stats import chi2
import numpy as np
from itertools import product
from scipy.optimize import minimize
import pandas as pd

'''Recursive partitioning is implemented as in Seok and Kang, Sci Rep, 2015 (https://www.nature.com/articles/srep10981). Legacy code, not included in final models for publication.'''

def _calculate_coincidence(X, mytuple):
        aa1, aa2 = mytuple
        return(X[(X[:,0] == aa1) & (X[:,1] == aa2), :].shape[0])

def CalculateNofCoincidence(X, AAs = None):
    '''
    This is N(xi, yj). It can take a list of things I expect to see in each column (e.g. a list of all AAs).
    This is pretty much the same code as ipa.pseudocounter_double, but I have removed the extra arguments and 
    the division to get frequency
    '''
    if AAs == None:
        AAs = list("ARNDCEQGHILKMFPSTWYV-")
    X = np.array(X)
    c = np.zeros(shape=(len(AAs), len(AAs)))

    for aa1, aa2 in product(AAs, AAs):
        x = AAs.index(aa1)
        y = AAs.index(aa2)
        c[x, y] = _calculate_coincidence(X, (aa1, aa2))

    assert c.sum() == X.shape[0]

    return(c)

def CalculateMarginals(coincidences, side = 'X'):
    '''
    side indicates either alpha or beta, coded as X or Y
    '''
    assert side in ['X', 'Y'], 'Side must be in [\'X\', \'Y\']'

    if side == 'X':
        axiskey = 1
    else:
        axiskey = 0

    marginals = coincidences.sum(axis = axiskey)
    return(marginals)

def SortCoincidenceArray(coincidences):
    mL = CalculateMarginals(coincidences, side = 'X')
    mR = CalculateMarginals(coincidences, side = 'Y')

    mLsorted = np.argsort(mL)
    mRsorted = np.argsort(mR)

    coincidences = coincidences[mLsorted, :]
    coincidences = coincidences[:, mRsorted]

    return(coincidences, [mLsorted, mRsorted])

def ChiSqSum(CA):
    '''
    Calculate chi squared sum to test for uniform distribution.
    '''
    m = CA.shape[0]
    n = CA.shape[1]
    Expected = CA.sum()/(m*n)
    chiM = (CA - Expected)**2
    chi = chiM.sum()/Expected    
    df = n*m - 1

    return chi, df

def ChiSqPVal(CoincidenceArray):
    chi, DF = ChiSqSum(CoincidenceArray)
    pVal = chi2.sf(chi, DF)
    return(pVal)

def residuals(M, c):
    return(((M - c)**2).sum())

def Tfunction(params, M, s):    
    '''
    This is the function we minimise to find the optimal partition
    '''
    c1, c2 = params
    M = list(sorted(M))
    M1 = np.array(M)[0:int(s+1)]
    M2 = np.array(M)[int(s+1):len(M)+1]
    T = residuals(M1, c1) + residuals(M2, c2)
    return(T)

def minimiseT(M, s):
    initial_guess = [0,0]
    myresult = minimize(Tfunction, initial_guess, args=(M, s))
    return(myresult.fun, myresult.x)

def findOptimalS(M):
    S = list(range(M.shape[0] - 1))
    Ts = {minimiseT(M, s)[0]:s for s in S}
    minT = min(Ts.keys())
    s = Ts[minT]

    return(minT, s)

def findPartitionPoint(coincidences):
    
    marginalsX = CalculateMarginals(coincidences, 'X')
    marginalsY = CalculateMarginals(coincidences, 'Y')

    if coincidences.shape[0] < 2:
        side = 'Y'
        T_R, sR = findOptimalS(marginalsY)
        s = sR
    elif coincidences.shape[1] < 2:
        side = 'X'
        T_L, sL = findOptimalS(marginalsX)
        s = sL
    else:
        T_L, sL = findOptimalS(marginalsX)
        T_R, sR = findOptimalS(marginalsY)
        if T_L > T_R:
            # print('Optimal split found on the left-hand side')
            side = 'X'
            s = sL
        else:
            # print('Optimal split found on the right-hand side')
            side = 'Y'
            s = sR
    
    return(side, s)

def TestOnePartition(coincidenceArray, pval, setSize):
    '''
    The coincidenceArray here will be the subset which was chosen
    in the partition before. 
    In the first instance, it will be the full coincidence matrix,
    then it will become smaller and smaller
    '''
    
    if (coincidenceArray.shape[0] == 1) and (coincidenceArray.shape[1] == 1):
        # Step 0: check shape. If shape == (1,1) then we proceed with norm
        final = 1
    elif coincidenceArray.sum().sum() == 0:
        # Step 0.1: check at least 1 obs is > 0
        final = 2 # we will return an array of 0s
    else:
        # Step 1: test for uniformity

        ChiPVal = ChiSqPVal(np.array(coincidenceArray))
        if ChiPVal > pval:
            # print('This partition is uniformly distributed')
            final = 1
            
        else:
            # print('This partition is not yet uniform')
            final = 0
    
    if final == 1:
        # Step2: calculate the total number of pairs in this partition
        CountPartition = np.array(coincidenceArray).sum()
        PrA = CountPartition/setSize
        # Step 3: each value in this coincidenceArray will be 
        #substituted with the average
        k = coincidenceArray.shape[0]
        l = coincidenceArray.shape[1]
        probEachPair = PrA/(k*l)
        cA = np.zeros(shape = (k, l)) + probEachPair
        cA = pd.DataFrame(cA, index=coincidenceArray.index, columns=coincidenceArray.columns)
        # Step 4: return new array
        return(cA, final)
    else:
        if final == 2:
            final = 1 # changing this because I want the table to go into the final_array (see fn below)
        return(coincidenceArray, final)

def _get_partitioned_array(DF):
    array = np.array(DF)
    PartitionSide, PartitionPoint = findPartitionPoint(array)
    # print('Partition: ', PartitionSide, PartitionPoint)
    sortedArray, (idxX, idxY) = SortCoincidenceArray(array)
    sortedArray = pd.DataFrame(sortedArray, index = DF.index[idxX], columns=DF.columns[idxY])

    if PartitionSide == 'Y':
        # split on the columns
        cols = sortedArray.columns.values
        new_cols0 = cols[0 : PartitionPoint + 1]
        new_cols1 = cols[PartitionPoint + 1 : len(cols)]
        new_array0 = sortedArray[new_cols0]
        new_array1 = sortedArray[new_cols1]
        assert new_array0.shape[1] + new_array1.shape[1] == array.shape[1]
    else:
        # split on the indices
        idxs = sortedArray.index.values
        new_idxs0 = idxs[0:PartitionPoint + 1]
        new_idxs1 = idxs[PartitionPoint + 1 : len(idxs)]
        new_array0 = sortedArray.loc[new_idxs0]
        new_array1 = sortedArray.loc[new_idxs1]
        assert new_array0.shape[0] + new_array1.shape[0] == array.shape[0]

    return(new_array0, new_array1)

def CoincidenceArrayForPositionPair(Ai, Bj, AAs = None, pval = 0.05, setSize = None):

    if AAs == None:
        AAs = list('ARNDCEQGHILKMFPSTWYV-')

    assert type(setSize) == int, 'setSize must be a defined integer'
    assert setSize == Ai.shape[0] == Bj.shape[0], 'I expect vector of aa to be the same at set size'
    X = np.stack([np.array(Ai), np.array(Bj)]).T

    InitialCoincidenceArray = CalculateNofCoincidence(X, AAs = AAs)
    InitialCoincidenceArray = pd.DataFrame(InitialCoincidenceArray, index=AAs, columns=AAs)
    # print('STARTING WITH:')
    # print(InitialCoincidenceArray)

    final = 0 
    it = 0
    not_final_arrays = {it:[InitialCoincidenceArray]}
    final_arrays = []

    while len(not_final_arrays[it]) > 0:
        # print('###################################################################')
        # print('Iteration: ', it, ', arrays to work on: ', len(not_final_arrays[it]))

        Target = len(not_final_arrays[it]) - 1 # I need this to make sure I don't exit this for loop too soon
        not_final_arrays[it + 1] = []
        for idx, array in enumerate(not_final_arrays[it]):
            # print(array)
            array, final = TestOnePartition(array, pval, setSize)
            # print('Is final? ', final)

            if final == 1:
                final_arrays.append(array)
            else:                
                new_array0, new_array1 = _get_partitioned_array(array)

                assert (new_array0.shape[0] > 0) and (new_array0.shape[1] > 0)
                not_final_arrays[it+1].append(new_array0)
                assert (new_array1.shape[1] > 0) and (new_array1.shape[1] > 0)
                not_final_arrays[it+1].append(new_array1)

            if idx == Target:
                it += 1
    

    assert len(set([i for x in final_arrays for i in x.index.values])) == InitialCoincidenceArray.shape[0]
    assert len(set([i for x in final_arrays for i in x.columns.values])) == InitialCoincidenceArray.shape[1]

    merged = pd.DataFrame()
    # once everything is in final arrays
    for i in range(len(final_arrays)):
        df = final_arrays[i]
        merged = merged.add(df, fill_value=0)

    merged = np.array(merged.reindex(AAs, axis = 0).reindex(AAs, axis=1)) # I need this to make sure it comes out in the correct order
    marginalsX = np.array(CalculateMarginals(merged, side = 'X'))
    marginalsY = np.array(CalculateMarginals(merged, side = 'Y'))

    try:
        assert round(merged.sum(),3) == 1., 'Sum of Pij is ' + str(round(merged.sum(), 3))
        assert round(marginalsX.sum(), 3) == 1.
        assert round(marginalsY.sum(), 3) == 1.
    except Exception as e:
        # print(merged)
        # print(marginalsX)
        # print(marginalsY)
        raise ValueError(e)

    return(merged, marginalsX, marginalsY)

