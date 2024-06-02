import functions.RecursivePartitioningFunctions as RP
import numpy as np
import pytest
import pandas as pd

@pytest.fixture
def Ai():
    Ai = np.array(
        ['A', 'A', 'A', 'B', 'C', 'C']
    )
    return(Ai)

@pytest.fixture
def Bj():
    Bj = np.array(
        ['A', 'A', 'A', 'B', 'C', 'C']
    )
    return(Bj)

@pytest.fixture
def coincidences():
    expected = np.array([
                        [3, 0, 0],
                        [0, 1, 0],
                        [0, 0, 2]
                        ])
    return expected

def test_CalculateNofCoincidence(Ai, Bj, coincidences):
    X = np.stack([np.array(Ai), np.array(Bj)]).T
    np.testing.assert_array_equal(RP.CalculateNofCoincidence(X, AAs = list('ABC')), coincidences)

def test_CalculateNofCoincidence_2(Ai, Bj):
    expectedA = 'ABCD'
    X = np.stack([np.array(Ai), np.array(Bj)]).T

    expected = np.array([
                        [3, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 2, 0],
                        [0, 0, 0, 0]
                        ])
    np.testing.assert_array_equal(RP.CalculateNofCoincidence(X, AAs = expectedA), expected)

def test_CalculateNofCoincidence_3():
    Ai = np.array(
        ['A', 'A', 'A', 'R', 'N', 'N']
    )
    Bj = np.array(
        ['A', 'A', 'A', 'R', 'N', 'N']
    )
    X = np.stack([np.array(Ai), np.array(Bj)]).T

    expected = np.array([
                        [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    ], dtype = float)
    np.testing.assert_array_equal(RP.CalculateNofCoincidence(X), expected)

def test_CalculateMarginals(coincidences):
    expected = np.array([3, 1, 2])

    np.testing.assert_array_equal(RP.CalculateMarginals(coincidences, 'X'), expected)
    np.testing.assert_array_equal(RP.CalculateMarginals(coincidences, 'Y'), expected)

def test_SortCoincidenceArray(coincidences):
    expected = np.array([
                        [1, 0, 0],
                        [0, 2, 0],
                        [0, 0, 3]
                        ])
    
    sortedArray = RP.SortCoincidenceArray(coincidences)[0]
    print(sortedArray)
    print(expected)

    np.testing.assert_array_equal(sortedArray, expected)

def test_ChiSqSum(coincidences):
    exp_chi = ((1-6/9)**2 + (2-6/9)**2 + (3-6/9)**2 + 6*(0-6/9)**2)/(6/9)
    DF = 8
    assert round(exp_chi,2) == round(RP.ChiSqSum(coincidences)[0],2)
    assert DF == RP.ChiSqSum(coincidences)[1]

def test_ChiSqPval(coincidences):
    assert 0.06 == round(RP.ChiSqPVal(coincidences), 2) # calculated manually

def test_residuals():
    v = np.array([0,3,5,6,8])
    m = 2
    assert RP.residuals(v, m) == 66

def test_Tfunction(coincidences):
    marginals = RP.CalculateMarginals(coincidences, 'X')
    print(marginals)
    s = 0
    c1 = 0
    c2 = 0
    params = [c1, c2]
    assert RP.Tfunction(params, marginals, s) == 14

    c1 = 1
    c2 = 2.5
    s = 0
    params = [c1, c2]
    assert RP.Tfunction(params, marginals, s) == .5

    c1 = 1.5
    c2 = 3
    s = 1
    params = [c1, c2]
    assert RP.Tfunction(params, marginals, s) == .5

    c1 = 1
    c2 = 2
    s = 1
    params = [c1, c2]
    assert RP.Tfunction(params, marginals, s) == 2

def test_minimiseT(coincidences):
    marginals = RP.CalculateMarginals(coincidences, 'X')

    s = 0
    assert round(RP.minimiseT(marginals, s)[0], 1) == 0.5
    np.testing.assert_array_equal(RP.minimiseT(marginals, s)[1].round(1), np.array([1.0, 2.5]))

    s = 2
    assert round(RP.minimiseT(marginals, s)[0], 1) == 2
    np.testing.assert_array_equal(RP.minimiseT(marginals, s)[1].round(1), np.array([2.0, 0.]))

def test_findOptimalS():
    marginals = np.array([3]*15 + list(range(10,15)))
    print(marginals)

    assert RP.findOptimalS(marginals)[0] == 10 # min T
    assert RP.findOptimalS(marginals)[1] == 14 # chosen s

@pytest.fixture
def MartinaAriannaExampleArray():
    Ai = np.array(list('MARTINA'))
    Bj = np.array(list('ARIANNA'))
    X = np.stack([np.array(Ai), np.array(Bj)]).T
    return(X)

@pytest.fixture
def MartinaAriannaExampleCoincidences(MartinaAriannaExampleArray):
    C = RP.CalculateNofCoincidence(MartinaAriannaExampleArray)
    return(C)

def test_findPartitionPoint(MartinaAriannaExampleCoincidences):
    side, s = RP.findPartitionPoint(MartinaAriannaExampleCoincidences)

    assert side == 'Y'
    assert s == 18

def test_findPartitionPoint_SmallerSet():
    # extra test, but helps in setting up final test
    Ai = np.array(list('MARTINAMARTINAMARTINAMARTINA'))
    Bj = np.array(list('ARIANNAARIANNAARIANNAARIANNA'))
    X = np.stack([np.array(Ai), np.array(Bj)]).T
    C = RP.CalculateNofCoincidence(X, AAs = list('ARINMTCSW'))
    side, s = RP.findPartitionPoint(C)

    assert side == 'Y'
    assert s == 6

def test_findPartitionPoint_2():
    Ai = np.array(list('MMMTINA'))
    Bj = np.array(list('ARIANNA'))
    X = np.stack([np.array(Ai), np.array(Bj)]).T
    C = RP.CalculateNofCoincidence(X)
    side, s = RP.findPartitionPoint(C)

    assert side == 'X'
    assert s == 19

def test_findPartitionPoint_3():
    C = np.array([
        [0],
        [1],
        [0],
        [1],
        [1]
    ])
    side, s = RP.findPartitionPoint(C)

    assert side == 'X'
    assert s == 1

def test_findPartitionPoint_4():
    C = np.array([
        [0, 1, 0, 1, 1]
    ])
    side, s = RP.findPartitionPoint(C)

    assert side == 'Y'
    assert s == 1

def test_TestOnePartition_Uniform(coincidences):
    FinalArray, final = RP.TestOnePartition(pd.DataFrame(coincidences), 0.05, 10)

    FinalArrayExpected = pd.DataFrame(np.array([
                        [6/90, 6/90, 6/90],
                        [6/90, 6/90, 6/90],
                        [6/90, 6/90, 6/90]
                        ]))
    
    pd.testing.assert_frame_equal(FinalArray, FinalArrayExpected)
    assert final == 1

def test_TestOnePartition_NotUniform(coincidences):

    coincidences = np.array([
                        [10, 0, 0],
                        [0, 11, 0],
                        [0, 0, 2]
                        ])

    FinalArray, final = RP.TestOnePartition(pd.DataFrame(coincidences), pval = 0.05, setSize = 10)

    FinalArrayExpected = pd.DataFrame(np.array([
                        [10, 0, 0],
                        [0, 11, 0],
                        [0, 0, 2]
                        ]))
    
    pd.testing.assert_frame_equal(FinalArray, FinalArrayExpected)
    assert final == 0

def test_TestOnePartition_TooSmall():
    FinalArray, final = RP.TestOnePartition(pd.DataFrame(np.array([[1]])), pval = 0.05, setSize = 10)

    FinalArrayExpected = pd.DataFrame(np.array([
                        [0.1]
                        ]))
    
    pd.testing.assert_frame_equal(FinalArray, FinalArrayExpected)
    assert final == 1

def test_TestOnePartition_AllZeros():
    FinalArray, final = RP.TestOnePartition(np.array([[0,0,0],[0,0,0]]), pval = 0.05, setSize = 10)

    FinalArrayExpected = np.array([
                        [0,0,0],
                        [0,0,0]
                        ])
    
    np.testing.assert_array_equal(FinalArray, FinalArrayExpected)
    assert final == 1

def test_partitioning():
    AAs = list('ARINMT')
    Ai = np.array(list('MARTINAMARTINAMARTINAMARTINA'))
    Bj = np.array(list('ARIANNAARIANNAARIANNAARIANNA'))
    X = np.stack([np.array(Ai), np.array(Bj)]).T
    CA = RP.CalculateNofCoincidence(X, AAs = AAs)
    CA = pd.DataFrame(CA, index=AAs, columns=AAs)

    new_array0, new_array1 = RP._get_partitioned_array(CA)

    expected0 = pd.DataFrame(index = ['R', 'I', 'N', 'M', 'T', 'A'], 
                            columns = ['M', 'T', 'R', 'I'],
                            data = [
                                [0,0,0,4],
                                [0,0,0,0],
                                [0,0,0,0],
                                [0,0,0,0],
                                [0,0,0,0],
                                [0,0,4,0]
                                ], dtype = float)
    
    expected1 = pd.DataFrame(index = ['R', 'I', 'N', 'M', 'T', 'A'], 
                            columns = ['N', 'A'],
                            data = [
                                [0,0],
                                [4,0],
                                [4,0],
                                [0,4],
                                [0,4],
                                [0,4]
                                ], dtype = float)

    pd.testing.assert_frame_equal(new_array0, expected0)
    pd.testing.assert_frame_equal(new_array1, expected1)

    new_array2, new_array3 = RP._get_partitioned_array(new_array1)

    expected2 = pd.DataFrame(index = ['R', 'I', 'N', 'M', 'T', 'A'], 
                            columns = ['N'],
                            data = [
                                [0],
                                [4],
                                [4],
                                [0],
                                [0],
                                [0]
                                ], dtype = float)
    
    expected3 = pd.DataFrame(index = ['R', 'I', 'N', 'M', 'T', 'A'], 
                            columns = ['A'],
                            data = [
                                [0],
                                [0],
                                [0],
                                [4],
                                [4],
                                [4]
                                ], dtype = float)
    
    pd.testing.assert_frame_equal(new_array2, expected2)
    pd.testing.assert_frame_equal(new_array3, expected3)

def test_CoincidenceArrayForPositionPair():
    Ai = np.array(list('MARTINAMARTINAMARTINAMARTINA'))
    Bj = np.array(list('ARIANNAARIANNAARIANNAARIANNA'))
    
    Fij, Fi, Fj = RP.CoincidenceArrayForPositionPair(Ai, Bj, AAs = list('ARINMT'), setSize=28)

    expected_Fij = np.array([
                    [1/7,1/7,0,0,0,0],
                    [0,0,1/7,0,0,0],
                    [0,0,0,1/7,0,0],
                    [0,0,0,1/7,0,0],
                    [1/7,0,0,0,0,0],
                    [1/7,0,0,0,0,0]
                    ])        
    np.testing.assert_array_equal(Fij, expected_Fij)
    np.testing.assert_array_equal(Fi, np.array([2/7, 1/7, 1/7, 1/7, 1/7, 1/7]))
    np.testing.assert_array_equal(Fj, np.array([3/7, 1/7, 1/7, 2/7, 0, 0]))