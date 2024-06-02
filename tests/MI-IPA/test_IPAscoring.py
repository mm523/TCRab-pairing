import numpy as np
import functions.IPAscoring as IPAscoring
import time, pytest

def test_greedy_algorithm():
    # expected output worked out manually in book 12.10.2022

    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    expected_output = np.array([
        [np.nan, np.nan, np.nan, 3],
        [np.nan, np.nan, 2, np.nan],
        [np.nan, 1/3, np.nan, np.nan],
        [1/3, np.nan, np.nan, np.nan]
    ])

    confScores = IPAscoring.greedy_scoring(myarray)

    assert ((confScores == expected_output) | (np.isnan(confScores) & np.isnan(expected_output))).all()

def test_get_ranked_pairs():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])

    expected_output = [(1,2), (0,3), (2,1), (3,0)] # this happens when no ties
    assert expected_output == IPAscoring.get_ranked_pairs(myarray, IDs=None)

    myarray1 = np.array([
        [0, 5, 3, 8],
        [7, 4, 8, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])

    expected_output1 = [(1,2), (0,3), (2,1), (3,0)] # with ties - this is one possible output
    expected_output2 = [(0,3), (1,2), (2,1), (3,0)]
    with pytest.warns(UserWarning, match='Tie detected'):
        output = IPAscoring.get_ranked_pairs(myarray1, IDs=None)
    try:
        assert expected_output1 == output
    except:
        assert expected_output2 == output

def test_hungarian_scoring():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    expected_output = np.array([
        [np.nan, 4, np.nan, np.nan],
        [4, np.nan, np.nan, np.nan],
        [np.nan, np.nan, 4, np.nan],
        [np.nan, np.nan, np.nan, 4]
    ])

    confScores = IPAscoring.hungarian_scoring(myarray)
    assert confScores[0,1] == expected_output[0,1]
    assert confScores[1,0] == expected_output[1,0]
    assert confScores[2,2] == expected_output[2,2]
    assert confScores[3,3] == expected_output[3,3]
    
def test_hungarian_pairs():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    IDs = [['1','1','2','2'],['1','1','2','2']]

    expected_pairs = [(0, 1), (1, 0), (2, 2), (3, 3)] # they will be these but in a random order

    confScores = IPAscoring.hungarian_scoring(myarray)
    print(confScores)
    with pytest.warns(UserWarning, match='Tie detected'):
        pairs = IPAscoring.get_ranked_pairs(confScores, IDs=None)
        print(pairs)
    
    assert sorted(expected_pairs) == sorted(pairs)

    expected_pairs = [(0, 1, '1'), (1, 0, '1'), (2, 2, '2'), (3, 3, '2')] # they will be these but in a random order

    confScores = IPAscoring.hungarian_scoring(myarray)
    print(confScores)
    with pytest.warns(UserWarning, match='Tie detected'):
        pairs = IPAscoring.get_ranked_pairs(confScores, IDs=IDs)
        print(pairs)
    
    assert sorted(expected_pairs) == sorted(pairs)

def test_hungarian_assignment():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    expected_pairs = [(1, 0), (0, 1), (2, 2), (3,3)]

    assigned = IPAscoring.hungarian_assignment(myarray)
    assigned_pairs = [(x, y) for x, y in zip(assigned[0], assigned[1])]

    assert np.array([x in expected_pairs for x in assigned_pairs]).all()

def test_hungarian_confidence():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    expected_output = np.array([
        [np.nan, 4, np.nan, np.nan],
        [4, np.nan, np.nan, np.nan],
        [np.nan, np.nan, 4, np.nan],
        [np.nan, np.nan, np.nan, 4]
    ])

    assignment = IPAscoring.hungarian_assignment(myarray)
    s = time.time()
    confScores = IPAscoring.hungarian_confidence(myarray, assignment)
    print('hungarian_confidence elapsed: ', time.time() - s)

    assert confScores[0,1] == expected_output[0,1]
    assert confScores[1,0] == expected_output[1,0]
    assert confScores[2,2] == expected_output[2,2]
    assert confScores[3,3] == expected_output[3,3]

def test_IPAscoring():
    myarray = np.array([
        [0, 5, 3, 8],
        [7, 4, 9, 5],
        [1, 2, 7, 4],
        [1, 0, 5, 7]
    ])
    expected_output_greedy = np.array([
        [np.nan, np.nan, np.nan, 3],
        [np.nan, np.nan, 2, np.nan],
        [np.nan, 1/3, np.nan, np.nan],
        [1/3, np.nan, np.nan, np.nan]
    ])
    expected_output_hungarian = np.array([
        [np.nan, 4, np.nan, np.nan],
        [4, np.nan, np.nan, np.nan],
        [np.nan, np.nan, 4, np.nan],
        [np.nan, np.nan, np.nan, 4]
    ])

    with pytest.raises(ValueError):
        IPAscoring.IPAscoring(myarray, '')

    np.testing.assert_array_equal(IPAscoring.IPAscoring(myarray, 'hungarian'), expected_output_hungarian)
    np.testing.assert_array_equal(IPAscoring.IPAscoring(myarray, 'greedy'), expected_output_greedy)

def test_IPAscoring_with_ID():
    myarray = np.array([
        [0, 5, 3, 8, np.nan, np.nan, np.nan, np.nan],
        [7, 4, 9, 5, np.nan, np.nan, np.nan, np.nan],
        [1, 2, 7, 4, np.nan, np.nan, np.nan, np.nan],
        [1, 0, 5, 7, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, 0, 5, 3, 8],
        [np.nan, np.nan, np.nan, np.nan, 7, 4, 9, 5],
        [np.nan, np.nan, np.nan, np.nan, 1, 2, 7, 4],
        [np.nan, np.nan, np.nan, np.nan, 1, 0, 5, 7]
    ])

    ids = [1,1,1,1,2,2,2,2]
    IDs = [ids, ids]

    expected_output_greedy = np.array([
        [np.nan, np.nan, np.nan, 3, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, 2, np.nan, np.nan, np.nan, np.nan, np.nan],
        [np.nan, 1/3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        [1/3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 3],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 2, np.nan],
        [np.nan, np.nan, np.nan, np.nan, np.nan, 1/3, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, 1/3, np.nan, np.nan, np.nan],
    ])
    expected_output_hungarian = np.array([
        [np.nan, 4, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        [4, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, 4, np.nan, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, 4, np.nan, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, np.nan, 4, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, 4, np.nan, np.nan, np.nan],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 4, np.nan],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 4]
    ])

    with pytest.raises(ValueError):
        IPAscoring.IPAscoring(myarray, '', IDs = IDs)

    s = time.time()
    np.testing.assert_array_equal(IPAscoring.IPAscoring(myarray, 'hungarian', IDs), expected_output_hungarian)
    np.testing.assert_array_equal(IPAscoring.IPAscoring(myarray, 'greedy', IDs), expected_output_greedy)
    print('IPA scoring ID elapsed: ', time.time()-s)