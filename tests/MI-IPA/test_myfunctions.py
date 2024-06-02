import functions.myfunctions as mf
import numpy as np
import pandas as pd
import pytest

def test_prepare_data():
    df =  {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"}
    
    df1 = {"cdr3a_IMGTgaps":"CAS---VTF", "cdr3b_IMGTgaps":"CAST---VTF", 
            "len_cdr3a_IMGTgaps":9, "len_cdr3b_IMGTgaps":10, 
            "cdr3a_IMGTgaps_padded":"C, A, S, -, -, -, V, T, F", 
            "cdr3b_IMGTgaps_padded":"C, A, S, T, -, -, -, V, T, F",}    

    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    output = pd.DataFrame.from_dict(df1, orient = "index").T
    vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps")

    for c in output.columns:
        print(output[c])
        print(vdj1[c])

    pd.testing.assert_frame_equal(vdj1.sort_index(axis=1), output.sort_index(axis=1), check_dtype=False)

def test_prepare_data_longer():
    df = {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"}
    
    df1 = {"cdr3a_IMGTgaps":"CAS---VTF", "cdr3b_IMGTgaps":"CAST---VTF", 
            "len_cdr3a_IMGTgaps":9, "len_cdr3b_IMGTgaps":10, 
            "cdr3a_IMGTgaps_padded":"C, A, S, -, -, -, -, -, -, -, V, T, F", 
            "cdr3b_IMGTgaps_padded":"C, A, S, T, -, -, -, -, -, -, V, T, F",}    
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    output = pd.DataFrame.from_dict(df1, orient = "index").T
    print(output)
    vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 13, max_len2 = 13)

    pd.testing.assert_frame_equal(vdj1.sort_index(axis=1), output.sort_index(axis=1), check_dtype=False)

def test_prepare_data_shorter_both():
    df = {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"} 
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    with pytest.raises(Exception) as e_info:
        vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 5, max_len2 = 5)
        assert e_info == "max_len1 too stringent, less than 2 records left"

def test_prepare_data_shorter_beta():
    df = {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"} 
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    with pytest.raises(Exception) as e_info:
        vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 9, max_len2 = 9)
        assert e_info == "max_len2 too stringent, less than 2 records left"


def test_prepare_data_tcr():
    s = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"
    df =  {"cdr3a_IMGTgaps":s,
            "cdr3b_IMGTgaps":s}
    
    df1 = {"cdr3a_IMGTgaps":s, "cdr3b_IMGTgaps":s, 
            "len_cdr3a_IMGTgaps":len(s), "len_cdr3b_IMGTgaps":len(s), 
            "cdr3a_IMGTgaps_padded":", ".join(list(s)), 
            "cdr3b_IMGTgaps_padded":", ".join(list(s)),}    

    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    output = pd.DataFrame.from_dict(df1, orient = "index").T
    vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", type = "tcr")

    for c in output.columns:
        print(output[c])
        print(vdj1[c])

    pd.testing.assert_frame_equal(vdj1.sort_index(axis=1), output.sort_index(axis=1), check_dtype=False)

def test_prepare_data_longer_tcr():
    s = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"
    s_out = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLA--GGRPEQYFGPGTRLTVT"

    df = {"cdr3a_IMGTgaps":s,
            "cdr3b_IMGTgaps":s}
    
    df1 = {"cdr3a_IMGTgaps":s, "cdr3b_IMGTgaps":s, 
            "len_cdr3a_IMGTgaps":len(s), "len_cdr3b_IMGTgaps":len(s), 
            "cdr3a_IMGTgaps_padded":", ".join(list(s_out)), 
            "cdr3b_IMGTgaps_padded":", ".join(list(s_out)),}    
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    output = pd.DataFrame.from_dict(df1, orient = "index").T
    print(output)
    vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = len(s) + 2, max_len2 = len(s) + 2, type = "tcr")

    pd.testing.assert_frame_equal(vdj1.sort_index(axis=1), output.sort_index(axis=1), check_dtype=False)

def test_prepare_data_longer_tcr():
    s = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"
    s1 = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLA-GGRPEQYFGPGTRLTVT"
    s_out = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLA--GGRPEQYFGPGTRLTVT"

    df = {"cdr3a_IMGTgaps":s,
            "cdr3b_IMGTgaps":s1}
    
    df1 = {"cdr3a_IMGTgaps":s, "cdr3b_IMGTgaps":s1, 
            "len_cdr3a_IMGTgaps":len(s), "len_cdr3b_IMGTgaps":len(s1), 
            "cdr3a_IMGTgaps_padded":", ".join(list(s_out)), 
            "cdr3b_IMGTgaps_padded":", ".join(list(s_out)),}    
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    output = pd.DataFrame.from_dict(df1, orient = "index").T
    print(output)
    vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = len(s) + 2, max_len2 = len(s) + 2, type = "tcr")

    pd.testing.assert_frame_equal(vdj1.sort_index(axis=1), output.sort_index(axis=1), check_dtype=False)

def test_prepare_data_shorter_both_tcr():
    df = {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"} 
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    with pytest.raises(Exception) as e_info:
        vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 5, max_len2 = 5, type = "tcr")
        assert e_info == "max_len1 too stringent, less than 2 records left"

def test_prepare_data_shorter_beta_tcr():
    df = {"cdr3a_IMGTgaps":"CAS---VTF",
            "cdr3b_IMGTgaps":"CAST---VTF"} 
            
    vdj = pd.DataFrame.from_dict(df, orient = "index").T
    with pytest.raises(Exception) as e_info:
        vdj1 = mf.prepare_data(vdj, "cdr3a_IMGTgaps", "cdr3b_IMGTgaps", max_len1 = 9, max_len2 = 9, type = "tcr")
        assert e_info == "max_len2 too stringent, less than 2 records left"


@pytest.fixture
def json_df():
    df = pd.DataFrame(
        columns=['Method', 'Meta', 'CDR3fix', 'Foo'],
        data=[
            ['{"foo": "0", "bar":""}','{"foo": "", "bar":""}','{"foo": "", "bar":""}','{"foo": "", "bar":"a"}'],
            ['{"foo": "", "bar":""}','{"foo": "1", "bar":""}','{"foo": "", "bar":"b"}','{"foo": "", "bar":""}'],
            ['{"foo": "", "bar":""}','{"foo": "", "bar":"c"}','{"foo": "2", "bar":""}','{"foo": "", "bar":""}'],
            ['{"foo": "", "bar":"d"}','{"foo": "", "bar":""}','{"foo": "", "bar":""}','{"foo": "3", "bar":""}'],
        ]
    )
    return df
class TestReadJsonColumn:
    @pytest.mark.parametrize(
        ('col', 'key', 'expected'),
        (
            ('Method', 'foo', pd.Series(data=['0',None,None,None],name='Method')),
            ('Method', 'bar', pd.Series(data=[None,None,None,'d'],name='Method')),
            ('Meta', 'foo', pd.Series(data=[None,'1',None,None],name='Meta')),
            ('Meta', 'bar', pd.Series(data=[None,None,'c',None],name='Meta')),
            ('CDR3fix', 'foo', pd.Series(data=[None,None,'2',None],name='CDR3fix')),
            ('CDR3fix', 'bar', pd.Series(data=[None,'b',None,None],name='CDR3fix')),
        )
    )
    def test_read_json_column(self, json_df, col, key, expected):
        result = mf.read_json_column(df=json_df, col=col, key=key)

        assert result.equals(expected)

    
    @pytest.mark.parametrize(
        'col', ('Foo', 'Bar')
    )
    def test_error_with_non_json_column(self, json_df, col):
        with pytest.raises(AssertionError):
            mf.read_json_column(df=json_df, col=col, key='foo')

    
#     def test_error_with_nonexistent_key(self, json_df):
#         with pytest.raises(ValueError):
#             mf.read_json_column(df=json_df, col='Method', key='baz')

def test_correlation():
    a = np.array([
        [1,1,2,2],
        [1,2,3,4],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    b = np.array([
        [1,1,2,2],
        [4,3,2,1],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    corr = np.array([
        [1, -0.89, 0, 0],
        [0.89,-1,0,0],
        [0,0,1,-1],
        [0,0,-1,1]
    ])

    np.testing.assert_array_equal(corr.round(2), (mf.correlation(a, b)).round(2))

def test_correlation_huge0():
    a = np.array([
        [1,1,2,2],
        [1,2,3,4],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    b = np.array([
        [1,1,2,2],
        [4,3,2,1],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    corr = np.array([
        [1, -0.89, 0, 0],
        [0.89,-1,0,0],
        [0,0,1,-1],
        [0,0,-1,1]
    ])

    abig = np.zeros(shape=(4098,4))
    bbig = np.zeros(shape=(5000,4))
    corrbig = np.zeros(shape=(4098,5000))
    corrbig[corrbig == 0] = np.nan

    abig[5:9] = a
    bbig[10:14] = b
    corrbig[5:9,10:14] = corr

    assert corrbig.shape == mf.correlation(abig, bbig).shape

    np.testing.assert_array_equal(corrbig.round(2), (mf.correlation(abig, bbig)).round(2))

def test_correlation_huge1():

    a = np.array([
        [1,1,2,2],
        [1,2,3,4],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    b = np.array([
        [1,1,2,2],
        [4,3,2,1],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    corr = np.array([
        [1, -0.89, 0, 0],
        [0.89,-1,0,0],
        [0,0,1,-1],
        [0,0,-1,1]
    ])

    abig = np.zeros(shape=(4098,4))
    bbig = np.zeros(shape=(5000,4))
    corrbig = np.zeros(shape=(4098,5000))
    corrbig[corrbig == 0] = np.nan

    abig[4094:4098] = a
    bbig[0:4] = b
    corrbig[4094:4098,0:4] = corr

    assert corrbig.shape == mf.correlation(abig, bbig).shape
    np.testing.assert_array_equal(corrbig.round(2), (mf.correlation(abig, bbig)).round(2))

def test_correlation_huge2():

    a = np.array([
        [1,1,2,2],
        [1,2,3,4],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    b = np.array([
        [1,1,2,2],
        [4,3,2,1],
        [1,np.nan, np.nan,1],
        [np.nan, 1,1,np.nan]
    ])

    corr = np.array([
        [1, -0.89, 0, 0],
        [0.89,-1,0,0],
        [0,0,1,-1],
        [0,0,-1,1]
    ])

    abig = np.zeros(shape=(4098,4))
    bbig = np.zeros(shape=(5000,4))
    corrbig = np.zeros(shape=(4098,5000))
    corrbig[corrbig == 0] = np.nan

    abig[0:4] = a
    bbig[4094:4098] = b
    corrbig[0:4, 4094:4098] = corr
    print(corrbig[0:10, 4090:4098])
    X = mf.correlation(abig, bbig)
    print(X[0:10, 4090:4098])

    assert corrbig.shape == X.shape
    np.testing.assert_array_equal(corrbig.round(2), X.round(2))

def test_covariance():
    df1 = pd.DataFrame({
    'a':[1,3,4,6,8],
    'b':[2,3,5,6,8],
    })

    df2 = pd.DataFrame({
        'c':[6,5,4,3,2],
        'd':[5,4,3,4,6]
    })

    x = np.array(df1.T)
    y = np.array(df2.T)

    np.testing.assert_array_equal(np.array([[-4.25,1.05],[-3.75,0.85]]), mf.covariance(x,y))

def test_slope():
    a = np.array([[15, 12, 8, 8, 7, 7, 7, 6, 5, 3],[15, 12, 8, 8, 7, 7, 7, 6, 5, 3]])
    b = np.array([[10, 25, 17, 11, 13, 17, 20, 13, 9, 15],[15, 12, 8, 8, 7, 7, 7, 6, 5, 3], [10, 25, 17, 11, 13, 17, 20, 13, 9, 15]])

    np.testing.assert_array_equal(np.array([[0.21, 1, 0.21], [0.21, 1, 0.21]]).round(2), mf.slope(a,b)[0].round(2))

def test_eucl_dist():
    assert mf.eucl_dist(np.array([[11,11,11]]), np.array([[1,1,1]])).round(2) == 2.89

def test_eucl():
    a = np.array([[15, 12, 8, 8, 7, 7, 7, 6, 5, 3],
                  [15, 12, 8, 8, 7, 7, 7, 6, 5, 3]])
    b = np.array([[10, 25, 17, 11, 13, 17, 20, 13, 9, 15],
                  [15, 12, 8, 8, 7, 7, 7, 6, 5, 3], 
                  [10, 25, 17, 11, 13, 17, 20, 13, 9, 15]])

    expected_output = np.array([[2.48,0.,2.48],
                                [2.48,0.,2.48]])

    np.testing.assert_array_equal(mf.eucl(a,b).round(2), expected_output)