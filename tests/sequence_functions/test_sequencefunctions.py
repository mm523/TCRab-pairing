import functions.sequencefunctions as sf
import numpy as np
import pandas as pd

a1ao7 = "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESS"
a1ao7_r = "-KEVEQNSGPLSVPEGAIASLNCTYSDRG------SQSFFWYRQYSGKSPELIMSIYS----NGDKED-----GRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDS--WGKLQFGAGTQVVVTP"
b1ao7 = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD"
b1ao7_r = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNH-------EYMSWYRQDPGMGLRLIHYSVG----AGITDQGEVP-NGYNVSRS-TTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"

def test_get_v_region_seq():
    a = sf.get_vregion_seq_with_gaps(a1ao7)
    a1 = sf.get_vregion_seq(a1ao7)

    assert a == a1ao7_r
    assert a1 == a1ao7_r.replace("-", "")

    b = sf.get_vregion_seq_with_gaps(b1ao7)
    b1 = sf.get_vregion_seq(b1ao7)

    assert b == b1ao7_r
    assert b1 == b1ao7_r.replace("-", "")

def test_get_cdrs_with_gaps():
    assert sf.get_cdr1_seq_with_gaps(a1ao7) == 'DRG------SQS'
    assert sf.get_cdr2_seq_with_gaps(a1ao7) == 'IYS----NGD'
    assert sf.get_cdr3_seq_with_gaps(a1ao7) == 'CAVTTDS--WGKLQF'

s1 = "AAAA"
s2 = "BAAA"
s3 = "BCAA"
s4 = "AAAA"

mylist = [s1, s2, s3, s4]
expected_array = np.array([[0, 1, 2, 0],
                            [1, 0, 1, 1],
                            [2, 1, 0, 2], 
                            [0, 1, 2, 0]])

def test_levd_distance_array():
    m = sf.array_levd(mylist)
    assert np.array_equal(m, expected_array)

myseqs = {"uno":{"my_aa_imgt":"AAAAABC"}, 
          "due":{"my_aa_imgt":"AAAAABCDEF"}, 
          "tre":{"my_aa_imgt":"AAAAA"},
          "quattro":{"my_aa_imgt":"BBB"}, 
          "cinque":{"my_aa_imgt":"BBBAA"}}

df = pd.DataFrame(myseqs).T.reset_index()
df = df.rename(columns = {"index":"pdb"})

myseqs1 = {"uno":{"my_aa_imgt":"AAAAABC", "my_aa_imgt_shorter":"AAAAA"}, 
          "due":{"my_aa_imgt":"AAAAABCDEF", "my_aa_imgt_shorter":"AAAAA"}, 
          "tre":{"my_aa_imgt":"AAAAA", "my_aa_imgt_shorter":"AAAAA"},
          "quattro":{"my_aa_imgt":"BBB", "my_aa_imgt_shorter":"BBB"}, 
          "cinque":{"my_aa_imgt":"BBBAA", "my_aa_imgt_shorter":"BBB"}}
df1 = pd.DataFrame(myseqs1).T.reset_index()
df1 = df1.rename(columns = {"index":"pdb"})

def test_find_more_duplicates():
    a, b = sf.find_more_duplicates(myseqs.keys(), myseqs,df, "my")
    assert b.equals(df1)

myseqs2 = {"uno":{"alpha_aa_imgt":"AAAAA", "beta_aa_imgt":"BBBBBBBB"},
          "unoB":{"alpha_aa_imgt":"AAAAA", "beta_aa_imgt":"BBBBBBBB"}, 
          "due":{"alpha_aa_imgt":"AAAAA", "beta_aa_imgt":"BABABABA"}, 
          "tre":{"alpha_aa_imgt":"AAAAA", "beta_aa_imgt":"BABACABA"},
          "quattro":{"alpha_aa_imgt":"BBB", "beta_aa_imgt":"AABABA"}, 
          "cinque":{"alpha_aa_imgt":"BBB", "beta_aa_imgt":"BABABABA"}}
df2 = pd.DataFrame(myseqs2).T
df2["pdb"] = df2.index.values
# df2 = df2.rename(columns = {"index":"pdb"})
mylist1=[myseqs2["uno"]["beta_aa_imgt"], myseqs2["unoB"]["beta_aa_imgt"], 
         myseqs2["due"]["beta_aa_imgt"], myseqs2["tre"]["beta_aa_imgt"]]
mylist2=[myseqs2["quattro"]["beta_aa_imgt"], myseqs2["cinque"]["beta_aa_imgt"]]
df3 = pd.DataFrame(sf.array_levd(mylist1))
df3.index = ["uno", "unoB", "due", "tre"]
df3.columns = ["uno", "unoB", "due", "tre"]
df4 = pd.DataFrame(sf.array_levd(mylist2))
df4.index = ["quattro", "cinque"]
df4.columns = ["quattro", "cinque"]
myseqs2_res = {"uno,unoB,due,tre":df3,
                "quattro,cinque":df4}

## I think the test below does not work because I do not have epitope in my dics above
# def test_find_groups_of_shared_chain():
#     groups = df2.reset_index().groupby(by = ["alpha_aa_imgt"])["pdb"].apply(list)
#     a, b = sf.find_groups_of_shared_chain(groups, "beta", df2)
#     # print(a)
#     assert a.keys() == myseqs2_res.keys()
#     assert a[list(a.keys())[0]].equals(myseqs2_res[list(a.keys())[0]])
            