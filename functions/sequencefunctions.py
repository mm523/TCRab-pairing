from Bio.SeqUtils import seq1
from Levenshtein import distance as levd
import numpy as np
from anarci import anarci
from itertools import permutations
import pandas as pd

def get_chain_seq(chain):
    ## get sequence from PDB chain object (structure[model][chain])
    r = [seq1(x.resname) for x in chain.get_residues()]
    return("".join(r).strip("X"))

def array_levd(mylist):
    # from list of sequences, obtain array of pairwise edit distances
    Matrix = np.zeros((len(mylist),len(mylist)),dtype=int)
    for i in range(0,len(mylist)):
        for j in range(0,len(mylist)):
            Matrix[i,j] = levd(mylist[i],mylist[j])
    return(Matrix)

def get_vregion_seq(s):
    # get sequence of only the V region
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering]).strip().replace("-", "")
    else:
        seq = None
    return(seq)

def get_vregion_seq_with_gaps(s):
    # get sequence of only the V region, with IMGT gaps
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering]).strip()
        # print(seq)
    else:
        seq = None
    return(seq)

def get_cdr3_seq(s):
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering if x[0][0] > 103 and x[0][0] <119]).strip().replace("-", "")
    else:
        seq = None
    return seq

def get_cdr3_seq_with_gaps(s):
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering if x[0][0] > 103 and x[0][0] <119]).strip()
    else:
        seq = None
    return seq

def get_cdr1_seq_with_gaps(s):
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering if x[0][0] > 26 and x[0][0] < 39]).strip()
    else:
        seq = None
    return seq

def get_cdr2_seq_with_gaps(s):
    imgt = anarci([("myseq", s)], "imgt")
    if imgt[0][0] != None:
        renumbering = imgt[0][0][0][0]
        seq = "".join([x[1] for x in renumbering if x[0][0] > 55 and x[0][0] < 66]).strip()
    else:
        seq = None
    return seq

def find_more_duplicates(all_pdbs, seqs, pdbs_shorter_seqs, chain):
    # this is to try and figure out if I have duplicates in my PDBs.
    # Because some chains have extra residues at beginning or end,
    # check if seq1 is inside seq2, if so replace seq2 with seq1
    pairlist = all_pdbs
    new_seqs = {}
    all_similarities = {}

    i = 0
    while len(pairlist) > 0:
        # iterating until there is nothing else that is somewhat similar
        i+=1
        print(i)
        substrings = {}
        
        for pdb1,pdb2 in permutations(all_pdbs,2):
            if i == 1:
                s1 = seqs[pdb1][chain + "_aa_imgt"].replace("-", "")
                s2 = seqs[pdb2][chain + "_aa_imgt"].replace("-", "")
            else:
                s1 = new_seqs[pdb1]
                s2 = new_seqs[pdb2]        
            
            assert s1!=np.nan, "sequence 1 is nan"
            assert s1 != "", "sequence 1 is empty"
            assert s2!=np.nan

            if s1 in s2:
                if s1 != s2:
                    # here s1 is inside s2, so I will update pdb2 sequence to have s1
                    assert len(s1) < len(s2)
                    substrings[(pdb1, pdb2)] = (s1, s2)
                    if i == 1:
                        all_similarities[(pdb1, pdb2)] = (s1, s2)
                    new_seqs[pdb1] = s1
                    new_seqs[pdb2] = s1
                    
                elif s1 == s2:
                    # these two sequences are identical, not interesting
                    if i == 1:
                        all_similarities[(pdb1, pdb2)] = (s1, s2)
                    new_seqs[pdb1] = s1
                    new_seqs[pdb2] = s2
                    
                else:
                    ValueError("This is a scenario I haven't considered")
            else:
                # these two sequences are not similar, keep both
                # print(pdb1, pdb2, "not similar")
                new_seqs[pdb1] = s1
                new_seqs[pdb2] = s2
                
        pairlist = substrings.keys()
    
    new_seqs_df = pd.DataFrame.from_dict(new_seqs, orient="index")
    new_seqs_df.columns = [chain + "_aa_imgt_shorter"]
    pdbs_shorter_seqs = pd.merge(pdbs_shorter_seqs, new_seqs_df,left_on = "pdb",  right_index=True)

    return all_similarities, pdbs_shorter_seqs

def find_groups_of_shared_chain(groups, paired_chain, sequences):
    # groups has a group of structure which share one chain (df.groupby(chain))
    # now, I want to see how different the paired chain is
    beta_results_levd = {}
    beta_results_levd_nodupl = {}
    # to_remove = []

    for group in groups:
        print("Group: ", group)
        alphas = []
        complexes = []
        for pdb in group:
            a = sequences.loc[sequences.pdb == pdb][paired_chain + "_aa_imgt"].values[0]
            p = sequences.loc[sequences.pdb == pdb]["epitope"].values[0]
            complexes.append(p)
            alphas.append(a)
        alpha_levd = pd.DataFrame(array_levd(alphas))
        alpha_levd.columns = group
        alpha_levd.index = group
        alpha_levd["epitope"] = complexes
        beta_results_levd[",".join(group)] = alpha_levd
        alpha_levd_nodupl = alpha_levd.drop_duplicates(keep="first")
        beta_results_levd_nodupl[", ".join(group)] = alpha_levd_nodupl
        # to_remove.append(list(set(alpha_levd_nodupl.columns.values) - set(alpha_levd_nodupl.index.values)))

        ## to_remove (as expected) returns the same structures as drop_duplicates on sequences

    return(beta_results_levd, beta_results_levd_nodupl)#, to_remove)
