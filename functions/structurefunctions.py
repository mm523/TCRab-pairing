import pandas as pd
from Bio import PDB
import pymol
from datetime import date

## the following functions try and parse the header information in many different ways, 
## to check what missing residues you have and decide how to deal with them. 

def GetTheSentences(infile):
    # get missing residues by parsing PDB as .txt file
    all_lines = []
    remark_found = 0
    with open(infile) as fp:
        for line in fp:
            all_lines.append(line)
            if "MISSING RESIDUES" in line:
                remark_num = line[0:10]
                remark_found = 1
                # print(remark_num)
    if remark_found == 1:
        lines = [line.replace(remark_num, "") for line in all_lines if remark_num in line]
    else:
        print("No missing residues")
        lines = []

    return lines

def get_missing_residues_from_header(dir, pdb):
    missing_header = GetTheSentences(dir + pdb + ".pdb")
    # print(missing_header)
    if len(missing_header) > 0:
        # print(missing_header)
        missing_header = [x.split("\n") for x in missing_header]
        # print(missing_header)
        missing_header = [[y.strip() for y in x] for x in missing_header][5:]
        # print(missing_header)
        missing_header = [y for x in missing_header for y in x if y !=""]
        # print(missing_header)
        m1 = [x.split(" ") for x in missing_header]
        # print(m1)
        assert missing_header[0] == "M RES C SSSEQI"
        m2 = [[y for y in x if y!=""]  for x in m1]
        df = pd.DataFrame(m2[1:])
        df.columns = ["RES", "C", "SSSEQI"]
        return df

def bioPDB_missing_residues(dir, pdb):
    h = PDB.parse_pdb_header(dir + pdb + ".pdb")
    missing = h["missing_residues"]
    df = pd.DataFrame([(x["res_name"], x["chain"], x["ssseq"]) for x in missing])
    df.columns = ["RES", "C", "SSSEQI"]
    return df

def pymol_missing_residues(pdb):
    basestr = "_pdbx_unobs_or_zero_occ_residues."
    attributes = ["auth_comp_id", "auth_asym_id", "auth_seq_id"]
    x = {}
    try:
        for i, name in enumerate(["RES", "C", "SSSEQI"]):
            x[name] = pymol.querying.cif_get_array(pdb, basestr + attributes[i])
        df = pd.DataFrame(x)
        df.columns = ["RES", "C", "SSSEQI"]
        return df
    except:
        print("No missing residue information from pymol")

## now some other structure-based functions

def clean_and_align_structures(pdblist, structures_dir, datalist, outf, coldict=None):
    """
    Takes the pdbs in pdblist, keeps only one set of chains as per datalist (which is a .csv, like you can extract from stcrdab),
    then colours them according to coldict and aligns them. If not provided,
    default colour is grey.
    Finally, it saves the pymol session in outf to be opened in pymol.
    """
    for pdb in pdblist:
        print(pdb)
        pymol.cmd.load(structures_dir + pdb + ".pdb")
        chains = [ch for ch in pymol.cmd.get_chains(pdb)]
        achain, bchain = datalist.loc[datalist.pdb == pdb, ["Achain", "Bchain"]].values[0]
        print(achain, bchain)
        chains_to_delete = [x for x in chains if x not in [achain, bchain]]
        if len(chains_to_delete)>0:
            mystring = "+".join(chains_to_delete)
            name = pdb + "_to_delete"
            mycommand = pdb + " & chain " + mystring
            print(mycommand)
            pymol.cmd.select(name, mycommand)
            pymol.cmd.remove(name)
        chains1 = [ch for ch in pymol.cmd.get_chains(pdb)]
        assert sorted(chains1) == sorted([achain, bchain])
        if coldict!= None and pdb in coldict.keys():
            pymol.cmd.color(coldict[pdb], pdb)
        else:
            pymol.cmd.color("grey", pdb)
        if pdb != pdblist[0]:
            pymol.cmd.align(pdb, pdblist[0])

    pymol.cmd.save(outf + str(date.today()) + "_pymol-"+"-".join(pdblist)+".pse")

def clean_and_align_structures_withep(pdblist, structures_dir, datalist, outf, name = None, coldict=None):
    """
    Takes the pdbs in pdblist, keeps only one set of chains as per datalist (which is a .csv, like you can extract from stcrdab),
    then colours them according to coldict and aligns them. If not provided,
    default colour is grey.
    Finally, it saves the pymol session in outf to be opened in pymol.
    Different from above because it keep pMHC.
    """
    for pdb in pdblist:
        print(pdb)
        pymol.cmd.load(structures_dir + pdb + ".pdb")
        chains = [ch for ch in pymol.cmd.get_chains(pdb)]
        achain, bchain, echain, mhc_chain1 = datalist.loc[datalist.pdb == pdb, ["Achain", "Bchain", "antigen_chain", "mhc_chain1"]].values[0]
        print(achain, bchain, echain, mhc_chain1)
        chains_to_delete = [x for x in chains if x not in [achain, bchain, echain, mhc_chain1]]
        if len(chains_to_delete)>0:
            mystring = "+".join(chains_to_delete)
            name = pdb + "_to_delete"
            mycommand = pdb + " & chain " + mystring
            print(mycommand)
            pymol.cmd.select(name, mycommand)
            pymol.cmd.remove(name)
        chains1 = [ch for ch in pymol.cmd.get_chains(pdb)]
        assert sorted(chains1) == sorted([achain, bchain, echain, mhc_chain1])
        if coldict!= None and pdb in coldict.keys():
            pymol.cmd.color(coldict[pdb], pdb)
        else:
            pymol.cmd.color("grey", pdb)
        if pdb != pdblist[0]:
            pymol.cmd.align(pdb, pdblist[0])
    
    if name != None:
        name = "-".join(pdblist)

    pymol.cmd.save(outf + str(date.today()) + "_pymol-"+name+".pse")