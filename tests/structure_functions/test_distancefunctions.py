from Bio import PDB
import functions.distancefunctions as df
import pymol
import numpy as np
import os
import pandas as pd

structures_dir = "test_data/tcrdb_raw/"
structures_dir1 = "test_data/pdb/"
structures_dir2 = "test_data/tcrdb_imgt/"
datalist = pd.read_csv("data/tcrdb_export_cleaned.csv").drop_duplicates(subset = ["pdb"], keep="first")
p = PDB.PDBParser()
pymol.cmd.set('fetch_path', structures_dir1, quiet=0)

pdblist = [x.split(".pdb")[0] for x in os.listdir(structures_dir)]
print(pdblist)

def test_BioPDBdistancefunction():
    for pdb in pdblist:
        print(pdb)
        s = pdb + ".pdb"
        pymol.cmd.fetch(pdb)
        struc = p.get_structure(pdb, structures_dir + s)
        struc_imgt = p.get_structure(pdb, structures_dir2 + s)
      
        chain_pairs = list(zip(datalist.loc[datalist.pdb == pdb]["Achain"], datalist.loc[datalist.pdb == pdb]["Bchain"]))

        for pair in chain_pairs:
            distances = {}
            print(pair)
            achain, bchain = pair[0], pair[1]
            d1 = df.CA_pairwise_distances(struc[0][achain],struc[0][bchain])
            d1pymol = df.pairwise_dist_pymol(pdb + " and chain " + achain, pdb + " and chain " + bchain, "3000")
            d1imgt = df.CA_pairwise_distances(struc_imgt[0][achain],struc[0][bchain])
            distances["raw"] = d1
            distances["pymol"] = d1pymol
            distances["imgt"] = d1imgt

            d_raw = pd.DataFrame.from_dict(distances["raw"], orient="index").reset_index(drop = "True")
            d_pymol = pd.DataFrame.from_dict(distances["pymol"], orient="index").reset_index(drop = "True")
            d_imgt = pd.DataFrame.from_dict(distances["imgt"], orient="index").reset_index(drop = "True")

            d = pd.concat([d_raw, d_pymol, d_imgt], axis = 1)
            d.columns = ["raw", "pymol", "imgt"]

            print(d)

            assert all(abs(i) < 0.15 for i in list(d.raw - d.pymol)), "max: " + str(max(list(d.raw - d.pymol))) + " and min: " + str(min(list(d.raw - d.pymol)))
            assert all(abs(i) < 0.15 for i in list(d.imgt - d.pymol)), "max: " + str(max(list(d.imgt - d.pymol))) + " and min: " + str(min(list(d.imgt - d.pymol)))
            assert all(abs(i) < 0.15 for i in list(d.raw - d.imgt)), "max: " + str(max(list(d.raw - d.imgt))) + " and min: " + str(min(list(d.raw - d.imgt)))
       