import functions.structurefunctions as stf
import numpy as np
import os
import pandas as pd
import pymol

structures_dir = "test_data/tcrdb_raw/"
structures_dir1 = "test_data/pdb/"
pdblist = [x.split(".pdb")[0] for x in os.listdir(structures_dir)]
pymol.cmd.set("cif_keepinmemory")
pymol.cmd.set('fetch_path', structures_dir1, quiet=0)

def test_header_retrieve():
    output = {}
    output1 = {}
    for pdb in pdblist:
        print(pdb)
        info = stf.get_missing_residues_from_header(structures_dir, pdb)
        info1 = stf.pymol_missing_residues(pdb)
        output[pdb] = info
        output1[pdb] = info1
        if info is not None:
            assert info.equals(info1)
        else:
            assert info == info1
    
    df_1ao7 = pd.DataFrame(np.array( [["GLU", "A", "275"], ["ASP", "D", "122"], 
                                ["PRO", "D", "123"], ["ALA", "D", "124"], 
                                ["VAL", "D", "125"], ["TYR", "D", "126"], 
                                ["GLN", "D", "127"], ["LEU", "D", "128"], 
                                ["ARG", "D", "129"], ["ASP", "D", "130"], 
                                ["SER", "D", "131"], ["LYS", "D", "132"], 
                                ["SER", "D", "133"], ["SER", "D", "134"], 
                                ["ASP", "D", "135"], ["LYS", "D", "136"], 
                                ["SER", "D", "137"], ["VAL", "D", "138"], 
                                ["CYS", "D", "139"], ["LEU", "D", "140"], 
                                ["PHE", "D", "141"], ["THR", "D", "142"], 
                                ["ASP", "D", "143"], ["PHE", "D", "144"], 
                                ["ASP", "D", "145"], ["SER", "D", "146"], 
                                ["GLN", "D", "147"], ["THR", "D", "148"], 
                                ["ASN", "D", "149"], ["VAL", "D", "150"], 
                                ["SER", "D", "151"], ["GLN", "D", "152"], 
                                ["SER", "D", "153"], ["LYS", "D", "154"], 
                                ["ASP", "D", "155"], ["SER", "D", "156"], 
                                ["ASP", "D", "157"], ["VAL", "D", "158"], 
                                ["TYR", "D", "159"], ["ILE", "D", "160"], 
                                ["THR", "D", "161"], ["ASP", "D", "162"], 
                                ["LYS", "D", "163"], ["THR", "D", "164"], 
                                ["VAL", "D", "165"], ["LEU", "D", "166"], 
                                ["ASP", "D", "167"], ["MET", "D", "168"], 
                                ["ARG", "D", "169"], ["SER", "D", "170"], 
                                ["MET", "D", "171"], ["ASP", "D", "172"], 
                                ["PHE", "D", "173"], ["LYS", "D", "174"], 
                                ["SER", "D", "175"], ["ASN", "D", "176"], 
                                ["SER", "D", "177"], ["ALA", "D", "178"], 
                                ["VAL", "D", "179"], ["ALA", "D", "180"], 
                                ["TRP", "D", "181"], ["SER", "D", "182"], 
                                ["ASN", "D", "183"], ["LYS", "D", "184"], 
                                ["SER", "D", "185"], ["ASP", "D", "186"], 
                                ["PHE", "D", "187"], ["ALA", "D", "188"], 
                                ["CYS", "D", "189"], ["ALA", "D", "190"], 
                                ["ASN", "D", "191"], ["ALA", "D", "192"], 
                                ["PHE", "D", "193"], ["ASN", "D", "194"], 
                                ["ASN", "D", "195"], ["SER", "D", "196"], 
                                ["ILE", "D", "197"], ["ILE", "D", "198"], 
                                ["PRO", "D", "199"], ["GLU", "D", "200"], 
                                ["ASP", "D", "201"], ["THR", "D", "202"], 
                                ["PHE", "D", "203"], ["PHE", "D", "204"], 
                                ["PRO", "D", "205"], ["SER", "D", "206"], 
                                ["PRO", "D", "207"], ["GLU", "D", "208"], 
                                ["SER", "D", "209"], ["SER", "D", "210"], 
                                ["ASN", "E", "  1"], ["ALA", "E", "  2"], 
                                ["GLU", "E", "131"], ["PRO", "E", "132"], 
                                ["SER", "E", "133"], ["GLU", "E", "134"], 
                                ["ALA", "E", "135"], ["GLU", "E", "136"], 
                                ["ILE", "E", "137"], ["SER", "E", "138"], 
                                ["HIS", "E", "139"], ["THR", "E", "140"], 
                                ["GLN", "E", "141"], ["LYS", "E", "142"], 
                                ["ALA", "E", "143"], ["THR", "E", "144"], 
                                ["LYS", "E", "180"], ["GLU", "E", "181"], 
                                ["GLN", "E", "182"], ["PRO", "E", "183"], 
                                ["ALA", "E", "184"], ["LEU", "E", "185"], 
                                ["ASN", "E", "186"], ["ASP", "E", "187"], 
                                ["SER", "E", "188"], ["ARG", "E", "189"], 
                                ["SER", "E", "220"], ["GLU", "E", "221"], 
                                ["ASN", "E", "222"], ["ASP", "E", "223"], 
                                ["GLU", "E", "224"], ["TRP", "E", "225"], 
                                ["THR", "E", "226"], ["GLN", "E", "227"], 
                                ["ASP", "E", "228"], ["ARG", "E", "229"]]))
    df_1ao7.columns = ["RES", "C", "SSSEQI"]
    # print(output["1ao7"])
    # print(df)
    output["1ao7"].SSSEQI = output["1ao7"].SSSEQI.astype("int")
    df_1ao7.SSSEQI = df_1ao7.SSSEQI.astype("int")
    assert output["1ao7"].equals(df_1ao7)
    assert output["2apf"] == None
    assert output["2cdf"] == None

    df_2ij0 = pd.DataFrame(np.array([["SER", "A", "1"], ["THR", "A", "2"], ["ASN", "A", "3"]]))
    df_2ij0.columns = ["RES", "C", "SSSEQI"]
    output["2ij0"].SSSEQI = output["2ij0"].SSSEQI.astype("int")
    df_2ij0.SSSEQI = df_2ij0.SSSEQI.astype("int")
    assert output["2ij0"].equals(df_2ij0)

    df_3qh3 = pd.DataFrame(np.array([["ALA", "D", "98"], ["GLY", "D", "99"], 
                                    ["GLY", "D", "100"], ["ARG", "D", "101"]]))
    df_3qh3.columns = ["RES", "C", "SSSEQI"]
    output["3qh3"].SSSEQI = output["3qh3"].SSSEQI.astype("int")
    df_3qh3.SSSEQI = df_3qh3.SSSEQI.astype("int")
    assert output["3qh3"].equals(df_3qh3)

