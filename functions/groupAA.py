import numpy as np

def translate_aa_into_groups(myarray):

    mapping = {
    # from https://www.researchgate.net/publication/308675542_Soy-based_adhesives_for_wood-bonding_a_review/figures?lo=1
    'G':'a', # simple, hydrophobic
    'P':'b', #cyclic, hydrophobic
    'F':'c', #aromatic, hydrophobic
    'A':'d', #aliphatic, hydrophobic
    'V':'d', #aliphatic, hydrophobic
    'I':'d', #aliphatic, hydrophobic
    'L':'d', #aliphatic, hydrophobic
    'S':'e', #hydroxylic, hydrophillic
    'T':'e', #hydroxylic, hydrophillic
    'D':'f', #acidic, hydrophillic
    'E':'f', #acidic, hydrophillic
    'C':'g', #sulphur, hydrophillic
    'N':'h', #amidic, hydrophillic
    'Q':'h', #amidic, hydrophillic
    'R':'i', #basic, hydrophillic
    'H':'i', #basic, hydrophillic
    'K':'j', #basic, amphipathic
    'W':'k', #aromatic, amphipathic
    'Y':'k', #aromatic, amphipathic
    'M':'l', #sulphur, amphipathic
    '-':'-'
}

    k = np.array(list(mapping.keys()))
    v = np.array(list(mapping.values()))

    out = np.zeros_like(myarray)
    for key,val in zip(k,v):
        out[myarray==key] = val
    
    return(out)