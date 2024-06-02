import os, pickle, json
from collections import Counter
from Levenshtein import distance as levd
import numpy as np
from sklearn.metrics import mutual_info_score, adjusted_mutual_info_score
from sklearn.metrics.pairwise import euclidean_distances
import pandas as pd
from itertools import product

def inplaceremove(x, c):
    for i in range(len(x)):
        if c in x:
            x.remove(c)
    return x

def save_file(f, n, path):
    if not os.path.exists(path):
        os.makedirs(path)
    with open(path + n, 'wb') as handle:
        pickle.dump(f, handle, protocol=-1)

def load_file(filename, path):
    with open(path + filename, "rb") as input_file:
        e = pickle.load(input_file)
    return e

def find_and_remove_multispecific_cdr3(combined, cdr_col = "CDR3", epitope_col="Epitope"):
    """
    This function is used to remove cdr3s from the set which are annotated to recognise
    more than one epitope.
    """
    ### Expects column with cdr3 to be called "CDR3"
    ### Expects columns with peptide to be called "Epitope"

    grouped_cdr3 = combined.groupby('CDR3')
    num_epitopes = grouped_cdr3.count()['Epitope']
    promiscuous = num_epitopes[num_epitopes > 1].index.tolist()
    true_promiscuous = []

    for cdr3 in promiscuous:
        epitopes = grouped_cdr3.get_group(cdr3)['Epitope']
        n = len(epitopes)
        has_centroid = False

        for epitope in epitopes:
            # Check for centroid
            num_edges = 0
            for other in epitopes:
                if other == epitope: continue
                dist = levd(epitope,other)
                if dist == 0:
                    raise RuntimeError('Levenshtein distance of 0 detected.')
                elif dist > 1: break
                else: num_edges += 1
            
            # If centroid, this is a pattern that is allowed
            # No need to continue search, so break
            if num_edges == n-1:
                has_centroid = True
                break
        
        if not has_centroid: true_promiscuous.append(cdr3)

    combined = combined[~combined['CDR3'].isin(true_promiscuous)].reset_index(drop=True)
    return combined

def mutual_information(X,Y):
    assert X.shape[0] == Y.shape[0], "Number of observations do not correspond"
    mutual_information = np.empty((X.shape[1], Y.shape[1]))
    for Xcol in range(X.shape[1]):
        for Ycol in range(Y.shape[1]):
            MI = mutual_info_score(X[:,Xcol], Y[:,Ycol])
            mutual_information[Xcol,Ycol] = MI
    return(mutual_information)

def adjusted_mutual_information(X,Y):
    assert X.shape[0] == Y.shape[0], "Number of observations do not correspond"
    mutual_information = np.empty((X.shape[1], Y.shape[1]))
    for Xcol in range(X.shape[1]):
        for Ycol in range(Y.shape[1]):
            MI = adjusted_mutual_info_score(X[:,Xcol], Y[:,Ycol])
            mutual_information[Xcol,Ycol] = MI
    return(mutual_information)

def calculate_entropy_of_residuals(alphas, betas): # not tested
    alpha_entropy = []
    beta_entropy = []

    for x in range(alphas.shape[1]):
        counts = pd.DataFrame.from_dict(Counter(alphas[x]), orient = "index")/alphas.shape[0]
        a = -sum(counts.values*np.log(counts.values))
        alpha_entropy.append(a)
    for x in range(betas.shape[1]):
        counts = pd.DataFrame.from_dict(Counter(betas[x]), orient = "index")/betas.shape[0]
        b = -sum(counts.values*np.log(counts.values))
        beta_entropy.append(b)
    
    return np.array(alpha_entropy), np.array(beta_entropy)

def prepare_data(epdf, col1, col2 = None, max_len1 = None, max_len2 = None, type = "cdr3"):
    '''
    This function takes a dataframe and 2 columns. The columns will be identifiers for either cdr3 or for TCR sequences
    which have been IMGT-renumbered. It then works by adding padding in the middle of the cdr3 so that the sequences 
    are all the same length.
    '''

    assert type in ["cdr3", "tcr"], "type must be one of [\"cdr3\", \"tcr\"]"

    epdf["len_" + col1] = [len(x) for x in epdf[col1]]    
    max_col1 = max(epdf["len_" + col1])

    if pd.notna(max_len1):
        print("set max length for " + col1 + ": cleaning up the dataset")
        epdf = epdf.loc[epdf["len_" + col1] <= max_len1]
        assert epdf.shape[0] >= 1, "max_len1 too stringent, less than 2 records left"
    else:
        max_len1 = max_col1
    
    assert max(epdf["len_" + col1])<=max_len1, "df max length cleaning not successful in " + col1
    
    if type == "cdr3":
        epdf[col1 + "_padded"] = [", ".join(list(x.split("-")[0]) + ["-"]*(max_len1 - len(x.replace("-", ""))) + list(x.split("-")[-1])) if len(x.split("-")) > 1 else ", ".join(list(x[0:round(len(x)/2)]) + ["-"]*(max_len1 - len(x)) + list(x[round(len(x)/2):])) for x in epdf[col1]]
    elif type == "tcr":
        # we assume TCRs can only have insertions that do not align them between positions 111 and 112
        epdf[col1 + "_padded"] = [", ".join(list(x[0:111]) + ["-"]*(max_len1 - len(x)) + list(x[111:len(x)])) for x in epdf[col1]]
    else:
        raise ValueError("type not recognised")
    
    # print(epdf[col1 + "_padded"])

    assert len(set([len(x.split(", ")) for x in epdf[col1 + "_padded"]])) == 1, col1 + " padding did not work. More than 1 lengths still present"
    assert list(set([len(x.split(", ")) for x in epdf[col1 + "_padded"]]))[0] == max_len1, col1 + " padding did not work. Length is not equal to max_len1"
    
    if type == "tcr":
        assert set([x.split(", ")[103] for x in epdf[col1 + "_padded"]]) == {"C"}, "position 104 on alpha is not all C"

    if not col2 is None:
        epdf["len_" + col2] = [len(x) for x in epdf[col2]]
        max_col2 = max(epdf["len_" + col2])

        if pd.notna(max_len2):
            print("set max length for " + col2 + ": cleaning up the dataset")
            epdf = epdf.loc[epdf["len_" + col2] <= max_len2]
            assert epdf.shape[0] >= 1, "max_len2 too stringent, less than 2 records left"
        else:
            max_len2 = max_col2
        
        assert max(epdf["len_" + col2])<=max_len2, "df max length cleaning not successful in " + col1

        if type == "cdr3":
            epdf[col2 + "_padded"] = [", ".join(list(x.split("-")[0]) + ["-"]*(max_len2 - len(x.replace("-", ""))) + list(x.split("-")[-1])) if len(x.split("-")) > 1 else ", ".join(list(x[0:round(len(x)/2)]) + ["-"]*(max_len2 - len(x)) + list(x[round(len(x)/2):])) for x in epdf[col2]]
        elif type == "tcr":
            # we assume TCRs can only have insertions that do not align them between positions 111 and 112
            epdf[col2 + "_padded"] = [", ".join(list(x[0:111]) + ["-"]*(max_len2 - len(x)) + list(x[111:len(x)])) for x in epdf[col2]]
        else:
            raise ValueError("type not recognised")
        
        assert len(set([len(x.split(", ")) for x in epdf[col2 + "_padded"]])) == 1, col2 + " padding did not work. More than 1 lengths still present"
        assert list(set([len(x.split(", ")) for x in epdf[col2 + "_padded"]]))[0] == max_len2, col2 + " padding did not work. Length is not equal to max_len2"

        if type == "tcr":
            assert set([x.split(", ")[103] for x in epdf[col2 + "_padded"]]) == {"C"}, "position 104 on beta is not all C"


    return epdf

def read_json_column(df: pd.DataFrame, col: str, key: str) -> pd.Series:  # function and associated tests from YN
    '''
    Given a copy of the a vdjdb dataframe slice, for each row, parse the json
    object from the specified column, and return the value assigned to the
    specified key.
    '''

    assert col in ('Method', 'Meta', 'CDR3fix')

    return df[col]\
        .transform(lambda x: json.loads(x)[key])\
        .transform(lambda x: None if x == "" else x)

def correlation(x1, y1):
    end = max(x1.shape[0], y1.shape[0])
    step = min(end, 1000)
    block1 = [x for x in range(0, end + step, step)]
    block2 = [x for x in range(step, end + step, step)]
    blocks = [(x,y) for x,y in zip(block1, block2)]

    corr1 = np.zeros(shape = (x1.shape[0], y1.shape[0]))
    corr1[corr1 == 0] = 5

    for xblock, yblock in product(blocks, blocks):
        xstart = max(0, xblock[0]-(step+1))
        ystart = max(0, yblock[0]-(step+1))
        xend = min(x1.shape[0],xblock[1])
        yend = min(y1.shape[0],yblock[1])
        x = x1[xstart:xend]
        y = y1[ystart:yend]
        x[np.isnan(x)] = 0
        y[np.isnan(y)] = 0

        x_m = x - np.mean(x,axis=1)[:,None]
        y_m = y - np.mean(y,axis=1)[:,None]

        X = np.sum(x_m**2,axis=1)
        Y = np.sum(y_m**2,axis=1)

        corr = np.dot(x_m, y_m.T)/np.sqrt(np.dot(X[:,None],Y[None]))
        assert corr.shape == (xend - xstart, yend - ystart)

        corr1[np.ix_(range(xstart, xend), range(ystart, yend))] = corr
    assert np.nanmax(abs(corr1.round(10))) <= 1.
    assert corr1.shape == (x1.shape[0], y1.shape[0])
    return(np.array(corr1))

def covariance(x1, y1):
    end = max(x1.shape[0], y1.shape[0])
    step = min(end, 1000)
    block1 = [x for x in range(0, end + step, step)]
    block2 = [x for x in range(step, end + step, step)]
    blocks = [(x,y) for x,y in zip(block1, block2)]

    cov1 = np.zeros(shape = (x1.shape[0], y1.shape[0]))
    cov1[cov1 == 0] = 5

    for xblock, yblock in product(blocks, blocks):
        xstart = max(0, xblock[0]-(step+1))
        ystart = max(0, yblock[0]-(step+1))
        xend = min(x1.shape[0],xblock[1])
        yend = min(y1.shape[0],yblock[1])
        x = x1[xstart:xend]
        y = y1[ystart:yend]
        x[np.isnan(x)] = 0
        y[np.isnan(y)] = 0

        x_m = x - np.mean(x,axis=1)[:,None]
        y_m = y - np.mean(y,axis=1)[:,None]

        cov = np.dot(x_m, y_m.T)/(x_m.shape[1] - 1)
        assert cov.shape == (xend - xstart, yend - ystart)

        cov1[np.ix_(range(xstart, xend), range(ystart, yend))] = cov

    assert cov1.shape == (x1.shape[0], y1.shape[0])
    return(np.array(cov1))

def slope(x1, y1):
    end = max(x1.shape[0], y1.shape[0])
    step = min(end, 1000)
    block1 = [x for x in range(0, end + step, step)]
    block2 = [x for x in range(step, end + step, step)]
    blocks = [(x,y) for x,y in zip(block1, block2)]

    slopes1 = np.zeros(shape = (x1.shape[0], y1.shape[0]))
    slopes1[slopes1 == 0] = 5
    intercepts1 = np.zeros(shape = (x1.shape[0], y1.shape[0]))

    for xblock, yblock in product(blocks, blocks):
        xstart = max(0, xblock[0]-(step+1))
        ystart = max(0, yblock[0]-(step+1))
        xend = min(x1.shape[0],xblock[1])
        yend = min(y1.shape[0],yblock[1])
        x = x1[xstart:xend]
        y = y1[ystart:yend]
        x[np.isnan(x)] = 0
        y[np.isnan(y)] = 0
        
        slopes = []
        intercepts = []
        for i in range(x.shape[0]):
            s, interc = np.polyfit(x[i], y.T, 1)[0:2]
            slopes.append(s)
            intercepts.append(interc)

        slopes = np.array(slopes)
        intercepts = np.array(intercepts)
        assert slopes.shape == (xend - xstart, yend - ystart)

        slopes1[np.ix_(range(xstart, xend), range(ystart, yend))] = slopes
        intercepts1[np.ix_(range(xstart, xend), range(ystart, yend))] = intercepts

    assert slopes1.shape == (x1.shape[0], y1.shape[0])
    return(np.array(slopes1), np.array(intercepts1))

def eucl_dist(x,y):
    return(np.linalg.norm(x-y)/np.mean([np.mean(x),np.mean(y)]))

def eucl(x1, y1):
    eucl1 = euclidean_distances(x1,y1)
    means = np.subtract.outer(np.mean(x1, axis=1),-np.mean(y1, axis=1))/2
    eucl1 = eucl1/means

    return(np.array(eucl1))