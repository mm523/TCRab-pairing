# Combining GA and IPA to infer paralog pairing

The scripts in this folder were taken as-is from https://github.com/carlosgandarilla/GA-IPA, but copied here for ease of use in the pipeline. Below the README from original repo.

In order to predict interaction partners among paralogs from two protein families, just from their sequences, we will use successively two distinct methods. The first one, Graph Alignment ([GA](https://doi.org/10.48550/arXiv.0905.1893)) algorithm, exploits phylogenetic relationships between interaction partners, via sequence similarity (i.e., it aims at identifying interologs). The second, Iterative Pairing Algorithm (IPA), relies on the inter-family coevolutionary signal as detectable by Direct Coupling Analysis ([DCA-IPA](https://doi.org/10.1073/pnas.1606762113)) and Mutual Information ([MI-IPA](https://doi.org/10.1371/journal.pcbi.1006401)).

Both signals can be used to computationally predict which proteins are specific interaction partners among the paralogs of two interacting protein families, starting just from their sequences. We show that combining them improves the performance of interaction partner inference, especially when the average number of potential partners is large and when the total data set size is modest.

Paper: "Combining phylogeny and coevolution improves the inference of iteration partners among paralogous proteins"(link to the paper) (Carlos A. Gandarilla-Perez, Sergio Pinilla, Anne-Florence Bitbol and Martin Weigt xxxxxxx , 2022)
<img width="1742" alt="Figure1" src="https://user-images.githubusercontent.com/43339953/191007155-de52cfd7-fdd1-4332-b891-f0743a5d755c.png">
We provide here the code to reproduce the key results and figures of the paper.

## Usage:

To run the code, you first need to install [julia](https://julialang.org/) and then run the following commands to test the demo:
```
include("./GA-IPA_main.jl")
n_replicates = 100
GA_IPA_DCA_robust(n_replicates) or GA_IPA_MI_robust(n_replicates)
```

where n_replicates is the number of realizations of GA experiment, each time taking different random matchings as starting points.

Here it is applied to a dataset of 5052 sequences of cognate histidine kinases and response regulators from the [P2CS](http://www.p2cs.org/) database. To run the code on the other datasets, unzip data from CoAlignments.zip, open the file GA_IPA_main.jl with your favorite editor, and replace protfile with the new protein file name. Make sure to change length of the two partners as well.

## Output - GA:

### Filename: 
* **"./HK-RR_GA_Orthology_TP_protB_M=$M._N_sweep=$n_sweep._N_replicates=$n_replicates.text"** *or* 
* **"./HK-RR_GA_", kNN, "-NN_TP_protB_M=$M._N_sweep=$n_sweep._N_replicates=$n_replicates.text"**

This file contains the results from the graph alignment. For M sequences where annealing is repeated n times, this will have shape (M+4, s).
Rows: 

(row:1) number of sequences in replicate (repeated in each column)

(row:2) mean number of pairs per species (repeated in each column)

Then each column is a repeat of n_repeated graph annealing runs. Each row within the column:

(row:3) true positives for that repeat

(row:4) energy of the final match for that repeat

(row:5 to M+4) final match for each sequence

### Filename:

* **"./HK-RR_Orthology_robust_pairs_GA_without_singletons_M=$M._N_replicates=$n_replicates.text"** *or* 
* **"./HK-RR_", kNN, "-NN_robust_pairs_GA_without_singletons_M=$M._N_replicates=$n_replicates.text"**

This file contains the pairs that are found to have robust pairing in the GA. 


(col:1) ID species

(col:2) ID sequence protA

(col:3) ID always repeated partner B

### Filename:

* **"./HK-RR_Orthology-NN_no_robust_pairs_GA_without_singletons_M=$M._N_replicates=$n_replicates.text"** *or* 
* **"./HK-RR_", kNN, "-NN_no_robust_pairs_GA_without_singletons_M=$M._N_replicates=$n_replicates.text"**

This file contains the pairs for which there is no robust partner found in GA. I *think* this is used as the file containing the right answers to feed to the DCA pairing

## Output - IPA:

**NOTE:** once the training set is aligned in the GA, it is treated as a golden set and therefore not aligned again. The only aligned sequences are the ones that did not get a robust aligned from the GA.

### Filename: 

* **"./HK-RR_DCA-IPA_TP_deltaE_data_GA_Orthology_robust_pairs_M=$M._N_replicates=$n_replicates._Ninc$Nincrement.txt"** *or*
* **"./HK-RR_DCA-IPA_deltaE_TP_data_GA_", kNN, "-NN_robust_pairs_M=$M._N_replicates=$n_replicates._Ninc$Nincrement.txt"** 

Output matrix: Each row corresponds to an iteration of the DCA-IPA.

(col:1) number of sequences NSeqs in concatenated alignment used as training set

(col:2) effective number of sequences Meff in concatenated alignment used as training set

(col:3) number of TP pairs

(col:4) number of FP pairs

(col:5) number of TP pairs in concatenated alignment used as training set

(col:6) number of FP pairs in concatenated alignment used as training set

### Filename: 

* **"./HK-RR_DCA-IPA_deltaE_Resf_GA_Orthology_robust_pairs_M=$M._N_replicates=$n_replicates._Ninc", Nincrement, "_NStart", M_start, "_round", Nrounds, ".txt** *or*
* **"./HK-RR_DCA-IPA_deltaE_Resf_GA_", kNN, "-NN_robust_pairs_M=$M._N_replicates=$n_replicates._Ninc", Nincrement, "_NStart", M_start, "_round", Nrounds, ".txt"**

This file includes the results for the last round for each predicted pair.

(col:1) species

(col:2) A index in initial alignment

(col:3) B index in initial alignment

(col:4) score by absolute energy of pairing

(col:5) gap wrt A

(col:6) gap wrt B

### Using julia functions in python

The following code runs: 

```
from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main

Main.include("./GraphAlignment/GA-IPA_main.jl")
Main.[function from imported file]
```

Note that the first time you use this, you might need to install a load of julia packages. Follow the instructions in the error message (as below):

```
from julia import Pkg
Pkg.add("SomePackageThatINeed")
```

Please also note that it's much easier to set up Julia if it is installed outside of the Python environment (because of the cert.pem file problem that lots of people have). If you want Julia's PyCall and PyPlot to point to a conda environment, you can then set it as follows: 

```
ENV["PYTHON"] = "/this/is/the/path/to/the/env"
```


