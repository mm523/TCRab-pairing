# TCRab-pairing
Scripts to reproduce the analyses in:

Milighetti M, Nagano Y, Henderson J, Hershberg U, Tiffeau-Mayer A, Bitbol AF, Chain B, **Intra- and inter-chain contacts determine TCR specificity: applying protein co-evolution methods to TCRαβ pairing**, bioRxiv 2024.05.24.595718; doi: https://doi.org/10.1101/2024.05.24.595718

#### NOTE: Graph alignment code

This code is from https://github.com/carlosgandarilla/GA-IPA, published as Gandarilla-Pérez CA, Pinilla S, Bitbol AF, Weigt M. Combining phylogeny and coevolution improves the inference of interaction partners among paralogous proteins. PLoS Comput Biol. 2023 Mar 30;19(3):e1011010. doi: 10.1371/journal.pcbi.1011010. PMID: 36996234; PMCID: PMC10089317.

Since the original implementation is in Julia, pyjulia was used to include the code in the pipeline. To make Julia work in Python, you need to install pyjulia: https://pyjulia.readthedocs.io/en/stable/index.html. PyJulia is already included in the environment, but you might need to run step 3 of installation. To make it work in conda, you also need the first hack of the troubleshooting guide.