# TCRab-pairing
Scripts to reproduce the analyses in:

Milighetti M, Nagano Y, Henderson J, Hershberg U, Tiffeau-Mayer A, Bitbol AF, Chain B, **Intra- and inter-chain contacts determine TCR specificity: applying protein co-evolution methods to TCRαβ pairing**, bioRxiv 2024.05.24.595718; doi: https://doi.org/10.1101/2024.05.24.595718

## List of scripts and figures generated

#### Data preprocessing

All the data preprocessing scripts are in data_preprocessing.

1. Structures:  ```clean_stcrdab_export.ipynb``` cleans up structures that I am not interested in, including edge cases. Then checks that the information I extract is what I want.
2. vdjdb: 
    1. ```preprocess_vdjdb.py``` does some initial cleaning of the export
    2. ```apply_stitchr.py``` takes the output from 1. and uses stitchr (https://github.com/JamieHeather/stitchr) to create full TCR sequences
    3. ```vdj_clean_for_MI.py``` applies the final cleaning steps and extracts cdrs for all MI applications.
3. Tanno dataset:
    1. ```Tanno_preprocessing.ipynb``` performs pre-processing as in Mayer and Callan, 2023
    2. ```tanno_clean_for_MI.py``` applies the final cleaning steps and extracts cdrs for all MI applications.
4. OTS dataset: ```OTS_preprocessing.ipynb``` selects the studies of interest and prepares data for MI calculation

#### Analysis of structures

The following scripts are in ```structure_analysis```

1. ```interloop_contacts_4structures.ipynb``` plots the circos plots for the contacts in the 4 example TCRs in Fig 1
2. ```extract_all_intrachain_pairwise_distance.ipynb``` and ```extract_all_interchain_pairwise_distance.ipynb``` extract all distances between and across chains and save them to file
3. ```intrachain_pairwise_distances.ipynb``` extracts and plots all pairwise distances between residues on the same chain
4. ```interchain_pairwise_distances.ipynb``` extracts and plots all pairwise distances between residues on different chains
5. ```threshold_analysis.ipynb``` shows the effect of choosing different thresholds on the number of detected contacts

These scripts reproduce all analysis in Figs 1-3, S1 and S9-S11

The functions needed to run these scripts are in the ```functions``` folder and are called directly within each script.

#### Analysis of VDJDb clustering

```cluster_alpha_beta.ipynb``` reproduces all the analysis in Figs 4, S12-S13.

#### MI calculation and analysis

1. ```totalMI_vdjdb.ipynb``` performs MI calculations for various subsamples from each epitope from the VDJDb dataset. Here padding for the MI calculations is added in the middle (as with IMGT numbering).
2. ```totalMI_tanno_oneNaiveSample.ipynb``` performs MI calculations for various subsamples from the naive sample from individual A1 from Tanno et al. Here padding for the MI calculations is added in the middle (as with IMGT numbering).
3. ```totalMI_vdjdb_changepadding.ipynb``` and ```totalMI_tanno_oneNaiveSample_changepadding.ipynb``` perform the same calculations as above, but by padding the sequences at the end instead of padding in the middle.
4. ```totalMI_analysis.ipynb``` performs the extrapolation with the subsamples and correction of MI by shuffle for each calculated MI category. Equivalent in ```totalMI_analysis_changepadding.ipynb``` for padding at the end.
5. ```V_to_cdr3_MI.ipynb``` plots the analyses in Figs 5-6 and S17 (with padding in the middle).
6. ```totalMI_compare_padding.ipynb``` compares the MI results depending on the padding, as in Fig S15.
7. ```MI_vs_gaps.ipynb``` evaluates the effects of gaps as in Fig S14
8. ```totalMI_ots.ipynb```, ```totalMI_analysis_ots.ipynb``` and ```OTS_plot_results.ipynb``` calculate MI for OTS set and analyse the results as in Fig S16

#### Pairing algorithms

##### MI-IPA

The functions to run the MI-IPA on the VDJDb dataset are contained in ```PairingVDJdb_MI.py```.

The script ```run_MIIPA_benchmark_eps``` runs the grid search on epitopes GLC and YLQ to find the best parameters.  The notebooks in ```VDJ_MI-IPA_results``` (```*_screen.ipynb```) can load the results from the grid search and plot the results as in Fig S5.

```run_MIIPA_all_eps.sh``` can be used to run all the epitopes from the command line. The parameters can be adjusted as desired. The notebook ```VDJ_MI-IPA_results/all_epitope_repeats.ipynb``` can load the results from the paring and plot the results as in Fig S18.

The notebook ```VDJ_MI-IPA_results/all_epitope_corr_with_MI.ipynb``` plots the correlations with MI as in Fig 9.

##### GA

The functions to run the MI-IPA on the VDJDb dataset are contained in ```PairingVDJdb_GA.py```.

> **NOTE: Graph alignment code.** 
This code is from https://github.com/carlosgandarilla/GA-IPA, published as Gandarilla-Pérez CA, Pinilla S, Bitbol AF, Weigt M. Combining phylogeny and coevolution improves the inference of interaction partners among paralogous proteins. PLoS Comput Biol. 2023 Mar 30;19(3):e1011010. doi: 10.1371/journal.pcbi.1011010. PMID: 36996234; PMCID: PMC10089317. 
Since the original implementation is in Julia, pyjulia was used to include the code in the pipeline. To make Julia work in Python, you need to install pyjulia: https://pyjulia.readthedocs.io/en/stable/index.html. PyJulia is already included in the environment, but you might need to run step 3 of installation. To make it work in conda, you also need the first hack of the troubleshooting guide.

The script ```run_GA_benchmark_eps.sh``` runs the grid search on epitopes GLC and YLQ to find the best parameters.  The notebooks ```VDJ_GA_results/Kscreen.ipynb``` can load the results from the grid search and plot the results as in Fig S6.

```run_GA_all_eps_lev_20.sh``` can be used to run all the epitopes from the command line. The parameters can be adjusted as desired. The notebook ```VDJ_GA_results/repeat_results.ipynb``` can load the results from the paring and plot the results as in Fig S19.

The notebook ```VDJ_GA_results/GA_perf_corr.ipynb``` plots the correlations with other clustering metrics as in Fig 9.

##### GA+MI-IPA

The functions to run the GA + MI-IPA on the VDJDb dataset are contained in ```PairingVDJdb_GAandMI.py```.

The script ```run_GAMI_benchmark_eps.sh``` runs the grid search on epitopes GLC and YLQ to find the best parameters.  The notebooks ```VDJ_GAMI_results/GLA_YLQ_repair_over_iterations.ipynb``` can load the results from the grid search and plot the results as in Fig S7.

```run_GAMI_all_eps.sh``` can be used to run all the epitopes from the command line. The parameters can be adjusted as desired. The notebook ```VDJ_GAMI_results/all_epitope_repeats.ipynb``` can load the results from the paring and plot the results as in Fig S20.

The notebook ```VDJ_MI-IPA_results/all_epitope_correlations.ipynb``` plots the correlations with other metrics as in Fig 9.

The summary result in Fig 8 can be generated with ```all_epitopes_pairing_GA_and_MI.ipynb```.

##### Benchmarking - TULIP and SCEPTR

TULIP and SCEPTR can be run using the scripts in ```pairing_prediction_other_software```. They assume that TULIP is ```git-cloned``` in a directory outside of the ```TCRab-pairing``` directory which the code accesses directly. SCEPTR can be ```pip install``` following instructions of the package.

1. The script ```make_dataset.ipynb``` creates the combinatorial pairs to be scored by TULIP. 
2. The script ```run_tulip.ipynb``` runs the prediction with TULIP and ```analyse_TULIP_results.ipynb``` outputs the pairing results
3. The script ```sceptr_prediction.ipynb``` runs prediction with SCEPTR and outputs the results
4. The script ```benchmarking_summary.ipynb``` produces the figures in Fig S21
