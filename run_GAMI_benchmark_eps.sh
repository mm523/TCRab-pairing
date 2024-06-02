: '
Script to run GA+MI-IPA on two benchmarking epitopes to find best k and best dist
'

declare -a EpitopeArray=('GLCTLVAML' 'YLQPRTFLL') 

for ep in "${EpitopeArray[@]}"
   do
      GA="data/output/pairing_GA/all_eps/GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --MI_repair yes --GA_thresh 0.95 --Vgene no --output lev20_repair_screen --GA_precomputed $GA --MI_repeats 10
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --MI_repair no --GA_thresh 0.95 --Vgene no --output lev20_repair_screen --GA_precomputed $GA --MI_repeats 10
      python PairingVDJdb_GAandMI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --MI_repair yes --GA_thresh 0 --Vgene no --output lev20_repair_screen --GA_precomputed $GA --MI_repeats 10
   done