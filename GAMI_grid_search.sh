declare -a EpitopeArray=('GLCTLVAML' 'YLQPRTFLL') 
GA_folder='/media/graphgang/DATA/PairingOutputNoSmallStudy/GA/all_eps_oct2025/'
input_data='data/vdj_cleaned_subset_for_MI_no-small-study.csv'
output_folder='data/output/GAMI/grid_search/'


for ep in "${EpitopeArray[@]}" #scan over epitope
   do
      GA="${GA_folder}GA-epitope-${ep}_ind-no_vgene-no_distance_type-lev_kNN-20_n_repeats-100.csv.gz"
      /home/graphgang/miniconda3/envs/MutualInformation/bin/python PairingVDJdb_GAandMI.py --epitope $ep --input $input_data --MI_repair yes --GA_thresh 0.95 --Vgene no --GA_precomputed $GA --MI_repeats 10 --output $output_folder
      /home/graphgang/miniconda3/envs/MutualInformation/bin/python PairingVDJdb_GAandMI.py --epitope $ep --input $input_data --MI_repair no --GA_thresh 0.95 --Vgene no --GA_precomputed $GA --MI_repeats 10 --output $output_folder
      /home/graphgang/miniconda3/envs/MutualInformation/bin/python PairingVDJdb_GAandMI.py --epitope $ep --input $input_data --MI_repair yes --GA_thresh 0 --Vgene no --GA_precomputed $GA --MI_repeats 10 --output $output_folder
   done
