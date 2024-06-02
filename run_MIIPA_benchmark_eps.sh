: '
Script to run MI-IPA on two benchmarking epitopes to find best lambda, theta, step and confidence
'

##### scan over lambda and w
declare -a arrL=(0.01 0.2 0.4 0.6 0.8 1)
declare -a arrW=(0.01 0.2 0.4 0.6 0.8 0.99)
ep='GLCTLVAML'
ep1='YLQPRTFLL'

for l in "${arrL[@]}" #scan over lambda
 do
   echo ${l}
   for w in "${arrW[@]}" # scan over weight condition
      do
         echo ${w}
         python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --output w_L_screen
         python PairingVDJdb_MI.py --epitope $ep1 --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --output w_L_screen
      done
 done

##### scan over step size
l=0.6
w=0.6
declare -a arrS=(1 2 3 4 5 7 8 9 15 20 25 30 50 100 150 200 250 300)
ep='GLCTLVAML'
ep1='YLQPRTFLL'

for step in "${arrS[@]}"
 do
   echo ${step}
   python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --output step_screen
   python PairingVDJdb_MI.py --epitope $ep1 --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --output step_screen
 done

### scan over confidence

l=0.6
l1=1
w=0.6
step=6
declare -a arrS=('greedy' 'none' 'hungarian')
ep='GLCTLVAML'
ep1='YLQPRTFLL'

for conf in "${arrS[@]}"
 do
   echo ${conf}
   python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --confidence $conf --output confidence_screen
   python PairingVDJdb_MI.py --epitope $ep --input data/vdj_cleaned_subset_for_MI.csv --L $l1 --weight $w --step $step --confidence $conf --output confidence_screen
   python PairingVDJdb_MI.py --epitope $ep1 --input data/vdj_cleaned_subset_for_MI.csv --L $l --weight $w --step $step --confidence $conf --output confidence_screen
   python PairingVDJdb_MI.py --epitope $ep1 --input data/vdj_cleaned_subset_for_MI.csv --L $l1 --weight $w --step $step --confidence $conf --output confidence_screen
 done