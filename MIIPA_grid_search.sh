: '
I use two epitopes because I want to check that it is general behaviour. GLC should have lots of signal, YLQ should have a 
bit less. They have 345 and 333 pairs each (so they should be equivalent in sample-size effect)

Default settings:

- vgene to no 
- method to none 
- prop_test to all 
- confidence to none 
- correlation to no 
- translation to no 
- small_ind to no 
- ones to keep 
- n_repeats to 10 
- step = 3  

Parameters I have swept:
- L
- w
- step
- vgene
- confidence

Other things I could check:
- corr
'

ep='GLCTLVAML'
ep1='YLQPRTFLL'
input='data/vdj_cleaned_subset_for_MI_no-small-study.csv'
output_folder='data/output/MIIPA/grid_search/'

#############################
##### scan over lambda and w
declare -a arrL=(0.01 0.2 0.4 0.6 0.8 1)
declare -a arrW=(0.01 0.2 0.4 0.6 0.8 0.99)

for l in "${arrL[@]}" #scan over lambda
 do
   echo ${l}
   for w in "${arrW[@]}" # scan over weight condition
      do
         echo ${w}
         python PairingVDJdb_MI.py --epitope $ep --input $input --L $l --weight $w --output "${output_folder}w_L_screen"
         python PairingVDJdb_MI.py --epitope $ep1 --input $input --L $l --weight $w --output "${output_folder}w_L_screen"
      done
 done

############################
#### scan over step size

l=0.6
w=0.6
declare -a arrS=(1 2 3 6 10 20 50 100 4 5 7 8 9 15 25 30 150 200 250 300)

for step in "${arrS[@]}"
 do
   echo ${step}
   python PairingVDJdb_MI.py --epitope $ep --input $input --L $l --weight $w --step $step --output "${output_folder}step_screen"
   python PairingVDJdb_MI.py --epitope $ep1 --input $input --L $l --weight $w --step $step --output "${output_folder}step_screen"
 done



## scan over confidence

l=0.6
l1=1
w=0.6
step=6
declare -a arrS=('greedy' 'none' 'hungarian')

for conf in "${arrS[@]}"
 do
   echo ${conf}
   python PairingVDJdb_MI.py --epitope $ep --input $input --L $l --weight $w --step $step --confidence $conf --output "${output_folder}confidence_screen"
   python PairingVDJdb_MI.py --epitope $ep --input $input --L $l1 --weight $w --step $step --confidence $conf --output "${output_folder}confidence_screen"
   python PairingVDJdb_MI.py --epitope $ep1 --input $input --L $l --weight $w --step $step --confidence $conf --output "${output_folder}confidence_screen"
   python PairingVDJdb_MI.py --epitope $ep1 --input $input --L $l1 --weight $w --step $step --confidence $conf --output "${output_folder}confidence_screen"
 done
