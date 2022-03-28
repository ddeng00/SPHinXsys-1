#!/usr/bin/bash
# generate_data.sh
# get directories
dir=$(pwd -P)
exec_path="$dir/heart"
out_dir="$dir/output"
data_dir="$dir/data"
# run initial relaxation
$exec_path -r &&
# ask for user input
read -e -p "Enter number of runs:" -i "200" run &&
read -e -p "Enter end time of each run:" -i "300" end &&
read -e -p "Enter number of stimuli:" -i "2" count &&
read -e -p "Enter radius of stimulus:" -i "10" radius &&
read -e -p "Enter sampling interval:" -i "10" interval &&
read -e -p "Enter time steps:" -i "3" step &&
# run multiple simulations in parallel
seq 0 $run | xargs -n 1 $exec_path --end_time $end --count $count --radius $radius --interval $interval --step $step --id &&
# run data processing
sphprocess $out_dir $data_dir $step