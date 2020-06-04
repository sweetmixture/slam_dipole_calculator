#!/bin/bash

main_input=$1
cur_path=$PWD

python_src_path="/home/uccawkj/src_tool/slam_toolkit/slam_dipole_package/calc_slam_dipole.py"

config_tmp="config"
mo_tmp="mo"
species_tmp="type"

if [ -z $main_input ];	then
	print "Error ... take the first argument with slam optimisation log file!"
	exit 1
fi

if [ ! -f "$cur_path"/"$main_input" ]; then
	print "Error ... main input file does not exist in the specified path!"
	exit 1
fi

POS_INTEGRAL=$( grep "Position Integral Reference :" $main_input | awk '{print $5}' )
MM_CNT=$(grep -A 1000 "Optimisation Meets Termination Condition, Final Configuration is" $main_input | grep -A 1 "CONFIGURATION_XYZ_SC_INFO" | tail -1 | awk '{print $1}')
QM_CNT=$(grep -A 1000 "Optimisation Meets Termination Condition, Final Configuration is" $main_input | grep -A 1 "CONFIGURATION_XYZ_SC_INFO" | tail -1 | awk '{print $2}')
CNT=$( echo "$MM_CNT + $QM_CNT" | bc )

readline_number=$( echo "$CNT + 1" | bc )
grep -A 1000 "Optimisation Meets Termination Condition, Final Configuration is" $main_input | grep -A $readline_number "CONFIGURATION_XYZ_SC_INFO" | tail -"$readline_number" > $config_tmp

echo "$QM_CNT" > $mo_tmp
readline_number=$( echo "$QM_CNT + 4" | bc )
grep -A 1000 "Optimisation Meets Termination Condition, Final Configuration is" $main_input | grep -A $readline_number "Lone Pair Molecular Orbital" | tail -"$QM_CNT" >> $mo_tmp

echo $MM_CNT $QM_CNT > $species_tmp
readline_number=$( echo "$MM_CNT + 4" | bc )
grep -A $readline_number "MM atoms/ions" $main_input | tail -"$MM_CNT" | awk '{ printf "%2s%2s%12.4f\n", $1, $2, $6}' >> $species_tmp
readline_number=$( echo "$QM_CNT + 4" | bc )
grep -A $readline_number "SP atoms/ions" $main_input | tail -"$QM_CNT" | awk '{ printf "%2s%20s\n", $1, $5}' >> $species_tmp
sed -i "s/\//   /g" $species_tmp


python $python_src_path $POS_INTEGRAL $config_tmp $mo_tmp $species_tmp $MM_CNT $QM_CNT > dip_slam.dat

rm $config_tmp $mo_tmp $species_tmp
