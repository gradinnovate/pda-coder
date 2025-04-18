#!/usr/bin/env bash
# gen_large_sample.sh: generate a large sample.txt for testing floorplanner
chipW=10000
chipH=10000
nsoft=50
nfixed=10
# Print header
echo "ChipSize $chipW $chipH"
echo "NumSoftModules $nsoft"
for ((i=1; i<=nsoft; i++)); do
    area=$((i * 1000))
    echo "SoftModule M$i $area"
done
echo "NumFixedModules $nfixed"
for ((i=1; i<=nfixed; i++)); do
    x=$(((i - 1) * 1000))
    echo "FixedModule FM$i $x 0 1000 1000"
done
# Generate nets: each soft module to a fixed module, plus chain of soft modules
declare -a nets_list
for ((i=1; i<=nsoft; i++)); do
    j=$(( (i % nfixed) + 1 ))
    w=$(( (i % 5) + 1 ))
    nets_list+=("Net M$i FM$j $w")
done
for ((i=1; i<nsoft; i++)); do
    nets_list+=("Net M$i M$((i+1)) 1")
done
echo "NumNets ${#nets_list[@]}"
for net in "${nets_list[@]}"; do
    echo "$net"
done