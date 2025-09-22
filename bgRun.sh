#!/bin/bash

file="run_cpu.pl"

Na=$(grep "^\$Na=" "$file" | awk -F '=' '{print $2}' | sed 's/[^0-9.-]*//g')
Nb=$(grep "^\$Nb=" "$file" | awk -F '=' '{print $2}' | sed 's/[^0-9.-]*//g')
tStop=$(grep "^\$tStop=" "$file" | awk -F '=' '{print $2}' | sed 's/[^0-9.-]*//g')
kba=$(grep "^\$kba=" "$file" | awk -F '=' '{print $2}' | sed 's/[^0-9.-]*//g')

# 修正bc命令的语法
N=$(echo "$Na + $Nb" | bc)
proportion=$(echo "scale=3; $Na / $Nb" | bc)

echo "Na=$Na, Nb=$Nb, N=$N, proportion=$proportion, tStop=$tStop, kba=$kba"

name=../conf-data/log/N${N}_kba${kba}_t${tStop}_p${proportion}.log
echo "name=$name"
mkdir -p ../conf-data
mkdir -p ../conf-data/log

nohup ./run_cpu.pl > $name 2>&1 &
echo "Simulation started, output is being logged to $name"