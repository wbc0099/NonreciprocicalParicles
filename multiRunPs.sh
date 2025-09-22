#!/bin/bash
N=1600

mkdir -p ../conf-data
mkdir -p ../conf-data/log

for i in $(seq 2 2 18); do
	for j in $(seq 1 2 17); do
	#for j in 1 5 9 13 17; do
		Nb=$((N/(i+10)*10))
		Na=$((N-Nb))
		kba=$(echo "scale=2; $j/100" | bc)
		sed -i "s!Na=*[0-9]*.[0-9]*;!Na=${Na};!" ./run_cpu.pl
		sed -i "s!Nb=[0-9]*.[0-9]*;!Nb=${Nb};!" ./run_cpu.pl
		sed -i "s!kba=-[0-9]*.[0-9]*;!kba=-${kba};!" ./run_cpu.pl
		echo "Running with Na=$Na, Nb=$Nb, kba=$kba"
		nohup ./run_cpu.pl > ../conf-data/log/{$Na}_{$Nb}_{$kba}_error.txt &
		sleep 1
	done
done