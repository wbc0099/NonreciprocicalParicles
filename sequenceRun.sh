#!/bin/bash 

N=1600

for i in 1 4 8 18; do
	for j in 1 5 10 20; do
		Na=$((N/(i+10)*i))
		Nb=$((N/(i+10)*10))
		kba=$(echo "scale=2; $j/100" | bc)
		sed -i "s!Na=*[0-9]*.[0-9]*;!Na=${Na};!" ./run_cpu.pl
		sed -i "s!Nb=[0-9]*.[0-9]*;!Nb=${Nb};!" ./run_cpu.pl
		sed -i "s!kba=-[0-9]*.[0-9]*;!kba=-${kba};!" ./run_cpu.pl
		nohup ./run_cpu.pl &
	done
done
