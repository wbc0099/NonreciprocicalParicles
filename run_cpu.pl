#!/usr/bin/perl -w

# system("g++ -O3 -lm -fopenmp nonreciprocal_self_assembly_MD.cpp -o nonreciprocal_self_assembly_MD");

# Simulation parameters
$InitCond    = 3;             # 0-import, 1-Uniform, 2-Band
$BounCond    = 2;             # 1-Rigid, 2-Periodic
$ExpoData    = 1;             # 1-only rho, -1-all fields
$CLorNL      = 1;             # 1-Cell List, 2-Neighbor List
$adapTStep   = 0;             # 
<<<<<<< HEAD
$Na=429;
$Nb=1071;
=======
$Na=457;
$Nb=1143;
>>>>>>> f33c8b4 (Initial commit)
$NNeighb     = 21;
$Lx          = 50;             # Boundary width.
$Ly          = 50;             # Boundary width.
$LCell       = 4;
$LSkin       = 0.4;
$tStart      = 0;
$tStop=10000;             #
$tStep       = 0.01;     #
$tExpo       = 10;         # Time resolution of output
$tResetCLNL    = 0.1;
$ReadSuccess = -8848;        # Used to check if input file is correctly read.

# Model parameters
$k0=1;
$kaa=0.01;
$kab=0.01;
$kba=-.05;
$kbb=0.01;
$r0=1;
$rd=2;
$gamma=1;
$kBT=0.0001;

open (FILE, '>input.dat');  #simulation parameters
print FILE "$InitCond\n";
print FILE "$BounCond\n";
print FILE "$ExpoData\n";
print FILE "$CLorNL\n";
print FILE "$adapTStep\n";
print FILE "$Na\n";
print FILE "$Nb\n";
print FILE "$NNeighb\n";
print FILE "$Lx\n";
print FILE "$Ly\n";
print FILE "$LCell\n";
print FILE "$LSkin\n";
print FILE "$tStart\n";  
print FILE "$tStop\n";  
print FILE "$tStep\n";
print FILE "$tExpo\n";
print FILE "$tResetCLNL\n";
print FILE "$k0\n";
print FILE "$kaa\n";
print FILE "$kab\n";
print FILE "$kba\n";
print FILE "$kbb\n";
print FILE "$r0\n";
print FILE "$rd\n";
print FILE "$gamma\n";
print FILE "$kBT\n";
print FILE "$ReadSuccess\n";
close (FILE);

$N=$Na+$Nb;
$proportion = sprintf("%.2f", $Na/$Nb);
$direData="N$N--proportion$proportion--kba$kba";

for ($sAddi = 1; $sAddi <1000; $sAddi=$sAddi+1) {
    if (-d "../data$sAddi") {
    } else {
	last;
    }
};

system("mkdir ../data$sAddi");
system("cp run_cpu.pl ../data$sAddi/");
system("mv input.dat ../data$sAddi/");
system("cp nonreciprocal_self_assembly_MD.cpp ../data$sAddi");
chdir("../data$sAddi");
system("g++ -O3 -lm nonreciprocal_self_assembly_MD.cpp -o nonreciprocal_self_assembly_MD$sAddi");
system("./nonreciprocal_self_assembly_MD$sAddi");

if (-d "../conf-data/$direData") {
    system("rm -f ../conf-data/$direData/*");
}
else { system("mkdir -p ../conf-data/$direData"); }

chdir("../");
system("mv data$sAddi/* conf-data/$direData");
system("rm -rf data$sAddi");
system("python3 source/plot.py conf-data/$direData");
<<<<<<< HEAD
system("vlc conf-data/$direData.avi")
=======
#system("vlc conf-data/$direData.avi")
>>>>>>> f33c8b4 (Initial commit)

#system("python plot.py $direData")
