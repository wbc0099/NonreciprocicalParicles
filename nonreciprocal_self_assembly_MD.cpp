#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <time.h>
#include <math.h>
#include <random> // The header for the generators.
#include <iomanip>
#include <omp.h>

//===============================================================
// Definitions
// Define the precision of real numbers, could be float/double.
#define real double
#define Pi 3.1415926535897932384626433832795
#define Zero 0

using namespace std;

// Initialize random number seed. To generate a random number, use:
// uniform_real_distribution<real> randUR; a=randUR(rng);
// or:
// uniform_int_distribution<int> randUI; a=randUI(rng);
using std::default_random_engine;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::normal_distribution;
int seed = time(0);
// int seed = 10032;
default_random_engine rng(seed);
uniform_real_distribution<real> randUR;
normal_distribution<real> randND;

// Declare variables --------------------------------------------
struct Configs {
    real * x;
    real * y;
    real * fx;
    real * fy;
    real * dx;
    real * dy;
    int * p1st;
    int * pNext;
    int * pNeighbor;
    real * dxExcess;
    real * dyExcess;
} C;

struct Parameters {
    int InitCond;
    int BounCond;
    int ExpoData;
    int CLorNL;
    int adapTStep;
    int N;
    int Na;
    int Nb;
    int NNeighb;
    real Lx;
    real Ly;
    real LCell;
    real LSkin;
    real tStep;
    real tStart;
    real tStop;
    real tExpo;
    real tResetCLNL;
    // Parameters for interactions
    real k0;
    real kaa;
    real kab;
    real kba;
    real kbb;
    real r0;
    real rd;
    real gamma;
    real kBT;
    // Other parameters
    int Ncx;
    int Ncy;
    int ReadSuccess;
} P;

int iStop;
real tNow;
real tStep;
real progress;

// Declare other global run-time variables.
clock_t tStart;

// Declare functions -------------------------------------------------------
// Host functions
void GetInput();
void InitConf();
void evolve();
void run();
void getCellList(Configs cc);
void getNeighborList(Configs cc);
void getForceCL(Configs cc);
void getForceNL(Configs cc);
void addForce(Configs cc, int i, int j, real dxe, real dye);
void update(Configs cc);
void BounPeri(Configs cc);
void BounRigi(Configs cc);
void bounMove(Configs cc);
void ShowProgress();
void ExpoConf(string str_t);
void MemAlloc();
void MemFree();
void MemAllocFields(Configs cc);
void MemFreeFields(Configs cc);

//===============================================================
int main() {
    // Get starting time of simulation.
    tStart = clock();

    // Get parameters from file.
    GetInput();
    if (P.ReadSuccess==-8848) {
        MemAlloc();
        InitConf();    
        evolve();
        MemFree();
    }
  
    clock_t tEnd = clock();
    cout<<"Simulation finished. ";
    cout<<"CPU time = ";
    cout<<double(tEnd - tStart)/CLOCKS_PER_SEC<<" sec"<<endl;
}

//===============================================================
void evolve() {
    std::string str_t;
    // This part is for testing -------------------------------  
    // iStop=1;
    // This part is for testing -------------------------------
  
    ShowProgress();
    if (P.tStart==0) {
        // getCellList(C);
        // getForce(C);
        ExpoConf("0");
    }

    while (iStop==0) {    
        
        // Evolution
        run();
    
        // Export.
        if (floor(tNow/P.tExpo)>floor((tNow-tStep)/P.tExpo)) {      
            int te=floor(tNow/P.tExpo);    
            str_t=to_string(te);
            ExpoConf(str_t);
            ShowProgress();      
        }

        // Check.
        if (tNow>P.tStop) {
            iStop=1;
        }
    
        // cout << tNow <<' ' <<iStop <<endl;
    }
}

//===============================================================
void run() {

    if (P.CLorNL==1) {        
        if (floor(tNow/P.tResetCLNL)>floor((tNow-tStep)/P.tResetCLNL)) {
        
            getCellList(C);
        }
        getForceCL(C);
    } else {
        if (floor(tNow/P.tResetCLNL)>floor((tNow-tStep)/P.tResetCLNL)) {
            getNeighborList(C);
        }
        getForceNL(C);
    };
        
    if (P.BounCond==1) {
        BounRigi(C);
    }
    update(C);
    if (P.BounCond==2) {
        BounPeri(C);
    }
  
    tNow=tNow+tStep;
}

// Get cell list ================================================
void getCellList(Configs cc) {
    int idx;
    int ic,jc;
    // Empty all cells.
    for (int j=0; j<P.Ncy; j++) {
        for (int i=0; i<P.Ncx; i++) {      
            idx=P.Ncx*j+i;
            cc.p1st[idx]=-1;
        }
    }
    for (int i=0; i<P.N; i++) {
        cc.pNext[i]=-1;
    }
    // Assign all particles.
    for (int i=0; i<P.N; i++) {
        ic=floor(cc.x[i]/P.LCell)+1;
        jc=floor(cc.y[i]/P.LCell)+1;
        if (ic<1) {
            ic=1;
        } else if (ic>P.Ncx-2) {
            ic=P.Ncx-2;
        }
        if (jc<1) {
            jc=1;
        } else if (jc>P.Ncy-2) {
            jc=P.Ncy-2;
        }
        idx=P.Ncx*jc+ic;
        cc.pNext[i]=cc.p1st[idx];
        cc.p1st[idx]=i;
    }
    if (P.BounCond==2) {
        int dj=P.Ncx;
        for (ic=0; ic<P.Ncx; ic++) {
            cc.p1st[ic]=cc.p1st[ic+dj*(P.Ncy-2)];
            cc.p1st[ic+dj*(P.Ncy-1)]=cc.p1st[ic+dj];
        }
        for (jc=0; jc<P.Ncy; jc++) {
            cc.p1st[jc*dj]=cc.p1st[jc*dj+dj-2];
            cc.p1st[jc*dj+dj-1]=cc.p1st[jc*dj+1];
        }
    }
}

// Get Neighbor list ============================================
void getNeighborList(Configs cc) {
    int idx,idxi0,idxj0,idx1,idx2;
    int i,j;
    real dr,dx,dy,dxe,dye;
    
    getCellList(cc);
    for (int i=0; i<P.N; i++) {
        idx=i*P.NNeighb;
        cc.pNeighbor[idx]=0;        
    }

        #pragma omp parallel for num_threads(8)
    for (int jc1=1; jc1<P.Ncy-1; jc1++) {
        for (int ic1=1; ic1<P.Ncx-1; ic1++) {    
            idx1=P.Ncx*jc1+ic1;
            i=cc.p1st[idx1];
            while (i>-1) {	
                for (int jc2=jc1-1; jc2<jc1+2; jc2++) {
                    for (int ic2=ic1-1; ic2<ic1+2; ic2++) {
                        idx2=P.Ncx*jc2+ic2;
                        j=cc.p1st[idx2];
                        while (j>-1) {
                            if (j!=i) {
                                dx=cc.x[i]-cc.x[j];
                                dy=cc.y[i]-cc.y[j];
                                dxe=0;
                                dye=0;
                                if (ic2==0) {
                                    dxe=P.Lx;
                                } else if (ic2==P.Ncx-1) {
                                    dxe=-P.Lx;
                                }
                                if (jc2==0) {
                                    dye=P.Ly;
                                } else if (jc2==P.Ncy-1) {
                                    dye=-P.Ly;
                                }
                                dr=sqrt((dx+dxe)*(dx+dxe)+(dy+dye)*(dy+dye));                     
                                if (dr<P.r0+2*P.rd+P.LSkin) {
                                    idxi0=i*P.NNeighb;
                                    cc.pNeighbor[idxi0]+=1;
                                    
                                    cc.pNeighbor[idxi0+cc.pNeighbor[idxi0]]=j;
                                    cc.dxExcess[idxi0+cc.pNeighbor[idxi0]]=dxe;
                                    cc.dyExcess[idxi0+cc.pNeighbor[idxi0]]=dye;
                                }
                            }
                            j=cc.pNext[j];
                        }
                    }
                }
                i=cc.pNext[i];
            }
        }
    }
}

//===============================================================
// Get inter-particle force using cell list
void getForceCL(Configs cc) {
    real dr,dxij,dyij;
    real dxe,dye;
    int idx1,idx2,i,j;
    // Initialize force
    for (int i=0; i<P.N; i++) {
        cc.fx[i]=0;
        cc.fy[i]=0;
    }

        #pragma omp parallel for
    for (int jc1=1; jc1<P.Ncy-1; jc1++) {
        for (int ic1=1; ic1<P.Ncx-1; ic1++) {    
            idx1=P.Ncx*jc1+ic1;
            i=cc.p1st[idx1];
            while (i>-1) {	
                for (int jc2=jc1-1; jc2<jc1+2; jc2++) {
                    for (int ic2=ic1-1; ic2<ic1+2; ic2++) {
                        idx2=P.Ncx*jc2+ic2;
                        j=cc.p1st[idx2];                        
                        while (j>-1) {
                            if (j!=i) {
                                dxe=0;
                                dye=0;
                                if (ic2==0) {
                                    dxe=P.Lx;
                                } else if (ic2==P.Ncx-1) {
                                    dxe=-P.Lx;
                                }
                                if (jc2==0) {
                                    dye=P.Ly;
                                } else if (jc2==P.Ncy-1) {
                                    dye=-P.Ly;
                                }
                                addForce(cc,i,j,dxe,dye);
                            }
                            j=cc.pNext[j];
                        }
                    }
                }
                i=cc.pNext[i];
            }
        }
    }
}


//===============================================================
// Get inter-particle force using neighbor list
void getForceNL(Configs cc) {
    real dr,dxij,dyij;
    real dxe,dye;
    int idxi0;
    // Initialize force
    for (int i=0; i<P.N; i++) {
        cc.fx[i]=0;
        cc.fy[i]=0;
    }

        #pragma omp parallel for num_threads(8)
    for (int i=0; i<P.N; i++) {
        idxi0=i*P.NNeighb;
        for (int j0=1; j0<cc.pNeighbor[idxi0]+1; j0++) {
            addForce(cc,i,cc.pNeighbor[idxi0+j0],cc.dxExcess[idxi0+j0],cc.dyExcess[idxi0+j0]);
        }
    }
}

//===============================================================
// Get inter-particle force
void addForce(Configs cc, int i, int j, real dxe, real dye) {
    real r0,rd,dr,ddr,dx12,dy12;    
    real f12,f21,k12,k21,ftemp;

    dx12=cc.x[i]-cc.x[j]+dxe;
    dy12=cc.y[i]-cc.y[j]+dye;
    dr=sqrt(dx12*dx12+dy12*dy12);
    r0=P.r0;
    rd=P.rd;
    
    if (dr<r0+rd) {
        // Repulsive forces
        // f12=P.k0/pow(r0,4)*pow(ddr,4)/dr;
        if (dr<r0) {
            f12=P.k0*(r0*r0*r0*r0/dr/dr/dr/dr-1)/dr;
            cc.fx[i] += f12*dx12;
            cc.fy[i] += f12*dy12;
            // cc.fx[j] += -f12*dx12;
            // cc.fy[j] += -f12*dy12;
        }
    // } else if (dr<r0) {
        // Active forces
        // Get prefactors
        if (i<P.Na) {
            if (j<P.Na) {
                k12=P.kaa;
                // k21=P.kaa;
            } else {
                k12=P.kab;
                // k21=P.kba;
            }
        } else {
            if (j<P.Na) {
                k12=P.kba;
                // k21=P.kab;
            } else {
                k12=P.kbb;
                // k21=P.kbb;
            }
        }

        // Get force amplitudes
        // ftemp=r0*r0/dr/dr/dr;
        ftemp=((r0+rd)*(r0+rd)/dr/dr-1)/dr;
        cc.fx[i] += k12*ftemp*dx12;
        cc.fy[i] += k12*ftemp*dy12;
        // cc.fx[j] += -k21*ftemp*dx12;
        // cc.fy[j] += -k21*ftemp*dy12;
    };
}

//===============================================================
// Update coordinates
void update(Configs cc) {
    real fT=sqrt(2*P.kBT*P.gamma*tStep);
    for (int i=0; i<P.N; i++) {
        real FRx=randND(rng);
        real FRy=randND(rng);
        cc.dx[i]=(cc.fx[i]*tStep+fT*FRx)/P.gamma;
        cc.dy[i]=(cc.fy[i]*tStep+fT*FRy)/P.gamma;
        cc.x[i]+=cc.dx[i];
        cc.y[i]+=cc.dy[i];
    }
}

//===============================================================
// Rigid wall boundary condition
void BounRigi(Configs cc) {
    real r0=P.r0;
    real dx,dy;
    real kf=P.k0/r0/r0;
    
    for (int i=0; i<P.N; i++) {
        if (cc.x[i]<r0) {
            dx=r0-cc.x[i];
            cc.fx[i]+=kf*dx*dx;
        } else if (cc.x[i]>P.Lx-r0) {
            dx=P.Lx-r0-cc.x[i];
            cc.fx[i]-=kf*dx*dx;
        }
        if (cc.y[i]<r0) {
            dy=r0-cc.y[i];
            cc.fy[i]+=kf*dy*dy;
        } else if (cc.y[i]>P.Ly-r0) {
            dy=P.Ly-r0-cc.y[i];
            cc.fy[i]-=kf*dy*dy;
        }
    }
}

//===============================================================
// Periodic boundary condition
void BounPeri(Configs cc) {
    for (int i=0; i<P.N; i++) {
        if (cc.x[i]<0) {
            cc.x[i]=cc.x[i]+P.Lx;
        } else if (cc.x[i]>P.Lx) {
            cc.x[i]=cc.x[i]-P.Lx;
        }
        if (cc.y[i]<0) {
            cc.y[i]=cc.y[i]+P.Ly;
        } else if (cc.y[i]>P.Ly) {
            cc.y[i]=cc.y[i]-P.Ly;
        }
    }
}

//===============================================================
void GetInput() {
    ifstream InputFile ("input.dat");

    P.ReadSuccess=0;
    // Simulation parameters  
    InputFile >> P.InitCond;
    InputFile >> P.BounCond;
    InputFile >> P.ExpoData;
    InputFile >> P.CLorNL;
    InputFile >> P.adapTStep;
    InputFile >> P.Na;
    InputFile >> P.Nb;
    InputFile >> P.NNeighb;
    InputFile >> P.Lx;
    InputFile >> P.Ly;
    InputFile >> P.LCell;
    InputFile >> P.LSkin;  
    InputFile >> P.tStart;  
    InputFile >> P.tStop;
    InputFile >> P.tStep;  
    InputFile >> P.tExpo;
    InputFile >> P.tResetCLNL;

    // Model parameters
    InputFile >> P.k0;
    InputFile >> P.kaa;
    InputFile >> P.kab;
    InputFile >> P.kba;
    InputFile >> P.kbb;  
    InputFile >> P.r0;
    InputFile >> P.rd;
    InputFile >> P.gamma;    
    InputFile >> P.kBT;
    InputFile >> P.ReadSuccess;
  
    P.N=P.Na+P.Nb;
    // Get size of cells
    // Assume we can fit integer number of cells in both x&y.
    // The additional two boxes in each dimension accounts for periodic boundary.
    P.Ncx=floor((P.Lx+0.001*P.LCell)/(P.LCell+0.0))+2;
    P.Ncy=floor((P.Ly+0.001*P.LCell)/(P.LCell+0.0))+2;
  
    InputFile.close();
    if (P.ReadSuccess!=-8848) {
        cout << "Error while reading the input file!" << endl;
    }
  
}

//===============================================================
void InitConf() {
    tNow=P.tStart;
    iStop=0;
    tStep=P.tStep;

    int te0=floor((P.tStart+0.1*P.tExpo)/P.tExpo);
    std::string InitFileName="data/conf_" +
        std::to_string(te0) + ".dat";  
    ifstream InputFile (InitFileName);

    for (int i=0; i<P.N; i++) {
        if (P.InitCond==0) {
            InputFile >> C.x[i] >> C.y[i];
        } else if (P.InitCond==1) {
            real ran1=randUR(rng);
            real ran2=randUR(rng);
            C.x[i]=P.Lx*ran1;
            C.y[i]=P.Ly*ran2;      
        } else if (P.InitCond==2) {
            real ran1=randUR(rng);
            real ran2=randUR(rng);
            C.y[i]=P.Ly*ran2;
            if (i<P.Na) {
                C.x[i]=P.Lx*(0.5+(ran1-0.5)*1*P.Na/(P.N+0.0));                
            } else {
                C.x[i]=P.Lx*((ran1-0.5)*1*P.Nb/(P.N+0.0));
                if (C.x[i]<0) {
                    C.x[i]=C.x[i]+P.Lx;
                }	                
            }

            if (i==0) {
                real ran1=randUR(rng);
                real ran2=randUR(rng);
                C.y[i]=P.Ly*ran2;
                C.x[i]=P.Lx*(0.5+(ran1-0.5)*1*P.Na/(P.N+0.0));                
            } else {
                real drmin=0;
                real drij2=0;
                while (drmin<0.9*P.r0*P.r0) {
                    real ran1=randUR(rng);
                    real ran2=randUR(rng);
                    C.y[i]=P.Ly*ran2;
                    if (i<P.Na) {
                        C.x[i]=P.Lx*(0.5+(ran1-0.5)*1*P.Na/(P.N+0.0));                
                    } else {
                        C.x[i]=P.Lx*((ran1-0.5)*1*P.Nb/(P.N+0.0));
                        if (C.x[i]<0) {
                            C.x[i]=C.x[i]+P.Lx;
                        }	                
                    }
                    drmin=2*P.r0*P.r0;
                    for (int j=0; j<i; j++) {
                        drij2=(C.x[i]-C.x[j])*(C.x[i]-C.x[j])+(C.y[i]-C.y[j])*(C.y[i]-C.y[j]);
                        if (drij2<drmin) {
                            drmin=drij2;
                        }
                    }                    
                }
            }
        } else if (P.InitCond==3) {
            if (i==0) {
            real ran1=randUR(rng);
            real ran2=randUR(rng);
            C.x[0]=P.Lx*ran1;
            C.y[0]=P.Ly*ran2;
            } else {
                real drmin=0;
                real drij2=0;
                while (drmin<0.9*P.r0*P.r0) {
                    real ran1=randUR(rng);
                    real ran2=randUR(rng);
                    C.x[i]=P.Lx*ran1;
                    C.y[i]=P.Ly*ran2;                    
                    drmin=2*P.r0*P.r0;
                    for (int j=0; j<i; j++) {
                        drij2=(C.x[i]-C.x[j])*(C.x[i]-C.x[j])+(C.y[i]-C.y[j])*(C.y[i]-C.y[j]);
                        if (drij2<drmin) {
                            drmin=drij2;
                        }
                    }                    
                }
            }
        }
    }

    if (P.InitCond==4) {    
        real dr0=1;
        int ii=0;
        real x0;
        real drm;
        real drm2;
        real drij2;
        for (int i=0; i<15; i++) {
            x0=-8.5*dr0+0.5*P.Lx;
            for (int j=0; j<18; j++) {
                C.x[ii]=x0+dr0*j;
                C.y[ii]=0.5+2*i*dr0+P.Ly/2;
                ii=ii+1;
            }
        }
        for (int i=0; i<15; i++) {
            x0=-8*dr0+0.5*P.Lx;
            for (int j=0; j<17; j++) {
                C.x[ii]=x0+dr0*j;
                C.y[ii]=1.5+2*i*dr0+P.Ly/2;
                ii=ii+1;
            }
        }
        for (int i=P.Na; i<P.Na+P.Nb; i++) {
            drm=0;
            while (drm<2*P.r0) {
                real ran1=randUR(rng);
                real ran2=randUR(rng);
                C.y[i]=P.Ly*ran2;
                C.x[i]=P.Lx*((ran1-0.5)*1.5*P.Nb/(P.N+0.0));
                if (C.x[i]<0) {
                    C.x[i]=C.x[i]+P.Lx;
                }
                drm2=10;
                for (int j=0; j<P.Na; j++) {
                    drij2=(C.x[i]-C.x[j])*(C.x[i]-C.x[j])+(C.y[i]-C.y[j])*(C.y[i]-C.y[j]);
                    if (drij2<drm2) {
                        drm2=drij2;
                    }
                }
                drm=pow(drm2,0.5);
            }
        }
    }
  
    if (P.InitCond==5) {    
        real dr0=1;
        int ii=0;
        real x0;
        real drm;
        real drm2;
        real drij2;
        for (int i=0; i<15; i++) {
            x0=-8.5*dr0+0.5*P.Lx;
            for (int j=0; j<18; j++) {
                C.x[ii]=x0+dr0*j;
                C.y[ii]=0.5+2*i+P.Ly/2;
                ii=ii+1;
            }
        }
        for (int i=0; i<15; i++) {
            x0=-8*dr0+0.5*P.Lx;
            for (int j=0; j<17; j++) {
                C.x[ii]=x0+dr0*j;
                C.y[ii]=1.5+2*i+P.Ly/2;
                ii=ii+1;
            }
        }
        for (int i=P.Na; i<P.Na+P.Nb; i++) {
            drm=0;
            while (drm<2*P.r0) {
                real ran1=randUR(rng);
                real ran2=randUR(rng);
                C.y[i]=P.Ly*ran2;
                C.x[i]=P.Lx*((ran1-0.5)*1.5*P.Nb/(P.N+0.0));
                if (C.x[i]<0) {
                    C.x[i]=C.x[i]+P.Lx;
                }
                drm2=10;
                for (int j=0; j<P.Na; j++) {
                    drij2=(C.x[i]-C.x[j])*(C.x[i]-C.x[j])+(C.y[i]-C.y[j])*(C.y[i]-C.y[j]);
                    if (drij2<drm2) {
                        drm2=drij2;
                    }
                }
                drm=pow(drm2,0.5);
            }
        }
    }
    InputFile.close();
}

//===============================================================
void ExpoConf(string str_t) {
    ofstream ConfFile;
    int PrecData=8;
  
    std::string ConfFileName="conf_" + str_t + ".dat";
    ConfFile.open ( ConfFileName.c_str() );

    for (int idx=0; idx<P.N; idx++) {
        if (P.ExpoData==1) {
            ConfFile<< fixed << setprecision(PrecData) <<C.x[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.y[idx]<<endl;
        } else if (P.ExpoData==-1) {
            ConfFile<< fixed << setprecision(PrecData) <<C.x[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.y[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.fx[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.fy[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.dx[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.dy[idx]<<' ';
            ConfFile<< fixed << setprecision(PrecData) <<C.pNeighbor[idx*P.NNeighb]<<endl;
        }
    }
    ConfFile.close();
}

//===============================================================
void MemAlloc() {
    // Allocate fields in host memory.
    C.x=new real[P.N];
    C.y=new real[P.N];
    C.fx=new real[P.N];
    C.fy=new real[P.N];
    C.dx=new real[P.N];
    C.dy=new real[P.N];
    C.pNext=new int[P.N];
    C.p1st=new int[P.Ncx*P.Ncy];
    C.pNeighbor=new int[P.N*P.NNeighb];
    C.dxExcess=new real[P.N*P.NNeighb];
    C.dyExcess=new real[P.N*P.NNeighb];
}

//===============================================================
void MemFree() {
    // Free host memory
    delete [] C.x;
    delete [] C.y;
    delete [] C.fx;
    delete [] C.fy;
    delete [] C.dx;
    delete [] C.dy;
    delete [] C.pNext;
    delete [] C.p1st;
    delete [] C.pNeighbor;
}

//===============================================================
void ShowProgress() {
    // Print progress.
    progress=(tNow-P.tStart)/(P.tStop-P.tStart);
    int barWidth = 50;
    clock_t tNow = clock();
    double tUsed=double(tNow-tStart)/CLOCKS_PER_SEC;
  
    std::cout << "Progress: ";  
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %";
    if (tNow==0) {
        std::cout <<"\r";
    } else {
        std::cout << ".  " << floor(tUsed/progress*(1-progress)) << "s remains.\r";
        // std::cout << ".  " << floor(tUsed/progress*(1-progress)) << "s, "<<tStep<<".\r";
    }

    std::cout.flush();
}

//===============================================================
