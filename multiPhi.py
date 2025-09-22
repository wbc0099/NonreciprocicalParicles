import tools
import os
import numpy as np
import matplotlib.pyplot as plt

path="../conf-data/"
items=os.listdir(path)
items=sorted(items,key=lambda x:
             (float(x.split("proportion")[1].split("--")[0]),
              float(x.split("kba-")[1])))
print(items)
N=[]
for item in items:
    data=os.path.join(path,item)+"/conf_1000.dat"
    data=tools.readFile(data)
    redParticleNumber=tools.readParameter(
        os.path.join(path,item)+"/run_cpu.pl",["$Na"])[0]
    bluePartilceNumber=tools.readParameter(
        os.path.join(path,item)+"/run_cpu.pl",["$Nb"])[0]
    redParticleNumber=int(redParticleNumber[:-1])
    bluePartilceNumber=int(bluePartilceNumber[:-1])
    sumParticleNumber=redParticleNumber + bluePartilceNumber
    print(f"redParticleNumber={redParticleNumber}")
    phiRed=tools.phi(data[:redParticleNumber], 60, 60, threshold=2, 
        filePath=os.path.join(path,item,"phi6.png"), order=6,
        drawPciture=False, returnSumPhi=True)
    phiBlue=tools.phi(data[redParticleNumber:], 60, 60, threshold=1.2, 
        filePath=os.path.join(path,item,"phi3.png"), order=3,
        drawPciture=False, returnSumPhi=True)
    print(f"phiRed={phiRed}, phiBlue={phiBlue}, phiSumDN={phiRed*redParticleNumber + phiBlue*(bluePartilceNumber)}")
    N.append((phiRed + phiBlue)/sumParticleNumber)
N=np.reshape(N,(9,5))
print("N=", N, "N.shape=", np.shape(N))
plt.imshow(N, cmap='hot', interpolation='nearest')
colorbar = plt.colorbar()
plt.savefig("./PS.png")


