import tools
import os
import numpy as np
import matplotlib.pyplot as plt

path="../conf-data/"
items=os.listdir(path)
items=[item for item in items if item.startswith("N")]
items=sorted(items,key=lambda x:
    (float(x.split("proportion")[1].split("--")[0]),
     float(x.split("kba-")[1])))
print(items)
items=sorted(items,key=lambda x:
    (float(x.split("proportion")[1].split("--")[0]),
     float(x.split("kba-")[1])))
N=[]
for item in items:
    print("item=", item)
    data=os.path.join(path,item)+"/conf_1000.dat"
    data=tools.readFile(data)
    # redParticleNumber,blueParticleNumber,xlim,ylim=tools.readParameter(
    #     os.path.join(path,item)+"/run_cpu.pl",["$Na","$Nb","$Lx","$Ly"])[0]
    # print(f"redParticleNumber={redParticleNumber}, blueParticleNumber={blueParticleNumber},xlim={xlim}, ylim={ylim}")
    redParticleNumber,blueParticleNumber,lxBox,lyBox=tools.readParameter(
        os.path.join(path,item)+"/run_cpu.pl",["$Na","$Nb","$Lx","$Ly"])
    print(f"redParticleNumber={redParticleNumber}, blueParticleNumber={blueParticleNumber}, lxBox={lxBox}, lyBox={lyBox}")
    data1=data[:int(redParticleNumber)]
    data2=data[int(redParticleNumber):]
    tools.plotParticles(data1[:,0], data1[:,1], float(lxBox), float(lyBox),
                        os.path.join(path,item,"plot.png"), x1=data2[:,0], y1=data2[:,1])
                                     
