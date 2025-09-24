import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import cv2
import sys

displacementX=5
displacementY=-10
blankAxis=0

path0 = './'

path1 = sys.argv[1]
path = path0 + path1

print(path)

if not os.path.exists(path+'/temp'):
    os.mkdir(path+'/temp')

with open(path + '/input.dat', "r", encoding='utf-8') as f:  #打开文本
    data = f.readlines()

Na = int(data[5])
Nb = int(data[6])
Lx = int(data[8])
Ly = int(data[9])
tStop = int(data[13])
tExport = float(data[15])
Nt = int(tStop/tExport)
text_split = path1[:len(path1)//2] + '\n' + path1[len(path1)//2:]

for i in range(Nt+1):
    data_i = np.loadtxt(path + '/conf_'+str(i)+'.dat')
    data_i[:,0] = (data_i[:,0]+displacementX)%Lx
    data_i[:,1] = (data_i[:,1]+displacementY)%Ly
    plt.figure(figsize=(10, 10))
    plt.scatter(data_i[0:Na-1,0],data_i[0:Na-1,1],s=40 ,c='r',marker='o')
    plt.scatter(data_i[Na:Na+Nb-1,0],data_i[Na:Na+Nb-1,1], s=40, c='b',marker='o')
    plt.xlim(0,Lx)
    plt.ylim(0,Ly)
    if blankAxis==1:
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlabel('x',fontsize=30)
        plt.ylabel('y',fontsize=30)
        title_fontdict = {
            'fontsize': 30,  # 字体大小
            'fontweight': 'bold'  # 字体粗细
            }
        plt.title("proportion={:.2f} T={}".format(Na/Nb, i*tStop), fontdict=title_fontdict)
    if blankAxis==0:
        plt.axis('off')
        plt.xticks([])
        plt.yticks([])
    plt.savefig(path + '/temp/fig_'+str(i)+'.png')
    plt.close() 

video_name = path.rstrip('/') +'.avi'  # 输出视频文件名
print(video_name)

# 获取第一张图片的尺寸
frame = cv2.imread(path+'/temp/fig_0.png')
height, width, layers = frame.shape

# 创建视频编写器
fourcc = cv2.VideoWriter_fourcc(*'XVID')
video = cv2.VideoWriter(video_name, fourcc, 10, (width, height))  # 帧率为 30

# 逐帧写入视频
for i in range(Nt+1):
    video.write(cv2.imread(path+'/temp/fig_'+str(i)+'.png'))

# 释放视频编写器
video.release()

os.system("vlc " + video_name)
