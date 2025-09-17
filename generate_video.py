import matplotlib
matplotlib.use('Agg')
import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import cv2
import shutil



def delete_folder_contents(folder_path):
    # 确保目录存在
    if os.path.exists(folder_path):
        # 遍历目录中的所有文件和子目录
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            # 如果是文件，则直接删除
            if os.path.isfile(item_path):
                os.remove(item_path)
            # 如果是目录，则递归删除目录中的内容
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)

    # 删除空目录
    if os.path.exists(folder_path):
        os.rmdir(folder_path)


# 指定要搜索的关键词
keyword1 = "proportion"  # 读取文件夹关键词
keyword2 = "conf"  # 位置坐标关键词
size = 8  # 绘图点大小

# 获取当前目录
current_dir = os.getcwd()

# 列出当前目录下包含关键词的所有文件夹
matching_folders = [folder for folder in os.listdir(current_dir) if
                    os.path.isdir(os.path.join(current_dir, folder)) and keyword1 in folder]

for folder in tqdm(matching_folders):
    next_dir_path = os.path.join(current_dir, folder)  # 下一级目录位置

    # 读取存储坐标数据的文件
    files = os.listdir(next_dir_path)
    dat_files = [file for file in files if file.endswith('.dat') and keyword2 in file]
    dat_files = sorted(dat_files, key=lambda x: int(x.split('_')[1].split('.')[0]))

    # 读取粒子数目的内容
    initial_parameter_path = 'input.dat'  # 初始信息存储的文件名
    file_path1 = os.path.join(next_dir_path, initial_parameter_path)
    initial_parameter = pd.read_csv(file_path1)  # initial_parameter 为所有参数
    Na = initial_parameter.values[4]
    Nb = initial_parameter.values[5]
    N = Na + Nb
    Lx = initial_parameter.values[7]
    Ly = initial_parameter.values[8]

    # 新建绘图的文件夹
    temp_path0 = os.path.join(next_dir_path, 'temp_png')
    os.makedirs(temp_path0, exist_ok=True)
    # 新建写入视频文件夹
    temp_path1 = os.path.join(current_dir, 'result_folder','video')
    os.makedirs(temp_path1, exist_ok=True)

    for ii in tqdm(range(0, len(dat_files))):
        temp_data = pd.read_csv(os.path.join(next_dir_path, dat_files[ii]), sep=' ', header=None)
        temp_data = temp_data.values
        r_a = temp_data[0:int(Na[0]), :]
        r_b = temp_data[int(Na[0]):, :]

        # 绘图
        plt.figure(figsize=(8, 8))
        plt.xlim(0, 100)
        plt.ylim(0, 100)
        plt.scatter(r_a[:, 0], r_a[:, 1], label='a', c='red', s=size)
        plt.scatter(r_b[:, 0], r_b[:, 1], label='b', c='blue', s=size)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.title("N:"+str(int(N[0])) + '_' + "Na:"+str(int(Na[0])) + '_' +"Nb:"+ str(int(Nb[0])) +"_Na/Nb:"+ "{:.2f}".format(float((Na / Nb)[0])),fontsize=30)
        plt.axis("equal")
        plt.legend(loc="upper right")
        plt.savefig(os.path.join(temp_path0, str(ii) + ".png"))
        plt.close()

    video_name_no_path = [str(int(N[0])), str(int(Na[0])), str(int(Nb[0])), "{:.2f}".format(float((Na / Nb)[0]))]
    video_name = os.path.join(temp_path1, "_".join(video_name_no_path) + ".avi")
    # 获取图片文件夹中的所有图片文件名并排序
    images = [img for img in os.listdir(temp_path0) if img.endswith(".png")]
    images = sorted(images, key=lambda x: int(x.split('.')[0]))

    # 获取第一张图片的宽度和高度
    frame = cv2.imread(os.path.join(temp_path0, images[0]))
    height, width, layers = frame.shape

    # 使用 VideoWriter 对象创建一个视频文件
    video = cv2.VideoWriter(video_name, cv2.VideoWriter_fourcc(*'DIVX'), 30, (width, height))

    # 将图片逐一写入视频文件
    for image in images:
        video.write(cv2.imread(os.path.join(temp_path0, image)))

    # 释放 VideoWriter 对象并销毁所有窗口
    cv2.destroyAllWindows()
    video.release()

    delete_folder_contents(temp_path0)
