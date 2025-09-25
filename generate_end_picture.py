import os
from tqdm import tqdm
import pandas as pd
from matplotlib import pyplot as plt

def plot_last_frame():
    # 关键词
    keyword1 = "proportion"  # 读取文件夹关键词
    keyword2 = "conf"  # 位置坐标关键词
    size = 40  # 绘图点大小

    # 获取当前目录
    current_dir = os.getcwd()
    
    # 列出当前目录下包含关键词的所有文件夹
    matching_folders = [folder for folder in os.listdir(current_dir) if
                        os.path.isdir(os.path.join(current_dir, folder)) and keyword1 in folder]

    for folder in tqdm(matching_folders):
        next_dir_path = os.path.join(current_dir, folder)  # 数据文件所在目录
        parent_dir = os.path.dirname(next_dir_path)  # 图片应存放的目录（与数据文件夹同级）

        # 读取存储坐标数据的文件
        files = os.listdir(next_dir_path)
        dat_files = [file for file in files if file.endswith('.dat') and keyword2 in file]
        dat_files = sorted(dat_files, key=lambda x: int(x.split('_')[1].split('.')[0]))
        
        # 如果没有数据文件，跳过
        if not dat_files:
            continue
        
        # 读取粒子数目的内容
        initial_parameter_path = 'input.dat'  # 初始信息存储的文件名
        file_path1 = os.path.join(next_dir_path, initial_parameter_path)
        initial_parameter = pd.read_csv(file_path1)  # initial_parameter 为所有参数
        Na = initial_parameter.values[4]
        Nb = initial_parameter.values[5]
        N = Na + Nb

        # 读取最后一个数据文件
        last_file = dat_files[-1]  # 取最后一个 .dat 文件
        temp_data = pd.read_csv(os.path.join(next_dir_path, last_file), sep=' ', header=None).values

        # 过滤超出范围 [0, 60] 的粒子
        r_a = temp_data[0:int(Na[0]), :]
        r_b = temp_data[int(Na[0]):, :]

        # 过滤 r_a 和 r_b 中超出范围的粒子
        r_a = r_a[(r_a[:, 0] >= 0) & (r_a[:, 0] <= 60) & (r_a[:, 1] >= 0) & (r_a[:, 1] <= 60)]
        r_b = r_b[(r_b[:, 0] >= 0) & (r_b[:, 0] <= 60) & (r_b[:, 1] >= 0) & (r_b[:, 1] <= 60)]

        # 绘图
        plt.figure(figsize=(10, 10))

        # 动态设置坐标轴范围
        x_min, x_max = min(r_a[:, 0].min(), r_b[:, 0].min()), max(r_a[:, 0].max(), r_b[:, 0].max())
        y_min, y_max = min(r_a[:, 1].min(), r_b[:, 1].min()), max(r_a[:, 1].max(), r_b[:, 1].max())
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

        plt.scatter(r_a[:, 0], r_a[:, 1], label='a', c='red', s=size)
        plt.scatter(r_b[:, 0], r_b[:, 1], label='b', c='blue', s=size)
        plt.axis('off')  # 去除坐标轴
        plt.tight_layout(pad=0)  # 自动调整布局以填充整个画布

        # 保存图片到数据文件夹的同级目录
        output_img_path = os.path.join(parent_dir, folder + "_last_frame.png")
        plt.savefig(output_img_path, bbox_inches='tight', pad_inches=0, dpi=600)
        plt.close()

if __name__ == "__main__":
    plot_last_frame()
