import subprocess
import os

# 获取当前文件夹路径
current_folder = '/home/undergraduate/wu_ben_chang/conf-data-all'

# 获取当前文件夹下的子文件和子文件夹的名称列表
item_names = os.listdir(current_folder)

# 存储未被命名的文件夹名称的列表
unnamed_folders = []

# 筛选未被命名的文件夹
for item_name in item_names:
    item_path = os.path.join(current_folder, item_name)
    if os.path.isdir(item_path):
        folder_name = item_name
        file_name = folder_name + ".avi"
        file_path = os.path.join(current_folder, file_name)
        if not os.path.exists(file_path):
            unnamed_folders.append(folder_name)

for item_name in unnamed_folders:
    subprocess.run("python ./plot.py "+item_name, shell=True, check=True)


