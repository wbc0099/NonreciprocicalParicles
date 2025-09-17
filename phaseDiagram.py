import re
from PIL import Image, ImageDraw, ImageFont
import os

# 获取当前文件夹中的所有 PNG 图片
image_files = [f for f in os.listdir() if f.endswith('.png')]

# 提取 proportion 和 kba 值
def extract_values(filename):
    match = re.match(r'N\d+--proportion([\d.]+)--kba-([\d.-]+)_last_frame\.png', filename)
    if match:
        proportion = float(match.group(1))
        kba = float(match.group(2))
        return proportion, kba
    return None, None

# 过滤并排序
images_with_values = []
for image_file in image_files:
    proportion, kba = extract_values(image_file)
    if proportion is not None and kba is not None:
        images_with_values.append((proportion, kba, image_file))

# 按 kba (x轴) 和 proportion (y轴) 排序
images_with_values.sort(key=lambda x: (x[1], -x[0]))  # 先按 kba，再按 proportion 降序（确保 0 刻度在最下面）

# 整合图片
# 假设每张图片的大小相同
sample_image = Image.open(images_with_values[0][2])
img_width, img_height = sample_image.size

# 计算大图的大小
unique_kba = sorted(set(kba for _, kba, _ in images_with_values))  # 唯一的 kba 值
unique_proportion = sorted(set(proportion for proportion, _, _ in images_with_values), reverse=True)  # 唯一的 proportion 值（从大到小）
cols = len(unique_kba)  # 列数由 kba 的唯一值决定
rows = len(unique_proportion)  # 行数由 proportion 的唯一值决定

# 设置边框参数
border_width = 5  # 每张图之间的边框宽度

# 创建大图
big_image_width = cols * img_width + (cols - 1) * border_width + 100  # 留出左边距和右边距
big_image_height = rows * img_height + (rows - 1) * border_width + 100  # 留出上边距和下边距
big_image = Image.new('RGB', (big_image_width, big_image_height), color=(255, 255, 255))
draw = ImageDraw.Draw(big_image)

# 设置字体（需要指定字体文件路径）
try:
    font = ImageFont.truetype("arial.ttf", 36)  # 使用 Arial 字体，字号调大
except IOError:
    font = ImageFont.load_default()  # 如果找不到字体，使用默认字体

# 将图片粘贴到大图上
for proportion, kba, image_file in images_with_values:
    img = Image.open(image_file)
    # 计算图片的位置
    col = unique_kba.index(kba)  # 根据 kba 确定列
    row = unique_proportion.index(proportion)  # 根据 proportion 确定行
    x_offset = col * (img_width + border_width) + 50  # 留出左边距
    y_offset = row * (img_height + border_width) + 50  # 留出上边距
    big_image.paste(img, (x_offset, y_offset))

    # 绘制每张图的边框
    draw.rectangle(
        [(x_offset - border_width, y_offset - border_width),
         (x_offset + img_width + border_width, y_offset + img_height + border_width)],
        outline="black", width=border_width
    )

# 保存大图
big_image.save('combined_image_with_borders.png')
print("整合完成，图片已保存为 combined_image_with_borders.png")