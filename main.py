### 数据预处理
# 1. 数据清洗
# 2. 数据填补
# 3. 数据归一化
## 注意数据质控、数据清洗、数据填补都需要对业务有一定的了解
## 数据归一化代谢大多根据QC归一化 转录组则有较多方式  FPKM TPM需要去了解
## 核心问题，当我们做了这些工作后 样本与样本之间是否可比？指标与指标之间是否可比？

### 样本基本趋势可视化
# 热图
# PCA图


### 分组差异分析与可视化
# 差异分析（两组 多组） t wilcox DEseq2 limma等等（需要学习底层原理）
# 差异代谢物 gene拿到之后可视化有火山图、热图、相关性热图、箱型图等
# 接着就是差异代谢物的生物学意义探索富集分析

### 富集分析
    ## ORA GESA GRSA
    ## https://github.com/Asa12138/ReporterScore
# 1. 富集分析背景集（基因和代谢物）下载
    ## R见这里 https://mp.weixin.qq.com/s/NjlRMqEmD5OqMrYsmKy6Iw
    ## 实际上是通过这个网页https://rest.kegg.jp/get/hsa05230/kgml/去处理的
    ## https://rest.kegg.jp/get/hsa05230 这个网页较为麻烦需要正则去处理了（也不好处理）





# 2. 富集分析通路下载（后续pathview画图，染色）
import requests
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt

# 下载并解析KGML文件
def get_kgml(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}/kgml"
    response = requests.get(url)
    return response.text

# 解析KGML文件中的节点坐标信息
def parse_kgml(kgml_content, target_ids):
    root = ET.fromstring(kgml_content)
    coords = {}
    for entry in root.findall(".//entry"):
        entry_id = entry.get('name')
        if entry_id and any(target in entry_id for target in target_ids):
            graphics = entry.find(".//graphics")
            if graphics is not None:
                x = float(graphics.get('x'))
                y = float(graphics.get('y'))
                coords[entry_id] = (x, y)
    return coords

# 下载通路图像
def get_image(pathway_id):
    url = f"https://rest.kegg.jp/get/{pathway_id}/image"
    response = requests.get(url)
    img_file = f"{pathway_id}.png"
    with open(img_file, 'wb') as f:
        f.write(response.content)
    return img_file

# 在图像上绘制标记
def draw_markers(img_file, coords, regulation_status):
    img = Image.open(img_file)
    draw = ImageDraw.Draw(img)
    
    # 设定标记的颜色
    color_map = {'upregulated': 'red', 'downregulated': 'blue'}
    
    for id, (x, y) in coords.items():
        base_id = id.replace('cpd:', '')  
        print(f"Processing ID: {id}, Base ID: {base_id}, Coordinates: ({x}, {y})")  # 调试信息  
       
        # 确保状态存在，默认为下调  
        status = regulation_status.get(base_id , 'downregulated')  # 默认状态为下调  
        color = color_map.get(status, 'blue')  
        draw.ellipse((x - 5, y - 5, x + 5, y + 5), outline='black', fill=color)  
    
    # modified_img_file = "modified_image.png"  
    # img.save(modified_img_file, dpi=(600, 600))
    #     draw.ellipse((x - 5, y - 5, x + 5, y + 5), outline='black', fill=color)
    
    modified_img_file = "modified_image.png"
    img.save(modified_img_file, dpi=(600, 600))  # 保存图像，设置dpi为600
    return modified_img_file

# 无法处理基因 可以创建一个矩型的辅助函数
# 示例：获取并修改mmu00010的图像，标记上调和下调的代谢物
pathway_id = "mmu00010"
target_ids = ["C01172", "C00022", "C00024"]  # 多个代谢物
regulation_status = {
    "C01172": "upregulated",  # 上调
    "C00022": "downregulated",  # 下调
    "C00024": "upregulated" # 上调
   
}

# 获取KGML内容并解析坐标
kgml_content = get_kgml(pathway_id)
coords = parse_kgml(kgml_content, target_ids)
print("Coordinates:", coords)  # 打印坐标

# 获取并修改图像
img_file = get_image(pathway_id)
modified_img_file = draw_markers(img_file, coords, regulation_status)

# 显示修改后的图像
modified_img = Image.open(modified_img_file)
plt.imshow(modified_img)
plt.axis('off')
plt.show()

# 3. 富集分析通路分类信息下载
    ## 主要借鉴：https://blog.csdn.net/qq_38774801/article/details/125981487
    ## 当然也可以爬虫保存成网页处理，自行完善吧
import re
import pandas as pd

# 读取文件
with open("a.txt", "r", encoding="utf-8") as file:
    lines = file.readlines()
print(len(lines))
# 提取包含有效信息的行
filtered_lines = [line for line in lines if re.search(r">([0-9])", line)]
print(len(filtered_lines))
# 初始化注释变量和最终存储列表
big_annotation = ""
small_annotation = ""
final_lines = []

# 遍历筛选后的行，提取信息
for line in filtered_lines:
    # 判断并更新大类注释
    if re.search(r">\d\. ", line):
        try:
            big_annotation = re.search(r">.*?<", line).group(0)
            big_annotation = re.split(r"\. ", big_annotation)[1].split("<")[0]
        except IndexError:
            continue  # 忽略可能的分割错误行

    # 判断并更新小类注释
    elif re.search(r">\d\.\d{1,2} ", line):
        try:
            small_annotation = re.search(r">.*?<", line).group(0)
            small_annotation = re.sub(r"^>\d\.\d{1,2} ", "", small_annotation).split("<")[0]
        except IndexError:
            continue

    # 提取路径信息
    elif re.search(r">\d{5}", line):
        try:
            # 提取路径 ID
            element1 = re.search(r">\d{5}", line).group(0)[1:]
            
            if "hsa+pathogen" not in line:
                # 提取 Pathway Identifier
                element2 = re.search(r"pathway/[a-zA-Z]{2,4}\d{5}", line).group(0).split("/")[1]
                
                # 提取 Pathway 名称
                element3_match = re.search(r'<a.*?>(.*?)</a>', line)
                element3 = element3_match.group(1) if element3_match else ""
            else:
                element2 = "organism:hsa+pathogen"
                element3_match = re.search(r'hsa\+pathogen.*?<', line)
                element3 = element3_match.group(0).split(">")[1].split("<")[0] if element3_match else ""

            # 组合信息并添加到最终列表
            final_line = [element1, element2, element3, big_annotation, small_annotation]
            final_lines.append(final_line)
        except AttributeError:
            continue  # 忽略无匹配项的行

# 创建 DataFrame 并保存为 Excel
df = pd.DataFrame(final_lines, columns=["ID", "Pathway Identifier", "Pathway", "Big Annotation", "Small Annotation"])
df.to_excel("kegg.xlsx", index=False)
