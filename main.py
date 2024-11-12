### 数据预处理
# 1. 数据清洗
# 2. 数据填补
# 3. 数据归一化


### 样本基本趋势可视化


### 分组差异分析与可视化


### 富集与通路分析
# 1. 富集分析背景集（基因和代谢物）下载
# 2. 富集分析通路基本信息下载
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
