import pandas as pd


def get_base_pair(prob_file):
    """
    读取氨基酸对应各种碱基对的概率，保存到字典中
    :param prob_file: 保存碱基对概率的文件绝对路径
    :return: 各种氨基酸对应各种碱基对概率的字典
    """
    base_pair_sheet = pd.read_excel(prob_file)
    print(base_pair_sheet)
    base_pair = {}
    tmp_base = None
    for i in range(len(list(base_pair_sheet['x']))):
        if not pd.isnull(list(base_pair_sheet['x'])[i]):
            tmp_base = list(base_pair_sheet['x'])[i]
            base_pair[list(base_pair_sheet['x'])[i]] = {}
            base_pair[list(base_pair_sheet['x'])[i]][list(base_pair_sheet['y'])[i]] = list(base_pair_sheet['z'])[i]
        else:
            base_pair[tmp_base][list(base_pair_sheet['y'])[i]] = list(base_pair_sheet['z'])[i]
    return base_pair

