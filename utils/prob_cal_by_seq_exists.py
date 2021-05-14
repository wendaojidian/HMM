import csv
import pandas as pd

from utils.prob_generate import prob_generate
from models.prob_aa_bp import *
from utils.prob_cal_util import *


def get_list_from_xlsx(file_path, line_name_amino_acid='AA', line_name_base_pair='Seq'):
    """
    从xlsx文件中读取氨基酸列表和碱基对列表
    :param file_path: 文件的绝对路径
    :param line_name_amino_acid: 氨基酸的列名
    :param line_name_base_pair: 碱基对的列名
    :return: 氨基酸序列和碱基对序列
    """
    sheet_data = pd.read_excel(file_path)
    aa_list = list(sheet_data[line_name_amino_acid])
    seq_list = list(sheet_data[line_name_base_pair])
    return aa_list, seq_list


def prob_cal_by_list(save_file_path, aa_list, seq_list, base_pair, word_prob, word_prob_2nd):
    with open(save_file_path, 'w', newline='') as f:
        f_csv = csv.writer(f)
        for i in range(len(aa_list)):
            aa = aa_list[i]
            seq = seq_list[i]

            prob_1st, prob_2nd = cal_prob(aa, seq, base_pair, word_prob), cal_prob_2(aa, seq, base_pair, word_prob_2nd)
            f_csv.writerow([aa, seq, prob_1st, prob_2nd])


if __name__ == '__main__':
    # 计算氨基酸之间概率的文件
    prob_between_aa_file = "/Users/liufan/program/PYTHON/张佳宁/HMM/Data/90_higher.csv"
    # 保存结果的csv文件
    save_file_path = "/Users/liufan/program/PYTHON/张佳宁/HMM/Data/E.coli-NGS-30bp-HMM-score.csv"
    # 计算氨基酸->碱基对的概率
    base_pair = get_base_pair("/Users/liufan/program/PYTHON/张佳宁/HMM/Data/base_pair_prob.xlsx")
    print("base_pair_prob", base_pair)
    # 计算一阶、二阶概率
    word_prob, word_prob_2 = prob_generate(prob_between_aa_file)
    # 获取待计算的序列
    aa_list, seq_list = get_list_from_xlsx("/Users/liufan/program/PYTHON/张佳宁/HMM/Data/E.coli-NGS-30bp-HMM.xlsx")
    # 计算概率并保存
    prob_cal_by_list(save_file_path, aa_list, seq_list, base_pair, word_prob, word_prob_2)
