import string
import pandas as pd

from models.seq_chance import *


def prob_generate(prob_cal_file_source):
    """
    读取原始氨基酸列表文件，生成一阶和二阶概率
    :return: 一阶概率、二阶概率
    """
    sheet = pd.read_csv(prob_cal_file_source)
    seq_list = list(sheet['aaseq'])

    seq_all = {}
    for word in string.ascii_uppercase:
        if word not in ['B', 'J', 'O', 'U', 'X', 'Z']:
            seq_all[word] = []

    for seq in seq_list:
        seq_all[seq[0]].append(seq)

    Seq_tmp = Seq_chance('All', seq_list)
    return Seq_tmp.word_prob, Seq_tmp.word_prob_2nd


if __name__ == '__main__':
    source_file = "/Users/liufan/program/PYTHON/张佳宁/HMM/Data/90_higher.csv"
    word_prob, word_prob_2 = prob_generate(source_file)
    print(word_prob[0][0])
