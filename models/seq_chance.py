import numpy as np
import pandas as pd

from models.word_index import *


class Seq_chance:
    """
    主要用于基于给定氨基酸序列统计概率图
    """

    def __init__(self, begin_word, seq_list):
        # 开头字母，不考虑时无用
        self.begin = begin_word
        # 氨基酸序列列表，包含多条氨基酸序列
        self.seq_list = seq_list
        # 初始化氨基酸概率图为20*20二维数组，下标通过宏变量INDEX_WORD查找
        self.word_prob = np.zeros((20, 20), dtype=int)
        # 计算一阶概率
        self.cal_prob()
        # 计算二阶概率
        self.word_prob_2nd = self.cal_prob_2nd()

    def cal_prob(self):
        for seq in self.seq_list:
            for i in range(len(seq) - 1):
                self.word_prob[INDEX_WORD[seq[i]]][INDEX_WORD[seq[i + 1]]] += 1
        self.word_prob = self.word_prob / np.sum(self.word_prob, axis=1).reshape((20, 1))

    def cal_prob_2nd(self):
        word_prob_2nd = {}
        for seq in self.seq_list:
            for i in range(len(seq) - 2):
                words_2nd = seq[i:i + 2]
                if words_2nd not in list(word_prob_2nd.keys()):
                    word_prob_2nd[words_2nd] = np.zeros(20, dtype=int)
                    word_prob_2nd[words_2nd][INDEX_WORD[seq[i + 2]]] += 1
                else:
                    word_prob_2nd[words_2nd][INDEX_WORD[seq[i + 2]]] += 1
        for seq in list(word_prob_2nd.keys()):
            word_prob_2nd[seq] = list(word_prob_2nd[seq] / np.sum(word_prob_2nd[seq]))

        df = pd.DataFrame(word_prob_2nd).T
        df.to_excel('prob_2nd.xlsx', index=True, header=False)
        return word_prob_2nd
