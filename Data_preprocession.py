import pandas as pd
import numpy as np
import string
import random
import math
import time
import xlwt
import csv


def init_index():
    """
    初始化氨基酸字母对应数组下标，类型为字典，例：{A:0,B:1}
    :return:
    """
    index_word = {}
    i = 0
    for word1 in string.ascii_uppercase:
        if word1 not in ['B', 'J', 'O', 'U', 'X', 'Z']:
            index_word[word1] = i
            i += 1
    return index_word


INDEX_WORD = init_index()
FILE_90_HIGHER = "./Data/90_higher.csv"
print(INDEX_WORD)


class Seq_chance():
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


def file_process():
    """
    读取原始氨基酸列表文件，生成一阶和二阶概率
    :return: 一阶概率、二阶概率
    """
    sheet = pd.read_csv(FILE_90_HIGHER)
    seq_list = list(sheet['aaseq'])

    seq_all = {}
    for word in string.ascii_uppercase:
        if word not in ['B', 'J', 'O', 'U', 'X', 'Z']:
            seq_all[word] = []

    for seq in seq_list:
        seq_all[seq[0]].append(seq)
    print(seq_all)

    Seq_tmp = Seq_chance('All', seq_list)

    # np.savetxt('prob_all.csv', Seq_tmp.word_prob, delimiter=',')
    # print(Seq_tmp.word_prob)

    # 按初始氨基酸分类计算
    # seq_prob_word={}
    # for begin_word in list(seq_all.keys()):
    #     Seq_tmp=Seq_chance(begin_word,seq_all[begin_word])
    #     seq_prob_word[begin_word]=Seq_tmp.word_prob
    # print(INDEX_WORD)
    # print(seq_prob_word)
    return Seq_tmp.word_prob, Seq_tmp.word_prob_2nd


def get_base_pair():
    """
    读取氨基酸对应各种碱基对的概率，保存到字典中
    :return: 各种氨基酸对应各种碱基对概率的字典
    """
    base_pair_sheet = pd.read_excel("./Data/base_pair_prob.xlsx")
    base_pair = {}
    tmp_base = None
    for i in range(len(list(base_pair_sheet['x']))):
        if not pd.isnull(list(base_pair_sheet['x'])[i]):
            tmp_base = list(base_pair_sheet['x'])[i]
            base_pair[list(base_pair_sheet['x'])[i]] = {}
            base_pair[list(base_pair_sheet['x'])[i]][list(base_pair_sheet['y'])[i]] = list(base_pair_sheet['z'])[i]
        else:
            base_pair[tmp_base][list(base_pair_sheet['y'])[i]] = list(base_pair_sheet['z'])[i]
    print(base_pair)
    return base_pair


def get_key(dct, value):
    """
    通过概率找碱基对
    :param dct: 某种氨基酸-碱基对字典
    :param value: 概率
    :return:
    """
    return [k for (k, v) in dct.items() if v == value]


def get_random_pair(base_pair_one):
    """
    按照概率随机选择碱基对
    :param base_pair_one:氨基酸-碱基对字典
    :return: 随机选择的碱基对及其概率
    """
    start_prob = 0
    index = 0
    prob_values = list(base_pair_one.values())
    randnum = random.random()
    for index, prob in enumerate(prob_values):
        start_prob += prob
        if randnum <= start_prob:
            break
    return list(base_pair_one.keys())[index], list(base_pair_one.values())[index]


def get_random_an(an_prob):
    """
    按照一阶概率图随机选择下一个氨基酸
    :param an_prob: 氨基酸的概率图
    :return: 随机选择的氨基酸及其概率
    """
    start_prob = 0
    index = 0
    randnum = random.random()

    for index, prob in enumerate(an_prob):
        start_prob += prob
        if randnum <= start_prob:
            break
    return get_key(INDEX_WORD, index)[0], an_prob[index]


def cal_prob_2(seq_an, seq_pair, base_pair_prob, word_prob_2nd):
    """
    计算二阶概率
    :param seq_an:氨基酸序列
    :param seq_pair: 碱基对序列
    :param base_pair_prob: 氨基酸-碱基对概率字典
    :param word_prob_2nd: 二阶氨基酸概率图
    :return:
    """
    # print(seq_an,seq_pair)
    prob_all = 1

    # 先计算前两个氨基酸对应碱基对的概率积
    an_st, an_nd = seq_an[0], seq_an[1]
    pair_st, pair_nd = seq_pair[0:3], seq_pair[3:6]

    prob_all = prob_all * base_pair_prob[an_st][pair_st] * base_pair_prob[an_nd][pair_nd]

    # 循环计算下一个氨基酸概率、氨基酸对应碱基对概率，累乘
    for i in range(2, 10):
        before_an = seq_an[i - 2:i]
        now_an = seq_an[i]

        now_pair = seq_pair[3 * i:3 * i + 3]

        prob_an_tmp = word_prob_2nd[before_an][INDEX_WORD[now_an]]
        prob_pair_tmp = base_pair_prob[now_an][now_pair]
        # print(before_an,now_an,now_pair,prob_an_tmp,prob_pair_tmp)

        prob_all = prob_all * prob_an_tmp * prob_pair_tmp
    # print(prob_all)
    # 如果在二阶图中的概率无限接近0，返回0
    if prob_all == 0.0:
        return 0
    else:
        return -math.log(prob_all)


def produce_seq(epochs, word_prob, base_pair_prob, word_prob_2nd):
    """
    按概率随机生成序列
    :param epochs: 迭代的次数
    :param word_prob: 氨基酸一阶概率图
    :param base_pair_prob: 各种氨基酸-碱基对概率字典
    :param word_prob_2nd: 氨基酸二阶概率图
    :return: None，直接写入到csv文件
    """
    time_start_tmp = time.time()
    with open("prob_concat_2.csv", 'w', newline='') as f:
        f_csv = csv.writer(f)
        # df=pd.DataFrame(columns=('aa_seq','gene_seq','rate_ln','rate_2nd'))
        for i in range(epochs):
            # 生成第一个：
            seq_an = get_key(INDEX_WORD, int(20 * random.random()))[0]
            seq_pair, pair_prob = get_random_pair(base_pair_prob[seq_an])

            before_an = seq_an

            # 循环14次，生成后十四个氨基酸及其对应的碱基对
            for j in range(9):
                an_tmp, an_prob_tmp = get_random_an(list(word_prob[INDEX_WORD[before_an]]))
                pair_tmp, pair_prob_tmp = get_random_pair(base_pair_prob[an_tmp])

                seq_an = seq_an + an_tmp
                seq_pair = seq_pair + pair_tmp

                pair_prob = pair_prob * an_prob_tmp * pair_prob_tmp

                before_an = an_tmp
            prob_2nd = cal_prob_2(seq_an, seq_pair, base_pair_prob, word_prob_2nd)
            f_csv.writerow([seq_an, seq_pair, -math.log(pair_prob), prob_2nd])
            # df=df.append({'aa_seq':seq_an,'gene_seq':seq_pair,'rate_ln':-math.log(pair_prob),'rate_2nd':prob_2nd},ignore_index=True)

            if i % (epochs / 100) == 0:
                time_end_tmp = time.time()
                print(i / epochs, ":", time_end_tmp - time_start_tmp)
                time_start_tmp = time.time()

        # df.to_csv('prob_concat.csv', index=True, header=False)


if __name__ == '__main__':
    # 计算氨基酸一阶和二阶概率
    word_prob, word_prob_2nd = file_process()
    # 计算氨基酸-碱基对概率洗点
    base_pair = get_base_pair()
    # 统计运行时间
    time_start = time.time()
    # 生成氨基酸及碱基对序列
    produce_seq(10000, word_prob, base_pair, word_prob_2nd)
    time_end = time.time()
    print('totally cost', time_end - time_start)
