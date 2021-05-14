import csv
import math
import random
import time

from models.prob_aa_bp import get_base_pair
from models.word_index import INDEX_WORD
from utils.prob_cal_by_seq_exists import get_list_from_xlsx, prob_cal_by_list
from utils.prob_cal_util import cal_prob_2
from utils.prob_generate import prob_generate




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


def produce_seq(epochs, word_prob, base_pair_prob, word_prob_2nd, save_file):
    """
    按概率随机生成序列
    :param epochs: 迭代的次数
    :param word_prob: 氨基酸一阶概率图
    :param base_pair_prob: 各种氨基酸-碱基对概率字典
    :param word_prob_2nd: 氨基酸二阶概率图
    :return: None，直接写入到csv文件
    """
    time_start_tmp = time.time()
    with open(save_file, 'w', newline='') as f:
        f_csv = csv.writer(f)
        for i in range(epochs):
            # 生成第一个：
            seq_an = get_key(INDEX_WORD, int(20 * random.random()))[0]
            seq_pair, pair_prob = get_random_pair(base_pair_prob[seq_an])

            before_an = seq_an

            # 循环9次，生成后九个氨基酸及其对应的碱基对
            for j in range(9):
                an_tmp, an_prob_tmp = get_random_an(list(word_prob[INDEX_WORD[before_an]]))
                pair_tmp, pair_prob_tmp = get_random_pair(base_pair_prob[an_tmp])

                seq_an = seq_an + an_tmp
                seq_pair = seq_pair + pair_tmp

                pair_prob = pair_prob * an_prob_tmp * pair_prob_tmp

                before_an = an_tmp
            prob_2nd = cal_prob_2(seq_an, seq_pair, base_pair_prob, word_prob_2nd)
            f_csv.writerow([seq_an, seq_pair, -math.log(pair_prob), prob_2nd])

            if i % (epochs / 100) == 0:
                time_end_tmp = time.time()
                print(i / epochs, ":", time_end_tmp - time_start_tmp)
                time_start_tmp = time.time()


if __name__ == '__main__':
    # 计算氨基酸之间概率的文件
    prob_between_aa_file = "/Users/liufan/program/PYTHON/张佳宁/HMM/Data/90_higher.csv"
    # 计算氨基酸->碱基对的概率
    base_pair = get_base_pair("/Users/liufan/program/PYTHON/张佳宁/HMM/Data/base_pair_prob.xlsx")
    # 计算一阶、二阶概率
    word_prob, word_prob_2 = prob_generate(prob_between_aa_file)
    # 计算概率并保存
    produce_seq(1000, word_prob, base_pair, word_prob_2, "/Users/liufan/program/PYTHON/张佳宁/HMM/Data/new_generate_2021.5.14.csv")
