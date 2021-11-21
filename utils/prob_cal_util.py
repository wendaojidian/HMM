import math
from models.word_index import *


def gen_aa_list_all(aa_list):
    """
    将氨基酸字符串转换成列表
    :param aa_list:
    :return:
    """
    return list(aa_list)


def gen_seq_list_all(seq_list):
    """
    将碱基对字符串转换成列表
    :param seq_list:
    :return:
    """
    seq_list_single = []
    i = 0
    while i < len(seq_list):
        seq_list_single.append(seq_list[i:i + 3])
        i += 3
    return seq_list_single


def cal_prob(seq_an, seq_pair, base_pair_prob, word_prob):
    aa_list = gen_aa_list_all(seq_an)
    seq_list = gen_seq_list_all(seq_pair)

    prob_all = base_pair_prob[aa_list[0]][seq_list[0]]
    for i in range(1, len(aa_list)):
        an_prob = word_prob[INDEX_WORD[aa_list[i-1]], INDEX_WORD[aa_list[i]]]
        pair_prob = base_pair_prob[aa_list[i]][seq_list[i]]
        prob_all = prob_all * an_prob * pair_prob
    return -math.log(prob_all)


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
