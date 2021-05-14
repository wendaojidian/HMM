from configparser import ConfigParser

from models.prob_aa_bp import get_base_pair
from utils.prob_cal_by_seq_exists import get_list_from_xlsx, prob_cal_by_list
from utils.prob_cal_random_gennerate import produce_seq
from utils.prob_generate import prob_generate


def main(i):
    """

    :param i: i=1:表示依据既有的序列计算分数，i=2:表示随机生成序列并计算分数。
    :return:
    """
    cp = ConfigParser()
    cp.read("config.cfg")
    prob_between_aa_file = cp.get('file_absolute_path', 'prob_between_aa_file')
    prob_aa_pair_file = cp.get('file_absolute_path', 'prob_aa_pair_file')
    list_source_file = cp.get('file_absolute_path', 'list_source_file')
    aa_line_name = cp.get('file_absolute_path', 'aa_line_name')
    seq_line_name = cp.get('file_absolute_path', 'seq_line_name')
    save_file_path_1 = cp.get('file_absolute_path', 'save_file_path_1')
    save_file_path_2 = cp.get('file_absolute_path', 'save_file_path_2')
    epochs = int(cp.get('file_absolute_path', 'epochs'))

    if i == 1:
        prob_aa_pair = get_base_pair(prob_aa_pair_file)
        word_prob, word_prob_2 = prob_generate(prob_between_aa_file)
        aa_list, seq_list = get_list_from_xlsx(list_source_file, aa_line_name, seq_line_name)
        prob_cal_by_list(save_file_path_1, aa_list, seq_list, prob_aa_pair, word_prob, word_prob_2)

    elif i == 2:
        prob_aa_pair = get_base_pair(prob_aa_pair_file)
        word_prob, word_prob_2 = prob_generate(prob_between_aa_file)
        produce_seq(epochs, word_prob, prob_aa_pair, word_prob_2,
                    save_file_path_2)


if __name__ == '__main__':
    main(1)
