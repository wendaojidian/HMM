import string


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
