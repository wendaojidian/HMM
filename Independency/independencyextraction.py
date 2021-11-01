import xlrd
from xlwt import Workbook


def load_dict(wk_bk, sheet_name):
    wk_dict = {}
    wk_all = wk_bk.sheet_by_name(sheet_name)

    for i in range(wk_all.nrows):
        list_value = wk_all.row_values(i)
        wk_dict[list_value[4]] = list_value
    return wk_dict, set(wk_all.col_values(4))


def write_to_xlsx(file_path, sheet_names, wk_uniques):
    ws = Workbook(encoding='UTF-8')
    for a in range(len(sheet_names)):
        sheet_name = sheet_names[a]
        w = ws.add_sheet(sheet_name)
        for i in range(len(wk_uniques[a])):
            for j in range(len(wk_uniques[a][i])):
                w.write(i, j, wk_uniques[a][i][j])
    ws.save(file_path)


if __name__ == '__main__':
    work_book = xlrd.open_workbook("../Data/All Filtered 2.0 Data.xlsx")
    print("step1 finished")

    wk1_dict, wk1_seqs = load_dict(work_book, 'Filtered 2.0-MG-WK-1')
    wk2_dict, wk2_seqs = load_dict(work_book, 'Filtered 2.0-MG-WK-2')
    wk3_dict, wk3_seqs = load_dict(work_book, 'Filtered 2.0-MG-WK-3')
    wk4_dict, wk4_seqs = load_dict(work_book, 'Filtered 2.0-MG-WK-4')
    print("step2 finished")

    wk1_unique = [wk1_dict[seq_unique] for seq_unique in
                  list(wk1_seqs.difference(wk2_seqs).difference(wk3_seqs).difference(wk4_seqs))]
    print("step3 finished")
    wk2_unique = [wk2_dict[seq_unique] for seq_unique in
                  list(wk2_seqs.difference(wk1_seqs).difference(wk3_seqs).difference(wk4_seqs))]
    print("step4 finished")
    wk3_unique = [wk3_dict[seq_unique] for seq_unique in
                  list(wk3_seqs.difference(wk2_seqs).difference(wk2_seqs).difference(wk4_seqs))]
    print("step5 finished")
    wk4_unique = [wk4_dict[seq_unique] for seq_unique in
                  list(wk4_seqs.difference(wk2_seqs).difference(wk3_seqs).difference(wk1_seqs))]
    print("step6 finished")

    write_to_xlsx("../Data/Unique Filtered 2.0 Data.xls",
                  [u'Filtered 2.0-MG-WK-1', u'Filtered 2.0-MG-WK-2', u'Filtered 2.0-MG-WK-3', u'Filtered 2.0-MG-WK-4'],
                  [wk1_unique, wk2_unique, wk3_unique, wk4_unique])
    print("step7 finished")
