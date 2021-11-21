from configparser import ConfigParser

from models.prob_aa_bp import get_base_pair

cp = ConfigParser()
cp.read("config.cfg", encoding='utf-8')
prob_between_aa_file = cp.get('file_absolute_path', 'prob_between_aa_file')
prob_aa_pair_file = cp.get('file_absolute_path', 'prob_aa_pair_file')

prob_aa_pair = get_base_pair(prob_aa_pair_file)
print(prob_aa_pair)
