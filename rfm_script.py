import os
config = snakemake.config
matrix = f"data/{config['run_id']}_combined_matrix.txt"
output = config["run_id"]
one_count = 0
all_count = -1
with open(matrix, 'r') as inp:
    for line in inp:
        all_count += 1
        if line[0] == '1':
            one_count += 1

remove_inds = snakemake.config['model_params']['columns_to_remove']
label_ind = snakemake.config['model_params']['class_label_index']
num_trees = snakemake.config['model_params']['num_trees']
max_depth = "-m " + str(snakemake.config['model_params']['max_depth']) if snakemake.config['model_params']['max_depth'] else ""
min_samples_split = "-s " + str(snakemake.config['model_params']['min_samples_split']) if snakemake.config['model_params']['min_samples_split'] else ""

cmd = f"python3 src/rf_model.py -i {matrix} -d {remove_inds} -t {label_ind} -n {num_trees} {max_depth} {min_samples_split} -c {one_count} -l {all_count} -o {output}"
os.system(cmd)
