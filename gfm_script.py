import os
control = snakemake.config["coordinates"]["control"]
disease = snakemake.config["coordinates"]["disease"]
control_matrix = f"{snakemake.config['output_matrix']['control']}_{snakemake.config['run_id']}"
disease_matrix = f"{snakemake.config['output_matrix']['disease']}_{snakemake.config['run_id']}"
coord_files = "-g " + " ".join(snakemake.config['matrix_params']['overlap_coordinate_files']) if snakemake.config['matrix_params']['overlap_coordinate_files'] else ""
normalized = "-z" if snakemake.config['matrix_params']['normalized'] else ""
variant_type = snakemake.config['matrix_params']['variant_type']
ref_genome = snakemake.config['matrix_params']['ref_genome']
randomized_num = snakemake.config['matrix_params']['randomized_num']
bigwig_files = "-b " + " ".join(snakemake.config['matrix_params']['bigwig_files']) if snakemake.config['matrix_params']['bigwig_files'] else ""
os.system(f"python3 src/generate_feature_matrix.py -c {control} {bigwig_files} {coord_files} -o {control_matrix} -t {variant_type} {normalized} -r {randomized_num} -rg {ref_genome}")
os.system(f"python3 src/generate_feature_matrix.py -c {disease} {bigwig_files} {coord_files} -o {disease_matrix} -t {variant_type} {normalized} -r {randomized_num} -rg {ref_genome} -f")
combined_matrix = f"data/{snakemake.config['run_id']}_combined_matrix.txt"
control_matrix = f"{snakemake.config['output_matrix']['control']}_{snakemake.config['run_id']}_normalized.txt" if snakemake.config['matrix_params']['normalized'] else f"{snakemake.config['output_matrix']['control']}_{snakemake.config['run_id']}.txt"
disease_matrix = f"{snakemake.config['output_matrix']['disease']}_{snakemake.config['run_id']}_normalized.txt" if snakemake.config['matrix_params']['normalized'] else f"{snakemake.config['output_matrix']['disease']}_{snakemake.config['run_id']}.txt"
os.system(f"cat {disease_matrix} > {combined_matrix} && tail -n +2 {control_matrix} >> {combined_matrix}")
