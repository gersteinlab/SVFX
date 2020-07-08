
######## Snakemake header ########
import sys; sys.path.extend(['/usr/local/lib/python3.6/site-packages', '/Users/jagathvytheeswaran/SVFX']); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X#\x00\x00\x00data/GRCh37.nr_deletions.common.bedq\x06X\'\x00\x00\x00data/GRCh37.nr_deletions.pathogenic.bedq\x07e}q\x08(X\x06\x00\x00\x00_namesq\t}q\n(X\x07\x00\x00\x00controlq\x0bK\x00N\x86q\x0cX\x07\x00\x00\x00diseaseq\rK\x01N\x86q\x0euh\x0bh\x06h\rh\x07ubX\x06\x00\x00\x00outputq\x0fcsnakemake.io\nOutputFiles\nq\x10)\x81q\x11X&\x00\x00\x00data/real_data_run_combined_matrix.txtq\x12a}q\x13h\t}q\x14sbX\x06\x00\x00\x00paramsq\x15csnakemake.io\nParams\nq\x16)\x81q\x17}q\x18h\t}q\x19sbX\t\x00\x00\x00wildcardsq\x1acsnakemake.io\nWildcards\nq\x1b)\x81q\x1c}q\x1dh\t}q\x1esbX\x07\x00\x00\x00threadsq\x1fK\x01X\t\x00\x00\x00resourcesq csnakemake.io\nResources\nq!)\x81q"(K\x01K\x01e}q#(h\t}q$(X\x06\x00\x00\x00_coresq%K\x00N\x86q&X\x06\x00\x00\x00_nodesq\'K\x01N\x86q(uh%K\x01h\'K\x01ubX\x03\x00\x00\x00logq)csnakemake.io\nLog\nq*)\x81q+}q,h\t}q-sbX\x06\x00\x00\x00configq.}q/(X\x06\x00\x00\x00run_idq0X\r\x00\x00\x00real_data_runq1X\x0b\x00\x00\x00coordinatesq2}q3(X\x07\x00\x00\x00controlq4X#\x00\x00\x00data/GRCh37.nr_deletions.common.bedq5X\x07\x00\x00\x00diseaseq6X\'\x00\x00\x00data/GRCh37.nr_deletions.pathogenic.bedq7uX\r\x00\x00\x00matrix_paramsq8}q9(X\x0c\x00\x00\x00bigwig_filesq:NX\x18\x00\x00\x00overlap_coordinate_filesq;]q<(X\x18\x00\x00\x00data/gc19_pc.prom.nr.bedq=X\x18\x00\x00\x00data/gc19_pc.3utr.nr.bedq>X\x18\x00\x00\x00data/gc19_pc.5utr.nr.bedq?X\x17\x00\x00\x00data/gc19_pc.cds.nr.bedq@X\x16\x00\x00\x00data/gc19_pc.ss.nr.bedqAX\x15\x00\x00\x00data/sensitive.nc.bedqBX\x1d\x00\x00\x00data/ultra.conserved.hg19.bedqCX/\x00\x00\x00data/wgEncodeBroadHmmGm12878HMM.Heterochrom.bedqDX"\x00\x00\x00data/H1-ESC_Dixon2015-raw_TADs.bedqEeX\n\x00\x00\x00normalizedqF\x88X\x0c\x00\x00\x00variant_typeqGX\x03\x00\x00\x00DELqHX\n\x00\x00\x00ref_genomeqIX\x10\x00\x00\x00data/hg19.fa.faiqJX\x0e\x00\x00\x00randomized_numqKK\x01uX\r\x00\x00\x00output_matrixqL}qM(X\x07\x00\x00\x00controlqNX\x13\x00\x00\x00data/control_matrixqOX\x07\x00\x00\x00diseaseqPX\x13\x00\x00\x00data/disease_matrixqQuX\x0c\x00\x00\x00model_paramsqR}qS(X\x11\x00\x00\x00columns_to_removeqTX\x0f\x00\x00\x00remove_inds.txtqUX\x11\x00\x00\x00class_label_indexqVK\x00X\t\x00\x00\x00num_treesqWK\x96X\t\x00\x00\x00max_depthqXKdX\x11\x00\x00\x00min_samples_splitqYK\x02uuX\x04\x00\x00\x00ruleqZX\x17\x00\x00\x00generate_feature_matrixq[X\x0f\x00\x00\x00bench_iterationq\\NX\t\x00\x00\x00scriptdirq]X\x1e\x00\x00\x00/Users/jagathvytheeswaran/SVFXq^ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/jagathvytheeswaran/SVFX/gfm_script.py';
######## Original script #########
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
control_matrix = f"{snakemake.config['output_matrix']['control']}_{snakemake.config['run_id']}_normalized.txt" if snakemake.config['matrix_params']['normalized'] else "{snakemake.config['output_matrix']['control']}_{snakemake.config['run_id']}.txt"
disease_matrix = f"{snakemake.config['output_matrix']['disease']}_{snakemake.config['run_id']}_normalized.txt" if snakemake.config['matrix_params']['normalized'] else "{snakemake.config['output_matrix']['disease']}_{snakemake.config['run_id']}.txt"
os.system(f"cat {disease_matrix} > {combined_matrix} && tail -n +2 {control_matrix} >> {combined_matrix}")