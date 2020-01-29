import argparse
from Bio import SeqIO
import csv
import sys
import pyBigWig as pbw
import tools
import find_overlap
import subprocess
import os

'''Find the mean absolute deviation (MAD) of a list, given the list and its
   average '''
def mad(lst, avg):
	if avg is None:
		return 0
	if len(lst) == 0:
		return 0
	ret = []
	for elem in lst:
		if elem == None:
			elem = 0
		ret.append(abs(float(elem) - avg))
	return sum(ret) / len(ret)

'''Convert a VCF file to a TSV file of coordinates '''
def vcf_to_coords(vcf_file, output_file):
	# Create a VCF Reader object for the VCF file
	reader = tools.get_reader(vcf_file)
	with open(output_file, 'w') as out:
		# Write the coordinates in TSV form
		for i in range(len(reader['variants/SVLEN'])):
			start = reader['variants/POS'][i]
			chrom = 'chr' + reader['variants/CHROM'][i]
			end = str(max(int(reader['variants/SVLEN'][i]) + int(start), int(reader['variants/END'][i])))
			af = reader['variants/AF'][i][0]
			out.write(str(chrom) + '\t' + str(start) + '\t' + str(end) + '\t' + str(af) + '\n')

'''Convert a CSV file to a TSV file of coordinates '''
def csv_to_coords(csv_filename):
	with open(csv_filename, 'rb') as csvfile:
		gene_reader = csv.reader(csvfile)
		out_name = csv_filename.split('.')[0] + '.txt'
		with open(out_name, 'w') as out:
			for row in gene_reader:
				coords = row[3]
				coords = coords.split(':')
				chrom = coords[0]
				if coords[1] != '-':
					coords = coords[1].split('-')
					start = coords[0]
					end = coords[1]
					out.write(chrom + '\t' + start + '\t' + end + '\n')
	return out_name

'''Add the elements of a list that may contain NoneType elements (Note: this
   differs from Python's built-in sum() function because it ignores NoneType elements,
   which may be generated in the processing of a BigWig file)'''
def add_lst(lst):
	ret = 0
	for num in lst:
		if num != None:
			ret += num
	return ret

''' Evaluate a BigWig file at a particular interval '''
def get_content(bw, chrom, start, end, metric):
	# The number of 10 BP bins in the intervali
	avg = bw.stats(chrom, start, end, type='mean')
	if metric == 'mean':
		return avg[0]
	else:
		lst = bw.stats(chrom, start, end, type='mean', nBins = (end - start) / 10)
		return mad(lst, avg[0])

''' Check if a cancer variant is intrachromosomal '''
def is_valid_variant(row, svtype):
	if row[0] == row[3] and svtype in row[10] and int(row[4]) - int(row[1]) <= 200000:
		return True
	return False


'''Convert a cancer bed file to a TSV coordinate file by filtering out invalid coordinates '''
def bed_to_coords(bed_files, sv_type):
	outputs = []
	for bed_file in bed_files:
		with open(bed_file, 'r') as input_file:
			with open(bed_file.split('.')[0] + svtype + '.coord', 'w') as output_file:
				rd = csv.reader(input_file, delimiter="\t")
				for row in rd:
					if is_valid_variant(row, sv_type):
						output_file.write('chr' + row[0] + '\t' + row[1] + '\t' + row[4]+ '\n')
				outputs.append(bed_file.split('.')[0] + svtype + '.coord')
	return outputs

''' Convert chromosome to single-character designation (e.g., 'chrX' to 'X')'''
def chrom_to_num(chrom):
	return chrom[3:]

''' Get GC content of an interval from the reference genome'''
def gc_cont(chrom, start, end, ref_dict):
	num_chrom = chrom_to_num(chrom)
	return tools.gc_content(ref_dict[num_chrom][start:end])


def get_max(bw, chrom, start, end, metric):
	# The number of 10 BP bins in the interval
	nBins = int((end - start) / 10)
	lst =  bw.stats(chrom, start, end, type='mean', nBins=nBins)
	return max(lst)

''' Add data to the matrix given a TSV coordinate file and feature files, along with a
    flag (indicating whether this is a gold-standard file or not) '''
def parse_coords(coord_file, data, bw_files, overlap_coords, flag):
	# Input variants
	variants = coord_file
	ref_dict = tools.dict_from_ref('hs37d5.fa')
	with open(variants, 'r') as variant_file:

		for filename in bw_files:
			#data[0].append(filename + '_mean')
			data[0].append(filename.split('/')[-1:][0] + '_mean')
			#data[0].append(filename + ' max')
			#data[0].append(filename + ' mad')
		# Constructing the header
		intersection_dicts = []
		for filename in overlap_coords:
			data[0].append(filename + '_overlap')
                        #s = 'bedtools intersect -a ' + coord_file + ' -b ' + filename + ' -wao > ' + filename + '_ints.tsv'
			s = 'bedtools intersect -a ' + coord_file + ' -b ' + filename + ' -wao > ' + '.'.join(coord_file.split('.')[:-1])+'_'+filename + '_ints.tsv'
			os.system(s)
                        #cmd = subprocess.Popen(s, shell=True, stdout=subprocess.PIPE)
                        #with open(filename + '_ints.tsv', 'r') as ints:
			with open('.'.join(coord_file.split('.')[:-1])+'_'+filename + '_ints.tsv', 'r') as ints:
                            d = {}
                            d['filename'] = filename
                            for line in ints:
                                line = line.split('\t')
                                chrom = line[0]
                                start = int(line[1])
                                end = int(line[2])
                                nucs = int(line[len(line)-1])
                                if (chrom, start, end) in d:
                                    d[(chrom, start, end)] += nucs
                                else:
                                    d[(chrom, start, end)] = nucs
                            intersection_dicts.append(d)

			cleanup = 'rm '+'.'.join(coord_file.split('.')[:-1])+'_'+filename + '_ints.tsv'
			os.system(cleanup)

		# Adding features for each SV to the feature matrix
		for line in variant_file:

			# Parsing the line for SV coordinates
			nums = line.split('\t')
			chrom = nums[0]
			start = int(nums[1])
			end = int(nums[2])

			# Skipping Y chromosome coordinates, as they may not be present
			# in feature data
			if chrom == 'chrY' or chrom == 'Y':
				continue

			# The feature vector for the current SV
			row = []


			# Gold-standard flag
			if flag:
				row.append(1)
			else:
				row.append(0)
			# SV coordinates
			row.append(chrom_to_num(chrom))
			row.append(start)
			row.append(end)

			# GC content
			#row.append(gc_cont(chrom, start, end, ref_dict))
			# Allele Frequency
			#if flag:
				#row.append(0.005)
			#else:
				#row.append(0.005)
			# BigWig data parsing

			for signal in bw_files:
				bw_sig = pbw.open(signal)
				avg = get_content(bw_sig, chrom, start, end, 'mean')
				row.append(avg)
				#row.append(get_max(bw_sig, chrom, start, end, 'max'))
				#row.append(get_content(bw_sig, chrom, start, end, 'mad'))

			for i, overlap_file in enumerate(overlap_coords):
				d = intersection_dicts[i]
				if d['filename'] != overlap_file:
					raise Exception('Dictionary Error')
				row.append(d[(chrom, start, end)] / (end - start))
				#row.append(find_overlap.bulk_overlap(chrom, start, end, overlap_file + '_ints.tsv', d))

			data.append(row)


if __name__ == '__main__':
	# Parse for input file types via flags
	parser = argparse.ArgumentParser(description='Generate a feature matrix from SV data')
	parser.add_argument('-v', '--vcf', nargs='+', help='SV files in VCF format')
	parser.add_argument('-c', '--coordinate', nargs='+', help='SV files in coordinate format')
	parser.add_argument('-b', '--bigwig', nargs='+', help='Feature files in BigWig format')
	parser.add_argument('-g', '--gene_overlap', nargs='+', help='Gene coordinates for overlap features')
	parser.add_argument('-o', '--output', help='Output file')
	parser.add_argument('-d', '--bed', nargs='+', help="SV files in BED format")
	parser.add_argument('-f', '--flag', action='store_true', help='Gold standard sample')
	parser.add_argument('-t', '--svtype', help='Type of SV for this model (e.g., DEL, DUP, INV)')
	# Convert the reference genome to a chromosome<=>sequence dictionary

	# Store cmd-parsed arguments into lists
	args = parser.parse_args()
	bw_files = args.bigwig
	vcf_files = args.vcf
	coord_files = args.coordinate
	output_file = args.output
	overlap_files = args.gene_overlap
	bed_files = args.bed
	bed_coords = []
	svtype = args.svtype
	if svtype not in ['DEL', 'DUP', 'INV']:
		raise ValueError('SVtype must be one of DEL, DUP, or INV')
	# Adding all input files in TSV coordinate format to the coord_files array
	if coord_files == None:
		coord_files = []
	if bed_files != None:
		bed_coords = bed_to_coords(bed_files, svtype)
	coord_files += bed_coords

	if vcf_files == None and coord_files == []:
		raise ValueError("Error: Must supply variants")
	if vcf_files != None:
		for vcf_file in vcf_files:
			output_name = vcf_file.split('.')
			output_name = output_name[0] + 'txt'
			vcf_coord = vcf_to_coords(vcf_file, output_name)
			coord_files.append(output_name)

	# Adding overlap interval coordinate files to the overlap_coords list
	overlap_coords = []
	if overlap_files != None:
		for filename in overlap_files:
			overlap_coords.append(filename)

	# Creating the feature matrix and its header
	data = [[]]
	data[0].append('Label (1 = gold standard, 0 = negative)')
	data[0].append('Chromosome')
	data[0].append('start')
	data[0].append('end')
	data[0].append('GC Content')
	#data[0].append('Allele Frequency')


	# Generate feature matrix data
	for coord_file in coord_files:
		print(parse_coords(coord_file, data, bw_files, overlap_coords, args.flag))


	# Output feature matrix to designated output file in TSV format
	with open(output_file, 'w') as tsvfile:
		writer = csv.writer(tsvfile, delimiter='\t')
		print(sum(1 for row in data))
		for row in data:
			writer.writerow(row)



