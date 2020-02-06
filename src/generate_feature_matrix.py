import argparse
from Bio import SeqIO
import csv
import sys
import pyBigWig as pbw
import tools
import find_overlap
import subprocess
import os
import pandas as pd
import numpy as np
from scipy import stats

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



def get_max(bw, chrom, start, end, metric):
    # The number of 10 BP bins in the interval
    nBins = int((end - start) / 10)
    lst =  bw.stats(chrom, start, end, type='mean', nBins=nBins)
    return max(lst)

''' Add data to the matrix given a TSV coordinate file and feature files, along with a
    flag (indicating whether this is a gold-standard file or not) '''
def parse_coords(coord_file, data, bw_files, overlap_coords, flag, length_flag):
    # Input variants
    variants = coord_file

    with open(variants, 'r') as variant_file:
        if not bw_files:
            bw_files = []
        for filename in bw_files:
            data[0].append(filename.split('/')[-1:][0] + '_mean')

        # Constructing the header
        intersection_dicts = []
        for filename in overlap_coords:
            data[0].append(filename + '_overlap')
            intersection_filename = '.'.join(coord_file.split('.')[:-1]) + '_' + filename.split('/')[-1] + '_ints.tsv'
            s = 'bedtools intersect -a ' + coord_file + ' -b ' + filename + ' -wao > ' + intersection_filename
            os.system(s)

            with open(intersection_filename, 'r') as ints:
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

            cleanup = 'rm ' + intersection_filename
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
            if length_flag:
                row.append(end-start)


            # BigWig data parsing

            for signal in bw_files:
                bw_sig = pbw.open(signal)
                avg = get_content(bw_sig, chrom, start, end, 'mean')
                if not avg:
                    avg = 0.0
                row.append(avg)


            for i, overlap_file in enumerate(overlap_coords):
                d = intersection_dicts[i]
                if d['filename'] != overlap_file:
                    raise Exception('Dictionary Error')
                row.append(d[(chrom, start, end)] / (end - start))


            data.append(row)

def generate_matrix(bw_files, coord_file, output_file, overlap_files,
                    svtype, length_flag):

    if svtype not in ['DEL', 'DUP', 'INV']:
        raise ValueError('SVtype must be one of DEL, DUP, or INV')

    # Adding all input files in TSV coordinate format to the coord_files array
    if coord_file == None:
        raise ValueError("Error: Must supply variants")


    # Adding overlap interval coordinate files to the overlap_coords list
    overlap_coords = []
    if overlap_files != None:
        for filename in overlap_files:
            overlap_coords.append(filename)

    # Creating the feature matrix and its header
    data = [[]]
    data[0].append('Label (1 = gold standard, 0 = negative)')
    data[0].append('Chromosome')
    data[0].append('Start')
    data[0].append('End')
    if length_flag:
        data[0].append('Length')


    # Generate feature matrix data
    parse_coords(coord_file, data, bw_files, overlap_coords, args.flag, length_flag)

    # Output feature matrix to designated output file in TSV format
    with open(output_file, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        for row in data:
            writer.writerow(row)

def randomize(root_file, randomized_num, ref_genome, blacklist_file):
    print("Generating random coordinates...")
    base = root_file.split('.')[0]
    for iter in range(randomized_num):
        if blacklist_file:
            randomized_sv = '''bedtools shuffle -i ''' + root_file + ''' -g ''' + reference_genome + ''' -chrom -excl ''' + blacklist_file + ''' | awk '{print $1\"\\t"$2\"\\t"$3}' > '''+ base + '''_randomized_''' + str(iter) + '.txt'
        else:
            randomized_sv = '''bedtools shuffle -i ''' + root_file + ''' -g ''' + reference_genome + ''' -chrom | awk '{print $1\"\\t"$2\"\\t"$3}' > '''+ base + '''_randomized_''' + str(iter) + '.txt'
        os.system(randomized_sv)


def generate_random_matrix(bw_files, coord_file, output_file, overlap_files,
                    svtype, randomized_num, reference_genome,
                    blacklist_file, length_flag):
    randomize(coord_file, randomized_num, reference_genome, blacklist_file)
    random_files = []
    for i in range(randomized_num):
        base = coord_file.split('.')[0]
        random_files.append(base + '_randomized_' + str(i) + '.txt')
    print("Generating matrices...")
    generate_matrix(bw_files=bw_files, coord_file=coord_file,
                    output_file=output_file, overlap_files=overlap_files,
                    svtype=svtype, length_flag=length_flag)

    for index, file in enumerate(random_files):
        generate_matrix(bw_files=bw_files, coord_file=file,
                        output_file=output_file + '_randomized_' + str(index) + '.txt',
                        overlap_files=overlap_files, svtype=svtype, length_flag=length_flag)

''' Z-Score Normalize features in the matrix. Ignores the length column '''
def normalize(output_file, randomized_num, length_flag):
    print("Normalizing features...")
    orig_matrix = pd.read_csv(output_file, sep="\t")
    columns = orig_matrix.columns;
    matrices = {output_file: orig_matrix}
    for index in range(randomized_num):
        matrix_name = output_file + '_randomized_' + str(index) + '.txt'
        matrices[matrix_name] = pd.read_csv(matrix_name, sep="\t")
    for column in columns:
        column_matrix = []
        if column in ['Chromosome', 'Start', 'End', 'Length'] or 'Label' in column:
            continue
        for matrix in matrices:
            vals = matrices[matrix][column].values[:]
            column_matrix.append(vals)

        if np.max(np.abs(np.array(column_matrix))) > 0:
            zs = stats.zscore(np.array(column_matrix).flatten())

            for i in range(len(matrices.keys())):
                skip = len(orig_matrix[column].values)

                matrices[list(matrices.keys())[i]][column] = zs[skip*i:skip*i+skip]

    matrices[output_file].to_csv(output_file.split('.')[0] + '_normalized.txt',
                                 index=False, sep="\t")





if __name__ == '__main__':
    # Parse for input file types via flags
    parser = argparse.ArgumentParser(description='Generate a feature matrix from SV data')
    parser.add_argument('-c', '--coordinate', help='SV files in coordinate format')
    parser.add_argument('-b', '--bigwig', nargs='+', help='Feature files in BigWig format')
    parser.add_argument('-g', '--gene_overlap', nargs='+', help='Gene coordinates for overlap features')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('-f', '--flag', action='store_true', help='Gold standard sample')
    parser.add_argument('-t', '--svtype', help='Type of SV for this model (e.g., DEL, DUP, INV)')
    parser.add_argument('-r', '--randomized_num', type=int, help='Number of sets of randomized coordinates to generate, if any')
    parser.add_argument('-rg', '--reference_genome', help='Reference genome file for randomized generation')
    parser.add_argument('-bf', '--blacklist_file', help='File of coordinates to blacklist')
    parser.add_argument('-z', '--z_normalize', action='store_true', help='Output Z-Score normalized matrices as well')
    parser.add_argument('-l', '--length_flag', action='store_true', help='Include SV length as a feature')

    print("Parsing arguments...")
    # Store cmd-parsed arguments into lists
    args = parser.parse_args()
    bw_files = args.bigwig
    coord_file = args.coordinate
    output_file = args.output
    overlap_files = args.gene_overlap
    svtype = args.svtype
    randomized_num = args.randomized_num
    reference_genome = args.reference_genome
    blacklist_file = args.blacklist_file
    norm = args.z_normalize
    length_flag = args.length_flag

    if randomized_num and not reference_genome:
        raise Exception('Must supply reference genome to generate randomized coordinates')

    if not randomized_num:
        generate_matrix(bw_files=bw_files, coord_file=coord_file,
                    output_file=output_file, overlap_files=overlap_files,
                    svtype=svtype, length_flag=length_flag)
        if norm:
            normalize(output_file=output_file, randomized_num=0, length_flag=length_flag)
    else:
        generate_random_matrix(randomized_num=randomized_num, bw_files=bw_files, coord_file=coord_file,
                    output_file=output_file, overlap_files=overlap_files,
                    svtype=svtype, reference_genome=reference_genome,
                    blacklist_file=blacklist_file, length_flag=length_flag)
        if norm:
            normalize(output_file=output_file, randomized_num=randomized_num, length_flag=length_flag)
    print('Cleaning up temporary files...')
    cleanup_cmd1 = 'rm ' + output_file + '_randomized*'
    cleanup_cmd2 = 'rm ' + coord_file.split('.')[0] + '_randomized*'
    os.system(cleanup_cmd1)
    os.system(cleanup_cmd2)
    print("Feature matrix generation complete!")
