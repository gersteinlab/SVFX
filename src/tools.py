from __future__ import print_function
import sys
from Bio import SeqIO
#import allel
import pyBigWig as pbw


''' Get GC content of a sequence '''
def gc_content(seq):
    gc = 0
    if len(seq) == 0:
        return 0
    for let in seq:
        if let == 'G' or let == 'C':
            gc += 1
    return float(gc) / float(len(seq))

def dict_from_ref(reference):
    return SeqIO.to_dict(SeqIO.parse(reference, 'fasta'))

#''' Initialize the VCF reader for a file '''
#def get_reader(file):
#    reader = allel.read_vcf(file, fields='*')
#    return reader

def sum_vals(vals):
	ret = 0
	for item in vals:
		if item == None:
			continue
		ret += item
	return ret

def histone_content(histone_file, chrom, start, end, metric):
	bins = int((end - start) / 10)
	bw = pbw.open(histone_file)
	print(bw.chroms(chrom))
	if metric != 'mean' and metric != 'max':
		raise ValueError('Histone content metric must be mean or max')
	print('chrom=', chrom)
	print('start=', start)
	print('end=', end)
	values = bw.stats(chrom, start, end, type='mean', nBins=bins)
	if metric == 'mean':
		return sum_vals(values)
	return max(values)


''' Create SIM instructions to emulate creation of an SV similar to given SV,
    given the ALT and CS fields of an SV'''
def sv_to_instruction(reader, i):
    s = ''
    alt = str(reader['variants/ALT'][i][0])
    cs = str(reader['variants/CS'][i])
    instr_dict = {
        'INS': 'INR ',
        'DUP': 'DUP ',
        'ALU': 'INR ',
        'DEL': 'DEL ',
        'INV': 'INV '
    }
    found = False
    for key in instr_dict.keys():
        if key in alt or key in cs:
            found = True
            s += instr_dict[key]
            break
    if not found:
        raise ValueError("VCF does not contain a metric in the\
                                            following strings: ", alt, ' ', cs)
    length = get_length(reader, i)
    s += str(length) + ' '
    s += str(reader['variants/CHROM'][i])
    return s


''' Get the length of an SV '''
def get_length(reader, i):
    pos = reader['variants/POS'][i]
    end = reader['variants/END'][i]
    length = reader['variants/SVLEN'][i]
    if length < 0:
        length = end - pos
    return length


''' Create a SIM file for SV's similar to given SV's '''
def create_sim(input_file, output_file):
    reader = get_reader(input_file)
    with open(output_file, 'w') as output:
        for i in range(len(reader['variants/POS'])):
            output.write(sv_to_instruction(reader, i) + '\n')

def get_histone_data(input_file, bw_file):
	data = []
	f = open(input_file, 'r')
	for line in f:
		nums = line.split(' ')
	data.append(histone_content(bw_file, nums[0], int(nums[1]), int(nums[2]), 'mean'))	



''' Create a histogram for GC or SV length data '''
def create_hist(input_file, reference, metric, bw_file=None):
    ref_dict = dict_from_ref(reference)
    reader = get_reader(input_file)
    if metric not in ['GC', 'histone', 'Histone', 'Length', 'gc', 'length']:
        raise ValueError("Type must be 'gc', 'length', or 'histone'")
    data = []
    length = 0
    for i in range(len(reader['variants/SVLEN'])):
        pos = reader['variants/POS'][i]
        chrom = reader['variants/CHROM'][i]
        length = get_length(reader, i)
        if metric == 'length' or metric == 'Length':
            data.append(length)
        elif metric == 'gc' or metric == 'GC':
            seq = ref_dict[chrom][pos:pos+length]
            data.append(gc_content(seq))
        elif metric == 'histone' or metric == 'Histone':
            if chrom != None:
                print('problematic chrom:', chrom)
                data.append(histone_content(bw_file, chrom, pos, pos+length, 'mean'))
	    
    if metric == 'length' or metric == 'Length':
        plt.hist(data, bins='auto', normed=True, range=(0, 7500), label="1 KG SV's", alpha=0.5) 

    elif metric == 'gc' or metric == 'GC':
        plt.hist(data, bins='auto', normed=True, label="1 KG SV's", alpha=0.5)

    elif metric == 'histone' or metric == 'Histone':	
        plt.hist(data, bins='auto', normed=True, label="1 KG Histone Marks", alpha=0.5)
