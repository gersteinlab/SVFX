import sys
import os
import argparse

for item in range(int(sys.argv[1]),int(sys.argv[2])):

    #randomizeSV = '''bedtools shuffle -i brca_somatic_del.bed -g hs37d5.fa.fai -chrom | awk '{print "chr"$1\"\\t"$2\"\\t"$3}' > ./randomized_del/brca_somatic_del.randomized.'''+str(item)+'.txt'
    randomizeSV = '''bedtools shuffle -i 1kg_del.baseline.bed -g hs37d5.fa.fai -chrom | awk '{print "chr"$1\"\\t"$2\"\\t"$3}' > ./1kg_randomized_del/1kg_del.randomized.'''+str(item)+'.txt'
    #randomizeSV = '''bedtools shuffle -i 1kg_del.bed -g hs37d5.fa.fai -chrom -excl EncodeBlackList.hs37d5.bed | awk '{print "chr"$1\"\\t"$2\"\\t"$3}' > ./1kg_randomized_del/1kg_del.randomized.'''+str(item)+'.txt'

    #print randomizeSV
    os.system(randomizeSV)
