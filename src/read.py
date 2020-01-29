import sys
import os

for item in range(int(sys.argv[1]),int(sys.argv[2])):


    #upatedFeatureMatrix = ' cat ./1kg_randomized_del/1kg_del_matrix.randomized.'+str(item)+'.updated.tsv | sort | uniq > ./1kg_randomized_del/1kg_del_matrix.randomized.'+str(item)+'.final.tsv'
    #upatedFeatureMatrix = '''awk 'BEGIN { FS = OFS = \"\\t\" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' ''' +' ./1kg_randomized_del/1kg_del_matrix.randomized.'+str(item)+'.tsv  > ./1kg_randomized_del/1kg_del_matrix.randomized.'+str(item)+'.updated.tsv'
    #removeMatrix = 'rm ./1kg_randomized_del/1kg_del_matrix.randomized.'+str(item)+'.tsv' 
 
    
    #upatedFeatureMatrix = '''awk 'BEGIN { FS = OFS = \"\\t\" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = 0 }; 1' ''' +' ./randomized_del/brca_del_matrix.randomized.'+str(item)+'.tsv  > ./randomized_del/brca_del_matrix.randomized.'+str(item)+'.updated.tsv'
    removeMatrix = 'rm ./randomized_del/brca_del_matrix.randomized.'+str(item)+'.tsv' 
    
    #print(upatedFeatureMatrix)
    #os.system(upatedFeatureMatrix)

    print(removeMatrix)
    os.system(removeMatrix)

