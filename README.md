# SVFX: Using Random Forests to Prioritize Structural Variants

SVFX is a machine learning based tool to assign pathogenic scores to large deletions and duplications.

## necessary python package installations:
pip install scikit-allel

pip install pyBigWig

pip install argparse

pip install sklearn

## `generate_feature_matrix.py`

Usage: `python3 generate_feature_matrix.py -c [SV files in tab-separated coordinate format (each line should start with chr)] -b [Feature files in BigWig format] -g [Coordinates for interval overlap-based features] -o [File root for output] -f [Present if class label 1, False otherwise] -t [SV type (e.g. DEL, DUP, INV)] -z [Present if matrix should be Z-Score normalized (length feature will not be normalized)] -r [Number of shuffled coordinates to generate for the randomization process] -rg [Reference genome index file - necessary for the -r flag] -l [If present, length will be included as a feature (and will be ignored during normalization)]`

Sample command (run from the root directory): `python3 src/generate_feature_matrix.py -c data/sample_input_SV.txt -g data/gc19_pc.prom.nr.bed  data/gc19_pc.3utr.nr.bed data/gc19_pc.5utr.nr.bed data/gc19_pc.cds.nr.bed data/gc19_pc.ss.nr.bed data/sensitive.nc.bed data/ultra.conserved.hg19.bed data/wgEncodeBroadHmmGm12878HMM.Heterochrom.bed data/H1-ESC_Dixon2015-raw_TADs.bed -o inputSVcoords -t DEL -z -r 1 -rg data/hs37d5.fa.gz.fai`

Given an input list of structural variants, generates a matrix of feature values calculated from the given BigWig and coordinate-overlap files for input to a model.

## `rf_model.py`

Usage: `python3 rf_model.py -i [input feature matrix] -d [File of indices of features to ignore (one per line)] -t [Index of the class label in the feature matrix] -n [Number of trees in the forest] -m [Maximum depth of a tree in the forest] -s [Minimum number of samples required to split an internal node] -c [Number of cancerous SV's in the matrix] -l [Total number of SV's in the matrix] -o [Output root filename]`

Given an input feature matrix, trains a set of 10 random forest classifiers on disjoint tenths of the dataset, then predicts class labels for each member of the dataset with the models that did not train on them. Saves the created models, as a Python list of trained sklearn models, to `[file root]_ten_models.pkl`. Saves the indices used to split the diseased data into tenths in `[file root]_cancer_indices.pkl`, and those used to split the control data in `[file root]_kg_indices.pkl`.

 Predictions on the training data (from the models that did not train on each specific SV) are outputted in dictionary format to `[file root]_predictions.pkl`. Finally, AU(ROC/PRC) on this data are calculated
 and outputted to `[file root]_ROC_PRC.txt`.

Note: The input matrix must be formatted in such a way that all control group rows follow all non-control rows (in other words, all rows with label 0 follow all rows with label 1). The easiest way to do this is to use `generate_feature_matrix.py` to generate the 0-label and 1-label matrices, then append the 0-label matrix (without title row) to the 1-label matrix.

## `find_overlap.py`

Contains code to calculate overlaps between sections of the genome; used by `generate_feature_matrix.py`.

## `roc_prc_gen.py`

Usage: `python3 roc_prc_gen.py [input file] [output file]`

Given a tab-separated two-column input file of predicted scores and true labels (with a header, as with the following example):

`Predicted	True`
`0.755	1`
`...`
`0.211	0`

Outputs AU(ROC/PRC) information, along with a plot of the ROC, saved to the specified output file.
