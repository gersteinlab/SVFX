# SVFX: Using Random Forests to Prioritize Structural Variants

## `generate_feature_matrix.py`

Usage: `python3 generate_feature_matrix.py -v [SV files in VCF format] -c [SV files in tab-separated coordinate format] -b [Feature files in BigWig format] -g [Gene coordinates for overlap features] -o [File root for output] -d [SV files in BED format] -f [True if class label 1, False otherwise] -t [SV type (e.g. DEL, DUP, INV)]`

Given an input list of structural variants (in VCF format with the `-v` flag or in tab-separated format with the `-c` flag), generates a matrix of feature values calculated from the given BigWig and coordinate-overlap files for input to a model.

## `rf_model.py`

Usage: `python3 rf_model.py -i [input feature matrix] -d [File of indices of features to ignore (one per line)] -t [Index of the class label in the feature matrix] -n [Number of trees in the forest] -m [Maximum depth of a tree in the forest] -s [Minimum number of samples required to split an internal node] -c [Number of cancerous SV's in the matrix] -l [Total number of SV's in the matrix] -o [Output root filename]`

Given an input feature matrix, trains a set of 10 random forest classifiers on disjoint tenths of the dataset, then predicts class labels for each member of the dataset with the models that did not train on them. Saves the created models, as a Python list of trained sklearn models, to `[file root]_ten_models.pkl`. Saves the indices used to split the diseased data into tenths in `[file root]_cancer_indices.pkl`, and those used to split the control data in `[file root]_kg_indices.pkl`.
 
 Predictions on the training data (from the models that did not train on each specific SV) are outputted in dictionary format to `[file root]_predictions.pkl`. Finally, AU(ROC/PRC) on this data are calculated 
 and outputted to `[file root]_ROC_PRC.txt`.

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
