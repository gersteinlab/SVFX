import numpy as np
import pandas as pd
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier
import matplotlib as plt
from sklearn.model_selection import cross_val_predict
import joblib
import sys
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description="Score a matrix of SV's with existing trained models")
parser.add_argument('-i', '--input_matrix', help='Input feature matrix')
parser.add_argument('-d', '--delete', help='File of feature matrix indices to ignore (one per line)')
parser.add_argument('-m', '--models', help='Serialized model list')
parser.add_argument('-t', '--target', help='Target index')
parser.add_argument('-o', '--output_file', help='Output filename, without extension')
args = parser.parse_args()
plt.use('Agg')


# Constructing the data for the forest from the file and list of indices to ignore
sv_file = args.input_matrix
ind_file = args.delete
sv_data = pd.read_table(sv_file)
remove_inds = []
with open(ind_file) as inds:
	try:
		for line in inds:
			remove_inds.append(int(line))
	except ValueError:
		print("Error: " + line + " is not a valid index")

with open(sv_file) as svs:
	total_length = sum(1 for line in svs) - 1

indices = list(range(sv_data.shape[1]))
for ind in remove_inds:
	indices.remove(ind)
data = sv_data.values[:, indices]

# Separately extracting target data
targets = sv_data.values[:, int(args.target)]

# Converting to pandas DataFrames, as they are what scikit-learn accepts
features = pd.DataFrame(data)
targets_frame = pd.DataFrame(targets)

# Scoring the SV's
lst = []
scores = {}
for i in range(total_length):
	scores[i] = []

lst = joblib.load(args.models)
for i in range(10):
	model = lst[i]
	model.n_jobs = 1
	for j in range(total_length):
		scores[j].append(model.predict_proba(data[j, :].reshape(1, -1)))

for key in scores:
	if len(scores[key]) > 0:
		scores[key] = sum(scores[key]) / len(scores[key])

# Outputting the scores
with open(args.output_file, 'w') as out:
	out.write("Score\tLabel\tChromosome\tStart\tEnd\n")
	for i in range(sum(1 for key in scores)):
		out.write(str(scores[i][0][1]) + '\t' + str(targets[i]) + '\t' + str(sv_data.values[i][1]) + '\t' + str(sv_data.values[i][2]) + '\t' + str(sv_data.values[i][3]) + '\n')
