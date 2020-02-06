# Relevant imports
import numpy as np
from sklearn.externals import joblib
import pandas as pd
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score
import argparse
from sklearn.externals import joblib
np.seterr(divide='ignore', invalid='ignore')

# Argument parsing
parser = argparse.ArgumentParser(description='Train a random forest classifier on a feature matrix and generate impact scores.')
parser.add_argument('-i', '--input_matrix', help='Input feature matrix')
parser.add_argument('-d', '--delete', help='File of feature matrix indices to ignore (one per line)')
parser.add_argument('-t', '--target', help='Target index')
parser.add_argument('-n', '--num_trees', help='Number of trees in the forest')
parser.add_argument('-m', '--max_depth', help='Maximum depth of a tree in the forest')
parser.add_argument('-s', '--min_samples_split', help='Minimum number of samples required to split an internal node')
parser.add_argument('-c', '--cancer_num', help='Number of Cancer SVs in the file')
parser.add_argument('-l', '--length', help='Total number of SVs in the file')
parser.add_argument('-o', '--output_file', help='Output filename, without extension')
args = parser.parse_args()

# Constructing the data for the forest from the file and list of indices to ignore
sv_data = pd.read_table(args.input_matrix)
num_cancer = int(args.cancer_num)
total_length = int(args.length)
cancer_indices = np.arange(0, num_cancer)
kg_indices = np.arange(num_cancer, total_length)
n_trees = int(args.num_trees)
m_depth = int(args.max_depth)
m_samples_split = int(args.min_samples_split)

print(cancer_indices)
print(kg_indices)

# Extracting feature data
np.random.shuffle(cancer_indices)
np.random.shuffle(kg_indices)
kg_trains = np.array_split(kg_indices, 10)
cancer_trains = np.array_split(cancer_indices, 10)
remove_inds = []
with open(args.delete) as inds:
	try:
		for line in inds:
			remove_inds.append(int(line))
	except ValueError:
		print("Error: " + line + " is not a valid index.")

# Extracting feature data
indices = list(range(sv_data.shape[1]))
for ind in remove_inds:
	indices.remove(ind)
data = sv_data.values[:, indices]

# Separately extracting target data
targets = sv_data.values[:, int(args.target)]

# Converting to pandas DataFrames, as they are what scikit-learn accepts
features = pd.DataFrame(data)
targets_frame = pd.DataFrame(targets)

# The array to hold the 10 models
model_list = []

# SVIS dictionary (scores[i] contains the 4 scores for the ith index)
scores = {}
for i in range(total_length):
	scores[i] = []

# The set of training indices for each of the 10 models
training_indices = [np.concatenate([cancer_trains[i], kg_trains[i]]) for i in range(10)]

# Training the ten models on disjoint tenths of the dataset
for i in range(10):
	print("Fitting model " + str(i))
	model_list.append(RandomForestClassifier(n_estimators=n_trees, max_depth=m_depth, min_samples_split=m_samples_split))
	X = data[training_indices[i], :]
	y = targets[training_indices[i]]
	y = y.astype("int")
	model_list[i].fit(X, np.ravel(y))

output = args.output_file.split('.')[0]

# Saving the objects through joblib so they can be accessed later if needed
joblib.dump(model_list, output + '_ten_models.pkl')
joblib.dump(cancer_trains, output + '_cancer_indices.pkl')
joblib.dump(kg_trains, output + '_kg_indices.pkl')

# Making predictions
num_pred = 9 * total_length
so_far = 1
for i in range(10):
	model = model_list[i]
	for j in range(total_length):
		if j not in training_indices[i]:
			print("Making prediction " + str(so_far) + " of " + str(num_pred))
			scores[j].append(model.predict_proba(data[j, :].reshape(1, -1)))
			so_far += 1

# Averaging predictions for each SV to get the final cross-validated score
for key in scores:
	if len(scores[key]) > 0:
		scores[key] = sum(scores[key]) / len(scores[key])

# Writing the AUROC to a file, so that we can compare

y_true = np.array(targets).astype(float)
y_score = [scores[i][0][1] for i in range(sum(1 for key in scores))]
with open(output + 'ROC_PRC.txt', 'w') as stats:
	stats.write("AUROC:\n")
	stats.write(str(roc_auc_score(y_true, y_score)))
	stats.write("AUPRC:\n")
	stats.write(str(average_precision_score(y_true, y_score)))

# Storing the prediction dictionary in a pkl file and the predictions in a text file
joblib.dump(scores, output +  '_predictions.pkl')
with open(output + '_predictions.txt', 'w') as out:
	out.write("Score\tChromosome\tStart\tEnd\n")
	for i in range(sum(1 for key in scores)):
		out.write(str(scores[i][0][1]) + '\t' + str(targets[i]) + '\t' + str(sv_data.values[i][1]) + '\t' + str(sv_data.values[i][2]) + '\t' + str(sv_data.values[i][3]) + '\n')
