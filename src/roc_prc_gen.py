from sklearn.externals import joblib
import numpy as np
from sklearn.metrics import roc_curve, auc, roc_auc_score, precision_recall_curve, average_precision_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
y_test = []
y_scores= []
with open(input_file, 'r') as inp:
	i = 0
	for line in inp:
		if i == 0:
			i += 1
			continue
		lin = line.split('\t')
		y_scores.append(float(lin[0]))
		y_test.append(float(lin[1]))



y_scores = np.array(y_scores)
y_test = np.array(y_test)
precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
fpr, tpr, threshold = roc_curve(y_test, y_scores, pos_label=1.0)
roc_auc = auc(fpr, tpr)


plt.title('Receiver Operating Characteristic')
plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.savefig(output_file)

print('Area under ROC: ', roc_auc_score(y_test, y_scores))
print('Area under PRC: ', auc(recall, precision))
print('AP score: ', average_precision_score(y_test, y_scores))
