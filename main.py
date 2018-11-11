import gzip
import pickle

from utils.svm_utils import build_vectors
from utils.utils import *
from utils.eval_utils import *
import glob

x_file = SVM_BLIND_VECTORS_DIR + 'x.csv'
y_file = SVM_BLIND_VECTORS_DIR + 'y.csv'
#build_vectors(glob.glob(PROFILES_BLIND_OUTPUT_CSV_DIR + '*'), x_file, y_file, BLIND_DSSP_DIR)
y_test = np.loadtxt(y_file, delimiter=",")
x_test = np.loadtxt(x_file, delimiter=",")

model = pickle.load(gzip.open('/home/vadik/src/lb2/data/svm/tests/test4/model_3.pkl.gz', 'r'))
y_pred = model.predict(x_test)

print(y_test)
print(y_pred)
print(get_accuracy(y_test, y_pred))
