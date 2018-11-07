import glob

from utils.gor_utils import *
from utils.eval_utils import *

np_result = gor_train(glob.glob(PROFILES_OUTPUT_CSV_DIR + '*'), INPUT_DSSP_FOLDER)
save_gor_as_csv(np_result)



acc_sum = 0
number = 0

for profile_file in glob.glob(PROFILES_BLIND_OUTPUT_CSV_DIR + '*'):
    y_pred = predict_profile(profile_file, GOR_TRAINED_FILE).casefold()
    y_test = get_dssp_from_file(get_relative_file(profile_file, BLIND_DSSP_DIR, 'dssp')).swapcase()
    if len(y_pred) == len(y_test):
        acc_sum += get_accuracy(y_test, y_pred)
        number += 1

acc = acc_sum/number
print(acc)
