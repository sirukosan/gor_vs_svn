from utils.svm_utils import *
import glob
import threading

for test_file in glob.glob(INPUT_CV_DIR + '*'):
    test_dir = SVM_TESTS_DIR + os.path.basename(test_file)
    os.mkdir(test_dir)
    with open(os.path.join(test_dir, 'log'), 'a') as log_file:
        log_file.write('computing vectors for ' + test_file + '\n')
        x_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'x.csv'
        y_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'y.csv'
        build_vectors(glob.glob(PROFILES_OUTPUT_CSV_DIR + '*'), x_file, y_file, log_file, INPUT_DSSP_FOLDER,
                      exclude=get_ids_from_list_file(test_file))

        log_file.write('loading data for ' + test_file + '\n')
        y_train = np.loadtxt(y_file, delimiter=",")
        log_file.write('y data loaded' + '\n')
        x_train = np.loadtxt(x_file, delimiter=",")
        log_file.write('x data loaded' + '\n')

        i = 0
        threads = []
        for par_tup in C_GAMMA_SET:
            process_thread = threading.Thread(target=loop, args=[i, test_dir, test_file, par_tup, x_train, y_train])
            threads.append(process_thread)
            process_thread.start()
            i += 1

        for thr in threads:
            thr.join()
