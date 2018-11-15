from utils.svm_utils import *
import glob
import threading


def make_models_for_cv():
    """
    Make models for cross-validation
    """
    for test_file in glob.glob(INPUT_CV_DIR + '*'):
        test_dir = SVM_TESTS_DIR + os.path.basename(test_file)
        os.mkdir(test_dir)
        with open(os.path.join(test_dir, 'log'), 'a') as log_file:
            log_file.write('computing vectors for ' + test_file + '\n')
            x_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'x.csv'
            y_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'y.csv'
            x_excepted_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'x_except.csv'
            y_excepted_file = os.path.join(test_dir, os.path.basename(test_file)) + '_' + 'y_except.csv'
            build_vectors(glob.glob(PROFILES_OUTPUT_CSV_DIR + '*'), x_file, y_file, dssp_dir=INPUT_DSSP_FOLDER,
                          log_out=log_file, excep=get_ids_from_list_file(test_file))
            build_vectors(glob.glob(PROFILES_OUTPUT_CSV_DIR + '*'), x_excepted_file, y_excepted_file,
                          dssp_dir=INPUT_DSSP_FOLDER, log_out=log_file, excep=get_ids_from_list_file(test_file),
                          exclude=False)

            log_file.write('loading data for ' + test_file + '\n')
            y_train = np.loadtxt(y_file, delimiter=",")
            log_file.write('y data loaded' + '\n')
            x_train = np.loadtxt(x_file, delimiter=",")
            log_file.write('x data loaded' + '\n')
            log_file.flush()

            i = 0
            threads = []
            for par_tup in C_GAMMA_SET:
                process_thread = threading.Thread(target=mock_loop,
                                                  args=[i, test_dir, test_file, par_tup, x_train, y_train])
                threads.append(process_thread)
                process_thread.start()
                i += 1

            for thr in threads:
                thr.join()


def make_final_model(c, gamma):
    """
    make model
    :param c: c param
    :param gamma: gamma param
    """
    test_dir = FINAL_TEST_DIR
    os.mkdir(test_dir)
    with open(os.path.join(test_dir, 'log'), 'a') as log_file:
        log_file.write('computing vectors for final test' + '\n')
        x_file = test_dir + 'x.csv'
        y_file = test_dir + 'y.csv'
        build_vectors(glob.glob(PROFILES_OUTPUT_CSV_DIR + '*'), x_file, y_file, dssp_dir=INPUT_DSSP_FOLDER,
                      log_out=log_file)

        log_file.write('loading data for final test' + '\n')
        y_train = np.loadtxt(y_file, delimiter=",")
        log_file.write('y data loaded' + '\n')
        x_train = np.loadtxt(x_file, delimiter=",")
        log_file.write('x data loaded' + '\n')
        log_file.flush()

        loop(0, test_dir, 'final', (c, gamma), x_train, y_train)


make_models_for_cv()
make_final_model(1.8, 0.2)
