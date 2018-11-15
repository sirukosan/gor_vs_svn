from sklearn.metrics import accuracy_score


def get_accuracy(y_test, y_pred):
    """
    Getting accuracy (Q3) between predicted and testing sequences
    :param y_test: test vector
    :param y_pred: predicted vector
    :return: q3 accuracy score
    """
    return accuracy_score(list(y_test), list(y_pred))
