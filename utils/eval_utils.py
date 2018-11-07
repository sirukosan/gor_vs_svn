from sklearn.metrics import accuracy_score


def get_accuracy(y_test, y_pred):
    return accuracy_score(list(y_test), list(y_pred))