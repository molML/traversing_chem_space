

class Model:
    def __init__(self, type: str = 'xgb'):
        self.type = type

    def optimize_hyperparameters(self, **kwargs):
        pass

    def train(self, x, y):
        pass

    def test(self):
        pass

    def predict(self):
        pass

    def __repr__(self):
        return f"{self.type} model wrapper"

    def __call__(self, x, *args, **kwargs):
        return self.predict(x)




class XGBoostEnsemble:
    def __init__(self):
        raise NotImplementedError


class GNNEnsemble:
    def __init__(self):
        raise NotImplementedError


class BayesianGNN:
    def __init__(self):
        raise NotImplementedError





