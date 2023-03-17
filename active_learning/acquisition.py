

class Acquisition:
    def __init__(self, type: str):
        self.type = type
        pass

    def acquire(self, n, y_hat, y_hat_sigma, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        pass
        # return idx


def random(data_size, n: int = 1):
    """ select n random samples """
    pass


def greedy_exploitative(y_hat, n: int = 1):
    """ Get the n highest predicted samples """
    pass


def greedy_explorative(y_hat_sigma, n: int = 1):
    """ Get the n most uncertain samples """
    pass


def clusters(n: int = 1):
    """ Select the n samples most representative of n clusters in the data """
    pass


def epsilon_greedy_exploitative(epsilon: float = 0.1, n: int = 1):
    """ greedy exploitative with a epsilon chance of selecting a random sample """
    pass


def epsilon_greedy_explorative(epsilon: float = 0.1, n: int = 1):
    """ greedy exploitative with a epsilon chance of selecting a random sample """
    pass

