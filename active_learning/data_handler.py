

class Handler:
    def __init__(self, n_train: int = 1000, balanced: float = None):
        self.acquired = {}
        self.train_data, self.screen_data = self.build_train_data(n_train=n_train, balanced=balanced)

    def build_train_data(self):
        pass

    def add(self):
        pass

    def get_data(self):
        pass

    def __call__(self, *args, **kwargs):
        return self.get_data()



