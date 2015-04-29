
import random

test = {'A': 0.01, 'B': 0.09, 'C':0.9}

class Sampler(object):
    def __init__(self, weights):
        self.weights=weights
        # sanity check
        sm=0.0
        for v in self.weights.itervalues():
            sm += v 
        if abs(1.0-sm) > 0.00001:
            raise Exception("Weights didn't sum to 1")

    def sample(self):
        r=random.random()
        sm=0.0
        for k, v in self.weights.iteritems():
            if r < v-sm: return k
            sm +=v
        return k

