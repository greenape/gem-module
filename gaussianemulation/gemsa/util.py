import os.path
import csv
from itertools import chain, product, permutations
import numpy as np

def read_main_effects(dir, params, iterations):
    """
    Read main effects file from GEMSA and translate it into
    something more usable. Expects params to be a list of
    names, in the same order as they are in GEMSA.
    """
    results = {'param':[], 'Main effect':[]}
    with open(os.path.join(dir, "maineffects.out")) as eff_in:
        raw = csv.reader(eff_in, delimiter=" ")
        for p in params:
            for i in range(iterations):
                try:
                    realisation = [float(elem) for elem in raw.next() if elem != ""]
                    results['Main effect'] += realisation
                    results['param'] += [p]*len(realisation)
                except:
                    pass
    results['realisation'] = list(chain.from_iterable([[i]*len(realisation) for i in range(iterations)]*len(params)))
    results['grid'] = list(chain.from_iterable([np.linspace(0, 1, len(realisation)) for i in range(iterations*len(params))]))
    return results

def read_joint_effects(dir, params, iterations, gridsize):
    """
    Read joint effects file from GEMSA and translate it into
    something more usable. Expects params to be a list of names,
    in the same order as they are in GEMSA.
    """
    results = {'x_param':[], 'y_param':[], 'Joint effect':[]}
    joints = reduce(lambda x, y: x + [y] if y[::-1] not in x else x, permutations(params, 2), [])
    with open(os.path.join(dir, "jointeffects.out")) as eff_in:
        raw = csv.reader(eff_in, delimiter=" ")
        for (a, b) in joints:
            for i in range(iterations):
                try:
                    # Line goes a1.b1, a1.b2 ... a1.bn, a2.b1, a2.b2 .. an
                    realisation = [float(elem) for elem in raw.next() if elem != ""]
                    results['Joint effect'] += realisation
                    results['x_param'] += [a]*len(realisation)
                    results['y_param'] += [b]*len(realisation)
                except:
                    pass
    results['realisation'] = list(chain.from_iterable([[i]*len(realisation) for i in range(iterations)]*len(joints)))
    results['x'], results['y'] = zip(*chain.from_iterable([product(np.linspace(0, 1, gridsize), np.linspace(0, 1, gridsize))
                                                for i in range(iterations*len(joints))]))
    return results
