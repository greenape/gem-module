import os.path
import csv
from itertools import chain
import numpy as np

def read_main_effects(dir, params, iterations):
	"""
	Read main effects file from GEMSA and translate it into
	something more usable. Expects params to be a list of
	names, in the same order as they are in GEMSA.
	"""
	results = {'param':[], 'output':[]}
	with open(os.path.join(dir, "maineffects.out")) as eff_in:
		raw = csv.reader(eff_in, delimiter=" ")
		for p in params:
			for i in range(iterations):
				try:
					realisation = [elem for elem in raw.next() if elem != ""]
					results['Main effect'] += realisation
					results['param'] += [p]*len(realisation)
				except:
					pass
	results['realisation'] = list(chain.from_iterable([[i]*len(realisation) for i in range(iterations*len(params))]))
	results['grid'] = list(chain.from_iterable([np.linspace(0, 1, len(realisation)) for i in range(iterations*len(params))]))
	return results
