# The module contains usefull functions for the pipeline

import pandas as pd
import re
from subprocess import check_output
from datetime import datetime

def read_resources(csv):
	df = pd.read_csv(csv, header = None, usecols = [0, 1], names = ["id", "path"])
	return df.set_index('id').T.to_dict('records')[0]

def ts():
	dt = datetime.now()
	return(dt.strftime("%d.%m.%Y %H:%M:%S"))

def check_log(log):
	pattern = 'error|warning|fail|invalid'
	out = []
	n = 0
	with open(log) as file:
		line = file.readline()
		while line:
			n += 1
			match = re.match(pattern, line, re.IGNORECASE)
			if match is not None:
				out.append([n, match.group()])
			line = file.readline()
	if out:
		for i in out: print(i)
	return out
