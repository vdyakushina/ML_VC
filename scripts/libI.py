import os
import re
import random
import string
import sys
from pprint import pprint
import re
import subprocess
import pathlib

sys.path.insert(1, os.environ['CLAUDIA'])
import AODDB
import Atlas
import Table
import FileSystem


def _unify_interpretation(string):
    if (string is None):
        return None
    string = str(string).lower()
    if (re.search(r"benign",string)):
        return "neutral"
    if (re.search(r"neutral", string)):
        return "neutral"
    if (re.search(r"suspected", string)):
        if (re.search(r"damaging", string)or(re.search(r"deleterious", string))):
            return "suspected damaging"
    if (re.search(r"pathogenic", string)):
        return "damaging"
    if (re.search(r"damaging", string)):
        return "damaging"
    if (re.search(r"deleterious", string)):
        return "damaging"
    if (re.search(r"vus", string)):
        return "vus"
    if (re.search(r"uncertain", string)):
        return "vus"
    if (re.search(r"unkown", string)):
        return "vus"
    if (re.search(r"Conflicting", string)):
        return "conflict"
    return None




















