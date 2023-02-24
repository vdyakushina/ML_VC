#!/usr/bin/env python3
import func
import sys

log = sys.argv[1]
out = func.check_log(log)
if out: sys.exit(f'{log} needs attention')


