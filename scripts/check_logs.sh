#!/bin/bash
# Parse all *.log files in the current folder
# for matches with error related messages
grep -iP "error|fail|warn|invalid" *.log
