#!/bin/bash
# A wrapper script for copy_fastq.py
# Initiate
filename=$1
declare -A path=(
	[scripts]=/pipeline/scripts
)

while read folder; do
	echo ${folder}
	${path[scripts]}/copy_fastq.py \
	${HOME}/data/${folder} \
	${HOME}/proc/${folder}
done < ${filename}
