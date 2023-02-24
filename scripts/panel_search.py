#!/usr/bin/env python3

from subprocess import run, PIPE
import sys

path='/home/onco-admin/ATLAS_software/aod-pipe/panel_info/'

bam_file=sys.argv[1]
panels_list=sys.argv[2]

def panel_search_f(bam, panels):
	panels=panels.split(',')
	d={i:int() for i in panels}
	for p in panels:
		panel=p+'.designed.bed'
		cov=run('samtools view -@6 -q 40 -b %s | bedtools coverage -a %s%s/%s -b stdin | rev | cut -f1,4 | rev' %(bam, path, p,panel), shell=True, stdout=PIPE).stdout.decode().strip().split('\n')
		cov_d=[float(i.split('\t')[0]) for i in cov]
		cov=[float(i.split('\t')[1]) for i in cov]
		fail=[i for i in range(len(cov_d)) if cov_d[i] <=9]
		for i in sorted(fail, reverse=True):
			cov[i]=0
		frac=sum(cov)/len(cov) #fraction of panel covered with bam
		off_target=float(run('samtools view -@6 -q 40 -b %s | bedtools intersect -v -bed -wa -a stdin -b %s%s/%s | wc -l' %(bam, path, p, panel), shell=True, stdout=PIPE).stdout.decode().strip())
		total=float(run('samtools view -c -@6 -q 40 %s' % bam, shell=True, stdout=PIPE).stdout.decode().strip())
		off_target_frac=off_target/total
		final=frac/(off_target_frac+0.00005)
		d[p]=final
	target=max(d, key=d.get)+'\t'+str(d[max(d, key=d.get)])
	return target

target_panel=panel_search_f(bam_file, panels_list)
print(target_panel)