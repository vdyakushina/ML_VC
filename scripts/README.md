## Scripts

* align.sh aligns the reads with bwa mem, sort and index them
	Usage: align.sh <fastq1> <fastq2> pref 
		fastq1, fastq2: path/to/filename_R{1,2}.fastq.gz
		pref: prefix of output bam file

* bqsr.py recalibrates the base quality of reads
	Usage: bqsr.py <bam> pref
		bam: path/to/file.bam
		pref: path/to/pref for output files

* call_freebayes.sh call variants with freebayes in parallel
	Usage: call_freebayes.sh <bam> pref
		bam: path/to/my.bam sorted and indexed
		pref: prefix for output file which is prefix.freebayes.raw.vcf 

* check_logs.sh parsers all \*.log files in the folder for error related 
	messages. 
	Usage: ./check_logs.sh 

* consolidate_json.R consolidates the content of all json files in the folder 
	Usage: Rscript consolidate_json.R <folder> project
		folder: path/to/folder with the resuls of processing having *.json files
		project: id of the project to save output as project.csv

* copy_fastq.py copies fastq.gz files to new location. 
	Identifies the sample ids from the basename.
	Puts the files in the raw subfolder.
	Usage: copy_fastq.py <input> <output>
		input: path/to/input folder 
		output: path/to/output folder

* cov.sh estimates the target coverage with mosdepth
	Usage: cov.sh <bam> target pref
		bam: path/to/filename.bam	
		target: target id of the file with targets
		pref: prefix of output files 

* depth.sh estimates the target coverage with samtools
	Usage: depth.sh <bam> target pref
		bam: path/to/filename.bam	
		target: target id of the file with targets
		pref: prefix of output files 

* depthnuc.py estimate average and target coverage with 'samtools depth',
  counts the metrics (eveness score, MADP, uniformity)
	Usage: python3 depthnuc.py <bam> target pref
		bam: path/to/filename.bam	
		target: target id of the file with targets
		pref: prefix of output files 

* depth_stat.py counts average coverage of the target and uniformity
	from the file generated with 'samtools coverage' command
	Usage: depth_stat.py <file>
		file: path/to/filename.txt 

* escore.R estimates the eveness score according to eq. 5b 
	http://www.nature.com/articles/jhg201621
	Usage: Rscript escore.R <file>
		file: path/to/filename.txt created with 'samtools coverage' command

* fastqc.sh runs QC with FastQC
	Usage: fastqc.sh <fastq1> <fastq2>
		fastq1, fastq2: path/to/filename_R{1,2}.fastq.gz

* multifolder_process.pl process each subfolder having the sample raw data with 
	process.py script
	Usage: folder_process.pl <folder>
		folder: /path/to/folder with subfolders to be processed

* mutect.py call for variants in tumor only mode
	Usage: mutect.py <bam> pref
		bam: path/to/sample.bam index file. The SM tag in bam file is required. 

* get_average_coverage.py reads sample.regions.bed.gz created 
	by mosdepth and estmate the average coverage of the target
	Usage: get_average_coverage.py <bed>
		bed: path/to/prefix.regions.bed.gz

* get_breadth.py parsers sample.mosdepth.region.dist.txt to extract 
	breadth 20 and 10
	Usage: get_breadth.py <file>
		file: path/tp/sample.mosdepth.region.dist.txt

* get_uniformity.py estimates the uniformity of coverage as the percent of targets
	covered higher than 20% of the mean coverage
	Usage: get_uniformity.py <bed> avgcov 
		bed: path/to/prefix.regions.bed.gz created with mosdepth
		avgcov: average coverage of the target estimated with get_average_coverage.py

* mapd.py estimates the MAPD metric as described [here](https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/interpret/metrics)
	Usage: python3 mapd.py <bed>
		bed: path/to/filename.bed having at least four columns (chr, start, end, cov) without header,
			where 'cov' contains the coverage of genomic regions.  

* ontarget.py estimates the coverage on target
	Usage: ontarget.py <bam> target
		bam: path/to/filename.bam	
		target: target id of the file with targets

* pipeline_ABC_ILLMN.py runs all steps of the pipeline.
	The output files are placed in the '/analysis' subfolder,
	which is created in the current folder. 
	Usage: pipeline_ABC_ILLMN.py <bam> target 
		bam: path/to/file.bam with aligned reads
		target: target id

* plot_chr_coverage.R plots the chromosomal coverage 
	estimated with samtools idxstats.
	Usage: Rscript plot_chr_coverage.R id
		id: sample id, assuming that id.bam and id.bam.bai are availalble.
			The plot is saved as id.chr.cov.png

* plot_comulative_coverage.R draws the graphics with comulative 
	coverage using the data obtained with mosdepth
	Usage: plot_comulative_coverage.R <file> cov pref
		file: path/tp/sample.mosdepth.region.dist.txt
		cov: average coverage
		pref: prefix for ouput file

* preprocess.py preprocess the data according to the pipeline
	Usage: preprocess.py <bam> <bed> <pref>
		bam: path/to/file.bam
		bed: path/to/file.bed with the genomic coordinates of the target
		pref: path/to/pref to save the results, where pref is a basename of the files 


* process.py runs the pipeline. 
	Usage: process.py <fastq1> <fastq2> pref 
		fastq1, fastq2: path/to/filename_R{1,2}.fastq.gz
		pref: prefix of output files, usually sample id

* run_ephagen.sh estimates the sensitivity with EphaGen solution
	Usage: run_ephagen <bam> pref
		bam: path/to/my.bam sorted and indexed
		pref: prefix for output files: prefix.{BRCA,CFTR}.{tsv,vcf}

* run_vd.sh runs variant calling with VarDictJava. 
	Usage: run_vd.sh <bam> target sample
		bam: path/to/sample.bam
		target: target id
		sample: sample id

* sinvict.py calls for variants 
	Usage: sinvict.py <bam> target <folder> sample
		bam: path/to/sample.bam indexed file
		target: target id
		folder: path/to/folder to save output
		sample: sample id

* sinvict_to_vcf.py convert the output file generated with sinvict.py into vcf format
	Usage: sinvict_to_vcf.py <infile> pref
		infile: path/to/filename.sinvict file
		pref: path/to/prefix of output vcf files (SNV and indels are saved in different files)
