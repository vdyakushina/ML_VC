#!/usr/bin/perl
# Script recalibrates the score of variant quality (VQSR) with GATK
# The output is written to the current folder.
# Author: Gennady Khvorykh

use strict; 
use warnings;
use Getopt::Long;
use Cwd;
use File::Temp;

my ($help);

GetOptions(
	help => \$help,
	h => \$help
);

help() if defined $help;

# Initiate
my ($input, $id) = @ARGV; 
my $res = "/pipeline/scripts/resources.csv";
help() unless defined $input;
help() unless defined $id;

my $ts1=time;

# Read resources
my %path = read_resources($res);

my $wd = getcwd();

# Check input file exists
die "$input doesn't exist" unless -f $input;

# Show input
print "Input: $input\n";
print "Sample: $id\n";

# === main ===

# Subset SNPs
my $output = "$id.snps.vcf";
system("
$path{gatk} SelectVariants \\
	-select-type SNP \\
	-V $input \\
	-O $output
") == 0 or die "error: Select SNPs failed: $!";

# Recalibrate the score of variant quality (VQSR)
$input=$output;
$output="$id.recall.vcf";
system("
$path{gatk} VariantRecalibrator \\
	-R $path{ref} \\
	--variant $input \\
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $path{bundle}/hapmap_3.3.b37.chr.vcf \\
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $path{bundle}/1000G_omni2.5.b37.chr.vcf \\
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $path{bundle}/1000G_phase1.snps.high_confidence.b37.chr.vcf \\
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $path{bundle}/dbsnp_138.b37.chr.vcf \\
	-an DP \\
	-an QD \\
	-an FS \\
	-an SOR \\
	-an MQ \\
	-an MQRankSum \\
 	-an ReadPosRankSum \\
	--mode SNP \\
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
	--output $output \\
	--tranches-file $id.tranches \\
	--rscript-file $id\_plots.R \\
") == 0 or die "error: VariantRecalibrator failed: $!";

# Apply the desired level of recalibration
my $recall=$output;
$output="$id.vqslod.vcf";
print "apply VQSR\n";
system("
$path{gatk} ApplyVQSR \\
	--reference $path{ref} \\
	--variant $input \\
	--mode SNP \\
	--truth-sensitivity-filter-level 99.0 \\
	--recal-file $recall \\
	--tranches-file $id.tranches \\
	--output $output
") == 0 or die "error: Apply VQSR failed: $!";

# Filter the variants based on the probabilities assigned by VQSR
$input = $output;
$output = "$id.recall.flt.vcf";
system("
$path{bcftools} view -i \'FILTER = \"PASS\"\' $input > $output
") == 0 or die "error: Filter variants failed: $!";

# Get the statistics from vcf obtained
$input = $output;
system("
$path{bcftools} stats $input > $id.stats
") == 0 or die "error: Get vcf file statistics failed: $!";
#clean();

my $ts2 = time;
my $dur = ($ts2 - $ts1)/60/60;
printf "Execution time, hours: %.1f\n", $dur;

# ==subroutines==

sub clean {
  system("rm -r tmp*");
}

sub read_resources {
	my $res = $_[0];
	my %out;
	open(my $fh, '<', $res) or die "Failed to open $res: $!";
	while(<$fh>){
		chomp;
		my @ar = split /,/;
		$out{$ar[0]} = $ar[1];
	}
	close($fh);
	return %out;
}

sub help {
	print <<EOT;

Script subsets SNP and runs VQSR  

Usage: 	vqsr.pl [options] <vcf> prefix 
	vcf: path/to/my.vcf file
	prefix: prefix of the output files 

Options: 
	--help, -h	Show this help.

Author: Gennady Khvorykh.

EOT
	exit 0;
}
