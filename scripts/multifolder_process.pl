#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Capture::Tiny 'capture_merged';

# Initiate
my $input = $ARGV[0];
my %path = (
	scripts => "/pipeline/scripts"
);
my %data;

# Check input
die "$input doesn't exist!" unless -d $input;

# Remove / at the end if any
chop $input if "/" eq substr $input, -1;

# Show input
print "The folder: $input\n";  

# Get list of subfolders
my @folders = glob("$input/*");

# Report the number of subfolders found
my $size = @folders;
print "Subfolders found: $size\n";

# Get list of samples
print "N\tSample\tFiles\n";
my $n = 1;
foreach (@folders){
	$_ =~ /\/([a-zA-Z0-9]*)$/;
	my $sample = $1;
	my @files = glob("$_/*");
	my $size = @files;
	print "$n\t$sample\t$size\n";
	push @{$data{$sample}}, @files if $size == 2;
	$n++;
}

foreach my $sample (keys %data){
	print "Processing $sample\n";
	my $log = "$sample.process.log";
	my @files = $data{$sample};
	open (my $fh, '>', $log) or die "Couldn't open $log";
	my $cmd = "$path{scripts}/process.py $files[0]->[0] $files[0]->[1] $sample";
	print $fh capture_merged {system("$cmd")};
	close($fh) ;
}
