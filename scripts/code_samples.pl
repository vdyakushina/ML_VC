#!/bin/perl
use strict; use warnings;
use File::Find;
use File::Basename;
use Data::Dumper;
use Data::Uniqid qw {uniqid};

my $input = $ARGV[0];
my %data;
my $meta = "/home/gkhvorykh/samples/meta.csv";
my @samples;

open (my $fh1, '<', "$input") or die "Failed to open $input: $!";
while(<$fh1>) {
	chomp;
	my @ar = split /,/;
	$data{$ar[0]} = $ar[1];
}

open(my $fh2, '>', "$meta") or die "Failed to open $meta: $1";

foreach my $folder (keys %data) {
	my $type = $data{$folder};
	my $project = basename($folder);
	# Find all fastq files in the folder
	find(\&proc, $folder);
	# Output the list of samples with meta
	foreach my $sample (@samples) {
		my $id = uniqid;
		print $fh2 "$id,$type,$project,$sample\n";
	}
	@samples = ();
}


# Parse sample id from each R1_001.fastq.gz file
sub proc {
	my $file  = $File::Find::name;
	my $id;
	return unless $file =~ /R1_001.fastq.gz$/;
	# Get sample id
	$file =~ /\/([a-zA-Z0-9_-]*)_L/;
	$id = $1;
	return if $id =~ /ndetermined/;
	push @samples, $id;
}

