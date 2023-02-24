use strict;
use warnings;
use XML::Simple;
use LWP::UserAgent;
use Data::Dumper;
use String::Util qw(trim);
my $ua = new LWP::UserAgent;


sub getSeq {
	my $build = shift;
	my $chr = shift;
	my $pos1 = shift;
	my $pos2 = shift;
	my $try = 0;
	FETCH:
#	print "http://genome.ucsc.edu/cgi-bin/das/$build/dna?segment=$chr:$pos1,$pos2\n";
	my $response = $ua->get("http://genome.ucsc.edu/cgi-bin/das/$build/dna?segment=$chr:$pos1,$pos2");
	unless ($response && $response->is_success) {
		++$try;
		print STDERR "try number $try\n";
		goto FETCH;
		}
	my $xmlString = $response->content;
	my @options = ();
	my $ref = XMLin($xmlString, @options);
	my $return = trim((((($ref)->{"SEQUENCE"})->{"DNA"})->{"content"}));
	$return =~ s/\n//g;
	return uc($return);
	}

open (READ, "<$ARGV[0]");

while (<READ>) {
	chomp;
	my @line  = split/\t/;
	my $chr = $line[0];
	my $pos = $line[1];
	my $ref = uc($line[3]);
	my $alt = uc($line[5]);
	my $dp	= $line[4];
	my $adf	= $line[8];
	my $adr	= $line[9];
	if ($adf =~ /(\d+)/) {$adf = $1}
	if ($adr =~ /(\d+)/) {$adr = $1}
	my $info = "DP=$dp;AD=$adf,$adr;AF=".($adf + $adr)/$dp;
	if ($ref =~ /\+/) {die};
	if ($ref =~ /\-/) {die};
	if ($alt =~ /^\+(\S+)$/) {
		print "$chr\t$pos\t.\t$ref\t$ref$1\t.\t.\t$info\n";
		} elsif ($alt =~ /^\-(\S+)$/) {
			my $refAct = uc(getSeq('hg19', $chr, $pos-1, $pos-1));
			print "$chr\t",$pos-1,"\t.\t$refAct$1\t$refAct\t.\t.\t$info\n";
			} else {
				print "$chr\t$pos\t.\t$ref\t$alt\t.\t.\t$info\n";
				}
	}

close READ;














