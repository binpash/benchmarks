#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long::Descriptive;
use GenOO::Data::File::FASTA;

my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Removes entries with the same sequence. The first occurence is retained"],
	[],
	['fasta=s', 'fasta file. Reads from STDIN if not provided'],
	['help|h', 'Print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $fp = GenOO::Data::File::FASTA->new(
	file => $opt->fasta
);

my %seqs;
while (my $r = $fp->next_record){
	my $seq = $r->sequence;
	if (exists $seqs{$seq}) {
		next;
	}
	$seqs{$seq} = 1;
	print ">" . $r->header . "\n" . $seq . "\n";
}
