#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long::Descriptive;
use GenOO::Data::File::FASTA;

# Define and read command line options
my ($opt, $usage) = describe_options(
	'%c %o',
	['fasta=s', 'FASTA file with sequences (default: STDIN)'],
	['name-suffix=s', 'suffix to add to name field'],
	['help|h', 'print usage message and exit'],
);
print($usage->text), exit if $opt->help;

my $fp = GenOO::Data::File::FASTA::->new(file => $opt->fasta);
while (my $r = $fp->next_record) {
	my $name = $r->header;
	if ($opt->name_suffix) {
		$name = $name . ':' . $opt->name_suffix;
	}
	print join("\t",
		$r->header,
		0,
		length($r->sequence),
		$name,
		length($r->sequence),
		'+'
	)."\n";
}

