#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Replaces whitespace characters in FASTA headers with given delimiter."],
	[],
	['input=s', 'input FASTA file; reads from STDIN if -', {required => 1}],
	['delim=s', 'delimiter that replaces whitespaces', {required => 1}],
	['keep=i@', 'index of the header fields to keep (optional)'],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $delim = $opt->delim;

my $FASTA = filehandle_for($opt->input);
while (my $l = $FASTA->getline) {
	if ($l !~ /^>/) {
		print $l;
		next;
	}
	chomp $l;
	my @fields = split(/\s+/, $l);
	if (defined $opt->keep) {
		@fields = @fields[@{$opt->keep}];
	}
	my $name = join($delim, @fields);
	print $name . "\n";
}

exit;

sub filehandle_for {
	my ($file) = @_;

	if ($file eq '-'){
		return IO::File->new("<-");
	}
	else {
		return IO::File->new($file, "<");
	}
}

