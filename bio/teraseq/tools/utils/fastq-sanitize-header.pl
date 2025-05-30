#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Replaces whitespace characters in FASTQ headers with given delimiter."],
	[],
	['input=s', 'input FASTQ file; reads from STDIN if -', {required => 1}],
	['delim=s', 'delimiter that replaces whitespaces', {required => 1}],
	['keep=i@', 'index of the header fields to keep (optional)'],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $delim = $opt->delim;

my $FASTQ = filehandle_for($opt->input);
while (my $name = $FASTQ->getline) {
	if ($name !~ /^\@/) {
		next;
	}
	chomp $name;
	my @fields = split(/\s+/, $name);
	if (defined $opt->keep) {
		@fields = @fields[@{$opt->keep}];
	}
	$name = join($delim, @fields);
	print $name . "\n";

	print $FASTQ->getline;

	my $second_name = $FASTQ->getline;
	chomp($second_name);
	$second_name =~ s/\s+/$delim/g;
	print $second_name . "\n";

	print $FASTQ->getline;
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

