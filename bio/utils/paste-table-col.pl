#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Print the input table adding a new column with a constant value."],
	[],
	['ifile=s',
		'Input file name. Reads from STDIN if -', {required => 1}],
	['col-name=s',
		'Name for new column. Assumes input file already has a header line'],
	['col-val=s',
		'Value for new column'],
	['help|h',
		'Print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $delim = "\t";

my $IN = filehandle_for($opt->ifile);

if (defined $opt->col_name) {
	my $header = $IN->getline();
	print $opt->col_name . $delim . $header;
}

while (my $line = $IN->getline) {
	print $opt->col_val . $delim . $line;
}
$IN->close();

exit;

sub filehandle_for {
	my ($file) = @_;

	if ($file eq '-'){
		open(my $IN, "<-");
		return $IN;
	}
	else {
		open(my $IN, "<", $file);
		return $IN
	}
}
