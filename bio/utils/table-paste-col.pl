#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use IO::File;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Print the input table adding a new column with a constant value. Assumes that the input file has a header line with column names"],
	[],
	['table=s', 'Input table file; reads from STDIN if "-"'],
	['ifile=s', 'Input file name; reads from STDIN if "-"', {hidden => 1}], # for backwards compatibility
	['col-name=s', 'Name for new column'],
	['col-val=s', 'Value for new column', {required => 1}],
	['sep=s', 'Column separator character [Default => "\t"]', {default => "\t"}],
	['at-end', 'Inserts new column as the right most one'],
	['verbose|v', 'Print progress'],
	['help|h', 'Print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

# For backwards compatibility: handles the deprecated ifile option.
if (defined $opt->ifile) {
	if (defined $opt->table) {
		die "Cannot specify both \'table\' and \'ifile\' parameters\n";
	}
	$opt->{'table'} = $opt->ifile
}
if (!defined $opt->table) {
	die "Mandatory parameter \'table\' missing\n";
}

my $sep = $opt->sep;

my $IN = filehandle_for($opt->table);

if (defined $opt->col_name) {
	my $header = $IN->getline();
	chomp($header);
	if ($opt->at_end) {
		print $header . $sep . $opt->col_name . "\n";
	} else {
		print $opt->col_name . $sep . $header . "\n";
	}
}

while (my $line = $IN->getline) {
	chomp($line);
	if ($opt->at_end) {
		print $line . $sep . $opt->col_val . "\n";
	} else {
		print $opt->col_val . $sep . $line . "\n";
	}
}
$IN->close();

exit;

sub filehandle_for {
	my ($file) = @_;

	my $f;
	if ($file eq '-'){
		$f = IO::File->new("<-") or die "cannot open stdin: $!\n";
	}
	else {
		$f = IO::File->new($file, "<") or die "cannot open file $file: $!\n";
	}
	return $f;
}
