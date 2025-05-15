#!/usr/bin/env perl

use Modern::Perl;
use Getopt::Long::Descriptive;
use IO::File;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o [FILE...]",
	["Appends one table after the other keeping only a single header."],
	[],
	['help|h', 'Print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my @files = @ARGV;

my $header;
foreach my $f (@files) {
	my $p = IO::File->new($f, "<") or die "error for $f: $!\n";
	my $h = $p->getline();
	if (!defined $header) {
		$header = $h;
		print $header;
	}
	if ($h ne $header) {
		die "error: different header for $f";
	}
	while (my $l = $p->getline()) {
		print $l;
	}
	close($p)
}
