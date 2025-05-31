#!/usr/bin/env perl

use Modern::Perl;
use Getopt::Long::Descriptive;
use IO::Interactive;
use GenOO::Data::File::SAM;


# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Add a new tag in every primary (non-duplicate) SAM record with the number of corresponding duplicate records plus one for the SAM record itself. This implies than primary records with no duplicates will have value of 1. The program assumes that each duplicate record is marked with the 0x0400 (1024) flag bit and that it contains a tag with the corresponding primary record qname."],
	[],
	['input=s', 'input SAM file; reads from STDIN if empty'],
	['count-tag=s', 'tag to output the duplicates count for primary records (Default: XC:i)', {default => 'XC:i'}],
	['primary-tag=s', 'the tag with the primary record qname for each duplicate', {required => 1}],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $primaryTag = $opt->primary_tag;
my $countTag = $opt->count_tag;

# abort if input option is empty and data are not coming from a pipe.
if (!defined $opt->input and IO::Interactive::is_interactive()) {
	say "Error: No input data.";
	print($usage->text);
	exit 1;
}

# open the SAM file.
my $p = GenOO::Data::File::SAM->new(file => $opt->input);

# print the header as is.
my $header = $p->header;
if (defined $header) {
	say $header;
}

# mark and count duplicates.
my @toCount;
my ($prev_start, $prev_rname) = ('', '');
while (my $r = $p->next_record) {
	if ($r->is_unmapped) {
		say $r->to_string();
		next;
	}

	# accumulate records with same rname and start.
	if ($r->rname eq $prev_rname and $r->start == $prev_start) {
		push @toCount, $r;
	} else {
		count_duplicates(\@toCount, $countTag, $primaryTag);
		map{say $_->to_string()} @toCount;

		@toCount = ($r);
		$prev_rname = $r->rname;
		$prev_start = $r->start;
	}
}
count_duplicates(\@toCount, $countTag, $primaryTag);
map{say $_->to_string()} @toCount;

exit;

###########################################################################
sub count_duplicates {
	my ($toCount, $collapsingTag, $primaryTag) = @_;

	my %count;
	foreach my $r (@toCount) {
		if ($r->is_secondary) {
			next;
		}
		if ($r->flag & 1024) {
			my $primaryTagValue = $r->tag($primaryTag);
			$count{$primaryTagValue}++;
		}
	}
	foreach my $r (@toCount) {
		if ($r->is_secondary) {
			next;
		}
		if ($r->flag & 1024) {
			next;
		}
		my $dupCount = $count{$r->qname} || 0;
		$r->add_tag($collapsingTag, $dupCount+1);
	}
}
