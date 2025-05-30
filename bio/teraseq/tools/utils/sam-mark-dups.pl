#!/usr/bin/env perl

use Modern::Perl;
use Getopt::Long::Descriptive;
use IO::Interactive;
use GenOO::Data::File::SAM;


# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Mark duplicate records in a SAM file. The input file must be sorted by coordinates (rname and pos). Secondary or unmapped records are not processed and are written in the output as is. Duplicate reads are defined as those that align on the same reference and at the same start and stop position. For paired end reads the start is defined as the leftmost position of the pair and the stop as the rightmost position. Duplicates are marked by setting the 0x0400 (1024) bit in the flag. If a barcode tag is provided then reads with different barcodes are treated independently."],
	[],
	['input=s', 'input SAM file; reads from STDIN if empty'],
	['barcode-tag=s', 'barcode tag; records with different barcodes are treated independently'],
	['primary-tag=s', 'primary tag; if defined, a tag is added to duplicates with the name of the corresponding primary record qname'],
	['paired', 'run in paired-end mode'],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $pairedMode = $opt->paired;
my $barcodeTag = $opt->barcode_tag;
my $primaryTag = $opt->primary_tag;

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
my @toMark;
my %markedDup;
my ($prev_start, $prev_rname) = ('', '');
while (my $r = $p->next_record) {
	if ($r->is_unmapped) {
		say $r->to_string();
		next;
	}

	# accumulate records with same rname and start.
	if ($r->rname eq $prev_rname and $r->start == $prev_start) {
		push @toMark, $r;
	} else {
		mark_duplicates(\@toMark, \%markedDup, $pairedMode, $barcodeTag, $primaryTag);
		map{say $_->to_string()} @toMark;

		@toMark = ($r);
		$prev_rname = $r->rname;
		$prev_start = $r->start;
	}
}
mark_duplicates(\@toMark, \%markedDup, $pairedMode, $barcodeTag, $primaryTag);
map{say $_->to_string()} @toMark;

exit;

###########################################################################
sub mark_duplicates {
	my ($toMark, $markedDup, $pairedMode, $barcodeTag, $primaryTag) = @_;

	my %group;
	foreach my $r (@toMark) {
		if ($r->is_secondary) {
			next;
		}
		my $key = extract_key($r, $pairedMode, $barcodeTag);
		push @{$group{$key}}, $r;
	}

	# mark duplicates in each group.
	foreach my $key (keys %group) {
		my @toResolve = ();
		my $records = $group{$key};
		foreach my $r (@$records) {
			# if we have seen a record with the same name before.
			if (exists $$markedDup{$r->qname}) {
				# and has been marked as duplicate
				if ($$markedDup{$r->qname}) {
					mark_as_duplicate($r);
					if($primaryTag) {
						$r->add_tag($primaryTag, $$markedDup{$r->qname});
					}
				}
				delete $$markedDup{$r->qname};
				next;
			}
			push @toResolve, $r;
		}
		my $random_r = $toResolve[int(rand(@toResolve))];

		foreach my $r (@toResolve) {
			if ($r == $random_r) {
				if ($pairedMode) {
					$$markedDup{$r->qname} = 0;
				}
			} else {
				mark_as_duplicate($r);
				if($primaryTag) {
					$r->add_tag($primaryTag, $random_r->qname);
				}
				if ($pairedMode) {
					$$markedDup{$r->qname} = $random_r->qname;
				}
			}
		}
	}
}

sub extract_key {
	my ($r, $pairedMode, $barcodeTag) = @_;

	my @keyCols = ($r->rname, $r->strand, $r->start, $r->stop);
	if ($pairedMode) {
		push @keyCols, $r->rnext, $r->pnext, $r->tlen;
	}
	if (defined $barcodeTag) {
		my $barcode = $r->tag($barcodeTag);
		if (!defined $barcode) {
			die "tag " . $barcodeTag . " not found for \"" . $r->to_string() . "\"\n";
		}
		push @keyCols, $barcode;
	}
	return join("|", @keyCols);
}

sub mark_as_duplicate {
	my ($r) = @_;
	
	my $flag = $r->flag;
	if (!($flag & 1024)) {
		$r->set_flag($flag + 1024);
	}
}
