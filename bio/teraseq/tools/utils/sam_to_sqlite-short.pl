#!/usr/bin/env perl
#
# use GenOO::Data::DB::DBIC::Species::Schema;
# my $schema = GenOO::Data::DB::DBIC::Species::Schema->connect(
# 	"dbi:SQLite:database=foo.db");
# $schema->_create_and_register_result_class_for(
# 	'sample', 'GenOOx::DBIC::ResultBase::SAM');
# $schema->deploy( { add_drop_table => 1 } );

=head1 NAME

CLIPSeqTools::PreprocessApp::sam_to_sqlite - Load a SAM file in an SQLite database.

=head1 SYNOPSIS

clipseqtools-preprocess sam_to_sqlite [options/parameters]

=head1 DESCRIPTION

Store alignments from a SAM file into and SQLite database. If SAM tag XC:i exists it will be used as the copy number of the record.

=head1 OPTIONS

  Input options.
    --sam_file <Str>       sam file to be stored in database. If not
                           specified STDIN is used.
    --records_class <Str>  type of records stored in SAM file. [Default:
                           GenOOx::Data::File::SAMstar::Record]

  Database options.
    --database <Str>       database name or path. Will be created.
    --table <Str>          database table. Will be created.
    --drop <Str>           drop table if it already exists.

  Other options.
    -v --verbose           print progress lines and extra information.
    -h -? --usage --help   print help message

=cut

package sam_to_sqlite;


# Make it an app command
use MooseX::App::Simple;
use MooseX::Getopt::Meta::Attribute::Trait::NoGetopt;


#######################################################################
#######################   Load External modules   #####################
#######################################################################
use Modern::Perl;
use autodie;
use DBI;


##############################################
# Import GenOO
use GenOO::Data::File::SAM;


#######################################################################
#######################   Command line options   ######################
#######################################################################
option 'sam_file' => (
	is            => 'rw',
	isa           => 'Str',
	documentation => 'sam file to be stored in database. If not specified STDIN is used.',
);

option 'records_class' => (
	is            => 'rw',
	isa           => 'Str',
	default       => 'GenOOx::Data::File::SAMstar::Record',
	documentation => 'type of records stored in SAM file.',
);

option 'database' => (
	is            => 'rw',
	isa           => 'Str',
	required      => 1,
	documentation => 'database name or path. Will be created.',
);

option 'table' => (
	is            => 'rw',
	isa           => 'Str',
	required      => 1,
	documentation => 'database table. Will be created.',
);

option 'drop' => (
	is            => 'rw',
	isa           => 'Bool',
	documentation => 'drop table if it already exists.',
);

option 'verbose' => (
	is            => 'rw',
	isa           => 'Bool',
	cmd_aliases   => 'v',
	default       => 0,
	documentation => 'print progress lines and extra information.',
);


#######################################################################
########################   Interface Methods   ########################
#######################################################################
sub validate_args {}

sub run {
	my ($self) = @_;

	warn "Validating arguments\n" if $self->verbose;
	$self->validate_args();

	# Load required classes
	eval 'require ' . $self->records_class;

	warn "Reading SAM\n" if $self->verbose;
	my $sam = GenOO::Data::File::SAM->new(
		file          => $self->sam_file,
		records_class => $self->records_class,
	);

	warn "Connecting to the database\n" if $self->verbose;
	my $dbh = DBI->connect('dbi:SQLite:database=' . $self->database) or die "Can't connect to database: $DBI::errstr\n";

	if ($self->drop) {
		warn "Dropping table " . $self->table . "\n" if $self->verbose;
		$dbh->do( q{DROP TABLE IF EXISTS } . $self->table );
	}

	warn "Creating table " . $self->table . "\n" if $self->verbose;
	{
		local $dbh->{PrintError} = 0; #temporarily suppress the warning in case table already exists

		$dbh->do(
			'CREATE TABLE '.$self->table.' ('.
				'id INTEGER PRIMARY KEY AUTOINCREMENT,'.
				'qname VARCHAR(250) NOT NULL,'.
				'flag UNSIGNED INT(10) NOT NULL,'.
				'rname VARCHAR(250) NOT NULL,'.
				'pos UNSIGNED INT(10) NOT NULL,'.
				'mapq UNSIGNED INT(10) NOT NULL,'.
				'cigar VARCHAR(250) NOT NULL,'.
				'rnext VARCHAR(250) NOT NULL,'.
				'pnext UNSIGNED INT(10) NOT NULL,'.
				'tlen UNSIGNED INT(4) NOT NULL,'.
				'seq VARCHAR(250) NOT NULL,'.
				'qual VARCHAR(250) NOT NULL,'.
				'tags TEXT,'.
				'strand INT(1) NOT NULL,'.
				'start UNSIGNED INT(10) NOT NULL,'.
				'stop UNSIGNED INT(10) NOT NULL,'.
				'copy_number UNSIGNED INT(6) NOT NULL DEFAULT 1,'.
				'sequence VARCHAR(250) NOT NULL,'.
				'mdz VARCHAR(250),'.
				'number_of_mappings UNSIGNED INT(5),'.
				'query_length UNSIGNED INT(4) NOT NULL,'.
				'alignment_length UNSIGNED INT(5) NOT NULL'.
			');'
		);

		if ($dbh->err) {
			die "Error: " . $dbh->errstr . "\n";
		}
	}


	warn "Loading data to table " . $self->table . "\n" if $self->verbose;
	$dbh->begin_work;
	my $insert_statement = $dbh->prepare(
		q{INSERT INTO } . $self->table . q{ (id, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags, strand, start, stop, copy_number, sequence, mdz, number_of_mappings, query_length, alignment_length) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)}
	);
	while (my $r = $sam->next_record) {
		my $copy_number = $r->copy_number;

		if ($r->is_unmapped) {
			next;
		}

		if (defined $r->tag('XC:i')) {
			$copy_number = $r->tag('XC:i');
		}

		my $number_of_mappings = $r->number_of_mappings;
		if (!defined $number_of_mappings) {
			$number_of_mappings = 1;
		}

		my $tags = join("\t", ($r->all_fields)[11..$r->count_fields-1]);
		$insert_statement->execute(undef, $r->qname, $r->flag, $r->rname,
			$r->pos, $r->mapq, $r->cigar, $r->rnext, $r->pnext, $r->tlen,
			$r->seq, $r->qual, $tags, $r->strand, $r->start,
			$r->stop, $copy_number, $r->query_seq, $r->mdz,
			$number_of_mappings, $r->query_length, $r->alignment_length);

		if ($sam->records_read_count % 100000 == 0) {
			$dbh->commit;
			$dbh->begin_work;
		}
	}
	$dbh->commit;

	warn "Building index on " . $self->table . "\n" if $self->verbose;
	$dbh->do(q{CREATE INDEX } . $self->table . q{_loc ON } . $self->table .q{ (rname, start);});

	warn "Disconnecting from the database\n" if $self->verbose;
	$dbh->disconnect;
}

#######################################################################
########################   Make it executable   #######################
#######################################################################
sam_to_sqlite->new_with_options->run();

1;
