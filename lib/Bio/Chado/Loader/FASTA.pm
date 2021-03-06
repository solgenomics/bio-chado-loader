package Bio::Chado::Loader::FASTA;

=head1 NAME

Bio::Chado::Loader::FASTA - Load FASTA sequences into a Chado schema

=head1 SYNOPSIS

mx-run Bio::Chado::Loader::FASTA [ options ] fasta_file fasta_file ...

=head1 DESCRIPTION

This is a MooseX::Runnable class that loads bioinformatics data given
in FASTA format into a Chado schema. You will have to run mx-run from 
the lib directory or include it in @INC using -I

=head2 Options

=cut

use Moose;
use Moose::Util::TypeConstraints;

use namespace::autoclean;

use POSIX;
use IO::Pipe;
use Digest::MD5 'md5_hex';
use Carp::Clan;

use autodie ':all';

with 'MooseX::Runnable';
with 'MooseX::Getopt';

has schema => (
    is         => 'rw',
    isa        => 'Bio::Chado::Schema',
    traits     => ['NoGetopt'],
    lazy_build => 1,
);

with 'Bio::Chado::Loader';

has '+organism_name' => (
	documentation =>
		  'exact species of organism, as stored in database, e.g. "Solanum lycopersicum" (required)',
	required => 1,
);

has '+db_attrs' => (
    traits => [ 'NoGetopt' ],
);

=head3 --create_features

Boolean flag.  If passed, autocreate features for each of the fasta sequences.

=cut

has 'create_features' => (
    documentation => <<'',
Default off.  If passed, autocreate features in the database.

    is => 'rw',
    isa => 'Bool',
   );

=head3 --type_name <string>

If creating features, the string SO term name to use for the type of these features.

=cut

has type_name => (
    documentation => <<'',
string Sequence Ontology (SO) term name to use for these sequences.  E.g. 'sequence_assembly' or 'EST'

    is       => 'ro',
    isa      => 'Str',

   );

=head3 --type_regex <string>

A regular expression that, when run against the ">ID defline" string,
will extract the feature's type.  Cannot specify both --type_name and
--type_regex.

=cut

{
  my $re = subtype as 'RegexpRef';
  coerce $re, from 'Str', via { eval "qr/$_/" or die "invalid regex '$_'\n" };
  MooseX::Getopt::OptionTypeMap->add_option_type_to_map( $re => '=s' );
  has type_regex => (
      documentation => <<'',
regular expression with a parentheses capture that will extract the feature's type name from the ">ID defline" string.

      is     => 'ro',
      isa    => $re,
      coerce => 1,
     );
}

sub BUILD  {
  my $self = shift;

  # check that either type_name or type_regex is set
  $self->type_name || $self->type_regex
      or die "must specify either a --type_name or --type_regex\n";
};

=head3 --analysis_name <string>

If creating features, pass this option to make an C<analysis> row with
C<programname> and optionally C<source> (see below), and connect the
features to it.

=cut

has analysis_name => (
    documentation => <<'',
connect features to this analysis (analysis.programname)

    is  => 'ro',
    isa => 'Str',

   );

=head3 --source <string>

If creating features, the string source name to use.

=cut

has source => (
    documentation => <<'',
string source name describing the source of these sequences. E.g. 'ITAG2'

    is  => 'ro',
    isa => 'Str',
   );

=head3 --large_residues <num>

Sequence length above which the residues will be stored as a
'large_residues' featureprop, instead of in the residues column.  Very
large residues values can cause performance issues when using
L<Bio::Chado::Schema>.  Default 4MB.

=cut

# NOTE: scientific notation OK for this
has large_residues => (
    documentation => <<'',
threshold (in bp) above which residues are stored in a 'large_residues' featureprop.  Default 4000000 (4M).

    is => 'ro',
    isa => 'Int',
    default => 4_000_000,
   );


has '_fixed_feature_type' => (
    is => 'ro',
    isa => 'Maybe[DBIx::Class::Row]',
    traits => ['NoGetopt'],
    lazy_build => 1,
   );
sub _build__fixed_feature_type {
    my ( $self ) = @_;

    return unless $self->type_name;

    my $cvt = $self->_get_so_term( $self->type_name )
      or die "SO term '".$self->type_name."' not found in database.\n";

    return $cvt;
}

sub _feature_type {
    my ( $self, $id, $defline ) = @_;
    return $self->_fixed_feature_type if $self->_fixed_feature_type;

    my ( $type_name ) = "$id $defline" =~ $self->type_regex
        or die "cannot parse type name from '$id $defline' using type regex "
               .$self->type_regex."\n";

    return $self->_get_so_term( $type_name );
}

use Memoize;
memoize('_get_so_term');
sub _get_so_term {
    my ( $self, $type_name ) = @_;
    return $self->schema->resultset('Cv::Cv')
                ->search({ 'me.name' => 'sequence' })
                ->search_related('cvterms', {
                    'cvterms.name'        => $type_name,
                    'cvterms.is_obsolete' => 0,
                   })
                ->single;
}

has organism => (
    is => 'ro',
    isa => 'DBIx::Class::Row',
    traits => ['NoGetopt'],
    lazy_build => 1,
   );
sub _build_organism {
    my ( $self ) = @_;
    my ( $genus, $species ) = split (/\s*_\s*/, $self->organism_name, 2);
    my $o = $self->schema
                 ->resultset('Organism::Organism')
                 ->find_or_create({ genus => $genus, species => $species })
        or die "Organism (genus: '$genus', species: '$species' not found, and could not create.\n";
    return $o;
}

sub _usage {

return <<USAGE;

At least one FASTA file must be given to load.

USAGE

}

sub run {
    my ( $self, @files ) = @_;

    die _usage() unless scalar(@files) > 0;

    print "Going to load " . @files . " FASTA files\n";
    my $stream = $self->fasta_stream( \@files );
    local $/ = "\n>";

    $self->schema->txn_do( sub {
        while( my $seq = <$stream> ) {
            my ( $id, $defline ) = $self->parse_fasta( \$seq );
            #print "Parsed $defline FASTA\n";
            $self->_find_or_create_feature( $id, $defline, \$seq );
        }
    });
}

sub _find_or_create_feature {
    my ( $self, $id, $defline, $seq ) = @_;

    my $seqlen = length $$seq;
    my $load_residues = $seqlen < $self->large_residues ? 1 : 0;

    my $feature_type = $self->_feature_type( $id, $defline );

    my @features =
        $feature_type
             ->search_related( 'features', {
                 'me.organism_id' => $self->organism->organism_id,
                 'me.name' => $id,
                });
    if( @features > 1 ) {
        die @features." features found for organism ".$self->organism->organism_id.", type_id ".$feature_type->cvterm_id.', name "'.$id."\", I cannot decide which to use.  Do you really want both of them in your database?\n";
    }

    my $feature = $features[0];

    if( $feature ) {
        print "Adding sequence to existing feature $id\n";
        print "Found $defline and updating\n";
        $feature->set_columns({
            seqlen      => $seqlen,
            md5checksum => md5_hex( $$seq ),
        });
        if( $load_residues ) {
            $feature->residues( $$seq );
        } else {
            $self->_load_large_residues_as_featureprop( $feature, $seq );
        }

        if( $self->analysis_name ) {
            $feature->find_or_create_related( 'analysisfeatures',{
                analysis => $self->_analysis,
            });
        }
    }
    elsif( $self->create_features ) {
        print "Creating feature $id\n";
        $feature =
            $feature_type
                 ->create_related( 'features', {
                     ( $load_residues ? ( residues => $$seq ) : () ),
                     organism    => $self->organism,
                     uniquename  => $id,
                     name        => $id,
                     seqlen      => $seqlen,
                     md5checksum => md5_hex( $$seq ),
                     is_analysis => $self->analysis_name ? 1 : 0,
                 });

        $self->_load_large_residues_as_featureprop( $feature, $seq )
            unless $load_residues;

        if( $self->analysis_name ) {
            $feature->add_to_analysisfeatures({
                analysis => $self->_analysis,
            });
        }
    }

    $feature or confess "Feature not found for organism '".$self->organism_name."' with id '$id' and type ".$feature_type->name." and could not create";
    $feature->update;
}

# if we did not load residues in the feature table, we will load
# them as a large_residues featureprop
sub _load_large_residues_as_featureprop {
    my ( $self, $feature, $seq ) = @_;

    my $large_seq_type =
        $self->schema->resultset('Cv::Cvterm')
             ->create_with({
                 name => 'large_residues',
                 cv   => 'feature_property',
                 db   => 'null',
                 dbxref => 'autocreated:large_residues',
             });

    $feature->delete_related('featureprops', {
        type_id  => $self->_large_seq_type->cvterm_id,
    });
    $feature->create_related('featureprops', {
        type_id => $self->_large_seq_type->cvterm_id,
        value   => $$seq,
    });

    $feature->residues( undef );
}

has '_analysis' => (
    is => 'ro',
    isa => 'DBIx::Class::Row',
    lazy_build => 1,
   ); sub _build__analysis {
       my ( $self ) = @_;

       $self->source or die "if setting --analysis_name, must provide a --source\n";

       $self->schema->resultset('Companalysis::Analysis')
            ->find_or_create({
                program        => $self->analysis_name,
                sourcename     => $self->source || $self->analysis_name,
                programversion => 'null',
              });
   }

has '_large_seq_type' => (
    documentation => <<'',
the cvterm for the 'large_residues' featureprop type

    is => 'ro',
    isa => 'DBIx::Class::Row',
    lazy_build => 1,
   ); sub _build__large_seq_type {
       my ( $self ) = @_;

       $self->schema->resultset('Cv::Cv')
            ->search({ 'me.name' => 'feature_property' })
            ->search_related('cvterms', {
                'cvterms.name' => 'large_residues',
                'cvterms.is_obsolete' => 0,
              })
            ->single

       ||

       $self->schema->resultset('Cv::Cvterm')
            ->create_with({
                name => 'large_residues',
                cv   => 'feature_property',
                db   => 'null',
                dbxref => 'autocreated:large_residues'
              });
   }

# in-place strips the id and defline off the seq and returns them
sub parse_fasta {
    my ( $self, $seq ) = @_;

    $$seq =~ s/^>?\s*(.+)\r?\n// or die "Error parsing fasta chunk '$$seq'";
    my ( $id, $defline ) = split /\s+/, $1, 2;
    $defline =~ s/^\s+|\s+$//g if $defline;
    $$seq =~ s/[\s>]//g;

    return ( $id, $defline );
}

sub fasta_stream {
    my ( $self, $files ) = @_;

    my $pipe = IO::Pipe->new;
    if( fork ) {

        $pipe->reader;
        return $pipe;

    } else {
        # fasta input spooler child

        $pipe->writer;

        # if we've got files, spool those
        if( $files && @$files ) {
            for my $file (@$files) {
                open my $f, '<', $file;
                $pipe->print( $_ ) while <$f>;
            }
        }
        # otherwise, read from stdin
        else {
            $pipe->print( $_ ) while <STDIN>;
        }

        $pipe->close;
        POSIX::_exit(0);
    }
}

1;
