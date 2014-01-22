package Bio::Chado::Loader::GFF3;
use strict;
use warnings;

=head1 NAME

Bio::Chado::Loader::GFF3 - Load GFF3 files into the Chado database schema

=head1 SYNOPSIS

  use Bio::Chado::Loader::GFF3;
  mx-run Bio::Chado::Loader::GFF3 [ options ] gff3_file gff3_file ...


=head1 DESCRIPTION

This is a MooseX::Runnable class that loads feature data given in GFF3 format into a Chado schema.
The GFF3 spec is available online at L<http://www.sequenceontology.org/gff3.shtml>

=over

=cut


use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
use Proc::ProcessTable;


has 'schema' => (
    is  => 'rw',
    isa => 'Bio::Chado::Schema',
    lazy => 1,
    builder => '_build_schema'
);

has 'file_name' => (
	isa       => 'Str',
	is        => 'rw',
	predicate => 'has_file_name',
	clearer   => 'clear_file_name'
);

has organism_name => (
    documentation => 'exact species of organism, as stored in database, e.g. "Solanum lycopersicum"',

    is       => 'rw',
    isa      => 'Str',
);

has organism_id => (
    documentation => 'id organism, as stored in database, e.g. 1 for Solanum lycopersicum',

    is       => 'rw',
    isa      => 'Str',
);

has 'pragma_lines' => (
	isa       => 'ArrayRef[Str]',
	is        => 'rw',
	predicate => 'has_pragma_lines',
	clearer   => 'clear_pragma_lines'
);

has 'features' => (
	documentation => 'hash of featureID or seqID => array of entire GFF record. TOFIX Partly redundant',
    is      => 'rw',
    #isa     => 'HashRef[Str]',
    isa     => 'HashRef[ArrayRef]',##Use constants to get arr member
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_feature    => 'set',
        count_features => 'count',
        feature_exists => 'exists',
    },
    );
    
has 'cvterms' => (
	documentation => 'hash of cvterm => 1. TOFIX Use counts instead of constant 1',
    is      => 'rw',
    isa     => 'HashRef[Str]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_cvterm    => 'set',
        count_cvterms => 'count',
    },
);

has 'is_analysis' => (
    documentation => <<'',
set true if this feature should be recorded as from an analysis

    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);

has 'cache' => (
	documentation => 'hash of cvterm.name => hash of feature.uniquename => feature.feature_id,featureloc.srcfeature_id,featureloc.locgroup,featureloc.rank',
    is      => 'rw',
    isa     => 'HashRef',
    default => sub { { } },
#    traits  => [ 'Hash' ],
    default => sub { { } },
#    handles => {
#        cache_exists=> 'exists',## Needed???
#        add_cache   => 'set',
#        count_cache => 'count',## Needed???
#    },
);

use constant SEQID         => 0;
use constant SOURCE        => 1;
use constant TYPE          => 2;
use constant FEATURE_START => 3;
use constant FEATURE_END   => 4;
use constant SCORE         => 5;
use constant STRAND        => 6;
use constant PHASE         => 7;
use constant ATTRIBUTES    => 8;

with 'Bio::Chado::Loader';

use autodie qw(:all);
use Data::Dumper;
use 5.010;

=item C<run ()>

TODO // See Fasta.pm

=cut

sub run {
    my ($self, %args) = @_;

    exit 0;
}

=item C<parse ()>

Parse a GFF3 file. Calls parse_line() for each line.

=cut

sub parse {
    my ($self, %args) = @_;
    open my $fh, "<", $self->file_name;
    while( my $line = <$fh> ) {
        next if $line =~ m/^\s*$/;
        chomp $line;
        next if $self->is_pragma_line($line);
        $self->parse_line($line);
    }
    #warn Dumper [ $self->features ];
}

=item C<parse_line ()>

Parse a GFF3 line to get feature details.

=cut

sub parse_line {
    my ($self, $line) = @_;
    my @fields = split (/\t/, $line);
    my $attribs = { map { split (/=/, $_)  }  split (/;/, $fields[ATTRIBUTES]) };

    $self->add_cvterm( $fields[TYPE]   => 1 );
    #$self->add_feature( $attribs->{ID} || $fields[SEQID] => 1 );
    $self->add_feature( $attribs->{ID} || $fields[SEQID] => \@fields );

    $self->_validate_parents($attribs);
}

=item C<_validate_parents ()>

Check if parent feature exists

=cut


sub _validate_parents {
    my ($self, $attribs) = @_;
    return unless $attribs->{Parent};
    my (@parents) = split (/,/, $attribs->{Parent});
    for my $p (@parents) {
        unless ( $self->feature_exists( $p ) ) {
            die "$p is an unknown Parent!";
        }
    }
}

=item C<is_pragma_line ()>

Check if line is a pragma

=cut

sub is_pragma_line {
    my ($self, $line) = @_;
    return $line =~ m/^##/ ? 1 : 0;
}

=item C<organism_exists ()>

Check if organism for GFF3 exists in the Chado database. Return organism::organism.organism_id 
if present otherwise 0. No organism should have id of 0.

=cut

sub organism_exists {
    my ($self) = @_;
    
    my $org_id = $self->schema->resultset('Organism::Organism')-> search (
    	{ species => $self->organism_name})->get_column('organism_id')->first();
    return $org_id ? $org_id : 0;
}

=item C<populate_cache ()>

Presuming organism_Populate cache hash with feature.uniquename=feature.feature_id,featureloc.srcfeature_id,
featureloc.locgroup,featureloc.rank,cvterm.name from cvterm,feature,featureloc relations. Return number of 
records added to cache. Count may be over estimation if some features have multiple locgroups(e.g. contigs).

=cut

sub populate_cache {
    my ($self) = @_;
    
    #setup DB dsn's
    my $fl_rs = $self->schema -> resultset('Sequence::Featureloc')->search(
    	{ 'organism_id'=> $self->organism_id },
    	{ join => [ 'feature' ] , prefetch=> [ 'feature']}
    	);

	my $ft_rs = $self->schema -> resultset('Sequence::Feature')->search(
		{ 'organism_id'=> $self->organism_id },
		{ join => [ 'type' ] , prefetch=> [ 'type']}
		);
		
    print STDERR "created resultsets\n\n";
    print STDERR "fl_rs count: ".$fl_rs->count()."\n";
    print STDERR "ft_rs count: ".$ft_rs->count()."\n";
    
    #create cache hash
    my $count = 0;
    while (my $fl_row = $fl_rs->next()){
		if($ft_rs->search(
			{'type_id' => $fl_row->feature->type_id(), 'uniquename' => $fl_row->feature->uniquename}
			)->count() == 1){
			my $ft_row = $ft_rs->search(
					{'type_id' => $fl_row->feature->type_id(), 'uniquename' => $fl_row->feature->uniquename}
					)->first();

#			$self->add_cache( {$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{feature_id} => $fl_row->feature->feature_id() );
#			$self->add_cache( {$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{locgroup} => $fl_row->locgroup() );
#			$self->add_cache( {$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{rank} => $fl_row->rank() );
#			$self->add_cache( {$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{scrcfeature_id} => $fl_row->srcfeature_id() );
			
			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{feature_id}=$fl_row->feature->feature_id();
			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{locgroup}=$fl_row->locgroup();
			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{rank}=$fl_row->rank();
			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{scrcfeature_id}=$fl_row->srcfeature_id();
			$count++;	
			if ($count % 10000 == 0) {
				print STDERR "processing $count\r";
				warn Dumper [ $self-> cache -> {$ft_row->type->name()}];
			}
		}
		else{
			die "Multiple features found for organism with same type_name and uniquename which violates Chado constraints. Exiting..";
		}
	}

    my ( $i, $t );
	$t = new Proc::ProcessTable;
	foreach my $got ( @{ $t->table } ) {
		next if not $got->pid eq $$;
		$i = $got->size;
	}
	print STDERR "Memory used: ".($i / 1024 / 1024)."\n";
	return ($count);
}

=item C<process_features ()>

Compare %features to %cache. Push content to write to disk for direct copy to DB. %features not found 
in %cache are written to exceptions.gff3. New featureloc records have locgroup=0 (primary location). Old 
locgroups are incremented by 1.

=cut

sub process_features {
    my ($self) = @_;
    

    
}

=item C<insert_features ($insert_data)>

Insert content from disk intermediate file into DB   

=cut

sub insert_features {
    my ($self) = @_;
    
    
    
}


###
1;   #do not remove
###

=pod

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Jonathan "Duke" Leto	<jonathan at leto.net>
    Surya Saha				<suryasaha at cornell.edu , @SahaSurya>   

=cut