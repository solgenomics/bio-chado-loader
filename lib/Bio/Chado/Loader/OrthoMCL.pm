package Bio::Chado::Loader::OrthoMCL;

=head1 NAME

Bio::Chado::Loader::OrthoMCL - Load OrthoMCL output into a Chado database

=head1 SYNOPSIS

mx-run Bio::Chado::Loader::OrthoMCL [ options ] orthomcl.out orthomcl.out ...

=head1 DESCRIPTION

This is a MooseX::Runnable class that loads OrthoMCL groups into a
Chado schema as features and feature_relationships.  Makes a feature
for the group, and then adds feature_relationships with type
C<member_of> between the group and its members.

=head2 Options

=cut

use Moose;

with (
    'MooseX::Runnable',
    'MooseX::Getopt',
    'MooseX::Role::DBIC',
);

sub schema { shift->dbic_schema }

has "+$_" => ( traits => ['NoGetopt'] )
    for qw( dbic_schema dbic_schema_options dbic_attrs );

has '+dbic_class' => (
    traits => ['NoGetopt'],
    default => 'Bio::Chado::Schema',
    );

has 'create_members' => (
    is => 'ro',
    isa => 'Bool',
    documentation => 'flag for creating new features (holding very little data) for group members that are not already in the database',
    default => 0,
   );

has 'group_type' => (
  is => 'ro',
  documentation => 'SO term name for the type of the group features to make.  default "gene_group"',
  default => 'gene_group',
 );

has '_group_term' => (
  is => 'ro',
  traits => ['NoGetopt'],
  lazy_build => 1,
);
sub _build__group_term {
    my ( $self ) = @_;
    $self->_find_cvterm( sequence => $self->group_type );
}



has 'member_type' => (
    is => 'ro',
    documentation => 'SO term name for the type of the member features.  default "gene"',
    default => 'gene',
  );

has '_member_term' => (
    is => 'ro',
    traits => ['NoGetopt'],
    lazy_build => 1,
    );

sub _build__member_term {
    my ( $self ) = @_;
    $self->_find_cvterm( sequence => $self->member_type );
}

sub _find_cvterm {
    my ( $self, $cv, $term ) = @_;
    my $t = $self->schema->resultset('Cv::Cv')
                 ->search( { 'me.name' => $cv })
                 ->search_related( 'cvterms', { 'cvterms.name' => $term } )
                 ->single;
    $t or die "Cannot find required '$term' term in $cv ontology.  Does this term exist?  Is the $cv ontology loaded?";
    return $t;
}

has '_member_of_term' => (
  is => 'ro',
  traits => ['NoGetopt'],
  lazy_build => 1,
);
sub _build__member_of_term {
    shift->_find_cvterm( sequence => 'member_of' );
}


sub run {
    my ( $self, @ortho_files ) = @_;

    # do everything in a single transaction
    $self->schema->txn_do( sub {

        for my $file ( @ortho_files ) {
            my $fh = IO::File->new( $file, '<' );
            while( <$fh> ) {
                my $group = $self->parse_orthomcl_group_line( $_ );
		print "loading $group->{group}{id} ($group->{group}{genes} genes)\n";
                $self->load_group( $group );
            }
        }

        die "dry run, rolling back.\n";
    });
}

sub load_group {
    my ( $self, $group ) = @_;

    # find or create features for all the group members
    my $features = $self->find_or_create_features( $group->{members} );

    # make a new feature for the group
    my $group_feature =
        $self->schema->resultset('Sequence::Feature')
                     ->create({
                         name       => $group->{group}{id}||die,
                         uniquename => $group->{group}{id}||die,
                         type       => $self->_group_term,
                         organism   => $self->_find_organism( 'any' ),
                       });

    my $member_of = $self->_member_of_term;
    # insert feature relationships for each of the group members
    for my $member ( values %$features ) {
        $self->schema->populate( 'Sequence::FeatureRelationship', [
            [  qw[ subject  type         object            ]],
            map {[ $_,      $member_of,  $group_feature    ]}
            values %$features
        ]);
    }
}

has '_org_cache' => ( is => 'ro', default => sub { {} } );
sub _find_organism {
    my ( $self, $name ) = @_;
    return $self->_org_cache->{$name} ||= $self->schema->resultset('Organism::Organism')->find({ species => $name });
}

sub find_or_create_features {
    my ( $self, $members ) = @_;

    # find existing features in a single query
    my %features =
        map { $_->name => $_ }
        $self->schema->resultset('Sequence::Feature')
                     ->search({
                         name => [ map $_->{id}, @$members ],
                         type_id => $self->_member_term->cvterm_id,
                     });

    # check that we have all of them, creating or dieing for any missing ones
    for my $member (@$members) {
        my $name = $member->{id};
        unless( $features{$name} ) {
            if( $self->create_members ) {
                $features{$name} =
                    $self->schema->resultset('Sequence::Feature')
                         ->create({
                             name       => $name,
                             uniquename => $name,
                             type       => $self->_member_term,
                             organism   => $self->_find_organism( $member->{species} ),
                         });
            } else {
                die "$.: no ".$self->member_type." feature found with name $name, and --create_members option not passed\n";
            }
        }
    }

    return \%features;
}

sub parse_orthomcl_group_line {
    my ( $self, $line ) = @_;

    my ( $cluster, @members ) =
        map {
            my ( $ident, $other ) = / ([^(]+) \( ([^\)]+) /x or die 'parse error';
            [ $ident, split / \s* , \s* /x, $other ]
        }
        split / (?<=\)) [:\s] \s* /x,
        $line;

    my %cluster = ( id => shift(@$cluster), map { reverse split /\s+/, $_ } @$cluster );
    @members = map +{ id => $_->[0], species => $_->[1] }, @members;

    return { group => \%cluster, members => \@members };
}

1;
