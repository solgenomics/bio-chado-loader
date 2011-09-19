package Test::Bio::Chado::Loader::FASTA;

use Test::Class::Most;
use Bio::Chado::Loader::FASTA;

my @test_fasta = qw( t/data/test1.fasta t/data/test2.fasta );

sub TEST_FASTA_SPOOLING : Test(1) {
    my $self = shift;
    my $loader = $self->_test_loader_1;
    my $s = $loader->fasta_stream( \@test_fasta );
    is( $self->_count_stream( $s ), 4, 'got 4 seqs from the seq stream' );
}

sub TEST_FASTA_PARSE : Test(6) {
    my $self = shift;
    my $loader = $self->_test_loader_1;
    my @test_seqs = (
        [ "foo234.12_X|27  bar bal bas[1]\r\nAATCGATCGATCG\nATCGTACGATCAG AGTAGCT\nAATCGATCGATCG\n",
          "foo234.12_X|27",
          "bar bal bas[1]",
          "AATCGATCGATCGATCGTACGATCAGAGTAGCTAATCGATCGATCG",
        ],
        [ "foo234.12_X|27\nAATCGATCGATCG\nATCGTACGATCAG AGTAGCT\n",
          "foo234.12_X|27",
          undef,
          "AATCGATCGATCGATCGTACGATCAGAGTAGCT",
         ],
       );

    for (@test_seqs) {
        my ( $seq, $tid, $tdefline, $tseq_after ) = @$_;
        my ( $id, $dl ) = $loader->parse_fasta( \$seq );
        is( $id,  $tid, 'right ID' );
        is( $dl,  $tdefline, 'right defline' );
        is( $seq, $tseq_after, 'right sequence' );
    }
}

sub _test_loader_1 {
    my $self = shift;
    Bio::Chado::Loader::FASTA->new(
        db_dsn   => 'dbi:SQLite:dbname=:memory:',
        type_name => 'sequence_assembly',
        create_features => 1,
        organism_name => 'Tyrannosaurus _ Tyrannosaurus rex',
        analysis_name => 'Test analysis',
        source => 'test_source',
        large_residues => 2000,
        @_
       );
}

sub TEST_LOAD : Test(6) {
    my $self = shift;
    my $loader = $self->_test_loader_1;
    $loader->schema->deploy;
    $loader->schema->txn_do( sub { $self->_populate( $loader->schema ) } );

    $loader->run( @test_fasta );

    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' )
        or diag explain( [ map { [ $_->feature_id, $_->type->name => $_->value ] } $loader->schema->resultset('Sequence::Featureprop')->all ] );

    # run the load again and test that we did not get duplicate features
    $loader->run( @test_fasta );
    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' )
        or diag explain( [ map { [ $_->feature_id, $_->type->name => $_->value ] } $loader->schema->resultset('Sequence::Featureprop')->all ] );

    # run it yet again with create_features off
    $loader->create_features( 0 );
    $loader->run( @test_fasta );
    my $schema = $loader->schema;
    is( $schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' )
        or diag explain( [ map { [ $_->feature_id, $_->type->name => $_->value ] } $loader->schema->resultset('Sequence::Featureprop')->all ] );

}

sub _test_loader_2 {
    my $self = shift;
    Bio::Chado::Loader::FASTA->new(
        db_dsn   => 'dbi:SQLite:dbname=:memory:',
        type_regex => 'type:(\S+)',
        create_features => 1,
        organism_name => 'Tyrannosaurus _ Tyrannosaurus rex',
        analysis_name => 'Test analysis',
        source => 'test_source',
        large_residues => 1000,
        @_
       );
}

sub TEST_LOAD_2 : Test(6) {
    my $self = shift;
    my $loader = $self->_test_loader_2;
    $loader->schema->deploy;
    $loader->schema->txn_do( sub { $self->_populate( $loader->schema ) } );

    $loader->run( @test_fasta );

    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 2, 'correct featureprop count' );

    for (
          [ 'test1.1' => 'beef_stew' ],
          [ 'test1.2' => 'sequence_assembly' ],
          [ 'test2.1' => 'chicken_soup' ],
          [ 'test2.2' => 'beef_stew' ],
        ) {

        is( $loader->schema->resultset('Sequence::Feature')->find({ name => $_->[0] })->type->name,
            $_->[1],
            "$_->[0] has type $_->[1]"
          );
    }
}

sub _populate {
    my ( $self, $schema) = @_;

    $schema->resultset('Cv::Cvterm')
           ->create_with({
               cv   => 'sequence',
               db   => 'null',
               name => $_,
           }) for qw/ sequence_assembly chicken_soup beef_stew /;

    $schema->resultset('Organism::Organism')
           ->create({
               genus => 'Tyrannosaurus',
               species => 'Tyrannosaurus rex',
           });
}

sub _count_stream {
    my ( $self, $s ) = @_;
    local $/ = "\n>";
    my $count = 0;
    $count++ while <$s>;
    return $count;
}

__PACKAGE__->runtests;

