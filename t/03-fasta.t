package Test::Bio::Chado::Loader::FASTA;

use Test::Class::Most;
use Bio::Chado::Loader::FASTA;

my @test_fasta = qw( t/data/test1.fasta t/data/test2.fasta );

sub TEST_FASTA_SPOOLING : Test(1) {
    my $self = shift;
    my $loader = Bio::Chado::Loader::FASTA->new;
    my $s = $loader->fasta_stream( \@test_fasta );
    is( $self->_count_stream( $s ), 4, 'got 4 seqs from the seq stream' );
}

sub TEST_FASTA_PARSE : Test(3) {
    my $self = shift;
    my $loader = Bio::Chado::Loader::FASTA->new;
    my $seq = "foo234.12_X|27  bar bal bas[1]\r\nAATCGATCGATCG\nATCGTACGATCAG AGTAGCT\n";
    my ( $id, $dl ) = $loader->parse_fasta( \$seq );
    is( $id,  'foo234.12_X|27', 'right ID' );
    is( $dl,  'bar bal bas[1]', 'right defline' );
    is( $seq, 'AATCGATCGATCGATCGTACGATCAGAGTAGCT', 'right sequence' );
}

sub TEST_LOAD : Test(6) {
    my $self = shift;
    my $loader = Bio::Chado::Loader::FASTA->new(
        db_dsn   => 'dbi:SQLite:dbname=:memory:',
        type_name => 'sequence_assembly',
        create_features => 1,
        organism_name => 'Tyrannosaurus _ Tyrannosaurus rex',
        analysis_name => 'Test analysis',
        source => 'test_source',
        large_residues => 800,
       );
    $loader->schema->deploy;
    $loader->schema->txn_do( sub { $self->_populate( $loader->schema ) } );

    $loader->run( @test_fasta );

    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' );

    # run the load again and test that we did not get duplicate features
    $loader->run( @test_fasta );
    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' );

    # run it yet again with create_features off
    $loader->create_features( 0 );
    $loader->run( @test_fasta );
    is( $loader->schema->resultset('Sequence::Feature')->count, 4, 'correct feature count' );
    is( $loader->schema->resultset('Sequence::Featureprop')->count, 1, 'correct featureprop count' );


}

sub _populate {
    my ( $self, $schema) = @_;

    $schema->resultset('Cv::Cvterm')
           ->create_with({
               cv   => 'sequence',
               db   => 'null',
               name => 'sequence_assembly',
           });

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

