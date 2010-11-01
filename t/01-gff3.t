#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;

sub TEST_IS_PRAGMA : Test(4) {
    my $loader = Bio::Chado::Loader::GFF3->new;

    ok(!$loader->is_pragma_line('Not a pragma!'),'is_pragma_line');
    ok($loader->is_pragma_line('## a pragma!'),'is_pragma_line');
    ok(!$loader->is_pragma_line('# not a pragma!'),'is_pragma_line');
    ok(!$loader->is_pragma_line(''),'is_pragma_line');
}

sub TEST_ITAG1_GENOMIC_REF : Test(3) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        filename => "t/data/ITAG1_genomic_ref_sample.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    $loader->parse();
    is($loader->count_cvterms, 2, 'found 2 unique cvterms');
    is($loader->count_features, 4, 'found 4 unique features');

}

__PACKAGE__->runtests;
