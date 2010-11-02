#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

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

sub TEST_CANONICAL_GENE : Tests {
    my $loader = Bio::Chado::Loader::GFF3->new(
        filename => "t/data/canonical_gene.gff3",
    );
    $loader->parse();
    is($loader->count_cvterms, 5, 'found 5 unique cvterms');
    is($loader->count_features, 14, 'found 14 unique features');

    cmp_set( [ keys %{$loader->features} ], [ qw/
        cds00001 cds00002 cds00003 cds00004
        exon00001 exon00002 exon00003 exon00004 exon00005
        gene00001 mRNA00001 mRNA00002 mRNA00003 tfbs00001/ ],
        'Found expected features',
    );
    cmp_set( [ keys %{$loader->cvterms} ], [ qw/
        CDS TF_binding_site exon gene mRNA/ ],
        'Found expected cvterms',
    );
}

sub TEST_UNKNOWN_PARENT : Tests {
    my $loader = Bio::Chado::Loader::GFF3->new(
        filename => "t/data/unknown_parent.gff3",
    );
    dies_ok sub { $loader->parse() }, qr/Bobby_Tables is an unknown Parent/;
}

__PACKAGE__->runtests;

