#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

#INVALID DUE TO USE OF  Bio::GFF3::LowLevel::Parser
#sub TEST_IS_PRAGMA : Test(4) {
#    my $loader = Bio::Chado::Loader::GFF3->new;
#
#    ok(!$loader->is_pragma_line('Not a pragma!'),'is_pragma_line');
#    ok($loader->is_pragma_line('## a pragma!'),'is_pragma_line');
#    ok(!$loader->is_pragma_line('# not a pragma!'),'is_pragma_line');
#    ok(!$loader->is_pragma_line(''),'is_pragma_line');
#}

sub TEST_ITAG1_GENOMIC_REF : Test(3) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/ITAG1_genomic_ref_sample.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    $loader->parse();
    is($loader->count_cvterms_gff, 2, 'found 2 unique cvterms');
    is($loader->count_features_gff, 4, 'found 4 unique features');
}

#FAILS SINCE >1 CDS WITH SAME NAME. ARE CDS STORED IN CHADO AT ALL??
#sub TEST_CANONICAL_GENE : Tests {
#    my $loader = Bio::Chado::Loader::GFF3->new(
#        file_name => "t/data/canonical_gene.gff3",
#    );
#    $loader->parse();
#    is($loader->count_cvterms_gff, 5, 'found 5 unique cvterms');
#    is($loader->count_features_gff, 14, 'found 14 unique features');
#
#    cmp_set( [ keys %{$loader->features_gff} ], [ qw/
#        cds00001 cds00002 cds00003 cds00004
#        exon00001 exon00002 exon00003 exon00004 exon00005
#        gene00001 mRNA00001 mRNA00002 mRNA00003 tfbs00001/ ],
#        'Found expected features',
#    );
#    cmp_set( [ keys %{$loader->cvterms_gff} ], [ qw/
#        CDS TF_binding_site exon gene mRNA/ ],
#        'Found expected cvterms',
#    );
#}

sub TEST_Solyc01g112300_2 : Test(4) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    $loader->parse();
    is($loader->count_cvterms_gff, 7, 'found 7 unique cvterms');
    is($loader->count_features_gff, 12, 'found 12 unique features');

    cmp_set( [ keys %{$loader->features_gff} ], [ qw/
        gene:Solyc01g112300.2 mRNA:Solyc01g112300.2.1
        exon:Solyc01g112300.2.1.1 five_prime_UTR:Solyc01g112300.2.1.0
        CDS:Solyc01g112300.2.1.1 intron:Solyc01g112300.2.1.1
        exon:Solyc01g112300.2.1.2 CDS:Solyc01g112300.2.1.2
        intron:Solyc01g112300.2.1.2 exon:Solyc01g112300.2.1.3
        CDS:Solyc01g112300.2.1.3 three_prime_UTR:Solyc01g112300.2.1.0/ ],
        'Found expected features',
    );
    cmp_set( [ keys %{$loader->cvterms_gff} ], [ qw/
        CDS five_prime_UTR three_prime_UTR exon intron gene mRNA/ ],
        'Found expected cvterms',
    );
}

#INVALID DUE TO USE OF  Bio::GFF3::LowLevel::Parser
#sub TEST_UNKNOWN_PARENT : Tests {
#    my $loader = Bio::Chado::Loader::GFF3->new(
#        file_name => "t/data/unknown_parent.gff3",
#    );
#    dies_ok sub { $loader->parse() }, qr/Bobby_Tables is an unknown Parent/;
#}

sub TEST_DUPLICATE_ID : Tests {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/duplicate_id.gff3",
    );
    dies_ok sub { $loader->parse() }, qr/Multiple features with same ID found. Exiting../;
}

sub TEST_MISSING_ID : Tests {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/missing_id.gff3",
    );
    dies_ok sub { $loader->parse() }, qr/No ID defined for feature in GFF. Exiting../;
}

__PACKAGE__->runtests;

