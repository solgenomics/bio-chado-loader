#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

sub TEST_ITAG1_GENOMIC_REF : Test(3) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/ITAG1_genomic_ref_sample.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    $loader->parse();
    is($loader->count_cvterms_gff, 2, 'found 2 unique cvterms');
    is($loader->count_features_gff, 4, 'found 4 unique features');
}

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

