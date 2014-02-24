#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

sub TEST_ITAG1_GENOMIC_REF : Test(6) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/ITAG1_genomic_ref_sample.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    is($loader->has_file_name(),1,'file name present');
    is($loader->file_name(),"t/data/ITAG1_genomic_ref_sample.gff3",'checked file name');
    is($loader->is_analysis(),0,'correct analysis flag');
    $loader->parse();
    is($loader->count_cvterms_gff, 2, 'found 2 unique cvterms');
    is($loader->count_features_gff, 4, 'found 4 unique features');
}

sub TEST_Solyc01g112300_2 : Test(5) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    $loader->parse();
    is($loader->count_cvterms_gff, 5, 'found 5 unique cvterms');
    is($loader->count_features_gff, 8, 'found 8 unique features');
    
    cmp_set( [ keys %{$loader->features_gff} ], [ qw/
        gene:Solyc01g112300.2 mRNA:Solyc01g112300.2.1
        exon:Solyc01g112300.2.1.1 intron:Solyc01g112300.2.1.1
        exon:Solyc01g112300.2.1.2 intron:Solyc01g112300.2.1.2 
        exon:Solyc01g112300.2.1.3 polypeptide:Solyc01g112300.2.1/ ],
        'Found expected features',
    );
    cmp_set( [ keys %{$loader->cvterms_gff} ], [ qw/
        polypeptide exon intron gene mRNA/ ],
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

