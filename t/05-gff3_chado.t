#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;


sub TEST_Solyc01g112300_2 : Test(4) {
    my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    $loader->parse();
    is($loader->count_cvterms, 7, 'found 7 unique cvterms');
    is($loader->count_features, 12, 'found 12 unique features');

    cmp_set( [ keys %{$loader->features} ], [ qw/
        gene:Solyc01g112300.2 mRNA:Solyc01g112300.2.1
        exon:Solyc01g112300.2.1.1 five_prime_UTR:Solyc01g112300.2.1.0
        CDS:Solyc01g112300.2.1.1 intron:Solyc01g112300.2.1.1
        exon:Solyc01g112300.2.1.2 CDS:Solyc01g112300.2.1.2
        intron:Solyc01g112300.2.1.2 exon:Solyc01g112300.2.1.3
        CDS:Solyc01g112300.2.1.3 three_prime_UTR:Solyc01g112300.2.1.0/ ],
        'Found expected features',
    );
    cmp_set( [ keys %{$loader->cvterms} ], [ qw/
        CDS five_prime_UTR three_prime_UTR exon intron gene mRNA/ ],
        'Found expected cvterms',
    );
}

sub TEST_DB_CONNECT : Test(4){
	my $loader = Bio::Chado::Loader::GFF3->new();
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
    is($loader->schema->resultset('Sequence::Feature')->count,2589666,'got correct nof rows in feature table');
}

sub TEST_DB_QUERY : Test(7){
	my $loader = Bio::Chado::Loader::GFF3->new();
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	is($loader->organism_exists() , 1, 'got correct org id from org table');
	ok($loader->organism_name('Platostoma coeruleum'), 'assigned org name');
	is($loader->organism_exists() , 2363,'got correct org id from org table');
}

__PACKAGE__->runtests;

