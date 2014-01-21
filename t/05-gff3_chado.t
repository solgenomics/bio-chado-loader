#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

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

sub TEST_DB_CACHE : Test(7){
	my $loader = Bio::Chado::Loader::GFF3->new();
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	is($loader->organism_exists() , 1, 'got correct org id from org table');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	ok($loader->populate_cache(), 'populated cache');
}
__PACKAGE__->runtests;

