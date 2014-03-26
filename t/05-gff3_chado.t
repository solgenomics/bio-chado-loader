#!perl
package Test::Bio::Chado::Loader::GFF3;

use Test::Class::Most;
use Bio::Chado::Loader::GFF3;
#use Carp::Always;

#util functions
sub mem_used {
	my ( $i, $t );
	$t = new Proc::ProcessTable;
	foreach my $got ( @{ $t->table } ) {
		next if not $got->pid eq $$;
		$i = $got->size;
	}
	print STDERR "Process id=", $$, "\n";
	print STDERR "Memory used(MB)=", $i / 1024 / 1024, "\n";
}

sub run_time {
	my ( $user_t, $system_t, $cuser_t, $csystem_t );
	( $user_t, $system_t, $cuser_t, $csystem_t ) = times;
	print STDERR "System time for process: $system_t\n";
	print STDERR "User time for process: $user_t\n\n";
}

#test functions
sub TEST_DB_CONNECT : Test(4){
	my $loader = Bio::Chado::Loader::GFF3->new();
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
    #is($loader->schema->resultset('Sequence::Feature')->count,2589670,'got correct nof rows in feature table');
}

sub TEST_DB_ORG : Test(7){
	my $loader = Bio::Chado::Loader::GFF3->new();
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	is($loader->organism_exists() , 1, 'got correct org id from org table');
	ok($loader->organism_name('Platostoma coeruleum'), 'assigned org name');
	is($loader->organism_exists() , 2363,'got correct org id from org table');
}

sub TEST_DB_CACHE_UNKNOWN_PARENT : Test(9){
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/db_unknown_parent.gff3",
    );
    ok($loader->debug(1), 'set debug verbosity flag');
    
    #Create cache first and then parse GFF 
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	is($loader->organism_exists() , 1, 'got correct org id from org table');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	run_time(); mem_used();
	ok($loader->populate_cache(), 'populated cache');
	run_time(); mem_used();
    dies_ok sub { $loader->parse() }, qr/dummy is an unknown Parent/;
}

#Change cache SQL search string in populate_cache() to run routine
#sub TEST_DB_CACHE_INSERT : Test(13){
#	my $loader = Bio::Chado::Loader::GFF3->new(
#        file_name => "t/data/insert_test.gff3",
#    );
#    
#    #Create cache first and then parse GFF 
#	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
#	ok($loader->db_user('test_usr'),'assigned user');
#	ok($loader->db_pass('test_usr'),'assigned pw');
#	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
#	is($loader->organism_exists() , 1, 'got correct org id from org table');
#	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
#	run_time(); mem_used();
#	ok(my $cnt = $loader->populate_cache(), 'populated cache');
#	print STDERR 'Cache has '.$cnt."\n";
#	run_time(); mem_used();
#
#    $loader->parse();
#    is($loader->count_cvterms_gff, 1, 'found 1 unique cvterms');
#    is($loader->count_features_gff, 3, 'found 3 unique features');
#    cmp_set( [ keys %{$loader->features_gff} ], [ qw/
#        gene:dummygene1 gene:dummygene2 gene:dummygene3/ ],
#        'Found expected features',
#    );
#    cmp_set( [ keys %{$loader->cvterms_gff} ], [ qw/ gene / ],
#        'Found expected cvterms',
#    );
#    run_time(); mem_used();
#    
#	ok($cnt=$loader->prepare_bulk_operation(), 'loaded data structures and wrote exception file');
#	print STDERR 'Prepped '.$cnt." recs for insertion\n";
#	ok($loader->bulk_upload(),'updated locgroups and inserted new rows into featureloc')
#}

sub TEST_DB_INSERT_Solyc01g112300 : Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    ok($loader->debug(1), 'set debug verbosity flag');
    
    #Create cache first and then parse GFF 
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	run_time(); mem_used();
	#Call populate_cache() before parse() in case parents are not in GFF but in DB already
	ok($loader->populate_cache(), 'populated cache');
	#Removed count test for public release since DB specific
	#is($loader->count_feature_uniquename_feature_id_cache(),1275149,'correct number of features read from DB into cache');
	run_time(); mem_used();

    $loader->parse();
	#only genes, mRNAs, exons, introns and  polypeptides 
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
    run_time(); mem_used();
    my $cnt;
	ok($cnt=$loader->prepare_bulk_operation(), 'loaded data structures and wrote exception file');
	is($loader->count_feature_ids_uniquenames_gff(),8,'recorded 8 pairs in feature_ids_uniquenames_gff hash');
	print STDERR 'Prepped '.$cnt." recs for insertion\n";
	ok($loader->bulk_upload(),'updated locgroups and inserted new rows into featureloc');
	run_time(); mem_used();
}

sub TEST_DB_DELETE_Solyc01g112300 : Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    ok($loader->debug(1), 'set debug verbosity flag');
    ok($loader->delete(1), 'set DELETION flag');
    
    #Connect 
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	#Call populate_cache() before parse() in case parents are not in GFF but in DB already
	ok($loader->populate_cache(), 'populated cache');

    $loader->parse(); #only genes, mRNAs, exons, introns and  polypeptides
	#$loader->dump_features_gff();
    my $cnt;
	ok($cnt=$loader->prepare_bulk_operation(), 'loaded data structures and wrote exception file');
	
	#NOTE: Do not create exons or polypeptide in mRNA_handler if features deleted manually already during testing
	print STDERR 'Prepped '.$cnt." recs for DELETION\n";
    ok($cnt=$loader->bulk_delete(),'DELETED rows from featureloc and decremented locgroups of remaining featureloc\'s');
    print STDERR 'bulk_delete deleted '.$cnt." records\n";
    run_time(); mem_used();
}


sub TEST_mRNA_HANDLER: Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    ok($loader->debug(1), 'set debug verbosity flag');
    $loader->parse();
    #only genes, mRNAs, exons, introns and  polypeptides 
    is($loader->count_cvterms_gff, 5, 'found 5 unique cvterms');
    is($loader->count_features_gff, 8, 'found 8 unique features');
    $loader->dump_features_gff();
}

sub TEST_mRNA_HANDLER_NOEXON: Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/Solyc01g112300.2.gff3_noexon",
    );
    isa_ok($loader, 'Bio::Chado::Loader::GFF3');
    $loader->parse();
    $loader->dump_features_gff();
}

sub TEST_DB_INSERT_CH01_100 : Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/ch01_100genes.gff3",
    );
    ok($loader->debug(1), 'set debug verbosity flag');
    
    #Create cache first and then parse GFF 
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	#Call populate_cache() before parse() in case parents are not in GFF but in DB already
	ok($loader->populate_cache(), 'populated cache');

    $loader->parse();
    my $cnt;
	ok($cnt=$loader->prepare_bulk_operation(), 'loaded data structures and wrote exception file');
	print STDERR 'Prepped '.$cnt." recs for insertion\n";
	ok($loader->bulk_upload(),'updated locgroups and inserted new rows into featureloc');
	run_time(); mem_used();
}

sub TEST_DB_DELETE_CH01_100 : Tests{
	my $loader = Bio::Chado::Loader::GFF3->new(
        file_name => "t/data/ch01_100genes.gff3",
    );
    ok($loader->debug(1), 'set debug verbosity flag');
    ok($loader->delete(1), 'set DELETION flag');
    
    #Connect 
	ok($loader->db_dsn("dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432"),'assigned dsn');
	ok($loader->db_user('test_usr'),'assigned user');
	ok($loader->db_pass('test_usr'),'assigned pw');
	ok($loader->organism_name('Solanum lycopersicum'), 'assigned org name');
	ok($loader->organism_id($loader->organism_exists()), 'assigned org id');
	#Call populate_cache() before parse() in case parents are not in GFF but in DB already
	ok($loader->populate_cache(), 'populated cache');

    $loader->parse(); #only genes, mRNAs, exons, introns and  polypeptides
	#$loader->dump_features_gff();
    my $cnt;
	ok($cnt=$loader->prepare_bulk_operation(), 'loaded data structures and wrote exception file');
	
	#NOTE: Do not create exons or polypeptide in mRNA_handler if features deleted manually already during testing
	print STDERR 'Prepped '.$cnt." recs for DELETION\n";
    ok($cnt=$loader->bulk_delete(),'DELETED rows from featureloc and decremented locgroups of remaining featureloc\'s');
    print STDERR 'bulk_delete deleted '.$cnt." records\n";
    run_time(); mem_used();
}

#$ENV{TEST_METHOD} = 'TEST_DB_CONNECT'; #set regex for routine(s) to run
__PACKAGE__->runtests;