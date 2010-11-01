#!perl

use Test::More tests => 3;

use lib './lib';

BEGIN {
    use_ok( 'Bio::Chado::Loader' );
    use_ok( 'Bio::Chado::Loader::GFF3' );
    use_ok( 'Bio::Chado::Loader::FASTA' );
}

diag( "Testing Bio::Chado::Loader $Bio::Chado::Loader::VERSION, Perl $], $^X" );
