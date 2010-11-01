#!perl

use Test::More tests => 2;

use lib './lib';

BEGIN {
    use_ok( 'Bio::Chado::Loader' );
    use_ok( 'Bio::Chado::Loader::GFF3' );
}

diag( "Testing Bio::Chado::Loader $Bio::Chado::Loader::VERSION, Perl $], $^X" );
