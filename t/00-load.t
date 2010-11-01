#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Bio::Chado::Loader' );
}

diag( "Testing Bio::Chado::Loader $Bio::Chado::Loader::VERSION, Perl $], $^X" );
