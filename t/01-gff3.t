#!perl

use Test::More tests => 4;
use Bio::Chado::Loader::GFF3;

diag( "Testing Bio::Chado::Loader $Bio::Chado::Loader::VERSION, Perl $], $^X" );

my $loader = Bio::Chado::Loader::GFF3->new;

ok(!$loader->is_pragma_line('Not a pragma!'),'is_pragma_line');
ok($loader->is_pragma_line('## a pragma!'),'is_pragma_line');
ok(!$loader->is_pragma_line('# not a pragma!'),'is_pragma_line');
ok(!$loader->is_pragma_line(''),'is_pragma_line');
