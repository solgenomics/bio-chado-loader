use strict;
use warnings;

use Test::More skip_all => 'Disabling since RB removed FixGMODBulkGFF3Polypeptides.pm' ;
#https://github.com/solgenomics/bio-chado-loader/commit/4580abb815a50b8cba6f1f7ab83266ac02abf232#diff-0

use_ok( 'Bio::Chado::Loader::FixGMODBulkGFF3Polypeptides' );

done_testing;
