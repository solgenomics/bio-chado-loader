use strict;
use warnings;

use Test::More;

use_ok( 'Bio::Chado::Loader::OrthoMCL' );

my $data = Bio::Chado::Loader::OrthoMCL->parse_orthomcl_group_line("ORTHOMCL24856(2 genes,1 taxa):   Solyc11g013710.1.1(tomato) Solyc11g013720.1.1(tomato)\n");
is_deeply(
    $data,
    {
        'group' => {
            'genes' => '2',
            'id' => 'ORTHOMCL24856',
            'taxa' => '1'
            },
        'members' => [
            {
                'id' => 'Solyc11g013710.1.1',
                'species' => 'tomato'
            },
            {
                'id' => 'Solyc11g013720.1.1',
                'species' => 'tomato'
            }
        ],
   },
   'parser works' );

done_testing;

