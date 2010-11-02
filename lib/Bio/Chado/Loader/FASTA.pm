package Bio::Chado::Loader::FASTA;
use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';

has schema => (
    is  => 'rw',
    isa => 'Bio::Chado::Schema',
);

has features => (
    is      => 'rw',
    isa     => 'HashRef[Str]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_feature    => 'set',
        count_features => 'count',
    },
);

with 'Bio::Chado::Loader';

sub run {
    exit 0;
}

use autodie qw(:all);
use Data::Dumper;
use 5.010;
1;
