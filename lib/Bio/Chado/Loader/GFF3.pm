package Bio::Chado::Loader::GFF3;

use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
use autodie qw(:all);
use 5.010;

has dbhost => (
    is      => 'ro',
    default => 'localhost',
    isa     => 'Str',
);
has dbuser => (
    is      => 'ro',
    default => 'postgres',
    isa     => 'Str',
);
has dbpass => (
    is      => 'ro',
    isa     => 'Str',
    default => '',
);
has dbname => (
    is      => 'ro',
    default => 'cxgn',
    isa     => 'Str',
);
has organism => (
    is      => 'ro',
    default => 'tomato',
    isa     => 'Str',
);

sub run {
    my ($self, %args) = @_;

    exit 0;
}

1;
