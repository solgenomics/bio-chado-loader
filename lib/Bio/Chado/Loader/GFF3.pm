package Bio::Chado::Loader::GFF3;

use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
with 'Bio::Chado::Loader';

use autodie qw(:all);
use 5.010;

sub run {
    my ($self, %args) = @_;

    exit 0;
}

1;
