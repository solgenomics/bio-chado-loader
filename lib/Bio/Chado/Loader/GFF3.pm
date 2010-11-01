package Bio::Chado::Loader::GFF3;

use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
with 'Bio::Chado::Loader';

use autodie qw(:all);
use Data::Dumper;
use 5.010;

sub run {
    my ($self, %args) = @_;

    exit 0;
}

sub parse {
    my ($self, %args) = @_;
    open my $fh, "<", $self->filename;
    while( my $line = <$fh> ) {
        next if $self->is_pragma_line($line);
        $self->parse_line($line);
    }
    #warn Dumper [ $self->cvterms ];
}

sub parse_line {
    my ($self, $line) = @_;
    my @fields = split /\t/, $line;

    $self->add_cvterm( $fields[2] => 1 );
}

sub is_pragma_line {
    my ($self, $line) = @_;
    return $line =~ m/^##/ ? 1 : 0;
}

1;
