package Bio::Chado::Loader::GFF3;

use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';

has schema => (
    is  => 'rw',
    isa => 'Bio::Chado::Schema',
);

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

=head1 NAME

Bio::Chado::Loader::GFF3 - Load GFF3 files into the Chado database schema

=head1 SYNOPSIS

  use Bio::Chado::Loader;


=head1 DESCRIPTION

The GFF3 spec is available online at L<http://www.sequenceontology.org/gff3.shtml>

=cut

1;
