package Bio::Chado::Loader::GFF3;

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

has cvterms => (
    is      => 'rw',
    isa     => 'HashRef[Str]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_cvterm    => 'set',
        count_cvterms => 'count',
    },
);

use constant SEQID         => 0;
use constant SOURCE        => 1;
use constant TYPE          => 2;
use constant FEATURE_START => 3;
use constant FEATURE_END   => 4;
use constant SCORE         => 5;
use constant STRAND        => 6;
use constant PHASE         => 7;
use constant ATTRIBUTES    => 8;

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
        next if $line =~ m/^\s*$/;
        chomp $line;
        next if $self->is_pragma_line($line);
        $self->parse_line($line);
    }
    #warn Dumper [ $self->features ];
}

sub parse_line {
    my ($self, $line) = @_;
    my @fields = split /\t/, $line;
    my $attribs = { map { split /=/, $_  }  split /;/, $fields[ATTRIBUTES] };

    $self->add_cvterm( $fields[TYPE]   => 1 );
    $self->add_feature( $attribs->{ID} => 1 );
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
