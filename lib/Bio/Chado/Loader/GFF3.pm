package Bio::Chado::Loader::GFF3;

=head1 NAME

Bio::Chado::Loader::GFF3 - Load GFF3 files into the Chado database schema

=head1 SYNOPSIS

  use Bio::Chado::Loader;


=head1 DESCRIPTION

The GFF3 spec is available online at L<http://www.sequenceontology.org/gff3.shtml>

=over

=cut


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
        feature_exists => 'exists',
    },
    );

has filename => (#Added SS, LM
    is      => 'rw',
    isa     => 'Str',
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

has is_analysis => (
    documentation => <<'',
set true if this feature should be recorded as from an analysis

    is      => 'ro',
    isa     => 'Bool',
    default => 0,
);

=item C<run ()>



=cut

sub run {
    my ($self, %args) = @_;

    exit 0;
}

=item C<parse ()>

Parse a GFF3 file

=cut

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

=item C<parse_line ()>

Parse a GFF3 line

=cut

sub parse_line {
    my ($self, $line) = @_;
    my @fields = split /\t/, $line;
    my $attribs = { map { split (/=/, $_)  }  split /;/, $fields[ATTRIBUTES] };

    $self->add_cvterm( $fields[TYPE]   => 1 );
    $self->add_feature( $attribs->{ID} || $fields[SEQID] => 1 );

    $self->_validate_parents($attribs);
}

=item C<_validate_parents ()>

Check if parent feature exists

=cut


sub _validate_parents {
    my ($self, $attribs) = @_;
    return unless $attribs->{Parent};
    my (@parents) = split (/,/, $attribs->{Parent});
    for my $p (@parents) {
        unless ( $self->feature_exists( $p ) ) {
            die "$p is an unknown Parent!";
        }
    }
}

=item C<is_pragma_line ()>

Check if line is a pragma

=cut

sub is_pragma_line {
    my ($self, $line) = @_;
    return $line =~ m/^##/ ? 1 : 0;
}

###
1;   #do not remove
###

=pod

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Jonathan "Duke" Leto	<duke@leto.net>
    Surya Saha			<suryasaha@cornell.edu , @SahaSurya>   

=cut