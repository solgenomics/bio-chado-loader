package Bio::Chado::Loader::GFF3;
use strict;
use warnings;

=head1 NAME

Bio::Chado::Loader::GFF3 - Load GFF3 files into the Chado database schema

=head1 SYNOPSIS

  use Bio::Chado::Loader::GFF3;
  mx-run Bio::Chado::Loader::GFF3 [ options ] gff3_file gff3_file ...


=head1 DESCRIPTION

This is a MooseX::Runnable class that loads feature data given in GFF3 format into a Chado schema.
The GFF3 spec is available online at L<http://www.sequenceontology.org/gff3.shtml>

=over

=cut


use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
#use Bio::GFF3::LowLevel  qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);
#use File::Basename;

has schema => (
    is  => 'rw',
    isa => 'Bio::Chado::Schema',
);

has 'file_name' => (
	isa       => 'Str',
	is        => 'rw',
	predicate => 'has_file_name',
	clearer   => 'clear_file_name'
);

has 'pragma_lines' => (
	isa       => 'ArrayRef[Str]',
	is        => 'rw',
	predicate => 'has_pragma_lines',
	clearer   => 'clear_pragma_lines'
);

has 'features' => (
    is      => 'rw',
    #isa     => 'HashRef[Str]',
    isa     => 'HashRef[ArrayRef]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_feature    => 'set',
        count_features => 'count',
        feature_exists => 'exists',
    },
    );
    
has 'cvterms' => (
    is      => 'rw',
    isa     => 'HashRef[Str]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_cvterm    => 'set',
        count_cvterms => 'count',
    },
);

has 'is_analysis' => (
    documentation => <<'',
set true if this feature should be recorded as from an analysis

    is      => 'ro',
    isa     => 'Bool',
    default => 0,
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

=item C<run ()>

TODO // See Fasta.pm

=cut

sub run {
    my ($self, %args) = @_;

    exit 0;
}

=item C<parse ()>

Parse a GFF3 file. Calls parse_line() for each line.

=cut

sub parse {
    my ($self, %args) = @_;
    open my $fh, "<", $self->file_name;
    while( my $line = <$fh> ) {
        next if $line =~ m/^\s*$/;
        chomp $line;
        next if $self->is_pragma_line($line);
        $self->parse_line($line);
    }
    warn Dumper [ $self->features ];
}

=item C<parse_line ()>

Parse a GFF3 line to get feature details.

=cut

sub parse_line {
    my ($self, $line) = @_;
    my @fields = split (/\t/, $line);
    my $attribs = { map { split (/=/, $_)  }  split (/;/, $fields[ATTRIBUTES]) };

    $self->add_cvterm( $fields[TYPE]   => 1 );
    #$self->add_feature( $attribs->{ID} || $fields[SEQID] => 1 );
    $self->add_feature( $attribs->{ID} || $fields[SEQID] => \@fields );

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

    Jonathan "Duke" Leto	<jonathan at leto.net>
    Surya Saha			<suryasaha at cornell.edu , @SahaSurya>   

=cut