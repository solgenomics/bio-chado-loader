package Bio::Chado::Loader;
use Moose::Role;

use Bio::Chado::Schema;

has db_dsn => (
    documentation => 'DBI dsn for our database',

    is      => 'rw',
    default => 'dbi:Pg:dbname=cxgn',
    isa     => 'Str',
);

has db_user => (
    documentation => 'database username',

    is      => 'rw',
    default => 'postgres',
    isa     => 'Str',
);

has db_pass => (
    documentation => 'database password',

    is      => 'rw',
    isa     => 'Str',
);

has db_attrs => (
    documentation => 'hashref of DBI database attributes',

    is      => 'rw',
    isa     => 'HashRef',
    default => sub { +{} },
);

requires 'schema';

# default _build_schema just connects with the db args.  consuming
# classes can use this if they just set lazy_build => 1 on their
# schema attrs
sub _build_schema {
    my ( $self ) = @_;

    return Bio::Chado::Schema->connect(
        $self->db_dsn,
        $self->db_user,
        $self->db_pass,
        $self->db_attrs,
       );
}

=head1 NAME

Bio::Chado::Loader - Base role for Chado Loaders

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Load data into chado.

=head1 AUTHOR

    Surya Saha				<suryasaha at cornell dot edu , @SahaSurya>   
    Lukas Mueller			<lam87 at cornell dot edu>
    Jonathan "Duke" Leto	<jonathan at leto dot net>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-chado-loader at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio::Chado::Loader>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Chado::Loader


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio::Chado::Loader>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio::Chado::Loader>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio::Chado::Loader>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio::Chado::Loader>

=back


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2010 Jonathan Duke Leto, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

1; # End of Bio::Chado::Loader
