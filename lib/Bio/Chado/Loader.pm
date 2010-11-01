package Bio::Chado::Loader;

use Bio::Chado::Schema;
use Moose::Role;

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


=head1 NAME

Bio::Chado::Loader - Base role for Chado Loaders

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Load data into chado.

=head1 AUTHOR

Jonathan Duke Leto, C<< <jonathan at leto.net> >>

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
