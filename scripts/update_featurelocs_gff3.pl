#!/usr/bin/perl

=head1 NAME

update_featurelocs_gff.pl

=head1 SYNOPSIS

update_featurelocs_gff.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password] -g [GFF3 file]

=head1 ARGUMENTS

 -g  GFF3 file for updating locations (required)
 -s  Database connection string e.g. "dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432" (required)
 -u  Database user name (required)
 -p  Database user password (required)
 -o  Organism name in quotes e.g. 'Solanum lycopersicum' (required)
 -d  Debug messages (0 or 1)
 -h  Help

=head1 AVAILABLITY

https://github.com/solgenomics/bio-chado-loader/blob/master/scripts/update_featurelocs_gff.pl

=cut

use strict;
use warnings;

use File::Slurp;
use Try::Tiny;
use Getopt::Std;
use Bio::Chado::Schema;

our ( $opt_g, $opt_s, $opt_u, $opt_p, $opt_o, $opt_d, $opt_h );
getopts('g:s:u:p:o:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_g || !$opt_s || !$opt_u || !$opt_p || !$opt_o ) {
	print "\nGFF3 file, database credentials and organism name are required. See help below\n\n\n";
	help();
}

#prep input data
my $organism = $opt_o; chomp $organism;
my $dsn = $opt_s; chomp $dsn;
my $user = $opt_u; chomp $user;
my $pass = $opt_p; chomp $pass;
my $gff_file = $opt_g;
my $gff      = read_file($gff_file)
  or die "Could not open GFF file: $gff_file\n";

if ($opt_d) {
	print STDERR "Params parsed..\n";
}

my $schema= Bio::Chado::Schema->connect($dsn, $user, $pass);
if (!$schema) { die "No schema is avaiable! \n"; }

if ($opt_d) {
	$schema->storage->debug(1);# print SQL statements
}






#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     The new chromosome backbone has to be loaded into the database using the GMOD bulk loader before new
     coordinates can be added. The source feature is required to be present before any new featurelocs
     can be placed on it.
     
    NOTE:
     Will break if no ID field in attributes in GFF record 

    Usage:
      update_featurelocs_gff.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password] -g [GFF3 file]
      
    Flags:

 -g  GFF3 file for updating locations (required)
 -s  Database connection string e.g. "dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432" (required)
 -u  Database user name (required)
 -p  Database user password (required)
 -o  Organism name in quotes e.g. 'Solanum lycopersicum' (required)
 -d  Debug messages (0 or 1)
 -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHOR

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__


