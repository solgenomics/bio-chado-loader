#!/usr/bin/perl

=head1 NAME

validate_gff.pl

=head1 SYNOPSIS

format_feature_names.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password]

=head1 COMMAND-LINE OPTIONS

 -s  Database connection string e.g. "dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432" (required)
 -u  Database user name (required)
 -p  Database user password (required)
 -o  Organism name in quotes e.g. 'Solanum lycopersicum' (required)
 -d  Debug messages (0 or 1)
 -h  Help

=cut

use strict;
use warnings;


use Try::Tiny;
use Getopt::Std;
use Bio::Chado::Schema;

our ( $opt_s, $opt_u, $opt_p, $opt_o, $opt_d, $opt_h );
getopts('s:u:p:o:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_s || !$opt_u || !$opt_p || !$opt_o ) {
	print "\nDatabse credentials and organism name are required. See help below\n\n\n";
	help();
}

#prep input data
my $organism = $opt_o;

my $gff_updated_file = $opt_u;
my $gff_updated      = read_file($gff_updated_file)
  or die "Could not open updated GFF file: $gff_updated_file\n";
my $fasta_updated_file;

if ($opt_d) {
	print STDERR "Params parsed..\n";
}













#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script corrects feature names if they are formatted like auto<feature_id> or <name>-<feature_id>. This typically happens when the GMOD bulk loader is used to add the same feature more than once.
     
    NOTE:
     You need to run this script since bio-chado-loader-gff updates featurelocs by comparing names from ID tag of attribute field from the GFF file and feature.uniquename from the CHADO database. Both identifiers need to be correctly formatted for the comparison to work. 

    Usage:
      format_feature_names.pl -o ["organism"] -s [DSN string] -u [DB user] -p [password]
      
    Flags:

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

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__





