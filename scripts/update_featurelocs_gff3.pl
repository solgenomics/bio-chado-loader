#!/usr/bin/perl

=head1 NAME

update_featurelocs_gff.pl

=head1 SYNOPSIS

update_featurelocs_gff.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password] -g [GFF3 file]

=head1 DESCRIPTION

 The new chromosome backbone has to be loaded into the database using the GMOD bulk loader before new
     coordinates can be added. The source feature is required to be present before any new featurelocs
     can be placed on it.

=head2 ARGUMENTS

 -g  GFF3 file for updating locations (required)
 -s  Database connection string e.g. "dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432" (required)
 -u  Database user name (required)
 -p  Database user password (required)
 -o  Organism name in quotes e.g. 'Solanum lycopersicum' (required)
 -r  Delete locations instead of adding
 -d  Debug messages (0 or 1)
 -h  Help

=head1 CAVEATS

 This script will break if no ID field is present in attributes in GFF record. This 
 may be a problem for non-gene/mRNA/exon/intron records although they may conform to 
 GFF3 standards.
 
=head1 AVAILABLITY

 https://github.com/solgenomics/bio-chado-loader/blob/master/scripts/update_featurelocs_gff.pl

=cut

use strict;
use warnings;

#util functions
sub mem_used {
	my ( $i, $t );
	$t = new Proc::ProcessTable;
	foreach my $got ( @{ $t->table } ) {
		next if not $got->pid eq $$;
		$i = $got->size;
	}
	print STDERR "Process id=", $$, "\n";
	print STDERR "Memory used(MB)=", $i / 1024 / 1024, "\n\n";
}

sub run_time {
	my ( $user_t, $system_t, $cuser_t, $csystem_t );
	( $user_t, $system_t, $cuser_t, $csystem_t ) = times;
	print STDERR "System time for process: $system_t\n";
	print STDERR "User time for process: $user_t\n\n";
}


use File::Slurp;
use Try::Tiny;
use Getopt::Std;
use Bio::Chado::Loader::GFF3;

our ( $opt_g, $opt_s, $opt_u, $opt_p, $opt_o, $opt_r, $opt_d, $opt_h );
getopts('g:s:u:p:o:r:d:h');
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
die "Could not open GFF file: $gff_file\n" if !(-e $gff_file);

if ($opt_d) {
	print STDERR "Params parsed..\n";
}

my $loader = Bio::Chado::Loader::GFF3->new(
       file_name => $gff_file,
       organism_name => $organism,
);

if ($opt_r) {
	$loader->delete(1);
}

if ($opt_d) {
	$loader->debug(1);# print SQL statements
}
  
#Create cache first and then parse GFF 
$loader->db_dsn($dsn);
$loader->db_user($user);
$loader->db_pass($pass);
$loader->organism_id($loader->organism_exists());

#Call populate_cache() before parse() in case parents are not in GFF but in DB already
$loader->populate__cache();
$loader->parse();
if ($opt_d) {
	run_time(); mem_used();
}

my $cnt;
$cnt=$loader->prepare_bulk_operation();
print STDERR 'Prepped '.$cnt." recs for insertion\n";

if ($opt_r) {
	$cnt = $loader->bulk_delete();
	print STDERR 'Deleted '.$cnt." records\n";
}
else{
	$loader->bulk_upload();
	print STDERR "\n\n\nupdated locgroups and inserted new rows into featureloc from $gff_file\n";
}


if ($opt_d) {
	run_time(); mem_used();
}

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     The new chromosome backbone has to be loaded into the database using the GMOD bulk loader or 
     Bio::Chado::Loader::FASTA before new coordinates can be added. The source feature is required 
     to be present before any new featurelocs can be placed on it.
     
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
 -r  Delete locations instead of adding
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


