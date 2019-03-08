#!/usr/bin/perl

=head1 NAME

format_polypeptide_names.pl

=head1 SYNOPSIS

 format_polypeptide_names.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password]

=head1 DESCRIPTION

 This script corrects polypeptide feature names and uniquenames when they are formatted like auto<feature_id> by the CHADO bulk loader. It get the name of the mRNA the polypeptide it derived from and uses it for the name and uniquename of the polypeptide record.
 You will need to run this script since bio-chado-loader-gff updates featurelocs by comparing names from ID tag of attribute field from the GFF file and feature.uniquename from the CHADO database. Both identifiers need to be correctly formatted for the comparison to work.


=head2 ARGUMENTS

 -s  Database connection string e.g. "dbi:Pg:dbname=ss_cxgn_uploadtest\;host=localhost\;port=5432" (required)
 -u  Database user name (required)
 -p  Database user password (required)
 -o  Organism name in quotes e.g. 'Solanum lycopersicum' (required)
 -d  Debug messages (0 or 1)
 -h  Help

=head1 CAVEATS

 This script *only* fixes the following
     1. feature.uniquenames with auto for polypeptide features.
     2. feature.names with polypeptide-auto names for polypeptide features.

 TODO Add -t option to test script. Compare uniquename pre and post fix.

=head1 AVAILABLITY

 https://github.com/solgenomics/bio-chado-loader/blob/master/scripts/format_feature_names.pl

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
	print "\nDatabase credentials and organism name are required. See help below\n\n\n";
	help();
}

#prep input data
my $organism = $opt_o; chomp $organism;
my $dsn = $opt_s; chomp $dsn;
my $user = $opt_u; chomp $user;
my $pass = $opt_p; chomp $pass;

if ($opt_d) {
	print STDERR "Params parsed..\n";
}

my $schema= Bio::Chado::Schema->connect($dsn, $user, $pass);
if (!$schema) { die "No schema is avaiable! \n"; }

if ($opt_d) {
	$schema->storage->debug(1);# print SQL statements
}

my $organism_id = $schema->resultset('Organism::Organism')->search(
							{'me.species' => $organism},
							)->single->organism_id();

# get cvterm for polypeptide
my $type_id = $schema->resultset('Cv::Cvterm')->search(
							{'me.name' => 'polypeptide'},
							)->single->cvterm_id();

# get polypeptides with auto uniquenames, not checking name field
my $ft_polypeptide_rs =  $schema->resultset('Sequence::Feature')
	  ->search( { 'organism_id' => $organism_id,
	  			  'type_id' => $type_id,
	  			  'uniquename' => {'like','auto%'}	},);

#get names for mRNA for selected polypeptides and store in hash
my %polypeptide_uniquename_mrna_name_hash;
while ( my $ft_polypeptide_row =  ($ft_polypeptide_rs->next ) ) {
	#find the parent mRNA feature, mRNA is the object, the protein is the subject in feature_relationship
	#the protein *should* have one parent object, which is the mRNA. Can add more checks in case the protein has more than one parent object feature
	my $ft_mrna_row = $ft_polypeptide_row->feature_relationship_subjects->single()->object();
	#print STDERR "class of \$ft_mrna_row: ",ref($ft_mrna_row),"\n";

	my $ft_mrna_row_name = $ft_mrna_row->name();
	$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->get_column('uniquename')} = $ft_mrna_row_name;
}

#update polypeptide name and uniquenames
$ft_polypeptide_rs->reset(); #reset counter to top
my $update_polypeptide_name_uniquename_sql = sub {
	my $counter = 0;
	while ( my $ft_row = $ft_polypeptide_rs->next() ) {
		$ft_row->set_column('name' => $polypeptide_uniquename_mrna_name_hash{ $ft_row->get_column('uniquename' ) } );
		$ft_row->set_column('uniquename' => ( 'polypeptide:'.$polypeptide_uniquename_mrna_name_hash{$ft_row->get_column('uniquename') } ) );
		$ft_row->update();
		$counter++;
	}
	print STDERR "$counter feature.name and uniquename rows updated for polypeptide records\n\n\n";
};

try {
	$schema->txn_do($update_polypeptide_name_uniquename_sql);
	}
	catch {# Transaction failed
		die "Could not update feature.uniquename for polypeptide records. Error: $!"
		if ( $_ =~ /Rollback failed/ );
		print STDERR "Error: " . $_;
	  };
if ($opt_d) {
	print STDERR "feature.name and uniquename rows updated for polypeptide records..\n";
}


#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script corrects feature names if they are formatted like auto<feature_id> or <name>-<feature_id>.
     This typically happens when the GMOD bulk loader is used to add the same feature more than once.
     You need to run this script since bio-chado-loader-gff updates featurelocs by comparing names from ID
     tag of attribute field from the GFF file and feature.uniquename from the CHADO database. Both identifiers
     need to be correctly formatted for the comparison to work.

    NOTE:
     This script *only* fixes the following
     1. feature.uniquenames with auto for polypeptide features.
     2. feature.uniquenames with feature_id suffix for gene,mRNA,exon,intron,ultracontig features.

     bio-chado-loader-gff may still fail if other features in CHADO have malformed uniquenames compared to your
     GFF file.


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

=head1 AUTHOR

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>
	Naama Menda <nm249@cornell.edu>

=cut

__END__
