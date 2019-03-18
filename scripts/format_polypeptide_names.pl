#!/usr/bin/perl

=head1 NAME

format_polypeptide_names.pl

=head1 SYNOPSIS

 format_polypeptide_names.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password]

=head1 DESCRIPTION

 This script corrects polypeptide feature names and uniquenames when they are formatted like auto<feature_id> by the CHADO bulk loader. It get the name of the mRNA the polypeptide it derived from and uses it for the name and uniquename of the polypeptide record.
 You will need to run this script since bio-chado-loader-gff updates featurelocs by comparing names from ID tag of attribute field from the GFF file and feature.uniquename from the CHADO database. Both identifiers need to be correctly formatted for the comparison to work.
 Sometimes the blk loader creates more than one polypeptide for one mRNA feature. This script only fixes the first one. You should delete the remaining auto named features manually (delete from feature where type_id = 75820 and uniquename like 'auto%';)


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

#get names for mRNA for selected polypeptides and store in hash, also a hash for fixed polypeptides
my %polypeptide_uniquename_mrna_name_hash;
my %mrna_fix_count_hash;

while ( my $ft_polypeptide_row =  ($ft_polypeptide_rs->next ) ) {
	#find the parent mRNA feature, mRNA is the object, the protein is the subject in feature_relationship
	#the protein *should* have one parent object, which is the mRNA. Can add more checks in case the protein has more than one parent object feature
	my $ft_mrna_row = $ft_polypeptide_row->feature_relationship_subjects->single()->object();
	#print STDERR "class of \$ft_mrna_row: ",ref($ft_mrna_row),"\n";

	my $ft_mrna_row_name = $ft_mrna_row->name();
	$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename()} = $ft_mrna_row_name;
	$mrna_fix_count_hash{$ft_mrna_row_name} = 0;
}

# print STDERR "Values in polypeptide_uniquename_mrna_name_hash: ".length(keys %polypeptide_uniquename_mrna_name_hash)."\n";
# print STDERR "Values in mrna_fix_count_hash: ".length(keys %mrna_fix_count_hash)."\n";

#update polypeptide name and uniquenames
$ft_polypeptide_rs->reset(); #reset counter to top
my $update_polypeptide_name_uniquename_sql = sub {
	my $counter = 0;
	my $multipopypeptides_counter = 0;
	while ( my $ft_polypeptide_row = $ft_polypeptide_rs->next() ) {
		#some polypeptides are derived from same mRNA so checking if the polypeptide for this mRNA has already been fixed
		#the remaining auto polypeptides for this organism can be safely deleted since they have same coords in featureloc
		#print STDERR "polypep uniqname before fix: ".$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename()}."\n";
		if ( $mrna_fix_count_hash{$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename()}} == 0 )
		{
			#print STDERR "polypep uniqname during fix: ".$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename()}."\n";
			#print STDERR "mRNA name during fix: ".$mrna_fix_count_hash{$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename()}}."\n";
			$mrna_fix_count_hash{ $polypeptide_uniquename_mrna_name_hash{ $ft_polypeptide_row->uniquename() } }++;
			$ft_polypeptide_row->set_column('name' => $polypeptide_uniquename_mrna_name_hash{ $ft_polypeptide_row->uniquename() } );
			$ft_polypeptide_row->set_column('uniquename' => ( 'polypeptide:'.$polypeptide_uniquename_mrna_name_hash{$ft_polypeptide_row->uniquename() } ) );
			$ft_polypeptide_row->update();
			$counter++;
		}
		else{
			$multipopypeptides_counter++;
		}
	}
	print STDERR "$multipopypeptides_counter mRNAs have more than 1 polypeptide. Please check and delete them in the feature and featureloc tables\n";
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

     This script corrects polypeptide feature names and uniquenames when they are formatted like auto<feature_id> by the CHADO bulk loader. It get the name of the mRNA the polypeptide it derived from and uses it for the name and uniquename of the polypeptide record.
	 You will need to run this script since bio-chado-loader-gff updates featurelocs by comparing names from ID tag of attribute field from the GFF file and feature.uniquename from the CHADO database. Both identifiers need to be correctly formatted for the comparison to work.
	 Sometimes the blk loader creates more than one polypeptide for one mRNA feature. This script only fixes the first one. You should delete the remaining auto named features manually (delete from feature where type_id = 75820 and uniquename like 'auto%';)


    Usage:
      format_polypeptide_names.pl -o ["organism"] -s ["DSN string"] -u [DB user] -p [password]

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
