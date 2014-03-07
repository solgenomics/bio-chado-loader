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


# auto, polypeptides
my $ft_polypeptide_rs =  $schema->resultset('Sequence::Feature')
	  ->search( { 'organism_id' => $organism_id,
	  			  'type_id' => 22027,
	  			  'uniquename' => {'like','auto%'}	},);

#update polypeptide uniquenames
my $update_polypeptide_uniquename_sql = sub {
	my $counter = 0;
	while ( my $ft_row = $ft_polypeptide_rs->next() ) {
		$ft_row->set_column('uniquename' => ( 'polypeptide:'.$ft_row->get_column('name') ) );
		$ft_row->update();
		$counter++;
	}
	print STDERR "$counter feature.uniquename rows updated for polypeptide records\n\n\n";
};

try {
	$schema->txn_do($update_polypeptide_uniquename_sql);
	}
	catch {# Transaction failed
		die "Could not update feature.uniquename for polypeptide records. Error: $!"
		if ( $_ =~ /Rollback failed/ );
		print STDERR "Error: " . $_;
	  };
	  

#fix name-feature_id for 'gene','mRNA','exon','intron'
my $name_condition = "in ('gene','mRNA','exon','intron')";
my @cv_rs_arr = $schema->resultset('Cv::Cvterm')->search(
{'name' => \$name_condition},
{'columns' => [qw/cvterm_id/]}
)->all();

my @type_ids;
foreach my $cv_rs (@cv_rs_arr){
	push @type_ids, $cv_rs->cvterm_id();
} 

my $uniquename_condition = "like \'%-\'||feature_id";
my $ft_uniquename_rs = $schema->resultset('Sequence::Feature')->search(
									 {
									   'organism_id' => $organism_id,
									   'uniquename'  => \$uniquename_condition,
									   'type_id'     => \@type_ids,
									 },);

#update 'gene','mRNA','exon','intron' uniquenames
my $update_malformed_uniquename_sql = sub {
	my $counter = 0;
	while ( my $ft_row = $ft_uniquename_rs->next() ) {
		#get uniquename, parse out feature_id, reassign
		my $feature_id = $ft_row->get_column('feature_id');
		my $uniquename = $ft_row->get_column('uniquename');
		$uniquename =~ s/\-$feature_id$//;
		$ft_row->set_column($uniquename);
		$ft_row->update();
		$counter++;
	}
	print STDERR "$counter feature.uniquename's updated for gene,mRNA,exon,intron records\n\n\n";
};

try {
	$schema->txn_do($update_malformed_uniquename_sql);
	}
	catch {# Transaction failed
		die "Could not update feature.uniquename for gene,mRNA,exon,intron records. Error: $!"
		if ( $_ =~ /Rollback failed/ );
		print STDERR "Error: " . $_;
	  };









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





