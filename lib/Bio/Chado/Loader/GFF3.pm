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

The new chromosome backbone has to be loaded into the database using the GMOD bulk loader before new 
coordinates can be added. You can also use Bio::Chado::Loader::Fasta and create new features during 
the process using --create_features. The source feature is required to be present before any new 
featurelocs can be placed on it.

CAVEAT
Will break if no ID field in attributes in GFF record

=over

=cut

require Bio::GFF3::LowLevel::Parser;

use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';

use namespace::autoclean;
use Proc::ProcessTable;
use Bio::GFF3::LowLevel::Parser;
use Bio::GFF3::LowLevel
  qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);
use Try::Tiny;

has schema => (
    is         => 'rw',
    isa        => 'Bio::Chado::Schema',
    traits     => ['NoGetopt'],
    lazy_build => 1,
);

has 'file_name' => (
	documentation =>
		  'name of gff3 file (required)',
	isa       => 'Str',
	is        => 'rw',
	predicate => 'has_file_name',
	clearer   => 'clear_file_name',
	required => 1,
);

with 'Bio::Chado::Loader';

has '+organism_name' => (
	documentation =>
		  'exact species of organism, as stored in database, e.g. "Solanum lycopersicum" (required)',
    required => 1,
);

has 'organism_id' => (
		documentation =>
		  'id organism, as stored in database, e.g. 1 for Solanum lycopersicum',
		is  => 'rw',
		isa => 'Str',
);

has 'debug' => (
				 documentation => 'Verbosity flag (0/1)',
				 is            => 'rw',
				 isa           => 'Int',
				 default       => 0,
);

has 'delete' => (
	  documentation => 'Delete features in GFF instead of inserting them (0/1)',
	  is            => 'rw',
	  isa           => 'Int',
	  default       => 0,
);

# hash of feature uniquename => hashref of entire GFF record. 
# It will break if no ID field in GFF record. uniquenames are used for 
# key since GFF3 mandates they be unique in a GFF3 file.
has '_features_gff' => (
	is      => 'ro',
	isa     => 'HashRef',
	traits  => ['Hash'],
	default => sub { {} },
	handles => {
				 add__features_gff    => 'set',
				 count__features_gff  => 'count',
				 _features_gff_exists => 'exists',
	},
);

=item C<dump__features_gff ()>

Prints the contents of the _features_gff to STDERR. For internal debugging use only. 

=cut

sub dump__features_gff {
	my ($self) = @_;
	my ( $feature_uniquename, $fields );

	while ( ( $feature_uniquename, $fields ) = each %{ $self->_features_gff } ) {
		print STDERR $feature_uniquename . " => "
		  . gff3_format_feature($fields);
	}
}

# hash of featureID => feature_uniquename for all features in GFF 
# that will be updated. Populated in prepare_bulk_operation
has '_feature_ids_uniquenames_gff' => (
	is      => 'ro',
	isa     => 'HashRef',
	traits  => ['Hash'],
	default => sub { {} },
	handles => {
				 add__feature_ids_uniquenames_gff   => 'set',
				 count__feature_ids_uniquenames_gff => 'count',
	},
);

# hash of cvterm => 1 from GFF. TOFIX Use counts instead of constant 1
has '_cvterms_gff' => (
	   is      => 'ro',
	   isa     => 'HashRef[Str]',
	   traits  => ['Hash'],
	   default => sub { {} },
	   handles => {
					add__cvterms_gff   => 'set',
					count__cvterms_gff => 'count',
	   },
);

has 'is_analysis' => (
			documentation =>
			  'set true if this feature should be recorded as from an analysis',
			is      => 'rw',
			isa     => 'Bool',
			default => 0,
);

# hash of cvterm.name => hash of feature.uniquename => feature.feature_id,featureloc.srcfeature_id. 
# Not recording featureloc.locgroup,featureloc.rank
has '_cache' => (
	is      => 'ro',
	isa     => 'HashRef',
	default => sub { {} },
);

# hash of feature uniquenames => feature_id for all features in _cache
has '_feature_uniquename_feature_id_cache' => (
		is      => 'ro',
		isa     => 'HashRef',
		traits  => ['Hash'],
		default => sub { {} },
		handles => {
					 add__feature_uniquename_feature_id_cache    => 'set',
					 count__feature_uniquename_feature_id_cache  => 'count',
					 _feature_uniquename_feature_id_cache_exists => 'exists'
		},
);

use autodie qw(:all);
use Data::Dumper;
#use 5.010;

sub _usage {

return <<USAGE;

At least one GFF3 file must be given to load.

USAGE

}

=item C<run ()>

TODO For mappable cmd line exec// See Fasta.pm

=cut

sub run {
	my ( $self, %args ) = @_;

	exit 0;
}

=item C<parse ()>

Parse a GFF3 file using Bio::GFF3::LowLevel::Parser.

=cut

sub parse {
	my ( $self, %args ) = @_;

	my $gff = Bio::GFF3::LowLevel::Parser->open( $self->file_name );
	while ( my $item = $gff->next_item ) {
		if ( ref $item eq 'ARRAY' ) {
			## $i is an arrayref of feature lines that have the same ID,
			## in the same format as returned by Bio::GFF3::LowLevel::gff3_parse_feature
			#warn Dumper [ $i ];
			for my $feature (@$item) {

				# for each location of this feature
				$self->parse_feature($feature);
			}
		}
		elsif ( $item->{directive} ) {
			if ( $item->{directive} eq 'FASTA' ) {
				my $fasta_filehandle = $item->{filehandle};
				## parse the FASTA in the filehandle with BioPerl or
				## however you want.  or ignore it.
				print "got fasta handle\n";

				#use bio-chado-loader-fasta to load fasta
			}
			elsif ( $item->{directive} eq 'gff-version' ) {
				print "it says it is GFF version $item->{value}\n";
			}
			elsif ( $item->{directive} eq 'sequence-region' ) {
				print( "found a sequence-region, sequence $item->{seq_id},",
					   " from $item->{start} to $item->{end}\n" );
			}
		}
		elsif ( $item->{comment} ) {
			## this is a comment in your GFF3 file, in case you want to do
			## something with it.
			print "that comment said: '$item->{comment}'\n";
		}
		else {
			die 'this should never happen!';
		}
	}
}

=item C<parse_feature ()>

Parse a feature HashRef to populate data structures for self and children/derived features. Calls itself recursively if there are child or derived features. Calls special mRNA handler if feature is mRNA.

CAVEAT: This will exit if no ID field is present

=cut

sub parse_feature {
	my $self         = shift;
	my $feature_hash = shift;
	
	#warn [Dumper $feature_hash];

	# add data
	if ($self->_features_gff_exists($feature_hash->{'attributes'}->{'ID'}->[0]))
	{
		#warn Dumper [ $self->_features_gff ];
		die "Multiple features with same ID as "
		  . $feature_hash->{'attributes'}->{'ID'}->[0]
		  . " found. Exiting..";
	}
	elsif ( !defined $feature_hash->{'attributes'}->{'ID'}->[0] )
	{    #exit if no ID field
		die "No ID defined for feature in GFF. Exiting..";
	}
	elsif ( $feature_hash->{'type'} eq 'mRNA' )
	{    #call mRNA handler for parsing exons, CDSs and UTRs
		$self->parse_mRNA_feature($feature_hash);
		return;
	}
	else {
		$self->add__features_gff($feature_hash->{'attributes'}->{'ID'}->[0] => $feature_hash );
		$self->add__cvterms_gff( $feature_hash->{'type'} => 1 );

		#warn Dumper [ $self->_features_gff ] if $feature_hash->{'attributes'}->{'ID'}->[0]=~ /exon/;
	}

	#recursively calling self for nested child features
	if ( $feature_hash->{'child_features'} ) {
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			$self->parse_feature( $feature_child->[0] );
		}
	}

	#recursively calling self for nested derived features
	if ( $feature_hash->{'derived_features'} ) {
		for my $feature_derived ( @{ $feature_hash->{'derived_features'} } ) {
			$self->parse_feature( $feature_derived->[0] );
		}
	}
}

=item C<parse_mRNA_feature ()>

Adds mRNA, exons and polypeptide to _features_gff if exons, CDSs and UTRs are present. Adds mRNA and creates exon & polypeptide features with CDS and UTR coordinates if both are present. Dies if only CDS or UTR is present.  

=cut

sub parse_mRNA_feature {
	my $self         = shift;
	my $feature_hash = shift;

	my %children = (
					 'exon'            => 0,
					 'intron'          => 0,
					 'CDS'             => 0,
					 'five_prime_UTR'  => 0,
					 'three_prime_UTR' => 0,
	);

	#test children presence/absence
	if ( $feature_hash->{'child_features'} ) {
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			die "mRNA has child feature other than five_prime_UTR, three_prime_UTR, CDS, exon, intron .. cannot handle "
			  if ( !exists $children{ $feature_child->[0]->{'type'} } );
			$children{ $feature_child->[0]->{'type'} }++;

			#$self->parse_feature($feature_child->[0]);
		}
	}
	else {
		if ( $self->debug ) {
			print STDERR "No child features saved for "
			  . $feature_hash->{'attributes'}->{'ID'}->[0] . "\n";
		}

		return;
	}

	#introns
	if ( $children{'intron'} ) {

		#save exons
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if ( $feature_child->[0]->{'type'} eq 'intron' ) {
				$self->parse_feature( $feature_child->[0] );
			}
		}
	}

	#polypeptide
	my ( @starts, @ends, $gene );
	$gene = $feature_hash->{'attributes'}->{'ID'}->[0];
	$gene =~ s/^mRNA\://;

	my $polypeptide_feature = {
		seq_id => $feature_hash->{'seq_id'},
		source => $feature_hash->{'source'},
		type   => 'polypeptide',
		score      => undef,
		strand     => $feature_hash->{'strand'},
		phase      => undef,
		attributes => {
				   ID     => [ 'polypeptide:' . $gene ],
				   Parent => [ $feature_hash->{'attributes'}->{'Parent'}->[0] ],
		},
		child_features   => [],
		derived_features => [],
	};

	if ( $children{'CDS'} ) {
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if ( $feature_child->[0]->{'type'} eq 'CDS' ) {
				push @starts, $feature_child->[0]->{'start'} - 1;
				push @ends,   $feature_child->[0]->{'end'};
			}
		}
	}
	elsif ( $children{'five_prime_UTR'} && $children{'three_prime_UTR'} ) {

		#Do we check if >1 3' or 5' UTRs???
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if ( $feature_child->[0]->{'type'} eq 'five_prime_UTR' ) {
				push @starts, $feature_child->[0]->{'start'} - 1;
			}
			elsif ( $feature_child->[0]->{'type'} eq 'three_prime_UTR' ) {
				push @ends, $feature_child->[0]->{'end'};
			}
		}
	}
	$polypeptide_feature->{'start'} = ( sort { $a <=> $b } @starts )[0];   #asc
	$polypeptide_feature->{'end'}   = ( sort { $b <=> $a } @ends )[0];     #desc
	#add polypeptide
	if ( $self->debug ) {
		print STDERR "\rCreating polypeptide feature for "
		  . $feature_hash->{'attributes'}->{'ID'}->[0];
	}

	$self->parse_feature($polypeptide_feature);

	#exons
	if ( $children{'exon'} ) {

		#save exons
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if ( $feature_child->[0]->{'type'} eq 'exon' ) {
				$self->parse_feature( $feature_child->[0] );
			}
		}
	}
	elsif ( !$children{'exon'}) {# this will not happen for tomato GFF3s
		if ( $self->debug ) {
			print STDERR "\rCreating exon features for "
			  . $feature_hash->{'attributes'}->{'ID'}->[0];
		}

		#check if only 1 3' and only 1 5' UTR
		die "More than one or no five_prime_UTR or three_prime_UTR for "
		  . $feature_hash->{'attributes'}->{'ID'}->[0]
		  if (    ( $children{'five_prime_UTR'} != 1 )
			   || ( $children{'three_prime_UTR'} != 1 ) );

		#get coords
		my ( %cds_starts, %cds_ends, $five_prime_UTR_start,
			 $three_prime_UTR_end );
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if ( $feature_child->[0]->{'type'} eq 'CDS' ) {
				$cds_starts{ $feature_child->[0]->{'attributes'}->{'ID'}->[0] }
				  = $feature_child->[0]->{'start'} - 1;
				$cds_ends{ $feature_child->[0]->{'attributes'}->{'ID'}->[0] } =
				  $feature_child->[0]->{'end'};
			}
			elsif ( $feature_child->[0]->{'type'} eq 'five_prime_UTR' ) {
				$five_prime_UTR_start = $feature_child->[0]->{'start'} - 1;
			}
			elsif ( $feature_child->[0]->{'type'} eq 'three_prime_UTR' ) {
				$three_prime_UTR_end = $feature_child->[0]->{'end'};
			}
		}

		#check if 5' before 1st CDS and 3' after last CDS
		die "five_prime_UTR starts at or after first CDS for "
		  . $feature_hash->{'attributes'}->{'ID'}->[0]
		  if ( $five_prime_UTR_start >=
			   ( ( sort { $a <=> $b } values %cds_starts )[0] ) );    #asc
		die "three_prime_UTR ends at or before last CDS for "
		  . $feature_hash->{'attributes'}->{'ID'}->[0]
		  if ( $three_prime_UTR_end <=
			   ( ( sort { $b <=> $a } values %cds_ends )[0] ) );      #desc

		#create exons
		my @cds_starts_values =
		  values %cds_starts
		  ; #need to do this outside each() loop since values uses same internal iterator
		while ( my ( $CDS_ID, $value ) = each(%cds_starts) ) {
			my $CDS_ID_trimmed = $CDS_ID;
			$CDS_ID_trimmed =~ s/^CDS\://;
			my $exon_feature = {
				seq_id => $feature_hash->{'seq_id'},
				source => $feature_hash->{'source'},
				type   => 'exon',
				score      => undef,
				strand     => $feature_hash->{'strand'},
				phase      => undef,
				attributes => {
					   ID     => [ 'exon:' . $CDS_ID_trimmed ],
					   Parent => [ $feature_hash->{'attributes'}->{'ID'}->[0] ],
				},
				child_features   => [],
				derived_features => [],
			};

			#if first CDS
			if ( $value == ( ( sort { $a <=> $b } @cds_starts_values )[0] ) ) {
				$exon_feature->{'start'} = $five_prime_UTR_start;
				$exon_feature->{'end'}   = $cds_ends{$CDS_ID};
			}
			elsif (
				  $value == ( ( sort { $b <=> $a } (@cds_starts_values) )[0] ) )
			{    #if last CDS
				$exon_feature->{'start'} = $value;
				$exon_feature->{'end'}   = $three_prime_UTR_end;
			}
			else {    #CDS in the middle
				$exon_feature->{'start'} = $value;
				$exon_feature->{'end'}   = $cds_ends{$CDS_ID};
			}

			#save exon
			#warn Dumper [ $exon_feature];
			$self->parse_feature($exon_feature);
		}
	}

	#remove mRNA sub features info to save memory
	$feature_hash->{'child_features'}   = [];
	$feature_hash->{'derived_features'} = [];

	#save mRNA
	$self->add__features_gff($feature_hash->{'attributes'}->{'ID'}->[0] => $feature_hash );
	$self->add__cvterms_gff( $feature_hash->{'type'} => 1 );

	if ( $self->debug ) {
		print STDERR "\n";
	}
}

=item C<organism_exists ()>

Check if organism for GFF3 exists in the Chado database. Return organism::organism. organism_id 
if present otherwise 0. No organism should have id of 0.

=cut

sub organism_exists {
	my ($self) = @_;

	my $org_id =
	  $self->schema->resultset('Organism::Organism')
	  ->search( { species => $self->organism_name } )->get_column('organism_id')
	  ->first();
	return $org_id ? $org_id : 0;
}

=item C<populate__cache ()>

Populate _cache hash with type_id->feature.uniquename->feature.feature_id,featureloc.srcfeature_id,
featureloc.locgroup,featureloc.rank,cvterm.name from cvterm,feature,featureloc relations. Call 
function before parse() in case parents are not in GFF but in DB already.  

CAVEAT: Count may be over estimation if some features have multiple locgroups(e.g. contigs).

=cut

sub populate__cache {
	my ($self) = @_;
	
	if ( $self->debug == 2 ) { #very high verbosity
		$self->schema->storage->debug(1);#print SQL statements
	}

	#setup DB dsn's
	#check names for only cvterms read from GFF file
	#will fail for remark records from assembly.gff3, are they needed in DB??
	my $ft_auto_err_rs = $self->schema->resultset('Sequence::Feature')->search(
									 {
									   'me.organism_id' => $self->organism_id,
									   'me.uniquename'  => { 'like', "%auto%" },
									   'type.name'      =>
										 [ keys %{ $self->_cvterms_gff } ],
									 },
									 { join => ['type'], prefetch => ['type'] });

	die "There are features in database with auto in feature.uniquename field. 
    	Please correct this in your database before running this script.
    	Use format_feature_names.pl in scripts dir to change the names.
    	This typically happens when the GMOD bulk loader is used to add the 
    	same feature more than once. Exiting..." if $ft_auto_err_rs->count() > 0;
    
    my $uniquename_condition = "like \'%-\'||feature_id";
    my $ft_uniquename_err_rs = $self->schema->resultset('Sequence::Feature')->search(
									 {
									   'me.organism_id' => $self->organism_id,
									   'me.uniquename'  => \$uniquename_condition,
									 },);

	die "There are features in database with feature_id in feature.uniquename field. 
    	Please correct this in your database before running this script.
    	Use format_feature_names.pl in scripts dir to change the names.
    	This typically happens when the GMOD bulk loader is used to add the 
    	same feature more than once. Exiting..." if $ft_uniquename_err_rs->count() > 0;
	

	my $fl_rs = $self->schema->resultset('Sequence::Featureloc')->search(
		#    	{ 'organism_id'=> $self->organism_id }, #for full DB
		{
		   'organism_id'        => $self->organism_id,
		   #'feature.uniquename' => { 'like', '%Solyc01g1123%' } #for testing, only 86 floc records
		   #'feature.uniquename' => { 'like', '%Solyc01g%' }, #50684 floc records
		   'feature.uniquename' => { 'like', '%Solyc01g0%' }, #?? floc records
		},    
		#{ 'organism_id'=> $self->organism_id , 'feature.uniquename' => { 'like', '%dummy%'}},#for testing, only few floc records
		{ join => ['feature'], prefetch => ['feature'] });

	my $ft_rs = $self->schema->resultset('Sequence::Feature')->search(
										{ 'organism_id' => $self->organism_id },
										{
										  join     => ['type'],
										  prefetch => ['type']
										});

	my $ft_fname_fid_rs = $self->schema->resultset('Sequence::Feature')->search(
									   { 'organism_id' => $self->organism_id },
									   {
										 columns => [qw/uniquename feature_id/]
									   },
												   );
	if ( $self->debug ) {
		print STDERR "fl_rs count: " . $fl_rs->count() . "\n";
		print STDERR "ft_rs count: " . $ft_rs->count() . "\n";
		print STDERR "ft_fname_fid_rs count: " . $ft_rs->count() . "\n";
	}

	my $count = 0;

	#create _feature_uniquename_feature_id_cache
	while ( my $ft_fname_fid_row = $ft_fname_fid_rs->next() ) {

		#add feature_id to _feature_uniquename_feature_id_cache
		$self->add__feature_uniquename_feature_id_cache(
			  $ft_fname_fid_row->uniquename() => $ft_fname_fid_row->feature_id() );
		$count++;
	}

	if ( $self->debug ) {
		print STDERR "Added "
		  . $self->count__feature_uniquename_feature_id_cache()
		  . " records to _feature_uniquename_feature_id_cache\n";
	}

	#create _cache hash
	$count = 0;
	while ( my $fl_row = $fl_rs->next() ) {
		if ( $ft_rs->search({
							   'type_id'    => $fl_row->feature->type_id(),
							   'uniquename' => $fl_row->feature->uniquename
							 })->count() == 1)
		{
			my $ft_row = $ft_rs->search(
								 {
								   'type_id'    => $fl_row->feature->type_id(),
								   'uniquename' => $fl_row->feature->uniquename
								 }
			)->first();

			$self->_cache->{ $ft_row->type->name() }
			  ->{ $fl_row->feature->uniquename() }->{feature_id} =
			  $fl_row->feature->feature_id();

#			$self->_cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{locgroup}=$fl_row->locgroup();
#			$self->_cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{rank}=$fl_row->rank();
			$self->_cache->{ $ft_row->type->name() }
			  ->{ $fl_row->feature->uniquename() }->{srcfeature_id} =
			  $fl_row->srcfeature_id();

#			#add feature_id to _feature_uniquename_feature_id_cache
#			$self->add__feature_uniquename_feature_id_cache( $fl_row->feature->uniquename() => $fl_row->feature->feature_id() );

			$count++;
			if ( $count % 100000 == 0 ) {
				print STDERR "\rprocessing $count";

				#warn Dumper [ $self-> _cache -> {$ft_row->type->name()}->{$fl_row->feature->uniquename()}];
				#warn Dumper [ $self-> _cache -> {$ft_row->type->name()}];
				my ( $i, $t );
				$t = new Proc::ProcessTable;
				foreach my $got ( @{ $t->table } ) {
					next if not $got->pid eq $$;
					$i = $got->size;
				}
				print STDERR "\nMemory used: " . ( $i / 1024 / 1024 ) . "MB \n";
				my ( $user_t, $system_t, $cuser_t, $csystem_t );
				( $user_t, $system_t, $cuser_t, $csystem_t ) = times;
				print STDERR "System time used: $system_t\n";
				print STDERR "User time used: $user_t\n";
			}
		}
		else {
			die
"Multiple features found for organism with same type_name and uniquename which violates Chado constraints. Exiting..";
		}
	}
	print STDERR "Added " . $count . " records to _cache\n";

	#return ($count);
}

=item C<prepare_bulk_operation ()>

Compare %features from GFF to %_cache. Prepare data structures to write to DB. GFF %features 
not found in %_cache are written to gff3.exceptions. 

CAVEAT: Only new featureloc's are handled right now.

=cut

sub prepare_bulk_operation {
	my ($self) = @_;

	#compare %features to %_cache
	my ( $feature_uniquename, $fields, $feature_gff, $disk_exception_str );
	my %counters = ( 'exceptions' => 0, 'inserts' => 0 );
	$disk_exception_str = '';

	#warn Dumper [ $self->_feature_uniquename_feature_id_cache ];
	while ( ( $feature_uniquename, $fields ) = each %{ $self->_features_gff } ) {

		#if feature and seq_id/scrfeature in _cache, try _feature_uniquename_feature_id_cache_exists??
		if (( $self->_cache->{ $fields->{'type'} }->{$feature_uniquename} )
			 && ( $self->_feature_uniquename_feature_id_cache
				  ->{ $fields->{'seq_id'} } ) )
		{

			#record feature_ids in class var
			$self->add__feature_ids_uniquenames_gff(
					 $self->_cache->{ $fields->{'type'} }->{$feature_uniquename}
					   ->{'feature_id'} => $fields->{'attributes'}->{'ID'}->[0] );

			$counters{'inserts'}++;
			if ( $counters{'inserts'} % 10000 == 0 ) {
				print STDERR "\r" . $counters{'inserts'} . " processed..\n";
			}
		}
		else {

			#create exception file str
			#warn Dumper [ $fields ];
			$disk_exception_str .= gff3_format_feature($fields);
			$counters{'exceptions'}++;
		}
	}
	die "No features to update." if $counters{'inserts'} == 0;

	if ( $counters{'exceptions'} > 0 ) {
		my ($exception_gff_fh);
		open( $exception_gff_fh, ">", $self->file_name . 'exceptions' )
		  or die("Could not create exception gff file: $!");
		print $exception_gff_fh $disk_exception_str;
		close($exception_gff_fh);
	}

	#warn Dumper [ $self->_feature_ids_uniquenames_gff];

	print STDERR $counters{'inserts'} . " records prepared for bulk operation\n";
	print STDERR $counters{'exceptions'} . " exception GFF records written to "
	  . $self->file_name. "\.exceptions \n";

	return $counters{'inserts'};
}

=item C<bulk_upload ()>

Insert content from in-memory data structures into DB using transactions. DB remains unmodified in case of update error.

=cut

sub bulk_upload {
	my ($self) = @_;

#   From current bulk loader - @tables array sets the order for which things will be inserted into the database
#	my @tables = (
#	   "organism",
#	   "analysis",
#	   "db",
#	   "dbxref",
#	   "cv",
#	   "cvterm",
#	   "feature",
#	   "featureloc",
#	   "feature_relationship",
#	   "featureprop",
#	   "feature_cvterm",
#	   "synonym",
#	   "feature_synonym",
#	   "feature_dbxref",
#	   "analysisfeature",
#	);

	$self->bulk_featureloc_upload();
}

=item C<bulk_featureloc_upload ()>

Insert content into featureloc from in-memory data structures into DB using transactions. 
DB remains unmodified in case of update error.

CAVEAT:
Gets srcfeature_id from _cache which is incorrect if parent was changed. Presuming that 
new location is the primary location of feature (locgroup=0, rank=0). Old locgroups for 
features being inserted are incremented by 1.

=cut

sub bulk_featureloc_upload {
	my ($self) = @_;

	#$self->schema->storage->debug(1);##SQL statement
	#warn Dumper[keys $self->_feature_ids_uniquenames_gff];
	my @feature_ids_to_update = keys( $self->_feature_ids_uniquenames_gff );
	my $fl_rs                 =
	  $self->schema->resultset('Sequence::Featureloc')
	  ->search( { 'feature_id' => \@feature_ids_to_update },
				{ 'order_by' => { -desc => [qw/feature_id locgroup/] } } );

	#create new rows
	my $create_sql = sub {
		my $counter = 0;
		while ( my ( $feature_id, $feature_uniquename ) =
				each $self->_feature_ids_uniquenames_gff )
		{
			my ( $strand, $phase );

			if ( $self->_features_gff->{$feature_uniquename}->{'strand'} eq '+' )
			{$strand = 1;}
			elsif ($self->_features_gff->{$feature_uniquename}->{'strand'} eq '-' )
			{$strand = -1;}
			else { $strand = 0; }

			if ( $self->_features_gff->{$feature_uniquename}->{'phase'} ) {
				$fl_rs->create(
					{
					   'feature_id' => $feature_id,

						#getting srcfeature_id from _cache which is incorrect if parent was changed
						#'srcfeature_id' => $self->_cache->{$self->_features_gff->{$feature_uniquename}->{'type'}}
						#	->{$feature_uniquename}->{srcfeature_id},
		
						#get feature_id of seq_id from _feature_uniquename_feature_id_cache
					   'srcfeature_id' =>	 $self->_feature_uniquename_feature_id_cache->{
										   $self->_features_gff->{$feature_uniquename}->{'seq_id'}
										 },
					   'fmin' => $self->_features_gff->{$feature_uniquename}->{'start'} - 1,
					   'fmax' => $self->_features_gff->{$feature_uniquename}->{'end'},
					   'strand' => $strand,
					   'phase'  => $self->_features_gff->{$feature_uniquename}->{'phase'},
					   'locgroup' => 0,
					   'rank'     => 0,
					}
				);
			}
			else {
				$fl_rs->create(
					{
					   'feature_id' => $feature_id,

				#getting srcfeature_id from _cache which is incorrect if parent was changed
				#'srcfeature_id' => $self->_cache->{$self->_features_gff->{$feature_uniquename}->{'type'}}
				#	->{$feature_uniquename}->{srcfeature_id},

				#get feature_id of seq_id from _feature_uniquename_feature_id_cache
					   'srcfeature_id' => $self->_feature_uniquename_feature_id_cache->{
										   $self->_features_gff->{$feature_uniquename}->{'seq_id'}
										 },
					   'fmin' => $self->_features_gff->{$feature_uniquename}->{'start'} - 1,
					   'fmax' => $self->_features_gff->{$feature_uniquename}->{'end'},
					   'strand'   => $strand,
					   'locgroup' => 0,
					   'rank'     => 0,
					}
				);
			}
			$counter++;
		}
		print STDERR "$counter rows added to featureloc\n";
	};

	#update locgroup=locgroup+1
	my $increment_locgroup_sql = sub {
		my $counter = 0;
		while ( my $fl_row = $fl_rs->next() ) {
			$fl_row->set_column('locgroup' => ( $fl_row->get_column('locgroup') + 1 ) );
			$fl_row->update();
			$counter++;
		}
		print STDERR "$counter featureloc rows locgroup fields updated\n";
		$self->schema->txn_do($create_sql);    # nested transaction for create
	};

	#update locgroups and insert rows into featureloc using transactions
	try {
		$self->schema->txn_do($increment_locgroup_sql);
	  }
	  catch {# Transaction failed
		die "Could not increment locgroups and/or create new rows. Error: $!"
		  if ( $_ =~ /Rollback failed/ );

		print STDERR "Error: " . $_;
	  };

}

=item C<bulk_delete ()>

Delete in database with content from in-memory data structures using transactions. DB remains unmodified in case of update error.

=cut

sub bulk_delete {
	my ($self) = @_;

	#check with user
	my $input;
	print STDERR "Are you sure you want to delete ".$self->count__feature_ids_uniquenames_gff()." rows (yes/no) : ";
	while ($input=<>) {
		last if ( ( $input eq "yes\n" ) or ( $input eq "no\n" ) );
		print STDERR "Please enter yes or no : ";
	}
	
	die "Exiting...." if $input eq "no\n";

#   From current bulk loader - @tables array sets the order for which things will be inserted into the database
#	my @tables = (
#	   "organism",
#	   "analysis",
#	   "db",
#	   "dbxref",
#	   "cv",
#	   "cvterm",
#	   "feature",
#	   "featureloc",
#	   "feature_relationship",
#	   "featureprop",
#	   "feature_cvterm",
#	   "synonym",
#	   "feature_synonym",
#	   "feature_dbxref",
#	   "analysisfeature",
#	);

	$self->bulk_featureloc_delete();
}

=item C<bulk_featureloc_delete ()>

DELETE content from featureloc in DB from in-memory data structures using transactions. 
DB remains unmodified in case of update error. Decrements locgroup of all remaining 
featureloc records for the feature with locgroups higher than the deleted record. 
Returns the number of floc records successfully deleted.

=cut

sub bulk_featureloc_delete {
	my ($self) = @_;

	#$self->schema->storage->debug(1);##SQL statement
	my $counter=0;
	while ( my ( $feature_id, $feature_uniquename ) = each $self->_feature_ids_uniquenames_gff ){

			my $srcfeature_id = $self-> schema -> resultset('Sequence::Feature')->search(
			{
			  'organism_id'=> $self->organism_id , 
			  'uniquename' => { 'like', $self->_features_gff->{$feature_uniquename}->{'seq_id'} },
			  },
			  {columns => [qw/feature_id/]},
			)->single()->feature_id();
			
			my $locgroup = $self-> schema -> resultset('Sequence::Featureloc')->search(
			{
				'feature.organism_id'=> $self->organism_id ,
				'feature.uniquename' => $feature_uniquename,
				'fmin' => $self->_features_gff->{$feature_uniquename}->{'start'} - 1,
				'fmax' => $self->_features_gff->{$feature_uniquename}->{'end'},
				'srcfeature_id' => $srcfeature_id,
			},
			  #{ join => [ 'feature' ] , prefetch=> [ 'feature']},
			  { join => [ 'feature' ]},
			)->single()->locgroup();
			
			my $decrement_locgroup_sql = sub {
				
				my $locgroup_condition = "> $locgroup";
				my $fl_rs = $self -> schema -> resultset('Sequence::Featureloc')->search(
				  {
				  'feature_id' => $feature_id,
				  'locgroup' => \$locgroup_condition,
				  },
				  { 'order_by' => { -asc => [qw/locgroup/]}}
				);
				
				while (my $fl_row=$fl_rs->next()){
					$fl_row->set_column('locgroup' => ($fl_row->get_column('locgroup') - 1));
					$fl_row->update();
				}
			};

			my $delete_sql = sub{
				my $retval = $self -> schema -> resultset('Sequence::Featureloc')->search(
				{
				  'feature_id' => $feature_id,
				  'fmin' => $self->_features_gff->{$feature_uniquename}->{'start'} - 1,
				  'fmax' => $self->_features_gff->{$feature_uniquename}->{'end'},
				  'srcfeature_id' => $srcfeature_id,
				  'locgroup' => $locgroup,
				  },
				)->delete_all();
				
				if($retval){print STDERR "DELETED featureloc record for $feature_uniquename\n";}
				
				$self->schema->txn_do($decrement_locgroup_sql);    # nested transaction for decrement				
				
			};
	
			
			#delete rows from featureloc and decrement locgroups using transactions
			try {
				$self->schema->txn_do($delete_sql);
			  }
			  catch {# Transaction failed
				die "Could not delete rows and/or decrement locgroups. Error: $!"
				  if ( $_ =~ /Rollback failed/ );
		
				print STDERR "Error: " . $_;
			  };
			  $counter++;
	}
	return $counter;
}

###
1;                                             #do not remove
###

=pod

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Surya Saha				<suryasaha at cornell dot edu , @SahaSurya>   
    Lukas Mueller			<lam87 at cornell dot edu>
    Jonathan "Duke" Leto	<jonathan at leto dot net>
    
=cut

