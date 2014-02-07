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

CAVEAT
Will break if no ID field in attributes in GFF record

=over

=cut

require Bio::GFF3::LowLevel::Parser;


use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
use Proc::ProcessTable;
use Bio::GFF3::LowLevel::Parser;
use Bio::GFF3::LowLevel qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);
use Try::Tiny;

has 'schema' => (
    is  => 'rw',
    isa => 'Bio::Chado::Schema',
    lazy => 1,
    builder => '_build_schema'
);

has 'file_name' => (
	isa       => 'Str',
	is        => 'rw',
	predicate => 'has_file_name',
	clearer   => 'clear_file_name'
);

has 'organism_name' => (
    documentation => 'exact species of organism, as stored in database, e.g. "Solanum lycopersicum"',
    is       => 'rw',
    isa      => 'Str',
);

has 'organism_id' => (
    documentation => 'id organism, as stored in database, e.g. 1 for Solanum lycopersicum',
    is       => 'rw',
    isa      => 'Str',
);

has 'features_gff' => (
	documentation => 'hash of feature uniquename => hashref of entire GFF record. It will break if no ID field in GFF record. uniquenames are used for key since GFF3 mandates they be unique in a GFF3 file.',
    is      => 'ro',
    isa     => 'HashRef',
    traits	=> ['Hash'],
    default => sub { { } },
    handles => {
        add_features_gff    => 'set',
        count_features_gff => 'count',
        features_gff_exists => 'exists',
    },
);

=item C<dump_features_gff ()>

Prints the contents of the features_gff to STDERR. For internal debugging use only. 

=cut

sub dump_features_gff {
    my ($self) = @_;

	my ($feature_uniquename,$fields);

	while (($feature_uniquename,$fields) = each %{$self->features_gff}){
		print STDERR $feature_uniquename." => ".gff3_format_feature($fields);
	}
}

has 'feature_ids_uniquenames_gff' => (
	documentation => 'hash of featureID => feature_uniquename for all features in GFF that will be updated. Populated in prepare_bulk_upload',
    is      => 'ro',
    isa     => 'HashRef',
    traits	=> ['Hash'],
    default => sub { { } },
    handles => {
        add_feature_ids_uniquenames_gff    => 'set',
        count_feature_ids_uniquenames_gff => 'count',
    },
);

has 'cvterms_gff' => (
	documentation => 'hash of cvterm => 1 from GFF. TOFIX Use counts instead of constant 1',
    is      => 'ro',
    isa     => 'HashRef[Str]',
    traits  => [ 'Hash' ],
    default => sub { { } },
    handles => {
        add_cvterms_gff    => 'set',
        count_cvterms_gff => 'count',
    },
);

has 'is_analysis' => (
    documentation => 'set true if this feature should be recorded as from an analysis',
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
);

has 'cache' => (
	documentation => 'hash of cvterm.name => hash of feature.uniquename => feature.feature_id,featureloc.srcfeature_id. Not recording featureloc.locgroup,featureloc.rank',
    is      => 'ro',
    isa     => 'HashRef',
    default => sub { { } },
);

has 'feature_uniquename_cache' => (
	documentation => 'hash of feature uniquenames => 1 for all features in cache',
    is      => 'ro',
    isa     => 'HashRef',
	traits	=> ['Hash'],
    default => sub { { } },
	handles => {
        add_feature_uniquename_cache    => 'set',
        count_feature_uniquename_cache => 'count',
        feature_uniquename_cache_exists => 'exists'
    },
);

with 'Bio::Chado::Loader';

use autodie qw(:all);
use Data::Dumper;
use 5.010;

=item C<run ()>

TODO For mappable cmd line exec// See Fasta.pm

=cut

sub run {
    my ($self, %args) = @_;

    exit 0;
}

=item C<parse ()>

Parse a GFF3 file using Bio::GFF3::LowLevel::Parser.

=cut

sub parse {
    my ($self, %args) = @_;
    
    my $gff = Bio::GFF3::LowLevel::Parser->open($self->file_name);
    while ( my $item = $gff->next_item ) {
		if ( ref $item eq 'ARRAY' ) {
			## $i is an arrayref of feature lines that have the same ID,
			## in the same format as returned by Bio::GFF3::LowLevel::gff3_parse_feature
			#warn Dumper [ $i ];
			for my $feature (@$item) {
				# for each location of this feature
				$self->parse_feature( $feature);
			}
		}
		elsif ( $item->{directive} ) {
			if ( $item->{directive} eq 'FASTA' ) {
				my $fasta_filehandle = $item->{filehandle};
				## parse the FASTA in the filehandle with BioPerl or
				## however you want.  or ignore it.
				print "got fasta handle\n";
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
    my $self = shift;
    my $feature_hash = shift;
    
    #warn Dumper [$self];
    #warn Dumper [ $feature_hash];
    
    # add data
    if ( $self->features_gff_exists($feature_hash->{'attributes'}->{'ID'}->[0]) ){
    	print STDERR "Duplicate found for ".$feature_hash->{'attributes'}->{'ID'}->[0]."\n";
    	#warn Dumper [ $self->features_gff ];
    	die "Multiple features with same ID as ".$feature_hash->{'attributes'}->{'ID'}->[0]." found. Exiting..";
    }
    elsif (!defined $feature_hash->{'attributes'}->{'ID'}->[0]){#exit if no ID field
    	die "No ID defined for feature in GFF. Exiting..";
    }
    elsif( $feature_hash->{'type'} eq 'mRNA' ){#call mRNA handler for parsing exons, CDSs and UTRs
    	$self->parse_mRNA_feature($feature_hash);
    	return;
    }
    else{
    	print STDERR "Adding feature ".$feature_hash->{'attributes'}->{'ID'}->[0]."\n";
	    $self->add_features_gff( $feature_hash->{'attributes'}->{'ID'}->[0] => $feature_hash);
	    $self->add_cvterms_gff( $feature_hash->{'type'} => 1 );
	    #warn Dumper [ $self->features_gff ] if $feature_hash->{'attributes'}->{'ID'}->[0]=~ /exon/;
    }
    
    #recursively calling self for nested child features
    if ( $feature_hash->{'child_features'} ) {
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			$self->parse_feature($feature_child->[0]);
		}
	}
	
	#recursively calling self for nested derived features
    if ( $feature_hash->{'derived_features'} ) {
		for my $feature_derived ( @{ $feature_hash->{'derived_features'} } ) {
			$self->parse_feature($feature_derived->[0]);
		}
	}
}

=item C<parse_mRNA_feature ()>

Saves mRNA, exons and polypeptide if exons, CDSs and UTRs are present. Saves mRNA and creates exon & polypeptide features with CDS and UTR coordinates if both are present. Dies if only CDS or UTR is present.  

=cut

sub parse_mRNA_feature {
    my $self = shift;
    my $feature_hash = shift;
    
    my %children = (
    	'exon' => 0,
    	'intron' => 0,
    	'CDS' => 0,
    	'five_prime_UTR' => 0,
    	'three_prime_UTR' => 0,
    );
     
    #test children presence/absence
    if ( $feature_hash->{'child_features'} ) {
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			die "mRNA has child feature other than five_prime_UTR, three_prime_UTR, CDS, exon .. cannot handle "
				if (!exists $children{$feature_child->[0]->{'type'}});
			$children{$feature_child->[0]->{'type'}}++;
			#$self->parse_feature($feature_child->[0]);
		}
	}
	else{
		print STDERR "No child features saved for ".$feature_hash->{'attributes'}->{'ID'}->[0]."\n";
		return;
	}
	
	#polypeptide
	my (@starts,@ends,$gene);
	$gene=$feature_hash->{'attributes'}->{'ID'}->[0];
	$gene=~ s/^mRNA\://;
	
	my $polypeptide_feature = {
		seq_id => $feature_hash->{'seq_id'},
		source => $feature_hash->{'source'},
		type   => 'polypeptide',
		#start  => '23486',
		#end    => '48209',
		score  => undef,
		strand => $feature_hash->{'strand'},
		phase  => undef,
		attributes => {
			ID => [
				'polypeptide:'.$gene
			],
			Parent => [
				$feature_hash->{'attributes'}->{'Parent'}->[0]
			],
		},
		child_features => [],
		derived_features => [],
	};

	if( $children{'CDS'}){
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if($feature_child->[0]->{'type'} eq 'CDS'){
				push @starts, $feature_child->[0]->{'start'} - 1 ;
				push @ends, $feature_child->[0]->{'end'};
			}
		}
	}
	elsif( $children{'five_prime_UTR'} && $children{'three_prime_UTR'}){
		#Do we check if >1 3' or 5' UTRs???
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if($feature_child->[0]->{'type'} eq 'five_prime_UTR'){
				push @starts, $feature_child->[0]->{'start'} -1;
			}
			elsif($feature_child->[0]->{'type'} eq 'three_prime_UTR'){
				push @ends, $feature_child->[0]->{'end'};
			}
		}
	}
	$polypeptide_feature->{'start'} = (sort { $a <=> $b } @starts)[0]; #asc
	$polypeptide_feature->{'end'} = (sort { $b <=> $a } @ends)[0]; #desc
	#add polypeptide
	$self->parse_feature($polypeptide_feature);
	
	#exons
	if( $children{'exon'}){
		#save exons
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if($feature_child->[0]->{'type'} eq 'exon'){
				print STDERR "Adding exon auto ".$feature_child->[0]->{'attributes'}->{'ID'}->[0]."\n";
				$self->parse_feature($feature_child->[0]);
			}
		}
	}
	elsif( !$children{'exon'} ){#this will not happen for tomato GFF3s
		print STDERR "Adding exons manually\n";
		#check if only 1 3' and only 1 5' UTR
		die "More than one or no five_prime_UTR or three_prime_UTR for ".$feature_hash->{'attributes'}->{'ID'}->[0]
			if (($children{'five_prime_UTR'} != 1) || ($children{'three_prime_UTR'} != 1) );
		
		#get coords
		my ( %cds_starts, %cds_ends, $five_prime_UTR_start, $three_prime_UTR_end);
		for my $feature_child ( @{ $feature_hash->{'child_features'} } ) {
			if($feature_child->[0]->{'type'} eq 'CDS'){
				print STDERR "Adding hash data for ".$feature_child->[0]->{'attributes'}->{'ID'}->[0]."\n";
				$cds_starts{$feature_child->[0]->{'attributes'}->{'ID'}->[0]} = $feature_child->[0]->{'start'} - 1; 
				$cds_ends{$feature_child->[0]->{'attributes'}->{'ID'}->[0]} = $feature_child->[0]->{'end'};
			}
			elsif($feature_child->[0]->{'type'} eq 'five_prime_UTR'){
				$five_prime_UTR_start = $feature_child->[0]->{'start'} - 1; 
			}
			elsif($feature_child->[0]->{'type'} eq 'three_prime_UTR'){
				$three_prime_UTR_end = $feature_child->[0]->{'end'}; 
			}
		}		
		
		#check if 5' before 1st CDS and 3' after last CDS
		die "five_prime_UTR starts at or after first CDS for ".$feature_hash->{'attributes'}->{'ID'}->[0]
			if ( $five_prime_UTR_start >= ((sort { $a <=> $b } values %cds_starts)[0])); #asc
		die "three_prime_UTR ends at or before last CDS for ".$feature_hash->{'attributes'}->{'ID'}->[0]
			if ( $three_prime_UTR_end <= ((sort { $b <=> $a } values %cds_ends)[0])); #desc
		
	
		#create exons
		#my ($CDS_ID,$value);
		
		print STDERR "cds_starts size ".(scalar keys %cds_starts)."\n";
		while ( my ($CDS_ID,$value) = each (%cds_starts)) {
			print STDERR "key: $CDS_ID\t val: $value\n";	
		}
		
		my @cds_starts_values = values %cds_starts;
		while ( my ($CDS_ID,$value) = each (%cds_starts)) {
			print STDERR "\tkey: $CDS_ID\t val: $value\n";
			my $CDS_ID_trimmed = $CDS_ID;
			$CDS_ID_trimmed =~ s/^CDS\://;
			my $exon_feature = {
				seq_id => $feature_hash->{'seq_id'},
				source => $feature_hash->{'source'},
				type   => 'exon',
				#start  => '23486',
				#end    => '48209',
				score  => undef,
				strand => $feature_hash->{'strand'},
				phase  => undef,
				attributes => {
					ID => [
						'exon:'.$CDS_ID_trimmed
					],
					Parent => [
						$feature_hash->{'attributes'}->{'ID'}->[0]
					],
				},
				child_features => [],
				derived_features => [],
			};			
			print STDERR "Created exon ".$exon_feature->{attributes}->{ID}->[0]."\n";

			#my @cds_starts_values = values %cds_starts;
			#if first CDS
			#if( $value == ((sort { $a <=> $b } (values %cds_starts))[0])){
			if( $value == ((sort { $a <=> $b } @cds_starts_values)[0])){
				$exon_feature->{'start'} = $five_prime_UTR_start;
				$exon_feature->{'end'} = $cds_ends{$CDS_ID};
			}	     
			#elsif( $value == ((sort { $b <=> $a } (values %cds_starts))[0])){#if last CDS
			elsif( $value == ((sort { $b <=> $a } (@cds_starts_values))[0])){#if last CDS
				$exon_feature->{'start'} = $value;
				$exon_feature->{'end'} = $three_prime_UTR_end;
			}
			else{#CDS in the middle
				$exon_feature->{'start'} = $value;
				$exon_feature->{'end'} = $cds_ends{$CDS_ID};
			}
			#save exon
			print STDERR "Saving exon ".$exon_feature->{attributes}->{ID}->[0]."\n";
			#warn Dumper [ $exon_feature];
			$self->parse_feature($exon_feature);
			print STDERR "Saved exon ".$exon_feature->{attributes}->{ID}->[0]."\n";
			#$CDS_ID=$CDS_ID_trimmed=undef;
		}
	}
	
	#remove mRNA sub features info to save memory
	$feature_hash->{'child_features'} = [];
	$feature_hash->{'derived_features'} = [];
	#save mRNA
	$self->add_features_gff( $feature_hash->{'attributes'}->{'ID'}->[0] => $feature_hash);
        $self->add_cvterms_gff( $feature_hash->{'type'} => 1 );
	
}

=item C<organism_exists ()>

Check if organism for GFF3 exists in the Chado database. Return organism::organism.organism_id 
if present otherwise 0. No organism should have id of 0.

=cut

sub organism_exists {
    my ($self) = @_;
    
    my $org_id = $self->schema->resultset('Organism::Organism')-> search (
    	{ species => $self->organism_name})->get_column('organism_id')->first();
    return $org_id ? $org_id : 0;
}

=item C<populate_cache ()>

Populate cache hash with feature.uniquename=feature.feature_id,featureloc.srcfeature_id,
featureloc.locgroup,featureloc.rank,cvterm.name from cvterm,feature,featureloc relations. 
Return number of records added to cache. Call function before parse() in case parents are
not in GFF but in DB already.  

CAVEAT: Count may be over estimation if some features have multiple locgroups(e.g. contigs).

=cut

sub populate_cache {
    my ($self) = @_;
    
    #setup DB dsn's
    my $fl_rs = $self->schema -> resultset('Sequence::Featureloc')->search(
#    	{ 'organism_id'=> $self->organism_id },
		{ 'organism_id'=> 1 , 'feature.uniquename' => { 'like', '%Solyc01g1123%'}},#for testing, only 69 floc records
#		{ 'organism_id'=> 1 , 'feature.uniquename' => { 'like', '%dummy%'}},#for testing, only 2 floc records
    	{ join => [ 'feature' ] , prefetch=> [ 'feature']}
    	);

	my $ft_rs = $self->schema -> resultset('Sequence::Feature')->search(
		{ 'organism_id'=> $self->organism_id },
		{ join => [ 'type' ] , prefetch=> [ 'type']}
		);
    print STDERR "fl_rs count: ".$fl_rs->count()."\n";
    print STDERR "ft_rs count: ".$ft_rs->count()."\n";
    
    #create cache hash
    my $count = 0;
    while (my $fl_row = $fl_rs->next()){
		if($ft_rs->search(
			{'type_id' => $fl_row->feature->type_id(), 'uniquename' => $fl_row->feature->uniquename}
			)->count() == 1){
			my $ft_row = $ft_rs->search(
					{'type_id' => $fl_row->feature->type_id(), 'uniquename' => $fl_row->feature->uniquename}
					)->first();

			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{feature_id}=$fl_row->feature->feature_id();
#			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{locgroup}=$fl_row->locgroup();
#			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{rank}=$fl_row->rank();
			$self->cache->{$ft_row->type->name()}->{$fl_row->feature->uniquename()}->{srcfeature_id}=$fl_row->srcfeature_id();
			
			$self->add_feature_uniquename_cache( $fl_row->feature->uniquename() => 1 );
			
			$count++;	
			if ($count % 100000 == 0) {
				print STDERR "processing $count\r";
				#warn Dumper [ $self-> cache -> {$ft_row->type->name()}->{$fl_row->feature->uniquename()}];
				#warn Dumper [ $self-> cache -> {$ft_row->type->name()}];
			    my ( $i, $t );
				$t = new Proc::ProcessTable;
				foreach my $got ( @{ $t->table } ) {
					next if not $got->pid eq $$;
					$i = $got->size;
				}
				print STDERR "Memory used: ".($i / 1024 / 1024)."MB \n";
			}
		}
		else{
			die "Multiple features found for organism with same type_name and uniquename which violates Chado constraints. Exiting..";
		}
	}
	return ($count);
}

=item C<prepare_bulk_upload ()>

Compare %features from GFF to %cache. Prepare data structures to write to DB. GFF %features 
not found in %cache are written to gff3.exceptions. New featureloc records have locgroup=0 and 
rank=0 (primary location). Old locgroups for features being inserted are incremented by 1.

CAVEAT: Only new featureloc's are handled right now.

=cut

sub prepare_bulk_upload {
    my ($self) = @_;

    #compare %features to %cache
    my ($feature_uniquename,$fields,$feature_gff,$disk_exception_str);
    my %counters=('exceptions' => 0, 'inserts' => 0);
    $disk_exception_str='';
    
    #warn Dumper [ $self->features ];
    while (($feature_uniquename,$fields) = each %{$self->features_gff}){
    	
    	if ($self->cache->{$fields->{'type'}}->{$feature_uniquename}){
    		#create intermediate file str
    		#print STDERR "Data string for $feature_uniquename\n"; 
    		
    		#record feature_ids in class var
    		$self->add_feature_ids_uniquenames_gff( $self->cache->{$fields->{'type'}}->{$feature_uniquename}->{feature_id} 
    			=> $fields->{'attributes'}->{'ID'}->[0]);

    		$counters{'inserts'}++;
    		if ($counters{'inserts'}%10000 == 0){
    			print STDERR "\r".$counters{'inserts'}." processed..\n";
    		}
    	}
    	else{
    		#create exception file str
    		#print STDERR "Exception GFF record for $feature_uniquename\n";
			#warn Dumper [ $fields ];
    		
    		#$disk_exception_str.=join("	", @$fields)."\n";
    		$disk_exception_str.=gff3_format_feature($fields);
    		$counters{'exceptions'}++;
    	}
    }
    die "No features to update." if $counters{'inserts'} == 0;
    
    if($counters{'exceptions'} > 0 ){
    	my ($exception_gff_fh);
    	open($exception_gff_fh,">",$self->file_name.'exceptions')
    		or die("Could not create exception gff file: $!");
    	print $exception_gff_fh $disk_exception_str;
    	close($exception_gff_fh);
    }
    
    #warn Dumper [ $self->feature_ids_uniquenames_gff];
    
    print STDERR $counters{'inserts'}." records prepared for insert\n";
    print STDERR $counters{'exceptions'}." exception GFF records written to ".$self->file_name."\.exceptions \n";
    
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
Gets srcfeature_id from cache which is incorrect if parent was changed
Presuming that new location is the primary location of feature (locgroup=0, rank=0)

=cut

sub bulk_featureloc_upload {
    my ($self) = @_;
    
    #$self->schema->storage->debug(1);##SQL statement
    #warn Dumper[keys $self->feature_ids_uniquenames_gff];
    my @feature_ids_to_update = keys ($self->feature_ids_uniquenames_gff);
    my $fl_rs = $self->schema -> resultset('Sequence::Featureloc')->search(
		{ 'feature_id'=>  \@feature_ids_to_update},
		{ 'order_by' => { -desc => [qw/feature_id locgroup/]}}
		);
	
	#create new rows
	my $create_sql = sub {
		my $counter=0;
		while ( my ($feature_id,$feature_uniquename) = each $self->feature_ids_uniquenames_gff){
			my ($strand,$phase);
			
			if( $self->features_gff->{$feature_uniquename}->{'strand'} eq '+'){ $strand =1;}
    		elsif( $self->features_gff->{$feature_uniquename}->{'strand'} eq '-'){ $strand =-1;}
    		else { $strand=0;}
    		
    		if ($self->features_gff->{$feature_uniquename}->{'phase'}){
				$fl_rs->create({
						'feature_id' => $feature_id,
						#getting srcfeature_id from cache which is incorrect if parent was changed
						'srcfeature_id' => $self->cache->{$self->features_gff->{$feature_uniquename}->{'type'}}
							->{$feature_uniquename}->{srcfeature_id},
						'fmin' => $self->features_gff->{$feature_uniquename}->{'start'} -1,
						'fmax' => $self->features_gff->{$feature_uniquename}->{'end'},
						'strand' => $strand,
						'phase' => $self->features_gff->{$feature_uniquename}->{'phase'},
						'locgroup' => 0,
						'rank' => 0,
					});	
    		}
    		else{
				$fl_rs->create({
						'feature_id' => $feature_id,
						#getting srcfeature_id from cache which is incorrect if parent was changed
						'srcfeature_id' => $self->cache->{$self->features_gff->{$feature_uniquename}->{'type'}}
							->{$feature_uniquename}->{srcfeature_id},
						'fmin' => $self->features_gff->{$feature_uniquename}->{'start'} -1,
						'fmax' => $self->features_gff->{$feature_uniquename}->{'end'},
						'strand' => $strand,
						'locgroup' => 0,
						'rank' => 0,
					});	
    		}
		$counter++;
		}
		print STDERR "$counter rows added to featureloc\n";
	};
	
	#update locgroup=locgroup+1
	my $increment_locgroup_sql = sub {
		my $counter=0;
		while (my $fl_row = $fl_rs->next()){
			$fl_row->set_column('locgroup' => ($fl_row->get_column('locgroup') + 1));
			$fl_row->update();
			$counter++;
		}
		print STDERR "$counter featureloc rows locgroup fields updated\n";
		$self->schema->txn_do($create_sql); # nested transaction for create
	};
    
    #update locgroups and insert rows into featureloc using transactions
    try {
		$self->schema->txn_do($increment_locgroup_sql);
	} catch {# Transaction failed
		die "Could not increment locgroups and/or create new rows. Error: $!" if ($_ =~ /Rollback failed/);
		#deal_with_failed_transaction();
		print STDERR "Error: ".$_;
	};
    
}

###
1;   #do not remove
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
