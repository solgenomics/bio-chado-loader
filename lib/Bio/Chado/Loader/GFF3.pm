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


use Moose;
with 'MooseX::Runnable';
with 'MooseX::Getopt';
use Proc::ProcessTable;
require Bio::GFF3::LowLevel;
use Bio::GFF3::LowLevel  qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);
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

has 'pragma_lines' => (
	isa       => 'ArrayRef[Str]',
	is        => 'rw',
	predicate => 'has_pragma_lines',
	clearer   => 'clear_pragma_lines'
);

has 'features_gff' => (
	documentation => 'hash of feature uniquename => hashref of entire GFF record via . It will break if no ID field in GFF record. uniquenames are used for key since GFF3 mandates they be unique in a GFF3 file.',
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

TODO // See Fasta.pm

=cut

sub run {
    my ($self, %args) = @_;

    exit 0;
}

=item C<parse ()>

Parse a GFF3 file. Calls parse_line() for each line.

=cut

sub parse {
    my ($self, %args) = @_;
    open my $fh, "<", $self->file_name;
    while( my $line = <$fh> ) {
        next if $line =~ m/^\s*$/;
        chomp $line;
        next if $self->is_pragma_line($line);
        $self->parse_line($line);
    }
}

=item C<parse_line ()>

Parse a GFF3 line to get feature details. 

CAVEAT: This will break if no ID field is present

=cut

sub parse_line {
    my ($self, $line) = @_;
    my $gff_line_hash = gff3_parse_feature($line);
    
    #warn Dumper [$gff_line_hash];
    
    #validate parents
    if (exists $gff_line_hash->{'attributes'}->{'Parent'}){
		$self->_validate_parents($gff_line_hash->{'attributes'}->{'Parent'}->[0]);    	
    }
    
    # add data, Will break if no ID field
    if ( $self->features_gff_exists($gff_line_hash->{'attributes'}->{'ID'}->[0]) ){
    	die "Multiple features with same ID found. Exiting..";
    }
    elsif (!defined $gff_line_hash->{'attributes'}->{'ID'}->[0]){
    	die "No ID defined for feature in GFF. Exiting..";
    }
    else{
	    $self->add_features_gff( $gff_line_hash->{'attributes'}->{'ID'}->[0] => $gff_line_hash);
	    $self->add_cvterms_gff( $gff_line_hash->{'type'} => 1 );	
    }
}

=item C<_validate_parents ()>

Check if parent feature exists

=cut


sub _validate_parents {
	my ($self, $parents) = @_;
    my (@parents) = split (/,/, $parents);
    for my $p (@parents) {
        unless ( $self->features_gff_exists( $p ) || $self->feature_uniquename_cache_exists( $p )) {
            die "$p is an unknown Parent!";
        }
    }
}

=item C<is_pragma_line ()>

Check if line is a pragma

=cut

sub is_pragma_line {
    my ($self, $line) = @_;
    return $line =~ m/^##/ ? 1 : 0;
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

CAVEAT: GFF values currently not handled are sequence. Only new featureloc's are handled right now.

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
#	   "organism", #dgg
#	   "analysis", #dgg
#	   "db", ## dgg
#	   "dbxref",
#	   "cv", ## dgg
#	   "cvterm", ## dgg
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
		die "Could not increment locgroups and/or create new rows. Error: $!" if ($_ =~ /Rollback failed/);# Rollback failed
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

    Jonathan "Duke" Leto	<jonathan at leto.net>
    Surya Saha				<suryasaha at cornell.edu , @SahaSurya>   

=cut