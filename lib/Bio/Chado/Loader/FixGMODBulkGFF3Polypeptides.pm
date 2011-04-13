package Bio::Chado::Loader::FixGMODBulkGFF3Polypeptides;
use Moose;
use namespace::autoclean;

with 'MooseX::Runnable';
with 'MooseX::Getopt';

has schema => (
    is         => 'rw',
    isa        => 'Bio::Chado::Schema',
    traits     => ['NoGetopt'],
    lazy_build => 1,
   );

has 'mrna_name_like' => (
    documentation => <<'',
LIKE clause for mRNA terms to process

    is => 'ro',
    isa => 'Str',
   );

with 'Bio::Chado::Loader';

sub run {
    my ( $self ) = @_;

    $self->schema->txn_do( sub {
        my $mrnas =
            $self->cvterm(qw( sequence mRNA ))
                ->search_related( 'features', {
                    $self->mrna_name_like ? ( 'me.name' => { like => $self->mrna_name_like } ) : (),
                  });

        while( my $m = $mrnas->next ) {
            my @polys =
                $m->search_related( 'feature_relationship_objects', {
                    'me.type_id' => $self->cvterm('relationship','derives_from')->cvterm_id,
                   })
                  ->search_related( 'subject', {
                      'subject.type_id' => $self->cvterm(qw( sequence polypeptide ))->cvterm_id,
                      'subject.name' => { like => 'polypeptide-auto%' },
                   });

            my $fix = shift @polys or die "no polys for ".$m->name;
            my ($fixloc) = my @fl = $fix->featureloc_features;
            die "invalid featurelocs for ".$fix->name unless @fl == 1;
            for my $p (@polys) {
                for ($p->featureloc_features) {
                    $fixloc->fmin == $_->fmin or die "mins not equal for ".$_->feature->name;
                    $fixloc->fmax == $_->fmax or die "maxs not equal for ".$_->feature->name;
                    $fixloc->srcfeature_id == $_->srcfeature_id or die "maxs not equal for ".$_->feature->name;
                }
                $p->delete;
            }

            $fix->name( $m->name );
            $fix->update;
        }
    });

    return 1;
}

my %cvt_cache;
sub cvterm {
    my ( $self, $cv, $cvterm ) = @_;

    $cvt_cache{$cv}{$cvterm} ||=
        $self->schema->resultset('Cv::Cv')
            ->search({ 'me.name' => $cv })
            ->search_related('cvterms', {
                'cvterms.name' => $cvterm,
                is_obsolete => 0,
              })
            ->single;

}

1;
