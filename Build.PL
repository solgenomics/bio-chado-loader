use Module::Build;
my $build = Module::Build->new(
    module_name          => 'Bio::Chado::Loader',
    dist_name            => 'Bio-Chado-Loader',
    dist_abstract        => 'Load data into the Chado database schema',
    license              => 'perl',
    recursive_test_files => 1,
    create_readme        => 1,
    add_to_cleanup       => [ 'Bio-Chado-Loader-*', 't/var' ],
    create_license       => 1,
    meta_merge           => {
        resources =>
          {
            repository => 'http://github.com/leto/Bio-Chado-Loader',
          },
    },
    build_requires => {
        'Test::Class::Most'   => 0,
    },
    requires       => {
        'autodie'             => 0,
        'IPC::System::Simple' => 0, # for autodie qw/:all/
	'Memoize'             => 0,
        'Moose'               => '1.0',
        'MooseX::Runnable'    => 0,
        'MooseX::Getopt'      => 0,
        'perl'                => '5.10.0',
        'Carp'                => 0,
        'Bio::Chado::Schema'  => '0.06001',
        'Try::Tiny'           => 0,
    },
);
$build->create_build_script;
