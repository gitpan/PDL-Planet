# Makefile.PL for a package defined by PP code.

use PDL::Core::Dev;            # Pick up development utilities
use ExtUtils::MakeMaker;
use Config;
use vars qw/$package/;

# -- Add new subroutines here! --

my @src = qw(util.c spa.c);
my @obj = @src;
map {s/\.[fc]/\.o/;} @obj; # swap .f, .c for .o

#
## Search for plplot library and include file
#

my $libdir = "./libimage";
my $incdir = "./libimage";

my %pmfiles = map { $_ => "\$(INST_LIBDIR)/Planet/$_" } glob ("images/*");
$pmfiles{'Planet.pm'} = '$(INST_LIBDIR)/Planet.pm';
#print "pmfiles = \n";
#for (keys %pmfiles) { print "$_ => $pmfiles{$_}\n"; }

$package = ["planet.pd",Planet,PDL::Planet];
WriteMakefile(
	      'NAME'  	     => 'PDL::Planet',
	      'AUTHOR'       => 'Douglas Hunt (dhunt@ucar.edu)',
	      'PREREQ_PM'    => {PDL         => 2.4.1, 
				 Image::Size => 0, 
				 PDL::Orbit  => 0},
	      'VERSION_FROM' => 'planet.pd',
	      'TYPEMAPS'     => [&PDL_TYPEMAP()],
	      'OBJECT'       => 'Planet.o ' . join (" ", @obj) . ' libimage/libimage.a',
	      'PM'           => \%pmfiles,
              'INC'          => &PDL_INCLUDE()." $plplot_include_path -I./libimage",
              'LIBS'         => ['-lgif -ljpeg -lpng -lm -lz -lpnm -ltiff -lm -lnetpbm'],
	      'clean'        => {'FILES'  => 
				   'Planet.pm Planet.xs Planet.o Planet.c Planet.bs libimage/*.o libimage/libimage.a libimage/config.* libimage/ephem libimage/Makefile'},
	      'dist'         => { COMPRESS => 'gzip', SUFFIX => 'gz' },	
	      );

sub MY::postamble { pdlpp_postamble($package); }	

# Add genpp rule; this will invoke PDL::PP on our PP file
# the argument is an array reference where the array has three string elements:
#   arg1: name of the source file that contains the PP code
#   arg2: basename of the xs and pm files to be generated  
#   arg3: name of the package that is to be generated
sub MY::postamble { 
  my $orig = pdlpp_postamble($package);

  my $cflags = '';
  # for AMD 64, use -fPIC
  if ($Config{'archname'} =~ /x86_64-linux/) {
    $cflags = '-fPIC';
  } 

  my $added = '';
  foreach my $s (@src) {
    my $o = $s;
    $o =~ s/\.[cf]/\.o/;
    
    if ($s =~ /\.c$/) { # C code
      $added .= "$o: $s\n\t$Config{'cc'} -c $cflags -o $o -I./libimage -O3 $s\n\n";
    }
    elsif ($s =~ /\.f$/) { # FORTRAN
      $added .= "$o: $s\n\tg77 -c -O3 -o $o $s\n\n";
    }
    
  }

  # compile libimage stuff
  $added .= "libimage/libimage.a: \n\tcd libimage; ./configure; make; cd ..\n\n";

  return $orig . $added;

}  
