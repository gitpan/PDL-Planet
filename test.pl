# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

use PDL;
use PDL::Orbit;
use Test;

eval "use PDL::Planet;";
unless ($@){
  plan tests => 13;
  print "ok 1\n";
} else {
  print "not ok 1\n";
}

sub testok ($$) {
  my $bool = shift;
  my $n    = shift;
  if ($bool) {
    print "ok $n\n";
  } else {    print "not ok $n\n";
  }
}


######################### End of black magic.

BEGIN { unlink glob ("testimg*.*") };

my $pl = PDL::Planet->new
                    ->read("./images/earth.jpg")
                    ->transform(TYPE => Orthographic, SIZE => [500,500], CENTER => [40,-105])
                    ->write("testimg2.jpg");
testok -s "testimg2.jpg" > 0, 2;

my $starttime = 728784000;
my $steptimes = (sequence(50)*180) + $starttime;
my ($lat, $lon, $height) = PDL::Orbit::sattrack ($steptimes, './champ_tle.txt');

my ($x, $y)   = $pl->llh2xy($lat, $lon, $height);

my $xok = pdl(qw[463 443 416 380 339 294 246 198 151 108 70 39 16 3 0 7 24 51 86 
                 -1 -1 -1 -1 -1 -1 -1 -1 -1 503 511 509 495 471 438 397 350 299 
                 246 193 142 95 55 22 0 -13 -14 -5 13 -1 -1]);
my $yok = pdl(qw[400 430 452 467 473 471 461 442 417 385 348 308 266 223 181 143 
                 108 79 56 -1 -1 -1 -1 -1 -1 -1 -1 -1 222 264 304 343 378 407 431 
                 449 458 460 455 441 421 394 362 326 286 245 205 166 -1 -1]);

testok (sum($x - $xok) == 0, 3);
testok (sum($y - $yok) == 0, 4);

my $glyph     = PDL::Planet->new->read("./images/satellite.gif")->write("testimg5.png");
testok -s "testimg5.png" > 0, 5;

$pl->paste($glyph, $x, $y)->rgb;

my $n = 6;
foreach my $type (qw(gif jpg bmp png pbm pgm ppm tif)) {
  $pl->write("./testimg$n.$type");
  testok -s "testimg$n.$type", $n++;
}











