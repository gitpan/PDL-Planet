use Test::More tests => 12;

use PDL;
use PDL::Orbit;
use PDL::Planet;

BEGIN { unlink glob ("testimg*.*") };

my $pl = PDL::Planet->new
                    ->read("./images/earth.jpg")
                    ->transform(TYPE => Orthographic, SIZE => [500,500], CENTER => [40,-105])
                    ->write("testimg2.jpg");
ok -s "testimg2.jpg", 'new/read/transform/write';

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

ok ( (sum($x - $xok) == 0) && (sum($y - $yok) == 0), "lat/lon/height -> X/Y");

my $glyph     = PDL::Planet->new->read("./images/satellite.gif")->write("testimg5.png");
ok -s "testimg5.png", "Read GIF glyph, write PNG glyph";

my $png       = PDL::Planet->new->read("./images/satellite.gif")->write_png;
my $file_png  = do { local( @ARGV, $/ ) = "testimg5.png"; <> } ; # slurp!
ok ( ($png eq $file_png), "Read GIF glyph, write PNG in-memory glyph");

$pl->paste($glyph, $x, $y)->rgb;

my $n = 6;
foreach my $type (qw(gif jpg bmp png tif)) {
  $pl->write("./testimg$n.$type");
  ok -s "testimg$n.$type", "Writing $type";
  $n++;
}

SKIP: {
  my $checkfuncs = do { local( @ARGV, $/ ) = "./libimage/checkfuncs.h"; <> } ; # slurp!
  skip "PNM support not enabled", 3 if ($checkfuncs !~ /\#define HAVE_LIBPNM/);

  foreach my $type (qw(pbm pgm ppm)) {
    $pl->write("./testimg$n.$type");
    ok -s "testimg$n.$type", "Writing $type";
    $n++;
  }
} 








