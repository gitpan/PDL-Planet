use Test::More tests => 4;

use PDL;
use PDL::Planet;

BEGIN { unlink glob ("testimg*.*") };

# This should show an image (testimg2.jpg) that is dark and the left hand
# side of the globe, depicting the sun shining over lat=0, lon=0.
my $pl = PDL::Planet->new
                    ->read("./images/earth.jpg")
                    ->transform(TYPE => Orthographic, SUNPOS => [0, 0], 
	                        SIZE => [500,500],    CENTER => [0,-90])
                    ->write("testimg2.jpg");
ok -s "testimg2.jpg", 'Set terminator via sun subpoint position';

# According to http://aa.usno.navy.mil/data/docs/AltAz.html
# on 2006, 3/21 at 12Z at lat=0, lon=0 the sun elevation will be 88.2 degrees
# (approximately zenith)
#
# This should show an image (testimg3.jpg) similar to testimg2.jpg.
#perldl> print TimeClass->new->set_ymdhms_gps(2006, 3, 21, 12, 0, 0)->get_julian;
#2453816 -> use this for test Julian date
$pl = PDL::Planet->new
                 ->read("./images/earth.jpg")
                 ->transform(TYPE => Orthographic, SUNTIME => 2453816,
                             SIZE => [500,500],    CENTER  => [0,-90])
                 ->write("testimg3.jpg");
ok -s "testimg3.jpg", 'Set terminator via time (vernal equinox)';

# According to http://aa.usno.navy.mil/data/docs/AltAz.html
# on 2006, 6/21 at 12Z at lat=0, lon=0 the sun elevation will be 66.6 degrees
# and the azimuth will be approximately due north.  So the sun should be on
# the right side of the image in the northern hemisphere.
#
#perldl> print TimeClass->new->set_ymdhms_gps(2006, 6, 21, 12, 0, 0)->get_julian;
#2453908 -> use this for test Julian date
$pl = PDL::Planet->new
                 ->read("./images/earth.jpg")
                 ->transform(TYPE => Orthographic, SUNTIME => 2453908,
                             SIZE => [500,500],    CENTER  => [0,-90])
                 ->write("testimg4.jpg");
ok -s "testimg4.jpg", 'Set terminator via time (summer solstice)';

$pl = PDL::Planet->new
                 ->read("./images/earth.jpg")
                 ->resize(1000,500)
                 ->terminator(2453908)
                 ->write("testimg5.jpg");
ok -s "testimg5.jpg", 'Set terminator via time (linear projection)';

### # This is commented out since gifsicle and imagemagick may not be installed 
### # (it takes a long time too)
### foreach my $hr (0..24) {
###   $hr = sprintf ("%02d", $hr);
###   $pl = PDL::Planet->new	
###     ->read("./images/earth.jpg")
###       ->transform(TYPE => Orthographic, SUNTIME => 2453908 + ($hr/24),
### 		  SIZE => [500,500],    CENTER  => [0,-90])
### 	->write("testimg5_.$hr.jpg");
###   system "convert testimg5_.$hr.jpg testimg5_.$hr.gif";
### }
### 
### system "gifsicle -O2 --colors=256 --delay=50 -o testimg5.gif testimg5_.*.gif";
### 
### ok -s "testimg5.gif", 'GIF animation of world turning';

1;



