/**----------------------------------------------------------------------
 * @name       RGB projection utilities
 *
 * @author     Doug Hunt
 * -----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**----------------------------------------------------------------------
 * @name       mag
 *
 * find the magnitude of a 3 vector
 *
 * @parameter
 * @ input:    v1 -- double array
 * @ output:   
 * @           mag -- double array
 *-----------------------------------------------------------------------*/
double mag (double v1[3]) {
  double mag;
  mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  return mag;
}

/* Earth radius -- ave of max and min Earth radii, km */
#define RE 6367.44395

struct mapdefaults {
  double minlat;
  double maxlat;
  double minlon;
  double maxlon;
  double centerlon;  /* center of projection  */
  double centerlat;  /* center of projection  */
  int    minx;       /* corresponds to minlon */
  int    maxx;
  int    miny;       /* corresponds to maxlat */
  int    maxy;
  double shrink;     /* factor to shrink the Earth to show more space around , 0-1) */
};

static struct mapdefaults mapset = {
 -90.0,
  90.0,
-180.0,
 180.0,

   0.0, // viewer lon (deg)
  40.0, // viewer lat (deg)

     0, 
 1999,
     0,
   999,
   1.0,
};

static double pi       = 3.1415926535897932384626;
static double ncxyz[3] = {-999.0, -999.0, -999.0};

/**----------------------------------------------------------------------
 * @name       setmap
 *
 * Set static map parameters
 *
 * @parameter
 * @ input:      
 * @ maxx
 * @ maxy
 *-----------------------------------------------------------------------*/
int setmap (int maxx, int maxy) {
  mapset.maxx   = maxx;
  mapset.maxy   = maxy;
}

/**----------------------------------------------------------------------
 * @name       setcenter
 *
 * Set static map parameters
 *
 * @parameter
 * @ input:      
 * @ lat
 * @ lon
 *-----------------------------------------------------------------------*/
int setcenter (double lat, double lon) {
  double d2r  = pi/180.0;
  double cxyz[3];

  mapset.centerlat = lat;
  mapset.centerlon = lon;

  // recompute projection center normal vector
  llh2xyz (lat * d2r, lon * d2r, 0.0,  cxyz);
  norm(cxyz, ncxyz);  

}

/**----------------------------------------------------------------------
 * @name       setshrink
 *
 * Set scale factor for Earth plotting (0 to 1, 1 means Earth fills whole
 * window.
 *
 * @parameter
 * @ shrink
 *-----------------------------------------------------------------------*/
int setshrink (double shrink) {
  mapset.shrink = shrink;
}


/**----------------------------------------------------------------------
 * @name       orthographic_llh2xy
 *
 * For an orthographic projection:  convert lat/lon/height to XY plot coords
 *
 * @parameter
 * @ input:    lat -- latitude (-90 to 90)
 * @           lon -- longitude (-180(E) to 180(W))
 * @           height -- km above the Earth's surface (spherical Earth)
 * @           
 * @ output:   x -- X value of this pixel
 * @           y -- Y value of this pixel
 *-----------------------------------------------------------------------*/
int orthographic_llh2xy (double lat, double lon, double height, int *xp, int *yp) {

  double x, y;

  double d2r  = pi/180.0;
  int   ysize = mapset.maxy - mapset.miny;
  int   xsize = mapset.maxx - mapset.minx;
  double asp  = (double)ysize/(double)xsize;
  double h    = (height + RE)/RE;  // height in Earth Radii
  double midx, midy;
  double cxyz[3], xyz[3], nxyz[3];
  double center_angle, cos_theta;
  double lon0 = mapset.centerlon * d2r;
  double lat0 = mapset.centerlat * d2r; 
  double viewer_radius;  // the apparent distance from the center of the earth of the
                         // current point, as seen from infinity

  *xp = -1;  /* error value */
  *yp = -1;  /* error value */

  lat *= d2r; // convert to radians
  lon *= d2r;

  // get rid of back side points

  // compute normalized vector to input lat/lon/height
  llh2xyz (lat, lon, height, xyz);
  norm(xyz, nxyz);  

  // compute angle between viewer location and input lat/lon/height
  // note:  ncxyz is a static variable computed when 'setcenter' is called
  dot(nxyz, ncxyz, &cos_theta);
  center_angle = acos(cos_theta);

  // check for satellites that are past the horizon but high
  // enough to be visible
  viewer_radius = cos (center_angle - pi/2.0) * mag(xyz);
  
  if ((center_angle > (pi/2.0)) && (viewer_radius < (RE))) return;

  // the projection itself -- From http://mathworld.wolfram.com/OrthographicProjection.html
  // X and Y range from -1 to 1, times h
  x  = h * cos(lat) * sin(lon - lon0);
  y  = h * (cos(lat0)* sin(lat) - sin(lat0)*cos(lat)*cos(lon - lon0));

  midx = ((double)xsize-1.0)/2.0;
  midy = ((double)ysize-1.0)/2.0;

  *xp = (int)(( x * midx * asp * mapset.shrink) + midx);  
  *yp = (int)(( y * midy * mapset.shrink)       + midy); 

  return;

}


/**----------------------------------------------------------------------
 * @name       xform_orthographic
 *
 * Convert a linear projection map into an orthographic projection
 * map of the same size.
 *
 * @parameter
 * @ input:    image -- pointer to input image
 * @           xsize -- x dimension of image
 * @           ysize -- y dimension of image
 * @ output:   
 * @           imageout -- pointer to output image (allocated by caller)
 *-----------------------------------------------------------------------*/
int xform_orthographic (unsigned char *image, int xsize, int ysize, unsigned char *imageout) {
  int i, j;
  int xi, yi, idx;
  double lat, lon;
  unsigned char this[3];

  setmap (xsize-1, ysize-1);
  for (j=0;j<ysize;j++) {
    for (i=0;i<xsize;i++) {	
      
      lat =  ((((double)j/(double)ysize) *-180.0) + 90);
      lon =  ((((double)i/(double)xsize) * 360.0) - 180);
      
      orthographic_llh2xy (lat, lon, 0.0, &xi, &yi);
      
      // bounds checks
      if (xi < 0) continue;
      if (yi < 0) continue;
      if (xi > xsize-1) continue;
      if (yi > ysize-1) continue;

      this[0] = image[3*xsize*j + 3*i    ];
      this[1] = image[3*xsize*j + 3*i + 1];
      this[2] = image[3*xsize*j + 3*i + 2];
      
      idx = 3*xsize*yi + 3*xi;
      imageout[idx  ] = this[0];
      imageout[idx+1] = this[1];
      imageout[idx+2] = this[2];

      idx = 3*xsize*yi + 3*(xi+1);
      imageout[idx  ] = this[0];
      imageout[idx+1] = this[1];
      imageout[idx+2] = this[2];

      idx = 3*xsize*(yi+1) + 3*xi;
      imageout[idx  ] = this[0];
      imageout[idx+1] = this[1];
      imageout[idx+2] = this[2];

      idx = 3*xsize*(yi+1) + 3*(xi+1);
      imageout[idx  ] = this[0];
      imageout[idx+1] = this[1];
      imageout[idx+2] = this[2];

    }
  }
}


/**----------------------------------------------------------------------
 * @name       xform_orthographic1
 *
 * Convert a linear projection map into an orthographic projection
 * map of the same size.  Uses a preferred reverse mapping method.
 *
 * @parameter
 * @ input:    image -- pointer to input image
 * @           xsize -- x dimension of image
 * @           ysize -- y dimension of image
 * @ output:   
 * @           imageout -- pointer to output image (allocated by caller)
 *-----------------------------------------------------------------------*/
int xform_orthographic1 (unsigned char *image, int xsize, int ysize, unsigned char *imageout) {
  int i, j;
  int xi, yi, idxout, idxin;
  double lat, lon, x, y, c, rho, phi, lam;
  double asp  = (double)ysize/(double)xsize;
  double d2r  = pi/180.0;
  double r2d  = 180.0/pi;
  double midx = ((double)xsize-1.0)/2.0;
  double midy = ((double)ysize-1.0)/2.0;
  double lon0 = mapset.centerlon * d2r;
  double lat0 = mapset.centerlat * d2r;

  for (j=0;j<ysize;j++) {
    for (i=0;i<xsize;i++) {	
      
      // scale i and j to range -1 to 1, adding in the shrink factor.
      x =      ((double)i - midx)/(midx * asp * mapset.shrink);
      y = -1.0*((double)j - midy)/(midy *       mapset.shrink);

      // printf ("y,x = %f, %f\n", y, x);

      // now do the inverse orthographic projection (from mathworld)
      rho = sqrt(x*x + y*y);
      if (rho > 1) continue; // check for other side of Earth

      c   = asin (rho);
      phi =        asin  (cos(c)*sin(lat0) + (y*sin(c)*cos(lat0))/rho);
      lam = lon0 + atan2 (x*sin(c), rho*cos(lat0)*cos(c) - y*sin(lat0)*sin(c));
      if (lam >  pi) lam -= 2.0*pi;
      if (lam < -pi) lam += 2.0*pi;

      // printf ("lat, lon = %f, %f\n", (r2d*phi), (r2d*lam));

      // convert from lat/lon radians to x/y coords 
      yi = (int)( (double)mapset.maxy * (((-r2d*phi) +  90.0)/180.0) );
      xi = (int)( (double)mapset.maxx * ((( r2d*lam) + 180.0)/360.0) );

      // printf ("yi, xi = %d, %d\n", yi, xi);

      // bounds checks
      if (xi < 0) continue;
      if (yi < 0) continue;
      if (xi > mapset.maxx) continue;
      if (yi > mapset.maxy) continue;

      idxin  = 3*mapset.maxx*yi + 3*xi;
      idxout = 3*xsize*j  + 3*i;
      imageout[idxout    ] = image[idxin    ];
      imageout[idxout + 1] = image[idxin + 1];
      imageout[idxout + 2] = image[idxin + 2];

    }
  }
}



/**----------------------------------------------------------------------
 * @name       llh2xyz
 *
 * Simple spherical to xyz coordinate conversion
 *
 * @parameter
 * @ input:    lat -- radians -pi/2 to pi/2
 * @           lon -- radians -pi to pi
 * @           height -- km
 * @ output:   
 * @           xyz -- double array, km
 *-----------------------------------------------------------------------*/
int llh2xyz (double lat, double lon, double height, double xyz[3]) {
  double colat = (pi/2.0 - lat);
  double r     = height + (RE);
  
  xyz[0] = r * cos(lon) * sin(colat);
  xyz[1] = r * sin(lon) * sin(colat);
  xyz[2] = r * cos(colat);
}

/**----------------------------------------------------------------------
 * @name       dot
 *
 * dot product of two 3 vectors
 *
 * @parameter
 * @ input:    v1
 * @           v2
 * @ output:   
 * @           dot -- double array
 *-----------------------------------------------------------------------*/
int dot (double v1[3], double v2[3], double *dot) {
  *dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/**----------------------------------------------------------------------
 * @name       norm
 *
 * normalize a 3 vector
 *
 * @parameter
 * @ input:    v1 -- double array
 * @ output:   
 * @           norm -- double array
 *-----------------------------------------------------------------------*/
int norm (double v1[3], double norm[3]) {
  double mag;
  mag = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  norm[0] = v1[0]/mag;
  norm[1] = v1[1]/mag;
  norm[2] = v1[2]/mag;
}


    



