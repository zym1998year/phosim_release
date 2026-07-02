// @brief CatalogObject
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "catalogobject.h"


CatalogObject::CatalogObject(CGParams* pCGPar)
{
  pCGParams =  pCGPar;
  //pRandom   =  prandom;
  // Need to set random seed  before we use random numbers
  fRandom.setSeed32( pCGParams->fRandomSeed); 

  // ************************************************
  // Load lsst chip center locations . Was retrieved from LSST's 
  // focalplanelayout.txt
  double lsstChipCenterXYmm[] = {-296.25, -254.0, -211.75, -169.25, -127.0, 
                                 -84.75, -42.25, 0.0, 42.25, 84.75, 127.0, 
                                 169.25, 211.75, 254.0, 296.25 };
  for (int i = 0; i < 15; i++) {
    //convert mm to deg: lsst focal plane = 20 arcsec/mm
    double chipXYDeg = lsstChipCenterXYmm[i]*20./3600.0;
    lsstChipCenterXYDeg.push_back(chipXYDeg);
  }
}
// ********************************************************************** 

void CatalogObject::GenerateSingleRandomRaDec()
{
  // ****************************************************************
  // We have to be a little clever in generating the random 
  // palcement. 
  // From:   // http://mathworld.wolfram.com/SpherePointPicking.html
  // 1: Pick randomly in an alt/az over a spherical cap centered at
  //  zenith. alt=90,az=0,0.
  // 2:Use the conversion program to convert the az/el to the 
  // correct ra/dec
  // ****************************************************************
  // Gen random az,el in cap with fZenithMaxRange. Set in Constructor
  // ***************************************
  //double azimuthRad=pRandom->uniform()  * 2 * PAL__DPI;  //Uniform in azimuth
  double azimuthRad=fRandom.uniform()  * 2 * PAL__DPI;  //Uniform in azimuth
                                          //zenith 0 to zenithMaxRange in cap
  //double zenithRad=acos( 1. - 2. * pRandom->uniform() * fZenithMaxRange);
  double zenithRad=acos( 1. - 2. * fRandom.uniform() * fZenithMaxRange);
  double elevationRad = PAL__DPI / 2. - zenithRad;
  
  double longRad = pCGParams->fCatalogFOVRADeg  * PAL__DPI / 180.0;
  double latRad  = pCGParams->fCatalogFOVDecDeg * PAL__DPI / 180.0;
  double haRad;
 
  // *********************************
  // Use generic conversion function from Slalib (From pallib actually for 
  // licensing issues)
  // ******************************** 

  palDh2e(azimuthRad,elevationRad ,latRad ,&haRad, &fDecRad );
  fRARad=longRad-haRad;
  if (fRARad<0.0) {
	fRARad= 2.0*PAL__DPI+fRARad;
  }
  return;
}
// ****************************************************************


void CatalogObject::GenerateGridRaDec(double RARad,  double DecRad, 
									  double FOVDeg, double StepSizeDeg)
// *****************************************************************
// Generates a square grid of focalplane locations(in deg) and 
// converts to equivalent ( tangent plane or gnomonic) projection to
// the celestial sphere in ra/dec centered at the specified ra/dec
// *****************************************************************
// Signs of directions will be determined later. Just go with 
// "whatever" for now.
// ****************************************************************
// Make up a vector of vectors grid of focal plane locations. 
// Outermost is "Y" which is along altitude coord.
// Results in Catalog members: "rightAscension" and "declination" vectors.
// ****************************************************************
{
  bool debugPrint = false;
  fRightAscension.clear();
  fDeclination.clear();


  // ****************************************************************
  // Determin xPintDeg,yPointDeg  (in focal plane);
  // *******************************
  bool makeLSSTChipCenterGrid = false;
  //bool makeLSSTChipCenterGrid = true;    //enable grid of  chip centers
  int  npoints;
  double halfWidth;
  if (makeLSSTChipCenterGrid ) {
    npoints = lsstChipCenterXYDeg.size();
    std::cout << " #Star Grid: LSST chip centers" << std::endl;
  }
  else{
    npoints =  FOVDeg / StepSizeDeg;
    // Since we want a point at the very center this should be odd
    if (npoints%2 == 0) {
      npoints=npoints+1;
    }
    halfWidth=( ( npoints - 1 ) /2 ) * StepSizeDeg;
    //std::cout << " #Star Grid spacing(deg):" << StepSizeDeg 
    //          << std::endl;
  }


  //std::cout << " #Number of points across grid: " << npoints << std::endl;
  //std::cout << " #Number of points in grid:     " << npoints*npoints 
  //		<< std::endl;
  
  if (npoints*npoints > kPrintFreq) {
	std::cout << " This may take a few minutes. Be patient!" << std::endl;
	std::cout <<" Grid Ra/Dec creation progress: each # = " << kPrintFreq 
			  <<" grid points completed" << std::endl;
	std::cout<< " " << std::flush; //Start off a little offset to lineup with
	                               // above 
  }
    
  double xPointDeg;
  double yPointDeg;

  for ( int j = 0; j < npoints; j++ ) {

    if (!makeLSSTChipCenterGrid) {
      yPointDeg= (- halfWidth) + j * StepSizeDeg;
    }
    else{
      yPointDeg = lsstChipCenterXYDeg.at(j);
    }

    for ( int i = 0 ;i < npoints ; i++ ) {

      if (!makeLSSTChipCenterGrid) {
        xPointDeg= (- halfWidth) + i * StepSizeDeg;
      }
      else{
        xPointDeg = lsstChipCenterXYDeg.at(i);
      }

      if (debugPrint) {
          std::cout << "xPointDeg, yPointDeg: "<<  xPointDeg << " " << yPointDeg
                    << std::endl;
      }

      // *******************************************************
      // Project xPointDeg and yPointDeg back to celestial sphere and 
      // get Ra/dec for this grid point.
      // ********************************************************
      double raGridRad;
      double decGridRad;

      // **********************************************************
      // Use palDtp2s. (Tangent Plane to Spherical)  gnomonic projection
      // function which we coppied to here. (Allowed under pallib/erfam
      // license)
      // **********************************************************

      palDtp2s(xPointDeg * PAL__DPI / 180., yPointDeg *  PAL__DPI / 180., 
               RARad,  DecRad, &raGridRad,  &decGridRad);

      // *****************************
      // Save the ra/dec of the grid point.
      // ****************************
      fRightAscension.push_back(raGridRad);
      fDeclination.push_back(decGridRad);

      if (debugPrint) {
          std::cout<<xPointDeg<< " " <<  yPointDeg << " " << raGridRad << " "
                   << decGridRad<<std::endl;
          std::cout << "at point: " <<  fRightAscension.size() 
                    << ": raGridRad,decGridRad: "<<  raGridRad << " " 
                    << decGridRad << std::endl;
      }

      if ( fRightAscension.size()%kPrintFreq == 0) {
          std::cout<< "#" << std::flush;
          if  ( fRightAscension.size()%(kPrintFreq*10) == 0) {
              std::cout<< "|" << std::flush;
          }
      }
    }
  }
  return;
}
// ****************************************************************


void CatalogObject::GetSphericalCellsInFOV(double FOVDiameterDeg) 
{
  // ******************************************************************
  // Find all the ra/dec cells that have any area in the FOV (centered at 
  // pCGParams->fCatalogFOVRADeg,pCGParams->fCatalogFOVDecDeg). As we find 
  // each cell place it in the vector "Cells". 
  // ******************************************************************
  // Basic approcch is to divide the celestial sphere (using theta,phi coords)
  // into cells. The height of each cell is ~kThetaHeightDeg. The width is 
  // ~kPhiWidthDeg (though that may change in the future to being a width 
  // inversly propotional to cos(Dec).  But, keep it simple for now)
  // Note: we make fine adjustments to the cell spacing so that they exactly 
  // cover all of the clestial sphere.
  // ******************************************************************
  // Use a brute force search. First increase the FOV radius by the larger of 
  // Phi angular Spacing or the Theta angular Spacing. Then, test each and 
  // every cell to see if it's center is within this expanded FOV. We can use 
  // the algorithum used for the SlaLib function SLA_DSEP to determine if the 
  // center of a cell is within the expanded FOV.
  // ********************************************************************

  //Theta,phi of FOV origin
  double originPhiDeg =  pCGParams->fCatalogFOVRADeg;
  double originThetaDeg = 90.-pCGParams->fCatalogFOVDecDeg;

  // make a grid on the sphere. Make sure we cover it all
  fNumThetaCells = 180/kThetaHeightDeg;  // note psossible round down
  fNumPhiCells   = 360/kPhiWidthDeg;
  
  // Adjust bin spacing to exactly fit(in case we had round down).
  // Note for phi the actutal spacing in distance on the sphere, the angular 
  // distance, is largest at dec=0 (theta=90.)
  double thetaCellSpacingDeg= 180./(double)fNumThetaCells;
  double phiCellSpacingDeg = 360./(double)fNumPhiCells;

  // ***********************************************************************
  // We need to find all those Theta,Phi cells that have any area in the FOV. 
  // To do this we add the larger of Theta angular distance spacing or phi 
  // angular distance spacing at dec=0 (theta=90.) to the FOV angular radius. 
  // This is to insure we find all possible cells, This is obviously way 
  // overkill (phi angular distance spacing is largest at the equator (dec=0)) 
  // I know,  but I think it results in us not missing anything.
  // ***********************************************************************
  double modifiedFOVRadiusDeg = FOVDiameterDeg/2;
  if (thetaCellSpacingDeg > phiCellSpacingDeg) {
    //Add another 10% for good luck
    modifiedFOVRadiusDeg += 1.1*thetaCellSpacingDeg; 
  }
  else{
    modifiedFOVRadiusDeg += 1.1*phiCellSpacingDeg;
  }   

  // Now iterate through ALL the cells, testing each one to see if it might
  // be in the FOV.
  // Following  could be more efficent(could make some obvious cuts on theta 
  // range for example) but this is clearer.
  for (int i=0; i < fNumThetaCells; i++ ) {
    double thetaDeg =   (i * thetaCellSpacingDeg ) + thetaCellSpacingDeg / 2.;
    for (int j=0; j <  fNumPhiCells; j++ ) {
      // Find Phi of the center of this cell
      double phiDeg = (j * phiCellSpacingDeg) + phiCellSpacingDeg / 2.;

      //See if cells center is within expanded FOV. Use Sofa eraSeps function
      double angDistDeg =  eraSeps(phiDeg * (PAL__DPI/180.), 
                                   (90.-thetaDeg) * (PAL__DPI/180.), 
                                   originPhiDeg *(PAL__DPI/180.), 
                                   (90.-originThetaDeg) *(PAL__DPI/180.) ) * 
                                                               (180./PAL__DPI);
      if (angDistDeg <= modifiedFOVRadiusDeg) {
        //This cell is in the field of view (FOV). Fill its values

        SphericalCell cell;
        cell.fRADeg      = phiDeg;            //Phi same as RA in deg
        cell.fDecDeg     = 90.-thetaDeg;     //Theta = 90-dec
        cell.fThetaIndex = i;
        cell.fPhiIndex   = j;
        // **********************
        // Now find area of this cell in deg**2
        // Get cell boundries
        // **********************
        cell.SetBoundries(thetaDeg, thetaCellSpacingDeg, phiDeg, 
                                                           phiCellSpacingDeg);
        cell.SetAreaDeg2();
        cell.SetGalacticFromRaDec();
 
        // And save this "hit" cell in the Catalog object vector fCells
        fCells.push_back(cell);
      }
    }
  }
  return;
}
// ***********************************************************************

int CatalogObject::ClosestElement(std::vector<int>* pVec, double value) 
// ************************************************************** 
// Retuns closest value in a sorted int vector to value (which is a double)
// pVec has no repeating elements.
// **************************************************************
{
  // First find lower value that bound us. The upper will be the next one 
  // since pVec is sorted.
  std::vector<int>::iterator low;
  std::vector<int>::iterator up;
  std::vector<int> vec = *pVec;  // I know this is inefficent.

  if ( value > vec.back() ) {
    return vec.back();
  }
  if (value < vec.front() ) {
    return vec.front();
  }
  up = std::upper_bound( vec.begin(), vec.end(), value); 
  low = up -1;
  
  // Now find which of *low or *up is closer to value   
  double lowDistance = fabs(value - *low);
  double hiDistance = fabs(value - *up);
  if ( lowDistance < hiDistance) {
    return *low;
  }
  else{
    return *up;
  }
}
// ******************************************************************
double CatalogObject::ClosestElement(std::vector<double>* pVec, double value) 
// ************************************************************** 
// Retuns closest value in a sorted double vector to value (which is a double)
// pVec has no repeating elements.
// **************************************************************
{
  // First find lower value that bound us. The upper will be the next one 
  // since pVec is sorted.
  std::vector<double>::iterator low;
  std::vector<double>::iterator up;
  std::vector<double> vec = *pVec;  // I know this is inefficent.

  if ( value > vec.back() ) {
    return vec.back();
  }
  if (value < vec.front() ) {
    return vec.front();
  }
  up = std::upper_bound( vec.begin(), vec.end(), value); 
  low = up -1;
  
  // Now find which of *low or *up is closer to value   
  double lowDistance = fabs(value - *low);
  double hiDistance = fabs(value - *up);
  if ( lowDistance < hiDistance) {
    return *low;
  }
  else{
    return *up;
  }
}
// ******************************************************************


// **************************************************************************
// All of the following pal* and era* functions are used under the following 
// copyright and license.
// Included explicitly here to minimize dependencies.
// ***************************************************************************

/*
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  All rights reserved.
**  
**  This library is derived, with permission, from the International
**  Astronomical Union's "Standards of Fundamental Astronomy" library,
**  available from http://www.iausofa.org.
**  
**  The ERFA version is intended to retain identical functionality to
**  the SOFA library, but made distinct through different function and
**  file names, as set out in the SOFA license conditions.  The SOFA
**  original has a role as a reference standard for the IAU and IERS,
**  and consequently redistribution is permitted only in its unaltered
**  state.  The ERFA version is not subject to this restriction and
**  therefore can be included in distributions which do not support the
**  concept of "read only" software.
**  
**  Although the intent is to replicate the SOFA API (other than
**  replacement of prefix names) and results (with the exception of
**  bugs;  any that are discovered will be fixed), SOFA is not
**  responsible for any errors found in this version of the library.
**  
**  If you wish to acknowledge the SOFA heritage, please acknowledge
**  that you are using a library derived from SOFA, rather than SOFA
**  itself.
**  
**  
**  TERMS AND CONDITIONS
**  
**  Redistribution and use in source and binary forms, with or without
**  modification, are permitted provided that the following conditions
**  are met:
**  
**  1 Redistributions of source code must retain the above copyright
**    notice, this list of conditions and the following disclaimer.
**  
**  2 Redistributions in binary form must reproduce the above copyright
**    notice, this list of conditions and the following disclaimer in
**    the documentation and/or other materials provided with the
**    distribution.
**  
**  3 Neither the name of the Standards Of Fundamental Astronomy Board,
**    the International Astronomical Union nor the names of its
**    contributors may be used to endorse or promote products derived
**    from this software without specific prior written permission.
**  
**  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
**  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
**  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
**  FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
**  COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
**  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
**  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
**  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
**  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
**  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
**  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
**  POSSIBILITY OF SUCH DAMAGE.
**  
*/
//#include "erfa.h"

double CatalogObject::eraAnp(double a)
/*
**  - - - - - - -
**   e r a A n p
**  - - - - - - -
**
**  Normalize angle into the range 0 <= a < 2pi.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range 0-2pi
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   double w;

   w = fmod(a, ERFA_D2PI);
   if (w < 0) w += ERFA_D2PI;

   return w;

}

double CatalogObject::eraPdp(double a[3], double b[3])
/*
**  - - - - - - -
**   e r a P d p
**  - - - - - - -
**
**  p-vector inner (=scalar=dot) product.
**
**  Given:
**     a      double[3]     first p-vector
**     b      double[3]     second p-vector
**
**  Returned (function value):
**            double        a . b
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double w;

   w  = a[0] * b[0]
      + a[1] * b[1]
      + a[2] * b[2];

   return w;

}
// ***************************************************************************

double CatalogObject::eraPm(double p[3])
/*
**  - - - - - -
**   e r a P m
**  - - - - - -
**
**  Modulus of p-vector.
**
**  Given:
**     p      double[3]     p-vector
**
**  Returned (function value):
**            double        modulus
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   return sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );

}
// ***************************************************************************

void CatalogObject::eraPxp(double a[3], double b[3], double axb[3])
/*
**  - - - - - - -
**   e r a P x p
**  - - - - - - -
**
**  p-vector outer (=vector=cross) product.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     axb      double[3]      a x b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double xa, ya, za, xb, yb, zb;

   xa = a[0];
   ya = a[1];
   za = a[2];
   xb = b[0];
   yb = b[1];
   zb = b[2];
   axb[0] = ya*zb - za*yb;
   axb[1] = za*xb - xa*zb;
   axb[2] = xa*yb - ya*xb;

   return;

}
// ***************************************************************************

void CatalogObject::eraS2c(double theta, double phi, double c[3])
/*
**  - - - - - - -
**   e r a S 2 c
**  - - - - - - -
**
**  Convert spherical coordinates to Cartesian.
**
**  Given:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**
**  Returned:
**     c        double[3]    direction cosines
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double cp;

   cp = cos(phi);
   c[0] = cos(theta) * cp;
   c[1] = sin(theta) * cp;
   c[2] = sin(phi);

   return;

}
// ***************************************************************************

double CatalogObject::eraSepp(double a[3], double b[3])
/*
**  - - - - - - - -
**   e r a S e p p
**  - - - - - - - -
**
**  Angular separation between two p-vectors.
**
**  Given:
**     a      double[3]    first p-vector (not necessarily unit length)
**     b      double[3]    second p-vector (not necessarily unit length)
**
**  Returned (function value):
**            double       angular separation (radians, always positive)
**
**  Notes:
**
**  1) If either vector is null, a zero result is returned.
**
**  2) The angular separation is most simply formulated in terms of
**     scalar product.  However, this gives poor accuracy for angles
**     near zero and pi.  The present algorithm uses both cross product
**     and dot product, to deliver full accuracy whatever the size of
**     the angle.
**
**  Called:
**     eraPxp       vector product of two p-vectors
**     eraPm        modulus of p-vector
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double axb[3], ss, cs, s;

/* Sine of angle between the vectors, multiplied by the two moduli. */
   eraPxp(a, b, axb);
   ss = eraPm(axb);

/* Cosine of the angle, multiplied by the two moduli. */
   cs = eraPdp(a, b);

/* The angle. */
   s = ((ss != 0.0) || (cs != 0.0)) ? atan2(ss, cs) : 0.0;

   return s;

}
// ****************************************************************************

#include "../raytrace/constants.h"
double CatalogObject::angularSep(double ra1, double dec1, double ra2, double dec2)

{
    double dra, cosdis;
    dra = fabs(ra1 - ra2);
    if (dra > PI) dra = 2*PI - dra;
    cosdis = sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(dra);
    if (cosdis > 1) cosdis = 1.0;
    if (cosdis < -1) cosdis = -1.0;
    return(acos(cosdis));
}

double CatalogObject::eraSeps(double al, double ap, double bl, double bp)
/*
**  - - - - - - - -
**   e r a S e p s
**  - - - - - - - -
**
**  Angular separation between two sets of spherical coordinates.
**
**  Given:
**     al     double       first longitude (radians)
**     ap     double       first latitude (radians)
**     bl     double       second longitude (radians)
**     bp     double       second latitude (radians)
**
**  Returned (function value):
**            double       angular separation (radians)
**
**  Called:
**     eraS2c       spherical coordinates to unit vector
**     eraSepp      angular separation between two p-vectors
**
**  Copyright (C) 2013-2015, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes above.
*/
{
   double ac[3], bc[3], s;

/* Spherical to Cartesian. */
   eraS2c(al, ap, ac);
   eraS2c(bl, bp, bc);

/* Angle between the vectors. */
   s = eraSepp(ac, bc);

   return s;

}
// ****************************************************************************

/*
*+
*  Name:
*     palDh2e

*  Purpose:
*     Horizon to equatorial coordinates: Az,El to HA,Dec

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDh2e( double az, double el, double phi, double * ha, double * dec );

*  Arguments:
*     az = double (Given)
*        Azimuth (radians)
*     el = double (Given)
*        Elevation (radians)
*     phi = double (Given)
*        Observatory latitude (radians)
*     ha = double * (Returned)
*        Hour angle (radians)
*     dec = double * (Returned)
*        Declination (radians)

*  Description:
*     Convert horizon to equatorial coordinates.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  Notes:
*     - All the arguments are angles in radians.
*     - The sign convention for azimuth is north zero, east +pi/2.
*     - HA is returned in the range +/-pi.  Declination is returned
*       in the range +/-pi/2.
*     - The latitude is (in principle) geodetic.  In critical
*       applications, corrections for polar motion should be applied.
*     - In some applications it will be important to specify the
*       correct type of elevation in order to produce the required
*       type of HA,Dec.  In particular, it may be important to
*       distinguish between the elevation as affected by refraction,
*       which will yield the "observed" HA,Dec, and the elevation
*       in vacuo, which will yield the "topocentric" HA,Dec.  If the
*       effects of diurnal aberration can be neglected, the
*       topocentric HA,Dec may be used as an approximation to the
*       "apparent" HA,Dec.
*     - No range checking of arguments is done.
*     - In applications which involve many such calculations, rather
*       than calling the present routine it will be more efficient to
*       use inline code, having previously computed fixed terms such
*       as sine and cosine of latitude.

*  History:
*     2012-02-08 (TIMJ):
*        Initial version with documentation taken from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1996 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

//#include "pal.h"
#include <math.h>

void CatalogObject::palDh2e ( double az, double el, double phi, double *ha, 
                              double *dec) 
{
  double sa;
  double ca;
  double se;
  double ce;
  double sp;
  double cp;

  double x;
  double y;
  double z;
  double r;

  /*  Useful trig functions */
  sa = sin(az);
  ca = cos(az);
  se = sin(el);
  ce = cos(el);
  sp = sin(phi);
  cp = cos(phi);

  /*  HA,Dec as x,y,z */
  x = -ca * ce * sp + se * cp;
  y = -sa * ce;
  z = ca * ce * cp + se * sp;

  /*  To HA,Dec */
  r = sqrt(x * x + y * y);
  if (r == 0.) {
    *ha = 0.;
  } else {
    *ha = atan2(y, x);
  }
  *dec = atan2(z, r);

  return;
}
// ************************************************************************
/*
*+
*  Name:
*     palDtp2s

*  Purpose:
*     Tangent plane to spherical coordinates

*  Language:
*     Starlink ANSI C

*  Type of Module:
*     Library routine

*  Invocation:
*     palDtp2s( double xi, double eta, double raz, double decz,
*               double *ra, double *dec);

*  Arguments:
*     xi = double (Given)
*        First rectangular coordinate on tangent plane (radians)
*     eta = double (Given)
*        Second rectangular coordinate on tangent plane (radians)
*     raz = double (Given)
*        RA spherical coordinate of tangent point (radians)
*     decz = double (Given)
*        Dec spherical coordinate of tangent point (radians)
*     ra = double * (Returned)
*        RA spherical coordinate of point to be projected (radians)
*     dec = double * (Returned)
*        Dec spherical coordinate of point to be projected (radians)

*  Description:
*     Transform tangent plane coordinates into spherical.

*  Authors:
*     PTW: Pat Wallace (STFC)
*     TIMJ: Tim Jenness (JAC, Hawaii)
*     {enter_new_authors_here}

*  History:
*     2012-02-08 (TIMJ):
*        Initial version with documentation taken from Fortran SLA
*        Adapted with permission from the Fortran SLALIB library.
*     {enter_further_changes_here}

*  Copyright:
*     Copyright (C) 1995 Rutherford Appleton Laboratory
*     Copyright (C) 2012 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Bugs:
*     {note_any_bugs_here}
*-
*/

//#include "pal.h"
//#include "pal1sofa.h"

#include <math.h>

void CatalogObject::palDtp2s ( double xi, double eta, double raz, double decz,
                               double *ra, double *dec ) 
{
  double cdecz;
  double denom;
  double sdecz;
  double d;

  sdecz = sin(decz);
  cdecz = cos(decz);
  denom = cdecz - eta * sdecz;
  d = atan2(xi, denom) + raz;
  *ra = eraAnp(d);
  *dec = atan2(sdecz + eta * cdecz, sqrt(xi * xi + denom * denom));

  return;
}
// ***************************************************************************

