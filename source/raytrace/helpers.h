///
/// @package phosim
/// @file helpers.h
/// @brief helpers
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <vector>
#include <stdlib.h>
#include <stdio.h>

inline void find_v (std::vector<double>  & y, double x, int *index) {
    int high, mid, low;
    int n = y.size();


    low = 0;
    high = n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x > y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }

inline void find_v_n (std::vector<double>  & y, int l, int n, double x, int *index) {
    int high, mid, low;

    low = l;
    high = l + n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x > y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }

inline void find_v_reverse (std::vector<double>  & y, double x, int *index) {
    int high, mid, low;
    int n = y.size();


    low = 0;
    high = n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x < y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }

inline void find (double *y, int n, double x, int *index) {
    int high, mid, low;

    low = 0;
    high = n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x > y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }

inline void find_reverse (double *y, int n, double x, int *index) {
    int high, mid, low;

    low = 0;
    high = n - 1;
    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (x < y[mid]) {
            low = mid;
        } else {
            high = mid;
        }
    }
     *index = low;
 }

inline int find_linear (double *xx, int n, double x, double *dj) {
    int j;
    *dj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    j = (int)(*dj);
    if (j > n - 2) j = n - 2;
    if (j < 0) j = 0;
    return j;
 }

inline int find_linearFloat (float *xx, int n, float x, float *dj) {
    int j;
    *dj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    j = (int)(*dj);
    if (j > n - 2) j = n - 2;
    if (j < 0) j = 0;
    return j;
 }

inline int find_linear_float_mc (float *xx, int n, float x, float random) {

    float rj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    int j = floorf(rj);
    rj -= j;
    if (random > rj) j++;
    if (j > n - 1) j = n - 1;
    if (j < 0) j = 0;
    return j;

}


inline int find_linear_v (std::vector<double> & xx, double x, double *dj) {
    int j;
    int n = xx.size();

    *dj = (x - xx[0])/(xx[n - 1] - xx[0])*(n - 1);
    j = (int)(*dj);
    if (j > n - 2) j = n - 2;
    if (j < 0) j = 0;

    return j;
}

inline void find_linear_wrap (double x, double dx, int n, int *i, int *j, double *di) {

    double xx, fxx;

    xx = x/dx + (static_cast<double>(n)/2 - 0.5);
    fxx = floor(xx);

    *i = static_cast<int>(fxx);
    *j = (*i) + 1;

    *i = (*i) % n;
    if ((*i) < 0) (*i) += n;

    *j = (*j) % n;
    if ((*j) < 0) (*j) += n;

    *di = xx - fxx;

}

inline void find_linear_wrapFloat (double x, double dx, int n, int *i, int *j, float *di) {

    double xx, fxx;

    xx = x/dx + (static_cast<double>(n)/2 - 0.5);
    fxx = floor(xx);

    *i = static_cast<int>(fxx);
    *j = (*i) + 1;

    *i = (*i) % n;
    if ((*i) < 0) (*i) += n;

    *j = (*j) % n;
    if ((*j) < 0) (*j) += n;

    *di = xx - fxx;

}


inline void normalize (double *vx, double *vy, double *vz) {

    double r;

    r = sqrt(((*vx)*(*vx)) + ((*vy)*(*vy)) + ((*vz)*(*vz)));

    *vx = *vx/r;
    *vy = *vy/r;
    *vz = *vz/r;

}

inline double interpolate_fraction (double *xx, double x, int index) {
    return((x-  xx[index])/(xx[index + 1] - xx[index]));
}

inline double interpolate_fraction_v (std::vector<double> & xx, double x, int index) {
    return((x-  xx[index])/(xx[index + 1] - xx[index]));
}

inline double interpolate (double *yy, double *xx, double x, int index) {

    return ((yy[index]*(xx[index + 1] - x) + yy[index + 1]*(x - xx[index]))/(xx[index + 1] - xx[index]));

}

inline double interpolate_float (float *yy, float *xx, float x, int index) {

    return ((yy[index]*(xx[index + 1] - x) + yy[index + 1]*(x - xx[index]))/(xx[index + 1] - xx[index]));

}

inline double interpolate_v (std::vector<double> & yy, std::vector<double> & xx, double x, int index) {

    return ((yy[index]*(xx[index + 1] - x) + yy[index + 1]*(x - xx[index]))/(xx[index + 1] - xx[index]));

}

inline double interpolate_linear (double *yy, int index,  double rindex) {

    double dr = rindex - index;
    return (yy[index]*(1.0 - dr) + yy[index + 1]*dr);

}

inline void rotation (double *xp, double *yp, double angle, double x, double y) {

    double sinAngle = sin(angle);
    double cosAngle = cos(angle);
    *xp = cosAngle*x + sinAngle*y;
    *yp = -sinAngle*x + cosAngle*y;

}

inline double interpolate_bilinear (double *zz, int nelem, int xindex, double rxindex, int yindex, double ryindex) {

    double dx = rxindex - xindex;
    double dy = ryindex - yindex;
    double dxdy = dx*dy;
    int location = xindex*nelem + yindex;
    int location2 = location + nelem;

    return (zz[location]*(1 - dx - dy + dxdy) +
            zz[location2]*(dx - dxdy)+
            zz[location + 1]*(dy - dxdy) +
            zz[location2 + 1]*dxdy);

}

inline float interpolate_bilinear_float(float *zz, int nelem, int xindex, float rxindex, int yindex, float ryindex) {

    float dx = rxindex - xindex;
    float dy = ryindex - yindex;
    float dxdy = dx*dy;
    int location = xindex*nelem + yindex;
    int location2 = location + nelem;

    return (zz[location]*(1 - dx - dy + dxdy) +
            zz[location2]*(dx - dxdy)+
            zz[location + 1]*(dy - dxdy) +
            zz[location2 + 1]*dxdy);


}

inline float interpolate_trilinear (float *zz, int nelemy, int nelemz, int xindex, float rxindex,
                                     int yindex, float ryindex, int zindex, float rzindex ) {

    float dx = rxindex - xindex;
    float dy = ryindex - yindex;
    float dz = rzindex - zindex;
    float dxdy = dx*dy;
    float dz1 = 1 - dz;
    int location = xindex*nelemy*nelemz + yindex*nelemz + zindex;
    int location2 = location + nelemy*nelemz;

    return (zz[location]*(1 - dx - dy + dxdy)*dz1 +
            zz[location2]*(dx - dxdy)*dz1 +
            zz[location + nelemz]*(dy - dxdy)*dz1 +
            zz[location2 + nelemz]*(dxdy)*dz1 +
            zz[location + 1]*(1 - dx - dy + dxdy)*dz +
            zz[location2 + 1]*(dx - dxdy)*dz +
            zz[location + nelemz + 1]*(dy - dxdy)*dz +
            zz[location2 + nelemz + 1]*(dxdy)*dz);

}

inline double interpolate_bilinear_float_wrap (float *zz, int nelem, int xindex, int xindex1, double dx, int yindex, int yindex1, double dy ) {

    double dxdy = dx*dy;

    return (zz[xindex*nelem  + yindex ]*(1 - dx - dy + dxdy) +
            zz[xindex1*nelem + yindex ]*(dx - dxdy) +
            zz[xindex*nelem  + yindex1]*(dy - dxdy) +
            zz[xindex1*nelem + yindex1]*(dxdy));

}



// Interpolation search
// asserts that: sorted_array[0] <= value
// returns position of first element less or equal to parameter value
inline int interpolationSearch(double sorted_array[], int n, double value) {
    int beg = 0;
    int end = n - 1;

    double begVal = sorted_array[beg];
    double endVal = sorted_array[end];

    if (endVal <= value) {
        return end;
    }
    //invariant for following loop: endVal>value

    while (begVal < endVal) {
        int    testPos = beg + int((value - begVal)*(end - beg)/(endVal - begVal));
        double testVal = sorted_array[testPos];

        if (value < testVal) {
            end = testPos;
            endVal = testVal;
        } else {
            if (testPos > beg){
                beg = testPos;
                begVal = testVal;
            } else {
                beg++;
                begVal = sorted_array[beg];
                if (begVal > value) {
                    return beg - 1;
                }
            }
        }
    }
    return beg;
}

inline double smallAnglePupilNormalize (double& x, double& y) {
  double xymag2 = x*x + y*y;
  //Prevent returning of nan for z component.
  if (xymag2 >= 1.0 ) {
    //At this point we will assume that z should be  0 and that we need 
    // to narmalize the x,y,z vector to have amagnitude of 1.0. Thus we will
    // normlize x and y
    double xymag = sqrt(xymag2);
    x = x/xymag;
    y = y/xymag;
    return 0.0;
  }
  else { 
    return(-sqrt(1.0 - xymag2));
  }
}


#ifdef ISO_LIBRARIES_ONLY
    inline bool isnan(double d) { return d != d; }
#endif
