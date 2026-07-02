///
/// @package phosim
/// @file mie.cpp
/// @brief mie physics
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///


#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <vector>

#include "mie.h"
#include "ancillary/readtext.h"
#include "helpers.h"
#include "constants.h"
#include "parameters.h"



void Mie::mieTyndall(double radius, double wavelength, double re, double im, double *extinction, double *scattering, int angles, std::vector<double>& phi, std::vector<double>& psf,
                     std::vector<std::complex<double> >& s1, std::vector<std::complex<double> >& s2,
                     std::vector<double>& mu, std::vector<double>& piN,
                     std::vector<double>& piNm2,std::vector<double>& piNm1,
                     std::vector<double>& tauN) {


    double x = 2*PI*radius/wavelength;
    std::complex<double> rindex = {re,im};

    //  std::vector<std::complex<double> > s1(angles), s2(angles);
    //std::vector<double> mu(angles), piN(angles), piNm2(angles), piNm1(angles), tauN(angles);

    double xn = x + 4.0*pow(x,1.0/3.0) + 2.0;
    int nm;
    if (xn > abs(x*rindex)) {
        nm = static_cast<int>(xn + 15);
    } else {
        nm = static_cast<int>(abs(x*rindex) + 15);
    }
    int ns = static_cast<int>(xn);

     std::vector<std::complex<double> > d(nm);

    for (int j=0; j < angles; j++) {
        mu[j] = cos(static_cast<double>(j)*PI/(static_cast<double>(angles - 1)));
        piNm2[j] = 0.0;
        piNm1[j] = 1.0;
        s1[j] = {0.0, 0.0};
        s2[j] = {0.0, 0.0};
    }

    d[nm - 1] = {0.0, 0.0};
    for (int n=0; n < (nm-1); n++) {
        d[nm - n - 2] = (static_cast<double>(nm - n)/(x*rindex)) - (1.0/ (d[nm-n-1]+static_cast<double>(nm - n)/(x*rindex)));
    }

    double psiNm2 = cos(x);
    double psiNm1 = sin(x);
    double chiNm2 = -sin(x);
    double chiNm1 = cos(x);
    std::complex<double> xiNm1 = {psiNm1,-chiNm1};
    std::complex<double> an;
    std::complex<double> bn;

    *scattering = 0.0;
    *extinction = 0.0;
    for (int n=0; n < ns; n++) {

        double np = static_cast<double>(n + 1);
        double fn = (2*np + 1)/(np*(np + 1));

        double psiN = (2*np - 1)*psiNm1/x - psiNm2;
        double chiN = (2*np - 1)*chiNm1/x - chiNm2;
        std::complex<double> xiN = {psiN,-chiN};

        an = ((d[n]/rindex + np/x)*psiN - psiNm1)/ ((d[n]/rindex + np/x)*xiN - xiNm1);
        bn = ((rindex*d[n] + np/x)*psiN - psiNm1)/ ((rindex*d[n] + np/x)*xiN - xiNm1);

        *scattering += (2*np + 1)*2/x/x* (norm(an) + norm(bn));
        *extinction += (2*np + 1)*2/x/x* (real(an) + real(bn));

        for (int j=0; j < angles; j++) {
            piN[j] = piNm1[j];
            tauN[j] = np*mu[j]*piN[j] - (np + 1.0)*piNm2[j];
            s1[j] += fn* (an*piN[j] + bn*tauN[j]);
            s2[j] += fn* (an*tauN[j] + bn*piN[j]);
        }

        psiNm2 = psiNm1;
        psiNm1 = psiN;
        chiNm2 = chiNm1;
        chiNm1 = chiN;
        xiNm1 = {psiNm1, -chiNm1};

        for (int j = 0;  j < angles; j++) {
            piNm1[j] = ((2*np + 1)*mu[j]*piN[j]- (np + 1)*piNm2[j])/np;
            piNm2[j] = piN[j];
        }

    }


    for (int i =0; i < angles; i++) {
        phi[i] = 180.0*(static_cast<double>(i))/(static_cast<double>(angles-1));
        psf[i] = (norm(s1[i]) + norm(s2[i]))/2.0;
    }

}


void Mie::setupMieGrid(int elementsWavelength, int elementsX, double minWavelength, double maxWavelength, double minX, double maxX, double *muValue, double *cumulativeMu, int elementsMu, double *absorption, double *scattering, std::vector<double> & refractionWavelength, std::vector<double> & refractionReal, std::vector<double> & refractionImag, double *muWavelength, double *muX, double *qmax, std::vector<double> & phi, std::vector<double> & psf,
                     std::vector<std::complex<double> > & s1, std::vector<std::complex<double> > & s2,
                     std::vector<double> & mu, std::vector<double> & piN,
                     std::vector<double> & piNm2,std::vector<double> & piNm1,
                     std::vector<double> & tauN) {

    double ext = 0.0;
    double sca = 0.0;
    double abs = 0.0;
    int nang=MIE_OUT_POINT;
    std::vector<double> phir, psfCumulative;
    //std::vector<double> phi, psf;
    //phi.reserve(nang);
    //psf.reserve(nang);
    phir.reserve(nang);
    psfCumulative.reserve(nang);
    double re = 0.0, im = 0.0;
    int index =0;
    *qmax = 0.0;
    //std::vector<std::complex<double> > s1, s2;
    //s1.reserve(nang);
    //s2.reserve(nang);
    //std::vector<double> mu, piN, piNm2, piNm1, tauN;
    //mu.reserve(nang);
    //piN.reserve(nang);
    //piNm2.reserve(nang);
    //piNm1.reserve(nang);
    //tauN.reserve(nang);
    for (int x = 0; x < elementsX; x++) {
        for (int w = 0; w < elementsWavelength; w++) {
            double wavelength = minWavelength + (maxWavelength - minWavelength)*static_cast<double>(w)/static_cast<double>(elementsWavelength - 1);
            double xValue = pow(10.0, minX + (maxX - minX)*(static_cast<double>(x))/static_cast<double>(elementsX - 1));
            double radius = xValue*wavelength/(2*PI);
            int cIndex = w*elementsX + x;
            muWavelength[w] = wavelength;
            muX[x] = xValue;

            find_v(refractionWavelength, wavelength, &index);
            re = interpolate_v(refractionReal, refractionWavelength, wavelength, index);
            im = interpolate_v(refractionImag, refractionWavelength, wavelength, index);

            mieTyndall(radius,wavelength,re,im,&ext,&sca,nang,phi,psf,s1,s2,mu,piN,piNm2,piNm1,tauN);
            abs = ext - sca;
            for (int i =0; i < nang; i++) phir[i] = phi[i]*PI/180.0;
            double total = 0.0;
            for (int i = 0; i < nang ; i++) total += psf[i]*sin(phir[i]);
            for (int i = 0; i < nang; i++) psfCumulative[i] = psf[i]*sin(phir[i])/total;
            for (int i = 1; i < nang; i++) psfCumulative[i] += psfCumulative[i-1];
            for (int i = 0; i < nang; i++) {
                cumulativeMu[cIndex*elementsMu + i] = psfCumulative[i];
                muValue[i] = cos(phir[i]);
            }
            absorption[cIndex] = abs;
            scattering[cIndex] = sca;
            if (ext > *qmax) *qmax = ext;
        }
    }

            // {
            // FILE *outdata;
            // char filename[4096];
            // sprintf(filename,"mu.txt");
            // outdata = fopen(filename, "w");
            // for (long kk=0; kk<MIE_OUT_POINT;kk++) {
            //     fprintf(outdata,"%lf %lf\n",muValue[kk],phir[kk]);
            // }
            // fclose(outdata);
            // }

}


void Mie::setupMieIndexRefraction (std::string dir, double minwavelength, double maxwavelength, double meanWater, double sigmaWater,double meanSeasalt, double sigmaSeasalt,double meanDust, double sigmaDust,double meanPollution, double sigmaPollution,double meanSmoke, double sigmaSmoke) {

    int elementsWavelength = MIE_IN_POINT;
    int elementsX = MIE_IN_POINT;
    int elementsMu = MIE_OUT_POINT;
    int elementsTau = MIE_TAU_POINT;
    double minDrop = 0.001;

    int nang = MIE_OUT_POINT;
    std::vector<double> phi, psf;
    phi.reserve(nang);
    psf.reserve(nang);
    std::vector<std::complex<double> > s1, s2;
    s1.reserve(nang);
    s2.reserve(nang);
    std::vector<double> mu, piN, piNm2, piNm1, tauN;
    mu.reserve(nang);
    piN.reserve(nang);
    piNm2.reserve(nang);
    piNm1.reserve(nang);
    tauN.reserve(nang);

    //WATER
    dropMeanWater = meanWater;
    dropSigmaWater = sigmaWater;
    dropMinWater = exp(log(dropMeanWater)-7.0*dropSigmaWater);
    dropMaxWater = exp(log(dropMeanWater)+7.0*dropSigmaWater);
    if (dropMinWater < minDrop) dropMinWater = minDrop;
    if (dropMaxWater > 100.0) dropMaxWater = 100.0;
    double maxx = log10(2*PI*dropMaxWater/(minwavelength/1000.0));
    double minx = log10(2*PI*dropMinWater/(maxwavelength/1000.0));
    readText::readCol(dir + "/optical/water.txt",refractionWavelengthWater, refractionRealWater, refractionImagWater);
    muValueWater =new double[elementsMu]();
    cumulativeMuWater = new double[elementsX*elementsWavelength*elementsMu]();
    absorptionWater = new double[elementsX*elementsWavelength]();
    scatteringWater = new double[elementsX*elementsWavelength]();
    totalExtinctionWater = new double[elementsWavelength]();
    totalProfileWater = new double[elementsWavelength*elementsMu]();
    muWavelengthWater = new double[elementsWavelength]();
    muXWater = new double[elementsX]();
    setupMieGrid(elementsWavelength, elementsX, minwavelength/1000.0, maxwavelength/1000.0, minx, maxx, muValueWater, cumulativeMuWater, elementsMu, absorptionWater, scatteringWater, refractionWavelengthWater, refractionRealWater, refractionImagWater, muWavelengthWater, muXWater, &qmaxWater, phi, psf, s1,  s2, mu,  piN, piNm2,piNm1, tauN);
    tauArrayWater = new double[elementsTau]();
    effsWater = new double[elementsTau*elementsWavelength]();
   setupMieScale(dropMinWater, dropMaxWater, dropMeanWater, dropSigmaWater, elementsTau,elementsWavelength, effsWater,tauArrayWater, muWavelengthWater, muXWater, absorptionWater, scatteringWater, totalExtinctionWater, cumulativeMuWater, totalProfileWater);
  
    //SEASALT
    dropMeanSeasalt = meanSeasalt;
    dropSigmaSeasalt = sigmaSeasalt;
    dropMinSeasalt = exp(log(dropMeanSeasalt)-7.0*dropSigmaSeasalt);
    dropMaxSeasalt = exp(log(dropMeanSeasalt)+7.0*dropSigmaSeasalt);
    if (dropMinSeasalt < minDrop) dropMinSeasalt = minDrop;
    if (dropMaxSeasalt > 100.0) dropMaxSeasalt = 100.0;
    maxx = log10(2*PI*dropMaxSeasalt/(minwavelength/1000.0));
    minx = log10(2*PI*dropMinSeasalt/(maxwavelength/1000.0));
    readText::readCol(dir + "/optical/nacl.txt",refractionWavelengthSeasalt, refractionRealSeasalt, refractionImagSeasalt);
    muValueSeasalt =new double[elementsMu]();
    cumulativeMuSeasalt = new double[elementsX*elementsWavelength*elementsMu]();
    absorptionSeasalt = new double[elementsX*elementsWavelength]();
    scatteringSeasalt = new double[elementsX*elementsWavelength]();
    totalExtinctionSeasalt = new double[elementsWavelength]();
    totalProfileSeasalt = new double[elementsWavelength*elementsMu]();
    muWavelengthSeasalt = new double[elementsWavelength]();
    muXSeasalt = new double[elementsX]();
    setupMieGrid(elementsWavelength, elementsX, minwavelength/1000.0, maxwavelength/1000.0, minx, maxx, muValueSeasalt, cumulativeMuSeasalt, elementsMu, absorptionSeasalt, scatteringSeasalt, refractionWavelengthSeasalt, refractionRealSeasalt, refractionImagSeasalt, muWavelengthSeasalt, muXSeasalt, &qmaxSeasalt,  phi, psf, s1,  s2, mu,  piN, piNm2,piNm1, tauN);
    tauArraySeasalt = new double[elementsTau]();
    effsSeasalt = new double[elementsTau*elementsWavelength]();
    setupMieScale(dropMinSeasalt, dropMaxSeasalt, dropMeanSeasalt, dropSigmaSeasalt, elementsTau,elementsWavelength, effsSeasalt,tauArraySeasalt, muWavelengthSeasalt, muXSeasalt, absorptionSeasalt, scatteringSeasalt, totalExtinctionSeasalt, cumulativeMuSeasalt, totalProfileSeasalt);

    //DUST
    dropMeanDust = meanDust;
    dropSigmaDust = sigmaDust;
    dropMinDust = exp(log(dropMeanDust)-7.0*dropSigmaDust);
    dropMaxDust = exp(log(dropMeanDust)+7.0*dropSigmaDust);
    if (dropMinDust < minDrop) dropMinDust = minDrop;
    if (dropMaxDust > 100.0) dropMaxDust = 100.0;
    maxx = log10(2*PI*dropMaxDust/(minwavelength/1000.0));
    minx = log10(2*PI*dropMinDust/(maxwavelength/1000.0));
    readText::readCol(dir + "/optical/si.txt",refractionWavelengthDust, refractionRealDust, refractionImagDust);
    muValueDust =new double[elementsMu]();
    cumulativeMuDust = new double[elementsX*elementsWavelength*elementsMu]();
    absorptionDust = new double[elementsX*elementsWavelength]();
    scatteringDust = new double[elementsX*elementsWavelength]();
    totalExtinctionDust = new double[elementsWavelength]();
    totalProfileDust = new double[elementsWavelength*elementsMu]();
    muWavelengthDust = new double[elementsWavelength]();
    muXDust = new double[elementsX]();
    setupMieGrid(elementsWavelength, elementsX, minwavelength/1000.0, maxwavelength/1000.0, minx, maxx, muValueDust, cumulativeMuDust, elementsMu, absorptionDust, scatteringDust, refractionWavelengthDust, refractionRealDust, refractionImagDust, muWavelengthDust, muXDust, &qmaxDust,  phi, psf, s1,  s2, mu,  piN, piNm2,piNm1, tauN);
    tauArrayDust = new double[elementsTau]();
    effsDust = new double[elementsTau*elementsWavelength]();
    setupMieScale(dropMinDust, dropMaxDust, dropMeanDust, dropSigmaDust, elementsTau,elementsWavelength, effsDust,tauArrayDust,muWavelengthDust, muXDust, absorptionDust, scatteringDust, totalExtinctionDust, cumulativeMuDust, totalProfileDust);

    //POLLUTION
    dropMeanPollution = meanPollution;
    dropSigmaPollution = sigmaPollution;
    dropMinPollution = exp(log(dropMeanPollution)-7.0*dropSigmaPollution);
    dropMaxPollution = exp(log(dropMeanPollution)+7.0*dropSigmaPollution);
    if (dropMinPollution < minDrop) dropMinPollution = minDrop;
    if (dropMaxPollution > 100.0) dropMaxPollution = 100.0;
    maxx = log10(2*PI*dropMaxPollution/(minwavelength/1000.0));
    minx = log10(2*PI*dropMinPollution/(maxwavelength/1000.0));
    readText::readCol(dir + "/optical/ammoniumsulfate.txt",refractionWavelengthPollution, refractionRealPollution, refractionImagPollution);
    muValuePollution =new double[elementsMu]();
    cumulativeMuPollution = new double[elementsX*elementsWavelength*elementsMu]();
    absorptionPollution = new double[elementsX*elementsWavelength]();
    scatteringPollution = new double[elementsX*elementsWavelength]();
    totalExtinctionPollution = new double[elementsWavelength]();
    totalProfilePollution = new double[elementsWavelength*elementsMu]();
    muWavelengthPollution = new double[elementsWavelength]();
    muXPollution = new double[elementsX]();
    setupMieGrid(elementsWavelength, elementsX, minwavelength/1000.0, maxwavelength/1000.0, minx, maxx, muValuePollution, cumulativeMuPollution, elementsMu, absorptionPollution, scatteringPollution, refractionWavelengthPollution, refractionRealPollution, refractionImagPollution, muWavelengthPollution, muXPollution, &qmaxPollution,  phi, psf, s1,  s2, mu,  piN, piNm2,piNm1, tauN);
    tauArrayPollution = new double[elementsTau]();
    effsPollution = new double[elementsTau*elementsWavelength]();
    setupMieScale(dropMinPollution, dropMaxPollution, dropMeanPollution, dropSigmaPollution, elementsTau,elementsWavelength, effsPollution,tauArrayPollution,muWavelengthPollution, muXPollution, absorptionPollution, scatteringPollution, totalExtinctionPollution, cumulativeMuPollution, totalProfilePollution);

    //SMOKE
    dropMeanSmoke = meanSmoke;
    dropSigmaSmoke = sigmaSmoke;
    dropMinSmoke = exp(log(dropMeanSmoke)-7.0*dropSigmaSmoke);
    dropMaxSmoke = exp(log(dropMeanSmoke)+7.0*dropSigmaSmoke);
    if (dropMinSmoke < minDrop) dropMinSmoke = minDrop;
    if (dropMaxSmoke > 100.0) dropMaxSmoke = 100.0;
    maxx = log10(2*PI*dropMaxSmoke/(minwavelength/1000.0));
    minx = log10(2*PI*dropMinSmoke/(maxwavelength/1000.0));
    readText::readCol(dir + "/optical/carbon.txt",refractionWavelengthSmoke, refractionRealSmoke, refractionImagSmoke);
    muValueSmoke =new double[elementsMu]();
    cumulativeMuSmoke = new double[elementsX*elementsWavelength*elementsMu]();
    absorptionSmoke = new double[elementsX*elementsWavelength]();
    scatteringSmoke = new double[elementsX*elementsWavelength]();
    totalExtinctionSmoke = new double[elementsWavelength]();
    totalProfileSmoke = new double[elementsWavelength*elementsMu]();
    muWavelengthSmoke = new double[elementsWavelength]();
    muXSmoke = new double[elementsX]();
    setupMieGrid(elementsWavelength, elementsX, minwavelength/1000.0, maxwavelength/1000.0, minx, maxx, muValueSmoke, cumulativeMuSmoke, elementsMu, absorptionSmoke, scatteringSmoke, refractionWavelengthSmoke, refractionRealSmoke, refractionImagSmoke, muWavelengthSmoke, muXSmoke, &qmaxSmoke,  phi, psf, s1,  s2, mu,  piN, piNm2,piNm1, tauN);
    tauArraySmoke = new double[elementsTau]();
    effsSmoke = new double[elementsTau*elementsWavelength]();
    setupMieScale(dropMinSmoke, dropMaxSmoke, dropMeanSmoke, dropSigmaSmoke, elementsTau, elementsWavelength, effsSmoke,tauArraySmoke, muWavelengthSmoke, muXSmoke, absorptionSmoke, scatteringSmoke, totalExtinctionSmoke, cumulativeMuSmoke, totalProfileSmoke);

}

void Mie::setupMieScale (double dropmin, double dropmax, double dropmean, double dropsigma, int elementsTau, int elementsWavelength, double *effs, double *tauArray, double *muWavelength, double *muX, double *absorption, double *scattering, double *totalExtinction, double *psfCumulative, double *totalProfile) {

    double avgCount = 0.0;
    double minchi = 1e30;
    double bestff = 0.0;
    //    double dropstep = (dropmax - dropmin)/(static_cast<double>(elementsTau));

    std::vector<double> ss(elementsTau);
    std::vector<double> avgExt(elementsWavelength);
    std::vector<double> extF(elementsTau*elementsWavelength);
    for (int i = 0; i < elementsTau; i++) {
        double vv = (-1.0 + static_cast<double>(i)/(static_cast<double>(elementsTau)-1)*2.0);
        double ww = sqrt(PI)*(0.5*vv + 1.0/24.0*PI*pow(vv, 3.0) + 7.0/960.0*pow(PI, 2.0)*pow(vv, 5.0) + 127.0/80640.0*pow(PI, 3.0)*pow(vv, 7.0) + 4369.0/11612160.0*pow(PI,4.0)*pow(vv, 9.0));
        ww=ww*sqrt(2.0);
        ss[i] = exp(log(dropmean) + ww*dropsigma);
        if (ss[i] > dropmax) ss[i] = dropmax;
        if (ss[i] < dropmin) ss[i] = dropmin;
    }
    int closestW = 0;
    double dist = 1e30;
    for (int i = 0 ;  i < elementsWavelength; i++) {
        totalExtinction[i] = 0.0;
        avgExt[i]=0.0;
        avgCount = 0.0;
        for (int j =0 ; j < elementsTau; j++) {
            int xx;
            int wx = i;
            double xq = 2*PI*ss[j]/muWavelength[i];
            find(muX,MIE_IN_POINT, xq, &xx);
            int cindex = wx*MIE_IN_POINT + xx;

            double abs = absorption[cindex];
            double sca = scattering[cindex];
            double ext = abs + sca;
            extF[i*elementsTau + j] = ext;
            totalExtinction[i] += ext*PI*ss[j]*ss[j]/static_cast<double>(elementsTau);
            avgExt[i] += PI*ss[j]*ss[j]/static_cast<double>(elementsTau);
            totalProfile[MIE_OUT_POINT*i + 0] = 0.0;
            for (int k=1; k<MIE_OUT_POINT; k++) {
                totalProfile[MIE_OUT_POINT*i + k] = totalProfile[MIE_OUT_POINT*i + k] + (psfCumulative[cindex*MIE_OUT_POINT + k] -
                                             psfCumulative[cindex*MIE_OUT_POINT + k - 1])*ext*PI*ss[j]*ss[j]*sca/ext;
            }
            avgCount += ext*PI*ss[j]*ss[j];
        }
        if (fabs(muWavelength[i] - 0.5) < dist) {
            dist = fabs(muWavelength[i] - 0.5);
            closestW = i;
        }
        for (int k=0; k<MIE_OUT_POINT; k++) {
            totalProfile[MIE_OUT_POINT*i + k] = totalProfile[MIE_OUT_POINT*i + k]/avgCount;
        }
        // {
        //     FILE *outdata;
        //     char filename[4096];
        //     sprintf(filename,"mie_%ld.txt",i);
        //     outdata = fopen(filename, "w");
        //     for (long kk=0; kk<MIE_OUT_POINT;kk++) {
        //         fprintf(outdata,"%e\n",totalProfile[MIE_OUT_POINT*i + kk]);
        //     }
        //     fclose(outdata);
        // }
    }
    double norm = totalExtinction[closestW];
    for (int i = 0 ;  i < elementsWavelength; i++) {
        totalExtinction[i] = totalExtinction[i] / norm;
    }


    minchi=1e30;
    bestff = 0.0;
    for (int j = 0 ;  j < elementsWavelength; j++) {
        double factor = totalExtinction[j]/avgExt[j]*norm;
        for (int i = 0 ; i < elementsTau; i++) {
            double tau = pow(10.0, -4.0 + 6.0*static_cast<double>(i)/static_cast<double>(elementsTau));
            double dropss=(dropmax*sqrt(factor)-dropmin*sqrt(factor))/(static_cast<double>(elementsTau));
            for (double ff = dropmin*sqrt(factor); ff < dropmax*sqrt(factor); ff += dropss) {
                double rr = 0.0;
                double taul = tau;
                if (tau > 10.0) taul = 10.0;
                //            for (int i = 0; i < elementsTau; i++) rr += (1 - exp(-tau)) - (1 - exp(-tau*ss[i]*ss[i]/ff/ff));
                double factor2 = extF[j*elementsTau + i];
                for (int i = 0; i < elementsTau; i++) rr += (1 - exp(-taul)) - (1 - exp(-taul*ss[i]*ss[i]*factor2/ff/ff));
                rr /= static_cast<double>(elementsTau);
                if (fabs(rr) < minchi) {
                    bestff = ff;
                    minchi  = fabs(rr);
                }
            }
            tauArray[i] = tau;
            effs[j*elementsTau+i] = bestff;
            //if (j==0 && i==0) printf("%ld %ld %lf %lf %lf %lf %lf %lf\n",i,j,tau,bestff,totalExtinction[j]/avgExt[j]*norm,dropmean,dropsigma,bestff/sqrt(factor));
        }
    }

}
