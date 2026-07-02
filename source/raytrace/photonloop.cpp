///
/// @package phosim
/// @file photonloop.cpp
/// @brief main photon loop
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


void* Image::threadFunction(void *voidArgs) {
    thread_args *args = (thread_args*)voidArgs;
    // int state, type;
    // pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,&state);
    // pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,&type);
    for (int i = 0;  i <args->runthread - 1; i++) {
        args->instance->photonLoop(args->ssource[i], args->thread, 0);
    }
    args->instance->photonLoop(args->ssource[args->runthread-1], args->thread, args->runthread);
    // pthread_setcancelstate(state,&state);
    // pthread_setcanceltype(type,&type);
    // pthread_testcancel();
    //pthread_detach(pthread_self());
    pthread_exit(NULL);
    //    return NULL;
}

void Image::photonLoop(int ssource, int thread, int finish) {

    Vector angle, largeAngle, largeAngle2, position;
    //positionDiffraction;
    Vector angleDiffraction, angleDeltaDiffraction, angleDeltaTotalDiffraction;
    Vector positionPrevious, anglePrevious, normal;
    float transmission, reflection;
    double distance;
    double extraLargeAngle;
    double shiftX = 0.0, shiftY = 0.0;
    long ray;
    long photonSource;
    long detRay;
    long backSourceOver;
    long sourceOver;
    long sourceLoop;
    double sourceSaturation;
    double sourceSA;
    double sourceLA;
    double sourceK;
    double sourceSAo;
    double sourceLAo;
    double sourceKo;
    int pmode;
    double satConst = 0.01;
    double backSigma = 0.0;
    int preGhost = 0;
    int newSurf = 0;
    int oldSurf;
    int pflag;
    int waveIndex;
    int waveSurfaceIndex;
    int straylightcurrent;
    int miss;
    int initFlag;
    int initGal;
    Photon photon;
    Clog localLog;
    double rate = 0.1;
    double photonSourceDouble;
    int lock4 = 0;
    int layer = 0;
    int sc = 0;
    int bindex=0;
    int backgroundProbOk=0;
    int nopsf=0;
    double twicevidotvn;

    //    SETUP PHOTONS
    photon.thread = thread;
    counterClear(&localLog);
    photon.maxcounter = 0;
    photonSourceDouble = nphot*sources.norm[ssource]/totalnorm;
    if (nsource == 1) photonSourceDouble = nphot;
    if ((poissonMode == 1) && (sources.type[ssource] >= 9)) {
        photonSource = random[thread].poisson(photonSourceDouble);
    } else {
        photonSource = static_cast<long>(photonSourceDouble);
    }
    if (telconfig != 0 && sources.type[ssource] != 0) photonSource = 0;

    if (opdfile == 1) {
        pthread_mutex_lock(&lock.lock4);
        lock4 = 1;
        for (int surfIdx = 0; surfIdx < nsurf; surfIdx++) {
            surface.innerRadius[surfIdx] = 0.0;
            surface.asphere(surfIdx, SURFACE_POINTS);
        }
    }

    ray = 0;
    initGal = 0;
    sourceOver = 1;
    sourceLA = 0.01;
    sourceSA = 0.01;
    sourceLAo = 0.01;
    sourceSAo = 0.01;
    sourceK = 0.01;
    sourceKo = 0.01;
    photon.sourceOver_m = 1;
    sourceLoop = 1;
    sourceSaturation = 1.0;
    photon.sourceSaturationRadius = -satbuffer;
    photon.sourceSaturationPreviousRadius = -satbuffer;
    photon.sourceSaturationFirst = 0;
    photon.sat_xpos=0;
    photon.sat_ypos=0;
    photon.sat_n=0;
    photon.sat_xx=0;
    photon.sat_yy=0;
    photon.sat_sigma = 0.0;
    photon.sat_maxr = 0.0;
    photon.sat_r = 0.0;
    obstruction.pupil = 0;
    photon.prtime = -1.0;
    if (sources.type[ssource] < 9) {
        rate = 0.05;
        if (sources.type[ssource]==8) rate=0.05;
        if (backGamma < 1.0) backSourceOver = (long)(backAlpha*sqrt(rate*photonSource));
        else backSourceOver = (long)(backAlpha/backGamma*sqrt(rate*photonSource));
        if (backAlpha <= 0.0) backSourceOver = static_cast<long>(backGamma);
        if (backSourceOver < 1) backSourceOver = 1;
        if (backSourceOver > photonSource) backSourceOver = photonSource;
    } else {
        backSourceOver = 1;
    }
    if (sources.type[ssource] >= 9) {
        backSigma = 0.0;
    } else {
        backSigma = backBeta*backRadius/3600.0*platescale;
    }
    if (sources.mag[ssource] < straylightcut && straylight == 1) {
        straylightcurrent = 1;
    } else {
        straylightcurrent = 0;
    }
    //   MAIN PHOTON LOOP
    photonloop:
    while (ray < photonSource) {

        //        if (diagnostic) diagnosticCount(ray,43);
        sourceOver = 1;
        photon.sourceOver_m = 1;
        pflag=0;
        pmode = 0;

        //   Get Wavelength and Time
        miss = getWavelengthTime(&photon, ssource);
        if (diagnostic) diagnosticCount(ray,0);
        photon.absoluteTime = photon.time - exptime/2.0 + timeoffset;


        //   Saturation Optimization
        if (saturation) {
            if (photon.sourceSaturationRadius > (-satbuffer + 2)  && sources.type[ssource] >= 9) {
                float sourceSaturationFrac = (static_cast<float>(ray)+1.0)/(static_cast<float>(photonSource)+1.0);
                if (satindex==1) {
                    if ((photon.sourceSaturationRadius > photon.sourceSaturationPreviousRadius)) {
                            photon.sourceSaturationFraction1 = photon.sourceSaturationFraction2;
                            photon.sourceSaturationRadius1 = photon.sourceSaturationRadius2;
                            photon.sourceSaturationFraction2 = photon.sourceSaturationFraction3;
                            photon.sourceSaturationRadius2 = photon.sourceSaturationRadius3;
                            photon.sourceSaturationFraction3 = photon.sourceSaturationFraction4;
                            photon.sourceSaturationRadius3 = photon.sourceSaturationRadius4;
                            photon.sourceSaturationFraction4 = sourceSaturationFrac;
                            photon.sourceSaturationRadius4 = photon.sourceSaturationRadius + static_cast<float>(satbuffer);
                            photon.sourceSaturationFirst += 1;
                            if (photon.sourceSaturationFirst >= 4) {
                                photon.sourceSaturationIndex = log(photon.sourceSaturationRadius4/photon.sourceSaturationRadius1)/
                                    log(photon.sourceSaturationFraction4/photon.sourceSaturationFraction1);
                                if (photon.sourceSaturationIndex < 0) photon.sourceSaturationIndex = 0.0;
                                if (photon.sourceSaturationIndex > 2) photon.sourceSaturationIndex = 2.0;
                            } else {
                                photon.sourceSaturationIndex = 0.0;
                            }
                            photon.sourceSaturationPreviousRadius = photon.sourceSaturationRadius;
                    }
                } else {
                    photon.sourceSaturationIndex = 0.0;
                }
                photon.sourceSaturationRadiusCorrected = ((static_cast<float>(photon.sourceSaturationRadius)+static_cast<float>(satbuffer))/(pow(sourceSaturationFrac, photon.sourceSaturationIndex)))*0.5 + 0.5*(static_cast<float>(photon.sourceSaturationRadius) + static_cast<float>(satbuffer)) - static_cast<float>(satbuffer);
            }
            if (photon.sourceSaturationRadiusCorrected > 0 && sources.type[ssource] >= 9) {


                sourceLoop = floor(sourceSaturation);
                if (random[thread].uniform() < (sourceSaturation - floor(sourceSaturation))) sourceLoop++;
                if (sourceLoop < 1) sourceLoop=1;

 
                if (random[thread].uniform() < 0.01) {
                    pmode = 0;
                    photon.sourceOver_m = 1;
                    sourceOver = 1;
                } else {
                    pmode = 1;
                    double ratio = (sourceLA/sourceSA)*(sourceSAo/sourceLAo);
                    photon.sourceOver_m = floor( ratio );
                    if (random[thread].uniform() < ratio  -
                        floor( ratio)) {
                        photon.sourceOver_m++;
                    }
                    double ratio2 = sourceKo/(sourceLAo+sourceSAo)*((sourceSA + sourceLA)/sourceK)*
                        (sourceLA/(sourceLA+sourceSA) + sourceSA/(sourceLA+sourceSA)*ratio);
                    //                    if ((thread==0) && ((ray%1000000)==0)) printf("%lf %lf %lf %lf %lf %lf\n",sourceKo/(sourceLAo+sourceSAo),sourceK/(sourceLA+sourceSA), (sourceLA/(sourceLA+sourceSA) + sourceSA/(sourceLA+sourceSA)*ratio),sourceKo/(sourceKo+sourceLAo+sourceSAo), sourceK/(sourceK+sourceLA+sourceSA),ratio2);
                    sourceOver = floor(ratio2);
                    if (random[thread].uniform() <  ratio2 - floor(ratio2)) sourceOver++;
                    if (static_cast<unsigned long>(sourceOver) > well_depth) sourceOver = well_depth;
                    if (sourceOver < 1) sourceOver = 1;

                    if (photon.sourceOver_m < 1) {
                        photon.sourceOver_m = 1;
                        sourceOver = 1;
                    }

                    
                }
                
            }
        }
        if (thread==0) if (((ray%10000000)==0) && (ray!=0)) {
                //printf("%ld %ld %d %lf %lf %lf %lf %lf %ld %ld %d %lf %lf\n",photon.sourceOver_m,sourceOver,photon.sourceSaturationRadius,sourceSaturation,sourceK/(sourceLA+sourceSA+sourceK),sourceKo/(sourceLAo+sourceSAo+sourceKo),sourceLA/(sourceLA+sourceSA),sourceLAo/(sourceLAo+sourceSAo),static_cast<long>((photonSource-ray)/1e6),static_cast<long>(ray/1e6),satbuffer,photon.sourceSaturationRadiusCorrected,photon.sourceSaturationIndex);


                 // FILE *outdafile;
                 // char tempstring[4096];


                 // snprintf(tempstring, sizeof(tempstring), "%s/dyntra_%s.txt", outputdir.c_str(), outputfilename.c_str());
                 // outdafile = fopen(tempstring, "w");
                 // for (int aa=0;aa<state.waveLimit;aa++) {
                 //     for (int bb=0;bb<state.surfaceLimit;bb++) {
                 //         for (int cc=0;cc<4;cc++) {
                 //             float tempa =state.dynamicTransmission[state.surfaceLimit*state.waveLimit*cc + aa*state.surfaceLimit+bb];
                 //             float tempb =state.dynamicTransmissionLow[state.surfaceLimit*state.waveLimit*cc + aa*state.surfaceLimit+bb];
                 //             int tempc =state.dynamicTransmissionInit[state.surfaceLimit*state.waveLimit*cc + aa*state.surfaceLimit+bb];
                 //             fprintf(outdafile, "%d %d %d %lf %lf %d \n", cc,aa,bb,tempa,tempb,tempc);
                 //         }
                 //     }
                 // }      
                 
                 // fclose(outdafile); 
            }

                
           
        if (sources.type[ssource] < 9 && photonSource > 100) nopsf=1;
        if (sources.type[ssource] < 9 && backGamma > 1.0) {
            sourceOver = static_cast<long>(backGamma);
            photon.sourceOver_m = sourceOver;
            if (backSourceOver*sourceOver > photonSource) {
                sourceOver = static_cast<long>(photonSource/backSourceOver);
                photon.sourceOver_m = static_cast<long>(photonSource/backSourceOver);
            }
            if (sourceOver < 1) {
                sourceOver = 1;
                photon.sourceOver_m = 1;
            }
        }
        if (miss) {
            if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
            countBad(&localLog, sourceOver*backSourceOver, &ray);
            goto photonloop;
        }

        // Get Angle
        getAngle(&angle, photon.time, ssource,&photon);
        if (diagnostic) diagnosticCount(ray,1);

        if (eventfile) {
            state.pEventLogging->logPhoton(angle.x, angle.y, photon.wavelength, 0);
            state.pEventLogging->logPhoton(photon.time, stof(sources.id[ssource]), 0, 1);
        }
        waveIndex = static_cast<int>((photon.wavelength*1000 - minwavelength));
        waveSurfaceIndex = waveIndex*state.surfaceLimit;

        initFlag = 0;
        if (throughputfile) addThroughput(&state.throughputLog, minwavelength, maxwavelength, -1, waveIndex, sourceOver*backSourceOver);

        // Dynamic Transmission Optimization
        int kmax;
        if (sources.type[ssource] < 9 && backGamma > 1.0) {
            kmax = 1;
        } else {
            kmax = sourceLoop;
            if (pmode == 0) kmax = 1;
        }
        for (int k = 0; k < kmax; k++) {
 
            if ((k == 0) || (k > 0 && straylightcurrent == 1)) {
                int lastSurface = -1;
                if ((k==0) && (straylightcurrent==0)) {
                    miss = dynamicTransmissionOptimizationZero(&lastSurface, &preGhost, waveSurfaceIndex, straylightcurrent,
                                                               &photon, sources.type[ssource]);
                    if (diagnostic) diagnosticCount(ray,2);

                } else {
                    miss = dynamicTransmissionOptimization(thread,k, &lastSurface, &preGhost, waveSurfaceIndex, straylightcurrent,
                                                           &photon, sources.type[ssource]);
                    if (diagnostic) diagnosticCount(ray,3);
                }
                if (miss == 1) {
                    if (throughputfile && lastSurface >= 0) {
                        for (int j = 0; j <= lastSurface; j++) {
                            addThroughput(&state.throughputLog, minwavelength, maxwavelength, j, waveIndex, sourceOver*backSourceOver);
                        }
                    }
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad_dt(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;

                }
            }
            vectorInit(&largeAngle);
            vectorInit(&largeAngle2);
            extraLargeAngle = 0.0;
            largeAngle.z = 1.0;
            shiftX = 0.0;
            shiftY = 0.0;

 
            miss = samplePupil(&position, ray, &photon);
            if (diagnostic) diagnosticCount(ray,4);
            if (miss) {
                if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                countBad(&localLog, sourceOver*backSourceOver, &ray);
                goto photonloop;
            }
            if (finiteDistance != 0.0) {
                double finiteDistanceR = sqrt(pow(finiteDistance, 2) + pow(position.x, 2) +
                                              pow(position.y, 2));
                angle.x += position.x/finiteDistanceR;
                angle.y += position.y/finiteDistanceR;
                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
            }
            if (opdfile && ray == 0) {
                angle.x = 0.0;
                angle.y = 1e-5;
                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
            }
            //            if (diffractionMode == 5) vectorCopy(position, &positionDiffraction);
            photon.xp = position.x;
            photon.yp = position.y;
            vectorInit(&largeAngle);
            if (initFlag == 0) {
                photon.shiftedAngle = spiderangle + photon.time*rotationrate*ARCSEC;
                photon.wavelengthFactor = pow(photon.wavelength, -0.2)/screen.wavelengthfactor_nom;
                initFlag = 1;
            }

            //  Diffraction
            if (diffractionMode >= 1 && spiderMode == 1) {
                miss = diffraction(&position, angle, &largeAngle, &photon);
                if (diagnostic) diagnosticCount(ray,5);
                if (miss) {
                    //                    if (k > 0) goto redodiff;
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }
            }
   
            //   Large Angle Scattering
            if ((nopsf==0) && (microroughness==1)) {
                for (int kq=0; kq<nsurf; kq++) {
                    if (surface.surfacetype[kq] == MIRROR) {
                        largeAngleScattering(&largeAngle, &photon, kq);
                    }
                }
            }
            //            printf("AB %e %e\n",largeAngle.x,largeAngle.y);
            if (diagnostic) diagnosticCount(ray,6);

            //   Second Kick
            if (nopsf==0 && diffractionMode == 1 && sources.type[ssource] !=  0) {
                secondKick(&largeAngle, &photon);
                //                printf("AC %e %e\n",largeAngle.x,largeAngle.y);
                if (diagnostic) diagnosticCount(ray,40);
            }

            //    Get Delta Angle
            miss = getDeltaAngle(&angle, &position, ssource, &shiftX, &shiftY, thread, &initGal, &photon);
            if (diagnostic) diagnosticCount(ray,7);
            if (miss) {
                if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                countBad(&localLog, sourceOver*backSourceOver, &ray);
                goto photonloop;
            }

            // Saturation Optimization Calculation
            if (photon.sourceSaturationRadiusCorrected > 0) {
                
                if ((sqrt(pow(largeAngle.x + shiftX, 2.0) + pow(largeAngle.y + shiftY, 2.0))/DEGREE*platescale/pixsize
                     > photon.sourceSaturationRadiusCorrected)  || (preGhost >= 2)) {
                    photon.sourceOver_m = 1;
                    //sourceSaturation -= ((1.0 - np)/np + 1.0)*satConst;
                    pflag=1;
                    //                                        sourceSaturation -= (exp(-np)/(1-exp(-np)) + 1.0)*satConst;
                    //                    if (sourceSaturation < 1) sourceSaturation = 1;
                    //                    if (thread==0) printf("A1 %lf\n",sourceSaturation);
                    break;
                }
            }
        }
        //if (photon.sourceSaturationRadius > 0) {
            //sourceSaturation += satConst;
            //            if (thread==0) printf("B %lf\n",sourceSaturation);
        //}
        if (np <= 0.0) sourceSaturation = 1;
        photon.counter = -1;
        if (opdfile) photon.op = 0.0;

        angle.x += largeAngle.x + shiftX;
        angle.y += largeAngle.y + shiftY;
        angle.z = smallAnglePupilNormalize(angle.x, angle.y);


        //BP
        if (backgroundProbability) {
            if (sources.type[ssource] < 9) {
                backgroundProbOk=1;
                double bpdx = (angle.x - bpVxMin)/bpStep;
                double bpdy = (angle.y - bpVyMin)/bpStep;
                int bpx = round(bpdx);
                int bpy = round(bpdy);
                if ((bpx > bpNx-1) || (bpy > bpNy-1) || (bpx < 0) || (bpy < 0)) {
                    backgroundProbOk=0;
                    goto jumpBP3;
                }
                bindex = bpx*bpNy + bpy;
                if ((state.bpRej[bindex] >= bpIndex) && (state.bpRej[bindex] >= bpIndex*state.bpAcc[bindex])) {
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad_dt(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }
            }
        }
    jumpBP3:;
     // ATMOSPHERE
        //int maxw=static_cast<int>((maxwavelength - minwavelength)*10);
        // for (int j=0; j< surfaceLimit;j++) {
        //     int d=0;
        //     for (int i =0; i < maxw; i++) {
        //         int q=state.dynamicTransmissionInit[i*surfaceLimit+j];
        //         if (q>0) d++;
        //     }
        //     printf("R %d %d\n",j,d);
        // }




        photon.airRefraction = 0.0;
        photon.airRefractionPrevious = 0.0;
        if (atmosphericdispcenter) {
            photon.airRefractionADC = 0.0;
            photon.airRefractionPreviousADC = 0.0;
        }
        photon.ncurr = 1.0;
        // Loop over Atmosphere Layers
        vectorInit(&angleDeltaTotalDiffraction);
        sc = 0;
        for (int l = -1; l < natmospherefile; l++) {
            if (spaceMode >= 0) {
                layer = l;
            } else {
                layer = natmospherefile - 2 - l;
            }

            // Atmospheric Dispersion
            photon.airRefractionPrevious = photon.airRefraction;
            if (atmosphericdispcenter) photon.airRefractionPreviousADC = photon.airRefractionADC;
            if (atmosphericdispcenter) photon.airRefraction = airIndexRefraction(&photon, layer, 0, &(photon.airRefractionADC));
            else photon.airRefraction = airIndexRefraction(&photon, layer, 1, &(photon.airRefractionADC));
            if (diagnostic) diagnosticCount(ray,39);
            photon.ncurr = 1.0 + photon.airRefraction;
            if (sources.type[ssource] !=  0) {
                vectorCopy(angle, &angleDiffraction);
                atmosphericDispersion(&angleDiffraction, &photon, layer);
                if (diagnostic) diagnosticCount(ray,38);
                vectorSubtract(angleDiffraction, angle, &angleDeltaDiffraction);
                vectorAdd(angleDeltaTotalDiffraction, angleDeltaDiffraction, &angleDeltaTotalDiffraction);
            }

             // Atmosphere Propagate
            atmospherePropagate(&position, angle, layer, diffractionMode, &photon);
            if (diagnostic) diagnosticCount(ray,8);
            if (sources.type[ssource] != 0) {
                if (atmospheremode == 2) {
                    if (layer >= 0) {

                        // Atmosphere Intercept
                        atmosphereIntercept(&position, layer, &photon);
                        if (diagnostic) diagnosticCount(ray,9);

                        // Atmosphere Refraction
                        if (nopsf==0) {
                            atmosphereRefraction(&angle, layer, diffractionMode, &photon);
                            if (diagnostic) diagnosticCount(ray,10);
                        }

                        // Clouds
                        if ((sources.type[ssource] == 4) && (sc==0)) {
                            transmission = cloudOpacity(layer, &photon);
                            if (diagnostic) diagnosticCount(ray,11);
                            if (transmissionCheck(transmission, 1+ (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                if (rejectScatterCloud(sources.backgroundAngle[ssource], &photon) == 0) {
                                    transmission = atmosphereOpacityMoon(angle, layer, &photon);
                                    if (diagnostic) diagnosticCount(ray,12);
                                    if (transmissionCheck(transmission, 4 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                        // nothing
                                    } else {
                                        transmission = cloudOpacityMoon(layer, &photon);
                                        if (diagnostic) diagnosticCount(ray,13);
                                        if (transmissionCheck(transmission, 3 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                            // nothing
                                        } else {
                                            sc += 1;
                                        }
                                    }
                                }
                            }
                        }

                        if ((sources.type[ssource] == 7) && (sc==0)) {
                            transmission = cloudOpacity(layer, &photon);
                            if (diagnostic) diagnosticCount(ray,11);
                            if (solarfactor != 1.0) transmission = 1.0 - (1.0 - transmission)*solarfactor;
                            if (transmissionCheck(transmission, 1 + (layer + 1)*6, waveSurfaceIndex, &photon))  {
                                if (rejectScatterCloud(sources.backgroundAngle[ssource], &photon) == 0) {
                                    transmission = atmosphereOpacitySun(angle, layer, &photon);
                                    if (diagnostic) diagnosticCount(ray,14);
                                    if (transmissionCheck(transmission, 6 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                        //nothing
                                    } else {
                                        transmission = cloudOpacitySun(layer, &photon);
                                        if (diagnostic) diagnosticCount(ray,15);
                                        if (transmissionCheck(transmission, 5 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                            //nothing
                                        } else {
                                            sc += 1;
                                        }
                                    }
                                }
                            }
                        }

                        if (((sources.type[ssource] != 4) && (sources.type[ssource] != 7)) || (sc > 0)) {
                            transmission = cloudOpacity(layer, &photon);
                            if (diagnostic) diagnosticCount(ray,11);
                            if (transmissionCheck(transmission, 1 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                int scattered = 1;
                                if ((sources.type[ssource] == 4) || (sources.type[ssource] ==7)) {

                                } else {
                                  int actScatter =1;
                                  if ((sources.type[ssource] == 0) || (sources.type[ssource] == 1) ||
                                      (sources.type[ssource] == 2) || (sources.type[ssource] == 5) || (sources.type[ssource] == 8)) actScatter = 0;
                                  scattered = mieScatter(-log(transmission),&largeAngle2,&photon,actScatter);
                                  if (diagnostic) diagnosticCount(ray,16);
                                }
                                if (scattered!=0) {
                                    if (sources.type[ssource] == 8) {
                                        sc++;
                                    } else {
                                        if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }
                                }
                            }
                        }
                    }

                       //Atmosphere Opacity
                        if ((sources.type[ssource]==4) && (sc==0)) {
                            transmission = atmosphereOpacity(angle, layer, &photon,waveSurfaceIndex, 2 + (layer + 1)*6);;
                            if (diagnostic) diagnosticCount(ray,17);
                            if (transmissionCheck(transmission, 2 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                if (rejectScatterAtmosphere(sources.backgroundAngle[ssource], &photon) == 0) {
                                    transmission = atmosphereOpacityMoon(angle, layer, &photon);
                                    if (diagnostic) diagnosticCount(ray,12);
                                    if (transmissionCheck(transmission, 4 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                        // nothing
                                    } else {
                                        transmission = cloudOpacityMoon(layer, &photon);
                                        if (diagnostic) diagnosticCount(ray,13);
                                        if (transmissionCheck(transmission, 3 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                            // nothing
                                        } else {
                                            sc += 1;
                                        }
                                    }
                                }
                            }
                        }

                        if ((sources.type[ssource] == 7) && (sc==0)) {
                            transmission = atmosphereOpacity(angle, layer, &photon, waveSurfaceIndex, 2 + (layer + 1)*6);
                            if (diagnostic) diagnosticCount(ray,17);
                            if (solarfactor != 1.0) transmission = 1.0 - (1.0 - transmission)*solarfactor;
                                if (transmissionCheck(transmission, 2 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                    if (rejectScatterAtmosphere(sources.backgroundAngle[ssource], &photon) == 0) {
                                        transmission = atmosphereOpacitySun(angle, layer, &photon);
                                        if (diagnostic) diagnosticCount(ray,14);
                                        if (transmissionCheck(transmission, 6 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                            //nothing
                                        } else {
                                            transmission = cloudOpacitySun(layer, &photon);
                                            if (diagnostic) diagnosticCount(ray,15);
                                            if (transmissionCheck(transmission, 5 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                                //nothing
                                            } else {
                                                sc += 1;
                                            }
                                        }
                                    }
                                }
                        }

                // Atmosphere Opacity
                    if (((sources.type[ssource] != 4) && (sources.type[ssource] != 7)) || (sc > 0)) {
                        transmission = atmosphereOpacity(angle, layer, &photon, waveSurfaceIndex, 2 + (layer + 1)*6);
                        if (diagnostic) diagnosticCount(ray,17);
                        if ((sources.type[ssource] == 0) || (sources.type[ssource] == 1) ||
                            (sources.type[ssource] == 2) || (sources.type[ssource] == 5) || (sources.type[ssource] == 8)) transmission = tauSuppression(transmission);
                            if (transmissionCheck(transmission, 2 + (layer + 1)*6, waveSurfaceIndex, &photon)) {
                                if (sources.type[ssource] == 8) {
                                    sc++;
                                } else {
                                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }
                    }
                }
            }

            if (eventfile) {
                state.pEventLogging->logPhoton(position.x, position.y, position.z, layer + 100);
            }

        }

        if (sources.type[ssource] == 3 || sources.type[ssource] == 4 || sources.type[ssource] == 6 ||
            sources.type[ssource] == 7 || sources.type[ssource]==8) {
            if (sc != 1) {
                if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                countBad(&localLog, sourceOver*backSourceOver, &ray);
                goto photonloop;
            }
        }

        vectorAdd(angle, angleDeltaTotalDiffraction, &angle);
        vectorAdd(angle, largeAngle2, &angle);
        extraLargeAngle = sqrt(largeAngle2.x*largeAngle2.x + largeAngle2.y*largeAngle2.y);



        // Atmosphere Diffraction
        if (diffractionMode == 2 && sources.type[ssource] != 0) {
            atmosphereDiffraction(&angle, &photon);
            if (diagnostic) diagnosticCount(ray,18);
        }

        // Dome Seeing
        if ((domeseeing > 0.0 || toypsf > 0.0) && (nopsf==0)) {
            domeSeeing(&angle, &photon);
            if (diagnostic) diagnosticCount(ray,19);
        }

        // Tracking
        if (trackingMode && (nopsf==0)) {
            tracking(&angle, photon.absoluteTime);
            if (diagnostic) diagnosticCount(ray,20);
        }

        if (telescopeMode == 0) {
            newSurf = nsurf - 2;
            double p0 = platescale/DEGREE/fabs(angle.z);
            angle.x -= rotatex;
            angle.y -= rotatey;
            normalize(&angle);
            if (nmirror % 2 == 0) {
                propagate(&position, angle, ((surface.height[nsurf - 1] + p0) - position.z)/angle.z);
                if (diagnostic) diagnosticCount(ray,21);
                angle.x -= (position.x/p0);
                angle.y -= (position.y/p0);
                angle.z = -1.0;
                normalize(&angle);
            } else {
                propagate(&position, angle, ((surface.height[nsurf - 1] - p0) - position.z)/angle.z);
                if (diagnostic) diagnosticCount(ray,21);
                angle.x -= (position.x/p0);
                angle.y -= (position.y/p0);
                angle.z = 1.0;
                normalize(&angle);
            }
            if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, 200);
        } else {
            newSurf = -1;
        }

        // OPTICS AND DETECTOR
        photon.direction = 1;
        photon.ghostFlag = 0;
        photon.saturationFlag = 0;
    surfaceloop: while (1) {
            oldSurf = newSurf;
            if (photon.direction == 1) {
                newSurf++;
            } else {
                newSurf--;
            }
        redostraylight:;

            // Find intercept
            if (newSurf>= 0 && newSurf<nsurf) {
                transform(&position, &angle, newSurf, 0, &photon);
                if (diagnostic) diagnosticCount(ray,22);
                if (surface.surfacetype[newSurf] == DETECTOR) {
                    transform(&position, &angle, newSurf + 1, focusFlag, &photon);
                    if (diagnostic) diagnosticCount(ray,22);
                }
                miss = findSurface(angle, position, &distance, newSurf, &photon);
                if (diagnostic) diagnosticCount(ray,23);
            } else {
                miss = 1;
            }

            //   Missed surface or wrong direction
            if (miss == 1 || (distance < 0)) {
                //                photon.ghostFlag += 1;
                if (straylightcurrent == 0) {
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;

                } else {
                    if (chooseSurface(&newSurf, &oldSurf, &photon) == 1) {
                        if (diagnostic) diagnosticCount(ray,24);
                        if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                        goto photonloop;

                    } else {
                        goto redostraylight;
                    }
                }
            }
            if (onlyadjacent == 1) {
                if ((newSurf > (oldSurf + 1)) || (newSurf < (oldSurf - 1))) {
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;

                }
            }
            propagate(&position, angle, distance);
            if (diagnostic) diagnosticCount(ray,21);
            if (opdfile) {

                if (newSurf == 0) {
                    photon.opdx = -position.z*angle.x/angle.z + position.x;
                    photon.opdy = -position.z*angle.y/angle.z + position.y;
                    distance = position.x*angle.x + position.y*angle.y + (position.z - 20000)*angle.z;
                }
                photon.op -= distance*photon.ncurr;
            }

            if (throughputfile == 1 && photon.direction == 1) {
                addThroughput(&state.throughputLog, minwavelength, maxwavelength, newSurf, waveIndex, sourceOver*backSourceOver);
            }

            //   DERIVATIVES
            interceptDerivatives(&normal, position, newSurf, &photon);
            if (diagnostic) diagnosticCount(ray,25);

            //   CONTAMINATION
            if (surface.surfacetype[newSurf] !=  DETECTOR && contaminationmode == 1) {
                miss = contaminationSurfaceCheck(position, &angle, newSurf, &photon);
                if (diagnostic) diagnosticCount(ray,26);
                if (miss) {
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;
                }

            }


            //   SURFACE COATINGS

            if ((surface.surfacetype[newSurf] == DETECTOR) && (coating.nonUniformity[newSurf]==1)) {
                if (surface.surfacecoating[newSurf]!=0 && coatingmode ==1) {
                        pixelPosition(&position, &photon);
                        if (diagnostic) diagnosticCount(ray,32);
                        if (photon.pos.x > maxx) photon.pos.x = maxx;
                        if (photon.pos.x < minx) photon.pos.x = minx;
                        if (photon.pos.y > maxy) photon.pos.y = maxy;
                        if (photon.pos.y < miny) photon.pos.y = miny;
                        int location = (floor((chip.nampx - 1)/DET_RESAMPLE)+1)*
                            floor((photon.pos.y - miny)/DET_RESAMPLE) +
                            floor((photon.pos.x - minx)/DET_RESAMPLE);
                        transmission = surfaceCoating(photon.wavelength/coating.thickness[newSurf][location]*coating.nominalThickness[newSurf], angle, normal, newSurf, &reflection, &photon);
                } else {
                    transmission = 1.0;
                }
            } else if ((surface.surfacetype[newSurf] != DETECTOR) && (coating.nonUniformity[newSurf]==1) ) {
                if (surface.surfacecoating[newSurf]!=0 && coatingmode ==1) {
#ifdef  PERTURBATION_COORDINATE
                    int eeint = photon.vvint;
#else
                    int eeint = photon.wwint;
#endif
                    int location = eeint*PERTURBATION_POINTS + photon.uuint;
                    transmission = surfaceCoating(photon.wavelength/coating.thickness[newSurf][location]*coating.nominalThickness[newSurf], angle, normal, newSurf, &reflection, &photon);
                } else {
                    transmission = 1.0;
                }
            } else {
                transmission = surfaceCoating(photon.wavelength, angle, normal, newSurf, &reflection, &photon);
                if (diagnostic) diagnosticCount(ray,27);
            }

            if (transmissionCheck(transmission, (natmospherefile + 1)*6 + 2 + newSurf*2, waveSurfaceIndex, &photon)) {
                if (straylightcurrent == 1 && ghost[newSurf] == 0) {
                    if (transmissionCheckNoRoll(reflection + transmission, (natmospherefile + 1)*6 + 2 + newSurf*2 + 1, waveSurfaceIndex, &photon)) {
                        if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                        countBad(&localLog, sourceOver*backSourceOver, &ray);
                        goto photonloop;

                    } else {
                        photon.direction = -photon.direction;
                        photon.ghostFlag += 1;
                        reflect(&angle, normal,&twicevidotvn);
                        if (diagnostic) diagnosticCount(ray,29);
                        transformInverse(&position, &angle, newSurf);
                        if (diagnostic) diagnosticCount(ray,28);
                        if (surface.surfacetype[newSurf] == DETECTOR) {
                            transformInverse(&position, &angle, newSurf + 1);
                            if (diagnostic) diagnosticCount(ray,28);
                        }
                        if (surface.surfacetype[newSurf] == MIRROR) {
                            if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                            countBad(&localLog, sourceOver*backSourceOver, &ray);
                            goto photonloop;

                        }
                        goto surfaceloop;
                    }
                } else {
                    if (pmode==1) sourceK+=sourceOver*backSourceOver; else sourceKo+=sourceOver*backSourceOver;
                    countBad(&localLog, sourceOver*backSourceOver, &ray);
                    goto photonloop;

                }
            }

            //   INTERACTIONS
            if (surface.surfacetype[newSurf] == MIRROR) {

                //   MIRROR
                reflect(&angle, normal,&twicevidotvn);
                surface.costotal[newSurf] += fabs(0.5*twicevidotvn);
                surface.coscount[newSurf] += 1.0;
                if (diagnostic) diagnosticCount(ray,29);
                transformInverse(&position, &angle, newSurf);
                if (diagnostic) diagnosticCount(ray,28);
                if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

            } else if (surface.surfacetype[newSurf] == LENS || surface.surfacetype[newSurf] == FILTER) {

                 //   LENS/FILTER
                newRefractionIndex(newSurf, &photon);
                if (diagnostic) diagnosticCount(ray,30);
                refract(&angle, normal, photon.nprev, photon.ncurr);
                if (diagnostic) diagnosticCount(ray,31);
                transformInverse(&position, &angle, newSurf);
                if (diagnostic) diagnosticCount(ray,28);
                if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

            } else if (surface.surfacetype[newSurf] == GRATING) {

                //   GRATING
                double wavelengthNm = photon.wavelength*1000.0;
                Vector angleOut;
                surface.ppGrating[newSurf]->diffract(
                                              angle.x, angle.y, angle.z, 
                                              angleOut.x, angleOut.y, 
                                              angleOut.z, wavelengthNm);
                vectorCopy(angleOut, &angle);
                transformInverse(&position, &angle, newSurf);
                if (eventfile) state.pEventLogging->logPhoton(position.x, 
                                                              position.y, 
                                                              position.z, 
                                                              newSurf + 200);
           } else if (surface.surfacetype[newSurf] == DETECTOR) {


                
                // check for widely scattered perturbations
                if (photon.sourceSaturationRadiusCorrected > 0) {
                    if (((sqrt(pow(largeAngle.x + largeAngle2.x + shiftX, 2.0) + pow(largeAngle.y + largeAngle2.y + shiftY, 2.0))/DEGREE*platescale/pixsize
                          > photon.sourceSaturationRadiusCorrected) && (photon.sourceOver_m > 1)) || (photon.ghostFlag>=2)) {
                            photon.sourceOver_m = 1;
                            pflag=1;
                    }
                }
                
                if (eventfile || opdfile) {
                    transformInverse(&position, &angle, newSurf + 1);
                    if (diagnostic) diagnosticCount(ray,28);
                    transformInverse(&position, &angle, newSurf);
                    if (diagnostic) diagnosticCount(ray,28);
                    if (eventfile) state.pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                    if (opdfile) {
                        if (ray == 0) {
                            state.epR[ssource] = -position.y*angle.z/angle.y;
                        } else if (ray == 1) {
                            state.cx[ssource] = position.x;
                            state.cy[ssource] = position.y;
                            state.cz[ssource] = position.z;
                            state.r0[ssource] = sqrt(pow(state.epR[ssource], 2) + pow(state.cx[ssource], 2) + pow(state.cy[ssource], 2));
                            pthread_mutex_unlock(&lock.lock4);
                            lock4 = 0;
                            pthread_mutex_lock(&lock.lock6);
                            remain--;
                            if (remain == 0) {
                                pthread_cond_broadcast(&lock.cond);
                            } else {
                                while (remain != 0) {
                                    pthread_cond_wait(&lock.cond, &lock.lock6);
                                }
                            }
                            pthread_mutex_unlock(&lock.lock6);
                            pthread_mutex_lock(&lock.lock7);
                            for (int surfIdx = 0; surfIdx < nsurf; surfIdx++) {
                                surface.innerRadius[surfIdx] = surface.innerRadius0[surfIdx];
                                surface.asphere(surfIdx, SURFACE_POINTS);
                            }
                            pthread_mutex_unlock(&lock.lock7);
                        }
                        //solve line-sphere intersection analytically
                        double ocx = position.x - state.cx[ssource];
                        double ocy = position.y - state.cy[ssource];
                        double ocz = position.z - state.cz[ssource];
                        double ocsqr = ocx*ocx + ocy*ocy + ocz*ocz;
                        double ocproj = ocx*angle.x + ocy*angle.y + ocz*angle.z;
                        distance = - ocproj + sqrt(ocproj*ocproj - ocsqr + state.r0[ssource]*state.r0[ssource]);
                        if (ray == 0) distance = state.epR[ssource];
                        photon.op -= distance;
                    }
                    transform(&position, &angle, newSurf, 0, &photon);
                    transform(&position, &angle, newSurf + 1, 0, &photon);
                }
                vectorCopy(position, &positionPrevious);
                vectorCopy(angle, &anglePrevious);

                int electrons = 1;
                if (photon.wavelength < silicon.maxWavelengthElectronPair) {
                    double meanelectrons = H_PLANCK*C_LIGHT/(photon.wavelength*1e-6)/E_CHARGE/silicon.pairenergy;
                    double ee= meanelectrons + random[thread].normal()*sqrt(meanelectrons*silicon.fano);
                    if (ee < 1) ee=1;
                    electrons = static_cast<int>(ee);
                }
                for (int e=0; e < electrons; e++) {

                    if (e > 0) {
                        position.x = positionPrevious.x;
                        position.y = positionPrevious.y;
                        position.z = positionPrevious.z;
                    }

                    detRay = 0;

            detectorloop: while (detRay < backSourceOver) {

                    if (diagnostic) diagnosticCount(ray,43);
                    if (detRay > 0) {
                        if (diagnostic) diagnosticCount(ray,44);
                        gaussianXY(&position,positionPrevious,backSigma,&photon);
                        vectorCopy(anglePrevious, &angle);
                    }

                    //   SENSOR
                    if (detectorMode) {

  

                        // Pixel position
                        pixelPosition(&position, &photon);
                        if (diagnostic) diagnosticCount(ray,32);
                        if (sources.type[ssource] < 9) {
                            if (diagnostic) diagnosticCount(ray,46);
                            if (!inBoundsBuffer(&(photon.pos),activeBuffer)) {
                                if (backgroundProbability) backgroundProbCount(bindex,backgroundProbOk,ssource);
                                if (pmode==1) sourceK+=sourceOver; else sourceKo+=sourceOver;
                                countBad(&localLog, sourceOver, &ray);
                                detRay++;
                                goto detectorloop;
                            }
                        }

                        // Relative Surface Height
                        double dh=0;
                        miss = getDeltaIntercept(position.x, position.y, &dh, newSurf, &photon);
                        if (diagnostic) diagnosticCount(ray,33);

                        // Detector Contamination
                        if (contaminationmode) {
                            miss = contaminationDetectorCheck(position, &angle, &photon);
                            if (diagnostic) diagnosticCount(ray,41);
                            if (miss) {
                                if (pmode==1) sourceK+=sourceOver; else sourceKo+=sourceOver;
                                countBad(&localLog, sourceOver, &ray);
                                detRay++;
                                goto detectorloop;
                            }
                        }

                        //Photon Silicon Propagate
                        miss = photonSiliconPropagate(&angle, &position, photon.wavelength, normal, dh, waveSurfaceIndex, &photon);
                        if (diagnostic) diagnosticCount(ray,34);
                        if (eventfile) {
                                transformInverse(&position, &angle, newSurf + 1);
                                transformInverse(&position, &angle, newSurf);
                                if (miss == 0) state.pEventLogging->logPhoton(position.x, position.y, position.z, 300);
                                if (miss != 0) state.pEventLogging->logPhoton(position.x, position.y, position.z, 303);
                                transform(&position, &angle, newSurf, 0, &photon);
                                transform(&position, &angle, newSurf + 1, 0, &photon);
                        }
                        if (miss) {
                            if (pmode==1) sourceK+=sourceOver; else sourceKo+=sourceOver;
                            countBad(&localLog, sourceOver, &ray);
                            detRay++;
                            goto detectorloop;
                        }

                        // Electron Silicon Propagate
                        int eCounter = 0;
                        do {

                            // pixelPosition(&position, &photon);
                            // if (diagnostic) diagnosticCount(ray,32);

                            miss = electronSiliconPropagate(&angle, &position, &photon, eCounter, nopsf);
                            if (diagnostic) diagnosticCount(ray,35);

                            eCounter++;
                            if (miss) eCounter=siliconSegments;
                            if (photon.zindex >= (SILICON_STEPS-siliconSubSteps)) eCounter=siliconSegments;

                            if (eventfile) {
                                position.z += sensorthickness/siliconSegments*photon.sensorDirection;
                                if (position.z*photon.sensorDirection >= photon.collect_z*photon.sensorDirection) {
                                    position.z = photon.collect_z;
                                }
                                if (eCounter < siliconSegments) {
                                    transformInverse(&position, &angle, newSurf + 1);
                                    transformInverse(&position, &angle, newSurf);
                                    state.pEventLogging->logPhoton(position.x, position.y, position.z, 302);
                                    transform(&position, &angle, newSurf, 0, &photon);
                                    transform(&position, &angle, newSurf + 1, 0, &photon);
                                } else {
                                    transformInverse(&position, &angle, newSurf + 1);
                                    transformInverse(&position, &angle, newSurf);
                                    state.pEventLogging->logPhoton(position.x, position.y, position.z, 301);
                                    transform(&position, &angle, newSurf, 0, &photon);
                                    transform(&position, &angle, newSurf + 1, 0, &photon);
                                }
                            }
                        } while (eCounter < siliconSegments);
                    }

                    pixelPositionDrift(&position, &photon);
                    if (diagnostic) diagnosticCount(ray,32);
                    if (eventfile) {
                        if (inBoundsDrift(&(photon.pos))) {
                            state.pEventLogging->logPhoton(static_cast<double>(photon.pos.x),static_cast<double>(photon.pos.y), 0.0, 304);
                        }
                    }
                    if (centroidfile) {
                        if (inBoundsDrift(&(photon.pos))) {
                            sourceXpos[ssource] += photon.pos.x*sourceOver;
                            sourceYpos[ssource] += photon.pos.y*sourceOver;
                            sourcePhoton[ssource] += sourceOver;
                        }
                    }
                    if (opdfile) {
                        if (ray > 0) {
                        int xx = floor(photon.opdx/maxr/2*(OPD_SCREEN_SIZE - 1.0) + OPD_SCREEN_SIZE/2.0);
                        int yy = floor(photon.opdy/maxr/2*(OPD_SCREEN_SIZE - 1.0) + OPD_SCREEN_SIZE/2.0);
                            if (xx >= 0 && xx < OPD_SCREEN_SIZE && yy >= 0 && yy < OPD_SCREEN_SIZE) {
                                pthread_mutex_lock(&lock.lock5);
                                *(state.opd + ssource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE + OPD_SCREEN_SIZE*yy + xx) += photon.op;
                                *(state.opdcount + ssource*OPD_SCREEN_SIZE*OPD_SCREEN_SIZE + OPD_SCREEN_SIZE*yy + xx) += 1;
                                pthread_mutex_unlock(&lock.lock5);
                            }
                        }
                    }

                    if (sources.type[ssource] < 9 && backGamma > 1.0) {
                        Vector newPosition;
                        vectorCopy(position, &newPosition);
                        int lmax;
                        lmax = photon.sourceOver_m;
                        photon.sourceOver_m = 1;
                        for (int l = 0; l < lmax; l++) {
                            if (l > 0) {
                                if (diagnostic) diagnosticCount(ray,44);
                                gaussianXY(&position,newPosition,backSigma/backDelta,&photon);
                            }
                            pixelPositionDrift(&position,&photon);
                            if (diagnostic) diagnosticCount(ray,32);
                            if (diagnostic) diagnosticCount(ray,45);
                            if (inBoundsDrift(&(photon.pos))) {
                                if (saturation) {
                                    saturateIn(ssource, &largeAngle, &photon, shiftX, shiftY, extraLargeAngle, &largeAngle2);
                                    if (diagnostic) diagnosticCount(ray,37);
                                } else {
                                    *(state.focalPlane + chip.extraX*(photon.pos.y - minyDrift) +
                                      (photon.pos.x - minxDrift)) += photon.sourceOver_m;
                                }
                                if (backgroundProbability) {
                                    if ((sources.type[ssource] < 9) && (backgroundProbOk)) {
                                        state.bpAcc[bindex] += 1;
                                    }
                                }
                                countGood(&localLog, photon.sourceOver_m, &ray);
                            } else {
                                if (backgroundProbability)  backgroundProbCount(bindex,backgroundProbOk,ssource);
                                if (pmode==1) sourceK+=photon.sourceOver_m; else sourceKo+=photon.sourceOver_m;
                                countBad(&localLog, photon.sourceOver_m, &ray);
                            }
                        }

                        photon.sourceOver_m = lmax;
                    } else {
                        if (diagnostic) diagnosticCount(ray,45);
                        if (inBoundsDrift(&(photon.pos))) {
                            if (backgroundProbability) {
                                if ((sources.type[ssource] < 9) && (backgroundProbOk)) {
                                    state.bpAcc[bindex] += 1;
                                }
                            }
                            if (saturation) {
                                 if (photon.sourceSaturationRadiusCorrected > 0) {
                                    if (sourceLA/(sourceLA+sourceSA) > (1-exp(-np))) {
                                        sourceSaturation -= satConst*satConst;
                                    } else {
                                        sourceSaturation += satConst*satConst;
                                    }
                                    //                                    sourceSaturation -= (exp(-np)/(1-exp(-np)))*satConst;
                                    if (pmode == 1) {
                                        if (pflag == 1) {
                                            sourceLA += 1.0;
                                        } else {
                                            sourceSA += 1.0;
                                        }
                                        sourceLAo *= 0.9999;
                                        sourceSAo *= 0.9999;
                                        sourceKo *= 0.9999;
                                        sourceLA *= 0.9999;
                                        sourceSA *= 0.9999;
                                        sourceK *= 0.9999;
                                    } else {
                                        if (pflag == 1) {
                                            sourceLAo += 1.0;
                                        } else {
                                            sourceSAo += 1.0;
                                        }
                                    }
                                
                                }
                                saturateIn(ssource, &largeAngle, &photon, shiftX, shiftY, extraLargeAngle, &largeAngle2);
                                if (diagnostic) diagnosticCount(ray,37);
                            } else {
                                *(state.focalPlane + chip.extraX*(photon.pos.y - minyDrift) + photon.pos.x - minxDrift) += photon.sourceOver_m;
                            }
                            countGood(&localLog, photon.sourceOver_m, &ray);
                            if (throughputfile) addThroughput(&state.throughputLog, minwavelength, maxwavelength, nsurf, waveIndex, photon.sourceOver_m);
                        } else {
                            if (backgroundProbability)  backgroundProbCount(bindex,backgroundProbOk,ssource);
                            if (saturation) {
                                saturateOut(ssource, &largeAngle, &photon, shiftX, shiftY, extraLargeAngle);
                                if (diagnostic) diagnosticCount(ray,42);
                            }
                            if (pmode==1) sourceK+=photon.sourceOver_m; else sourceKo+=photon.sourceOver_m;
                            countBad(&localLog, photon.sourceOver_m, &ray);
                        }
                    }
                    detRay++;
                    }
                }
                break;
            }
        }
    }

    counterAdd(&localLog, &state.counterLog);
    counterAdd(&localLog, &state.globalLog);

    if (lock4 == 1) pthread_mutex_unlock(&lock.lock4);
    if (finish != 0) {
        pthread_mutex_lock(&lock.lock8);
        openthread[thread] = 0;
        openthreads -= finish;
        pthread_mutex_unlock(&lock.lock8);
    }

}
