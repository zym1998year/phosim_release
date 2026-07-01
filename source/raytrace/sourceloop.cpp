///
/// @package phosim
/// @file sourceloop.cpp
/// @brief source loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Colin Burke (Purdue)
///
/// @warning Please treat results using PhoSim
/// with caution as this project is in active development.
///
/// @warning This material is copyrighted and
/// subject to a open source license with
/// restrictions.  See COPYING for details.
///

#include <thread>
#include <chrono>

int Image::sourceLoop() {

    char tempstring[4096];
    int sourceCounter;
    std::vector<std::thread> thread(numthread);
    thread_args *args;

    //    OPTIMIZATION INITIALIZATION
    args = (thread_args*)malloc(numthread*sizeof(thread_args));
    openthread = static_cast<int*>(malloc(numthread*sizeof(int)));
    state.surfaceLimit = (natmospherefile + 1)*6 + nsurf*2 + 3;
    state.waveLimit = (static_cast<int>(maxwavelength) - static_cast<int>(minwavelength))+ 1;
    if (backgroundProbability > 0) {
        state.bpAcc = static_cast<std::atomic<unsigned int>*>(calloc((bpNx)*(bpNy),sizeof(std::atomic<unsigned int>)));
        state.bpRej = static_cast<std::atomic<unsigned int>*>(calloc((bpNx)*(bpNy),sizeof(std::atomic<unsigned int>)));
    }
    if (diagnostic) state.diagnosticCounter = static_cast<std::atomic<long>*>(calloc(1000,sizeof(std::atomic<long>)));
    state.dynamicTransmission = static_cast<std::atomic<float>*>
        (malloc(state.surfaceLimit*state.waveLimit*4*sizeof(std::atomic<float>)));
    state.dynamicTransmissionLow = static_cast<std::atomic<float>*>
        (malloc(state.surfaceLimit*state.waveLimit*4*sizeof(std::atomic<float>)));
    state.dynamicTransmissionInit = static_cast<std::atomic<int>*>
        (malloc(state.surfaceLimit*state.waveLimit*4*sizeof(std::atomic<int>)));
    for (int i = 0; i < state.waveLimit; i++) {
        for (int j = 0; j < state.surfaceLimit; j++) {
            for (int k = 0; k < 4; k++) {
                state.dynamicTransmission[state.surfaceLimit*state.waveLimit*k + i*state.surfaceLimit + j] = 0.0;
                state.dynamicTransmissionLow[state.surfaceLimit*state.waveLimit*k + i*state.surfaceLimit + j] = 1.0;
                state.dynamicTransmissionInit[state.surfaceLimit*state.waveLimit*k + i*state.surfaceLimit + j] = 0;
            }
        }
    }
    if (checkpointcount != 0) readCheckpoint(checkpointcount);
    counterInit(&state.counterLog);
    counterClear(&state.globalLog);

    // std::mutex / std::condition_variable members are default-constructed (no init needed)
    if (opdfile == 1) {
        remain = 0;
        sourceperthread = 1;
    }
    // integration sequence bookkeeping
    int nsamples = nframes + nskip; // total number of samples in sequence
    int nframei = 0; // current frame number within group
    int ngroupi = 0; // current group number
    std::string outputfilenameold = outputfilename;
    int tframe = 1;
    if (ngroups > 1 || nsamples > 1) tframe = round(vistime/nsnap/devvalue);
    
    //    INTEGRATION SEQUENCE LOOP
    for (int nframe = 0; nframe < tframe; nframe++) {

        for (int i = 0; i < numthread; i++) {
            openthread[i] = 0;
        }
        openthreads = 0;
        int subsource = 0;

        std::ostringstream outfile;
        if (ngroups > 1 || nsamples > 1) {
            ngroupi = floor ( nframe/nsamples);
            std::cout << "------------------------------------------------------------------------------------------" << std::endl;
            std::cout << "Integration Sequence: Group " << (ngroupi + 1) << " Frame " << (nframei + 1) << std::endl;
            std::cout << "------------------------------------------------------------------------------------------" << std::endl;
            outfile << outputfilenameold << "_R" << std::setfill('0') << std::setw(3) << nframe;
        } else {
            outfile << outputfilenameold;
        }
        outputfilename = outfile.str();

        nframei++;
        if (nframe>0) {
            tai += devvalue/86400.0;
            timeoffset += devvalue/86400.0;
        }
            //time offset

        //    LOGGING INITIALIZATION
        if (eventfile) {
            state.pEventLogging = new EventFile((long)(nphot*20), workdir, eventFitsFileName + ".fits");
        } else {
            state.pEventLogging = NULL;
        }

        if (throughputfile) initThroughput(&state.throughputLog, minwavelength, maxwavelength, nsurf);
        float pauseTime = 10000.0*sourceperthread/numthread;

        //    SOURCE TYPE LOOP
        for (int sourceType = 966; sourceType >= 0; sourceType--) {

            sourceCounter = 0;

            if (static_cast<int>(round(sourceType*checkpointtotal/966)) == checkpointcount) {

                //    SOURCE LOOP
                // long subsource = 0;

                for (int source = 0; source < nsource; source++) {

                    if ((sourceType >= 0   && sourceType <=  99 && sources.type[source] == 0 &&
                         (static_cast<int>(source/backgroundTypes)%100) == sourceType) ||
                        (sourceType >= 100 && sourceType <= 199 && sources.type[source] == 1 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 100)) ||
                        (sourceType >= 200 && sourceType <= 299 && sources.type[source] == 2 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 200)) ||
                        (sourceType >= 300 && sourceType <= 399 && sources.type[source] == 3 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 300)) ||
                        (sourceType >= 400 && sourceType <= 499 && sources.type[source] == 4 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 400)) ||
                        (sourceType >= 500 && sourceType <= 599 && sources.type[source] == 5 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 500)) ||
                        (sourceType >= 600 && sourceType <= 699 && sources.type[source] == 6 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 600)) ||
                        (sourceType >= 700 && sourceType <= 799 && sources.type[source] == 7 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 700)) ||
                        (sourceType >= 800 && sourceType <= 899 && sources.type[source] == 8 &&
                         (static_cast<int>(source/backgroundTypes)%100) == (sourceType - 800)) ||
                        (sourceType == 900 && sources.type[source] == 9 && sources.mag[source] > 32.25) ||
                        (sourceType >= 901 && sourceType <= 965 && sources.type[source] == 9 &&
                         sources.mag[source] >  (482.25 - static_cast<double>(sourceType)/2.0) &&
                         sources.mag[source] <= (482.75 - static_cast<double>(sourceType)/2.0))
                        || (sourceType == 966 && sources.type[source] == 9 && sources.mag[source] <= -0.25)) {

                        int count = 0;
                        while  (openthreads >= numthread*sourceperthread) {
                            // int state, type, status=0;
                            // status=pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,&state);
                            // status=pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED,&type);
                            std::this_thread::sleep_for(std::chrono::nanoseconds(static_cast<long long>(pauseTime)));
                            count++;
                            // pthread_setcancelstate(state,&state);
                            // pthread_setcanceltype(type,&type);
                        }
                        if (count > 0) {
                            pauseTime = pauseTime * pow(static_cast<double>(count)/100.0, 0.5/static_cast<double>(numthread));
                            if (pauseTime < 1.0) pauseTime = 1.0;
                        }
                        //                        printf("%d %lf\n",count,pauseTime/1000000000);

                        lock.lock8.lock();
                        int bestthread = sourceperthread + 1;
                        for (int i = 0; i < numthread; i++) {
                            if ((openthread[i] < sourceperthread) && (openthread[i] < bestthread)) {
                                subsource = i;
                                bestthread = openthread[i];
                            }
                        }
                        openthread[subsource] += 1;
                        openthreads++;
                        sourceCounter++;
                        args[subsource].ssource[openthread[subsource] - 1] = source;
                        int ok = 0;
                        if (openthread[subsource] == sourceperthread) ok=1;
                        lock.lock8.unlock();

                        if (ok) {
                            args[subsource].instance = this;
                            args[subsource].thread = subsource;
                            args[subsource].runthread = sourceperthread;
                            if (opdfile == 1) {
                                if (((sourceCounter - 1) % numthread) == 0) {
                                    if (source < static_cast<int>(floor(nsource/numthread))*numthread) {
                                        remain = numthread;
                                    } else {
                                        remain = nsource - static_cast<int>(floor(nsource/numthread))*numthread;
                                    }
                                }
                            }
                            if (opdfile == 0) {
                                thread[subsource] = std::thread(&Image::threadFunction, &args[subsource]);
                                thread[subsource].detach();

                            } else {
                                thread[subsource] = std::thread(&Image::threadFunction, &args[subsource]);
                            }
                            //if (opdfile == 0) pthread_detach(thread[subsource]);
                        }
                        if (opdfile == 1) {
                            if ((((sourceCounter - 1) % numthread) == (numthread - 1)) ||
                                (source == nsource - 1)) {
                                for (int ss = 0; ss < numthread; ss++) {
                                    if (openthread[ss] > 0 && thread[ss].joinable()) {
                                        thread[ss].join();
                                    }
                                }
                            }
                        }

                    }

                }
            }
            if (sourceType >= 0   && sourceType < 100) snprintf(tempstring, sizeof(tempstring), "Dome Light     %3d%% ", (100 - sourceType));
            if (sourceType >= 100 && sourceType < 200) snprintf(tempstring, sizeof(tempstring), "Airglow-Coll   %3d%% ", (200 - sourceType));
            if (sourceType >= 200 && sourceType < 300) snprintf(tempstring, sizeof(tempstring), "Airglow-Phot   %3d%% ", (300 - sourceType));
            if (sourceType >= 300 && sourceType < 400) snprintf(tempstring, sizeof(tempstring), "Scattered Moon %3d%% ", (400 - sourceType));
            if (sourceType >= 400 && sourceType < 500) snprintf(tempstring, sizeof(tempstring), "Scattered Moon %3d%% ", (500 - sourceType));
            if (sourceType >= 500 && sourceType < 600) snprintf(tempstring, sizeof(tempstring), "Zodiacal       %3d%% ", (600 - sourceType));
            if (sourceType >= 600 && sourceType < 700) snprintf(tempstring, sizeof(tempstring), "Scattered Sun  %3d%% ", (700 - sourceType));
            if (sourceType >= 700 && sourceType < 800) snprintf(tempstring, sizeof(tempstring), "Scattered Sun  %3d%% ", (800 - sourceType));
            if (sourceType >= 800 && sourceType < 900) snprintf(tempstring, sizeof(tempstring), "Artificial     %3d%% ", (900 - sourceType));
            if (sourceType == 900)                     snprintf(tempstring, sizeof(tempstring), "Astro Object m>32.0 ");
            if (sourceType >= 901 && sourceType < 966) snprintf(tempstring, sizeof(tempstring), "Astro Object m=%4.1f ", 482.5 - static_cast<double>(sourceType)/2.0);
            if (sourceType == 966)                     snprintf(tempstring, sizeof(tempstring), "Astro Object m< 0.0 ");
            if (sourceCounter > 0) counterCheck(&state.counterLog, sourceCounter, tempstring);
        }
        for (int i=0; i<numthread; i++) {
            if (openthread[i] != sourceperthread && openthread[i] != 0) {
                args[i].instance = this;
                args[i].thread = i;
                args[i].runthread = openthread[i];
                thread[i] = std::thread(&Image::threadFunction, &args[i]);
                thread[i].detach();
            }
        }
        while (openthreads != 0) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(static_cast<long long>(pauseTime)));
        }

        // COSMIC RAYS
        if ((!opdfile) && (backgroundMode>0) && (raydensity>0.0)) {
            long rays;
            double rate;
            double depth;
            double cosz=cos(zenith);
            depth = 1100.0*exp(-groundlevel/7000.0);
            if (spaceMode > 0) depth = 0.0;
            if (spaceMode > 0) cosz=1.0;

            rate = raydensity*17000.0*exp(-depth/121.0/cosz);
            if (cosz < 0.0 && spaceMode <=0) rate = 0.0;
            particleLoop(rate, 1e9, 100e9, -2.7, 0.0, 3, 1, &rays);
            snprintf(tempstring, sizeof(tempstring), "Cosmic Protons      ");
            counterCheck(&state.counterLog, rays, tempstring);

            rate = raydensity*1700.0*exp(-depth/121.0/cosz);
            if (cosz < 0.0 && spaceMode <= 0) rate = 0.0;
            particleLoop(rate, 1e9, 100e9, -2.7, 0.0, 4, 1, &rays);
            snprintf(tempstring, sizeof(tempstring), "Cosmic Alphas       ");
            counterCheck(&state.counterLog, rays, tempstring);

            rate = raydensity*200.0*exp(-depth/50.0/cosz);
            if (cosz < 0.0 && spaceMode <= 0) rate = 0.0;
            particleLoop(rate, 1e9, 100e9, -2.7, 0.0, 1, 1, &rays);
            snprintf(tempstring, sizeof(tempstring), "Cosmic Electrons    ");
            counterCheck(&state.counterLog, rays, tempstring);

            rate = raydensity*1600.0*exp(-depth/350.0)*(1 - exp(-depth/100.0));
            particleLoop(rate, 1e9, 100e9, -2.7, 2.0, 2, 1, &rays);
            snprintf(tempstring, sizeof(tempstring), "Secondary Muons     ");
            counterCheck(&state.counterLog, rays, tempstring);

            rate = raydensity*100.0*exp(-depth/180.0)*(1 - exp(-depth/100.0));
            particleLoop(rate, 1e9, 100e9, -2.7, 2.0, 2, 1, &rays);
            snprintf(tempstring, sizeof(tempstring), "Secondary Electrons ");
            counterCheck(&state.counterLog, rays, tempstring);

            if (spaceMode == 0) {
                particleLoop(100000.0*raydensity, 1e5, 2e6, 0.0, 0.0, 0, 1, &rays);
                snprintf(tempstring, sizeof(tempstring), "Terrestrial Gammas  ");
                counterCheck(&state.counterLog, rays, tempstring);
            }
        }

        counterFinish(&state.counterLog);
        // OUTPUT DATA
        if ((nframe % nsamples) < nframes) {
            if (checkpointcount == checkpointtotal) writeImageFile();
            if (opdfile) writeOPD();
            if (checkpointcount !=  checkpointtotal) writeCheckpoint(checkpointcount);
            if (centroidfile) writeCentroidFile(workdir, outputfilename, sourcePhoton, sourceXpos, sourceYpos, sources.id, nsource);
            if (throughputfile) writeThroughputFile(workdir, outputfilename, &state.throughputLog, minwavelength, maxwavelength, nsurf);
            if (eventfile) state.pEventLogging->eventFileClose();
        }

            if (diagnostic) {
                //        if (((ssource % 11) ==0) && (thread==0)) {
            printf("----------\n");
            //            double total = 0;
            for (int ii=0; ii < 1000; ii++) {
                //total+= state.diagnosticCounter[ii];
            }
            int used[1000];
            for (int ii=0; ii< 1000; ii++) used[ii]=0;
            for (int ii=0; ii < 1000; ii++) {
                long highvalue=0;
                int bv=0;
                for (int jj=0; jj < 1000; jj++) {
                    long value = state.diagnosticCounter[jj]*1000;
                    if (value >= highvalue && used[jj]==0) {
                        highvalue=value;
                        bv=jj;
                    }
                }
                used[bv]=1;
                if (highvalue>0) {
                    long value = state.diagnosticCounter[bv]*1000;
                    if (bv==0)  printf("$$%40s %22e %22.6lf\n","getWavelengthTime",(double)value,30.0); //
                    if (bv==1)  printf("$$%40s %22e %22.6lf\n","getAngle",(double)value,11.0); //
                    if (bv==2)  printf("$$%40s %22e %22.6lf\n","dynamicTransmissionOptimizationZero",(double)value,97.0); //
                    if (bv==3)  printf("$$%40s %22e %22.6lf\n","dynamicTransmissionOptimization",(double)value,215.0); //
                    if (bv==4)  printf("$$%40s %22e %22.6lf\n","samplePupil",(double)value,45.0); //
                    if (bv==5)  printf("$$%40s %22e %22.6lf\n","diffraction",(double)value,146.0);
                    if (bv==6)  printf("$$%40s %22e %22.6lf\n","largeAngleScattering",(double)value,12.0); //
                    if (bv==7)  printf("$$%40s %22e %22.6lf\n","getDeltaAngle",(double)value,157.0);
                    if (bv==8)  printf("$$%40s %22e %22.6lf\n","atmospherePropagate",(double)value,19.0);//
                    if (bv==9)  printf("$$%40s %22e %22.6lf\n","atmosphereIntercept",(double)value,19.0);//
                    if (bv==10) printf("$$%40s %22e %22.6lf\n","atmosphereRefraction",(double)value,174.0);
                    if (bv==11) printf("$$%40s %22e %22.6lf\n","cloudOpacity",(double)value,48.0);
                    if (bv==12) printf("$$%40s %22e %22.6lf\n","atmosphereOpacityMoon",(double)value,26.0);
                    if (bv==13) printf("$$%40s %22e %22.6lf\n","cloudOpacityMoon",(double)value,17.0);
                    if (bv==14) printf("$$%40s %22e %22.6lf\n","atmosphereOpacitySun",(double)value,26.0);
                    if (bv==15) printf("$$%40s %22e %22.6lf\n","cloudOpacitySun",(double)value,17.0);
                    if (bv==16) printf("$$%40s %22e %22.6lf\n","mieScatter",(double)value,228.0);
                    if (bv==17) printf("$$%40s %22e %22.6lf\n","atmosphereOpacity",(double)value,52.0);
                    if (bv==18) printf("$$%40s %22e %22.6lf\n","atmosphereDiffraction",(double)value,119.0);
                    if (bv==19) printf("$$%40s %22e %22.6lf\n","domeSeeing",(double)value,13.0);
                    if (bv==20) printf("$$%40s %22e %22.6lf\n","tracking",(double)value,14.0); //
                    if (bv==21) printf("$$%40s %22e %22.6lf\n","propagate",(double)value,3.0); //
                    if (bv==22) printf("$$%40s %22e %22.6lf\n","transform",(double)value,55.0);
                    if (bv==23) printf("$$%40s %22e %22.6lf\n","findSurface",(double)value,99.0);
                    if (bv==24) printf("$$%40s %22e %22.6lf\n","chooseSurface",(double)value,34.0);
                    if (bv==25) printf("$$%40s %22e %22.6lf\n","interceptDerivatives",(double)value,60.0);
                    if (bv==26) printf("$$%40s %22e %22.6lf\n","contaminationSurfaceCheck",(double)value,39.0);
                    if (bv==27) printf("$$%40s %22e %22.6lf\n","surfaceCoating",(double)value,32.0);
                    if (bv==28) printf("$$%40s %22e %22.6lf\n","transformInverse",(double)value,47.0);
                    if (bv==29) printf("$$%40s %22e %22.6lf\n","reflect",(double)value,44.0);//
                    if (bv==30) printf("$$%40s %22e %22.6lf\n","newRefractionIndex",(double)value,24.0);
                    if (bv==31) printf("$$%40s %22e %22.6lf\n","refract",(double)value,41.0);
                    if (bv==32) printf("$$%40s %22e %22.6lf\n","pixelPosition",(double)value,2.0);//
                    if (bv==33) printf("$$%40s %22e %22.6lf\n","getDeltaIntercept",(double)value,42.0);
                    if (bv==34) printf("$$%40s %22e %22.6lf\n","photonSiliconPropagate",(double)value,96.0);
                    if (bv==35) printf("$$%40s %22e %22.6lf\n","electronSiliconPropagate",(double)value,124.0);
                    if (bv==36) printf("$$%40s %22e %22.6lf\n","pixelFraction",(double)value,2.0);
                    if (bv==37) printf("$$%40s %22e %22.6lf\n","saturateIn",(double)value,75.0);
                    if (bv==38) printf("$$%40s %22e %22.6lf\n","atmosphericDispersion",(double)value,126.0);
                    if (bv==39) printf("$$%40s %22e %22.6lf\n","airIndexRefraction",(double)value,44.0);
                    if (bv==40) printf("$$%40s %22e %22.6lf\n","secondKick",(double)value,52.0);
                    if (bv==41) printf("$$%40s %22e %22.6lf\n","contaminationDetectorCheck",(double)value,30.0);
                    if (bv==42) printf("$$%40s %22e %22.6lf\n","saturateOut",(double)value,28.0); //
                    if (bv==43) printf("$$%40s %22e %22.6lf\n","photonLoop",(double)value,187.0);
                    if (bv==44) printf("$$%40s %22e %22.6lf\n","gaussianXY",(double)value,3.0);//
                    if (bv==45) printf("$$%40s %22e %22.6lf\n","inBounds",(double)value,1.0);//
                    if (bv==46) printf("$$%40s %22e %22.6lf\n","inBoundsBuffer",(double)value,1.0);//
                }
            }
           printf("----------\n");
        }
            //    }

    }

    return(0);

}
