///
/// @package phosim
/// @file particleloop.cpp
/// @brief particle loop
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

void Image::particleLoop(double flux, double energyMin, double energyMax, double energyIndex, double angularIndex, int particleType, int downward, long *nparticles) {

    long long particles;
    double particlesDouble;
    double energy = 0.0;
    double theta = 0.0;
    double mass = 0.0;
    double phi = 0.0;
    double step = 1e-3;
    int miss = 0;
    int charge = 0;
    int spin = 0;
    Photon photon;
    Vector position;
    Vector angle;
    Vector tempAngle;
    Vector largeAngle;
    Vector largeAngle2;
    Vector origPosition;
    Vector deltaInitAngle[1000];
    Vector deltaPosition[1000];
    double deltaAngle[1000];
    double deltaPhiAngle[1000];
    double deltaEnergy[1000];
    int deltaFlag = 0;
    int particle = 0;
    double frac1, frac2;
    double rutherfordAngle[1000];
    double rutherford[1000];
    long ray;
    double kleinnishinaAngle[1000];
    double kleinnishina[1000];
    double over = 1.0;
    Clog localLog;

    counterClear(&localLog);
    particle=particleType;
    vectorInit(&largeAngle);
    vectorInit(&largeAngle2);
    vectorInit(&position);
    vectorInit(&origPosition);
    vectorInit(&angle);
    ray = 0;
    double thetacut=0.5*PI/180.0;
    for (int i = 0; i < 1000; i++) {
        rutherfordAngle[i] = PI*i/(1000-1);
        if (rutherfordAngle[i] >= thetacut) {
            rutherford[i]=1.0/(pow(sin(rutherfordAngle[i]/2.0),4.0))*(sin(rutherfordAngle[i]));
        } else {
            rutherford[i]=0.0;
        }
    }
    double total = 0.0;
    for (int i = 0; i < 1000; i++) total+=rutherford[i];
    for (int i = 0; i < 1000; i++) rutherford[i]=rutherford[i]/total;
    for (int i = 1; i < 1000; i++) rutherford[i]=rutherford[i]+rutherford[i-1];
 
    double area = (pixelsx+2*activeBuffer)*pixsize*1e3*(pixelsy+2*activeBuffer)*pixsize*1e3*2.0 +
        (pixelsx+2*activeBuffer)*pixsize*1e3*sensorthickness*1e3*2.0 +
        (pixelsy+2*activeBuffer)*pixsize*1e3*sensorthickness*1e3*2.0;
    frac1 = (pixelsx+2*activeBuffer)*pixsize*1e3*sensorthickness*1e3*2.0/area;
    frac2 = (pixelsy+2*activeBuffer)*pixsize*1e3*sensorthickness*1e3*2.0/area;
    particlesDouble = flux * 2 * PI /(angularIndex + 1) * (area*1e-6*1e-6)*exptime;
    if (particleType == 0) particlesDouble = particlesDouble/over;
    if (poissonMode==1) {
        particles = random[0].poisson(particlesDouble);
    } else {
        particles = static_cast<long>(particlesDouble);
    }
    //printf("Particles %ld\n",particles);
    if (particleType !=0) *nparticles = particles;
    else *nparticles = 0;
    for (long long i = 0; i < particles ; i++ ) {

        // choose mass
        if (particleType == 0) {mass = 0; charge = 0; spin=2;}
        if (particleType == 1) {mass = M_ELECTRON; charge = -1; spin=1;}
        if (particleType == 2) {mass = M_MUON; charge = -1; spin=1;}
        if (particleType == 3) {mass = M_PROTON; charge = 1; spin=1;}
        if (particleType == 4) {mass = M_ALPHA; charge = 2; spin=0;}

        particle=particleType;

        // choose Energy
        double r = random[0].uniform();
        if (energyIndex == -1) {
            energy = exp((log(energyMax) - log(energyMin))*r + log(energyMin));
        } else {
            energy = pow((pow(energyMax, energyIndex + 1) - pow(energyMin, energyIndex + 1))*r +
                         pow(energyMin, energyIndex + 1), 1.0/(energyIndex + 1));
        }
        energy += mass*C_LIGHT*C_LIGHT/E_CHARGE;

        // choose Angles
        theta = acos(pow(1-random[0].uniform(),1.0/(angularIndex + 1)));
        if (downward == 1) theta = PI - theta;
        //need tilt for vertical detectors and zenith
        theta = theta - zenith - body[nsurf-1][2];
        phi = 2 * PI * random[0].uniform();
        angle.x = sin(theta)*cos(phi);
        angle.y = sin(theta)*sin(phi);
        angle.z = cos(theta);

        // choose x, y, z
        double rr=random[0].uniform();
        double ss=random[0].uniform();
        if (rr < frac1) {
            position.x = (random[0].uniform() * (pixelsx + 2*activeBuffer) - pixelsx/2 - activeBuffer) * pixsize;
            if (ss<= 0.5) position.y = (pixelsy/2 + activeBuffer) * pixsize; else position.y = -(pixelsy/2 + activeBuffer) * pixsize;
            position.z = -(random[0].uniform())*sensorthickness;
        } else if (rr < frac1+frac2) {
            if (ss<= 0.5) position.x = (pixelsx/2 + activeBuffer) * pixsize; else position.x = -(pixelsx/2 + activeBuffer) * pixsize;
            position.y = (random[0].uniform() * (pixelsy + 2*activeBuffer) - pixelsy/2 - activeBuffer) * pixsize;
            position.z = -(random[0].uniform())*sensorthickness;
        } else {
            position.x = (random[0].uniform() * (pixelsx + 2*activeBuffer) - pixelsx/2 - activeBuffer) * pixsize;
            position.y = (random[0].uniform() * (pixelsy + 2*activeBuffer) - pixelsy/2 - activeBuffer) * pixsize;
            if (ss <= 0.5) position.z = 0.0; else position.z = -sensorthickness;
        }


        if (deltaFlag >= 1) {
            if (deltaFlag >= 1000) {
                deltaFlag--;
            } else {
                mass=M_ELECTRON;
                charge = -1;
                spin = 1;
                particle=1;
                vectorCopy(deltaInitAngle[deltaFlag-1],&angle);
                vectorCopy(deltaPosition[deltaFlag-1],&position);
                shift_mu(&angle,deltaAngle[deltaFlag-1],deltaPhiAngle[deltaFlag-1]);
                energy=deltaEnergy[deltaFlag-1] + mass*C_LIGHT*C_LIGHT/E_CHARGE;
                if (energy < mass) energy=mass;
                //printf("beta %lf %lf\n",deltaEnergy[deltaFlag-1],sqrt(1.0-pow(mass*C_LIGHT*C_LIGHT/E_CHARGE/energy, 2.0)));
                deltaFlag--;
            }
        }

        mass=mass*C_LIGHT*C_LIGHT/E_CHARGE;

        double scattering=0.0;
        if (particle == 0) {
            for (int k = 0; k < 1000; k++) {
                kleinnishinaAngle[k] = PI*k/(1000-1);
                double pp = 1.0/(1.0 + energy/(M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE)*(1.0 - cos(kleinnishinaAngle[k])));
                double alpha=1.0/137.04;
                double rc=3.8616e-13;
                kleinnishina[k] = 0.5*alpha*alpha*rc*rc*pp*pp*(pp+1.0/pp-pow(sin(kleinnishinaAngle[k]),2.0))*sin(kleinnishinaAngle[k])*PI*1.0/(1000-1)*2*PI;
                if (kleinnishina[k] < 0) printf("Error: %e %e %e %d\n",pp,kleinnishinaAngle[k],energy,k);
            }
            double total = 0.0;
            for (int k = 0; k < 1000; k++) {
                total+=kleinnishina[k];
            }
            for (int k = 0; k < 1000; k++) kleinnishina[k]=kleinnishina[k]/total;
            for (int k = 1; k < 1000; k++) {
                kleinnishina[k]=kleinnishina[k]+kleinnishina[k-1];
            }
            scattering = total*over;
        }

        //        long counter =0;
        while (1) {
        // propagate
        // make step size of 1 micron
        // still doing this in mm
        position.x += angle.x * step;
        position.y += angle.y * step;
        position.z += angle.z * step;

        // check if in silicon
        if ((abs(position.x/pixsize) < pixelsx) &&
            (abs(position.y/pixsize) < pixelsy) &&
            ((position.z < 0) && (position.z > -sensorthickness))) {

            // if not photon produce ionization & do delta rays
            if (particle != 0) {

                double K=K_DEDX/E_CHARGE*1e4;
                double Z=silicon.atomicZ;
                double A=silicon.atomicA;
                double excitation=silicon.excitation;
                double density=silicon.density;
                double pairenergy=silicon.pairenergy;
                double fano=silicon.fano;

                double mecc=M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE;
                double beta=sqrt(1.0-pow(mass/energy, 2.0));
                double gamma=1.0/sqrt(1.0-beta*beta);
                double plasmaenergy=28.816*sqrt(density*Z/A);
                double densityeffect=log(plasmaenergy/excitation)+log(beta*gamma)-0.5;
                double Tmax = 2.0*mecc*pow(beta, 2.0)*pow(gamma, 2.0)/(1.0+2*gamma*mecc/mass + pow(mecc/mass, 2.0));
                double Tcut = 100.0*excitation;
                double de=-density*K*pow(charge, 2.0)*Z/A/beta/beta*(0.5*log(2.0*mecc*pow(beta, 2.0)*pow(gamma, 2.0)*Tcut/excitation/excitation)-beta*beta/2.0*(1+Tcut/Tmax)-densityeffect)*(1e-4*step*1e3);
                energy += de;
                if (energy <= mass) {
                    de = -(energy+de-mass);
                    energy = mass;
                }
                double meanelectrons=-de/pairenergy;
                double ee=meanelectrons+random[0].uniform()*sqrt(meanelectrons*fano);
                if (ee < 0) ee=0;
                long electrons=static_cast<long>(ee);
                //                counter += electrons;

                // Delta rays
                double stepe=0.1;
                for (double aa=0.0; aa<9.0; aa+=stepe) {
                    double T = pow(10.0, aa);
                    if ((T > Tcut) && (T < Tmax)) {
                        double dT = pow(10.0, aa+stepe) - pow(10.0, aa);
                        double dx = step*1e3*1e-4;
                        double F = 0.0;
                        if (spin == 1) F=(1-beta*beta*T/Tmax + T*T/2.0/energy/energy);
                        if (spin == 0) F=(1-beta*beta*T/Tmax);
                        double dN=0.5*K*charge*charge*density*Z/A/beta/beta*F/T/T*dT*dx;
                        if (random[0].uniform() < dN) {
                            if (deltaFlag >= 1000) {
                                printf("Warning:  Too many delta rays\n");
                            } else {
                                double pe=sqrt(T*T+2*mecc*T);
                                double pmax =sqrt(Tmax*Tmax + 2*mecc*Tmax);
                                deltaAngle[deltaFlag] = T/pe*pmax/Tmax;
                                deltaPhiAngle[deltaFlag] = 2*PI*random[0].uniform();
                                vectorCopy(position,&deltaPosition[deltaFlag]);
                                vectorCopy(angle,&deltaInitAngle[deltaFlag]);
                                deltaEnergy[deltaFlag] = T;
                                //printf("Delta Ray of %lf keV and cos angle of %lf\n",T/1e3,deltaAngle[deltaFlag]);
                                energy=energy - T;
                                if (energy <= mass) energy=mass;
                                deltaFlag++;
                                i--;
                            }
                        }
                    }
                }


                double epsilon0=8.854e-12;
                beta=sqrt(1.0-pow(mass/energy,2.0));
                gamma=1.0/sqrt(1.0-beta*beta);
                double massm=mass*1.782e-36;
                double momentum=massm*beta*gamma*C_LIGHT;
                double n = (density*1e3/(A*M_PROTON));
                double sigma = (2*PI*Z*Z*pow(E_CHARGE/epsilon0,2.0)/64.0/PI/PI);
                sigma=sigma/C_LIGHT/C_LIGHT;
                sigma = sigma/beta/beta/momentum/momentum*n*E_CHARGE*E_CHARGE;
                sigma = sigma * ((2.0/pow(sin(thetacut/2.0),2.0)-2.0) - beta*beta*(-4.0*log(sin(thetacut/2.0))));
                double mfp = 1.0/(sigma);
                double prob=(step*1e-3)/mfp;
                if (random[0].uniform() < prob) {
                    int index;
                    find(rutherford,1000,random[0].uniform(),&index);
                    //printf("Rutherford scatter %e\n",rutherfordAngle[index]*180/PI);
                    shift_mu(&angle,cos(rutherfordAngle[index]),random[0].uniform()*2*PI);
                }

                vectorCopy(position,&origPosition);

                for (long j = 0; j < electrons; j++) {
                    vectorCopy(origPosition,&position);
                    //                    int notFinished = 1;
                    int eCounter = 0;


                    if (detectorMode == 1) {
                        //while (notFinished == 1) {
                            do {


                            pixelPosition(&position,&photon);
                            pixelFraction(&position,&photon);

                            photon.z0=0.0;
                            photon.thread=0;
                            photon.xindex = find_linear_float_mc(silicon.temperatureGrid, silicon.numTemperature, sensorTempNominal + sensorTempDelta, random[0].uniform());
                            photon.location = floor((chip.nampx + 2*activeBuffer)/DET_RESAMPLE)*floor((photon.pos.y - miny + activeBuffer)/DET_RESAMPLE) + floor((photon.pos.x - minx + activeBuffer)/DET_RESAMPLE);

                            if (photon.pos.x < minx - activeBuffer || photon.pos.x > maxx + activeBuffer ||
                                photon.pos.y < miny - activeBuffer || photon.pos.y > maxy + activeBuffer) {
                                long xL, yL;
                                xL = photon.pos.x;
                                yL = photon.pos.y;

                                if (photon.pos.x < minx - activeBuffer) xL = minx - activeBuffer;
                                if (photon.pos.x > maxx + activeBuffer) xL = maxx + activeBuffer;
                                if (photon.pos.y < miny - activeBuffer) yL = miny - activeBuffer;
                                if (photon.pos.y > maxy + activeBuffer) yL = maxy + activeBuffer;
                                photon.location = floor((chip.nampx + activeBuffer)/DET_RESAMPLE)*floor((yL - miny + activeBuffer)/DET_RESAMPLE) + floor((xL - minx + activeBuffer)/DET_RESAMPLE);
                            }
                            photon.collect_z = photon.z0 - sensorthickness;

                            miss = electronSiliconPropagate(&angle, &position, &photon, eCounter, 0);
                            
                            //                            if (miss == 1) notFinished = 0;
                            if (photon.z0 > photon.collect_z) {
                                if (position.z <= photon.collect_z) {
                                    //                                    notFinished = 0;
                                    position.z = photon.collect_z;
                                }
                            }
                            if (photon.z0 <= photon.collect_z) {
                                if (position.z >= photon.collect_z) {
                                    //                                    notFinished = 0;
                                    position.z = photon.collect_z;
                                }
                            }
                            eCounter++;
                            if (miss) eCounter=siliconSegments;
                            if (photon.zindex >= (SILICON_STEPS-siliconSubSteps)) eCounter=siliconSegments;

                            } while (eCounter < siliconSegments);
                            //if (eCounter > (static_cast<double>(SILICON_STEPS)/static_cast<double>(siliconSubSteps))) notFinished = 0;

                            //                    }
                    }

                    pixelPosition(&position,&photon);
                    photon.sourceOver_m = 1;
                    photon.pos.x = static_cast<long>(floor(position.x/pixsize + pixelsx/2));
                    photon.pos.y = static_cast<long>(floor(position.y/pixsize + pixelsy/2));
                    photon.sourceOver_m=1;
                    if (photon.pos.x >= minx && photon.pos.x <= maxx && photon.pos.y >= miny && photon.pos.y <= maxy) {

                        if (saturation) {
                            double shiftX = 0.0;
                            double shiftY = 0.0;
                            double extraLargeAngle = 0.0;
                            sources.spatialtype.push_back(0);
                            sources.type.push_back(6);

                            saturateIn(0, &largeAngle, &photon, shiftX, shiftY, extraLargeAngle, &largeAngle2);
                        } else {
                            *(state.focalPlane + chip.nampx*(photon.pos.y - miny) +
                              (photon.pos.x - minx)) += photon.sourceOver_m;
                        }
                        photon.sourceOver_m = 1;
                        countGood(&localLog, photon.sourceOver_m, &ray);
                    } else {
                        photon.sourceOver_m = 1;
                        countBad(&localLog, photon.sourceOver_m, &ray);
                    }

                }
                vectorCopy(origPosition,&position);
            }


            if (particle == 0) {
            // if photon consider compton scattering
                double A=28.055;
                double density=2.32;
                double prob = step*1e-3*(scattering*(density*1e3/(A*M_PROTON)));
                if (random[0].uniform() < prob) {
                    *nparticles = *nparticles + 1;
                    //printf("Compton Scatter of ray %ld!\n",i);
                    int index;
                    find(kleinnishina,1000,random[0].uniform(),&index);
                    double phi = random[0].uniform()*2*PI;
                    double pp = 1.0/(1.0 + energy/(M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE)*(1.0 - cos(kleinnishinaAngle[index])));
                    double eprime = pp*energy;
                    double ee=energy - eprime + M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE;
                    //printf("Energies: %e %e %e %e %e\n",energy,eprime,ee,pp,cos(kleinnishinaAngle[index]));
                    double pe=sqrt(ee*ee - pow(M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE, 2.0));
                    double mue = (energy - eprime*cos(kleinnishinaAngle[index]))/(pe*C_LIGHT);
                    double phie = phi + PI;
                    // now shift
                    energy = eprime;
                    vectorCopy(angle,&tempAngle);
                    shift_mu(&angle,cos(kleinnishinaAngle[index]),phi);
                    if (deltaFlag >= 1000) {
                        printf("Warning:  Too many compton scatters %d\n",deltaFlag);
                    } else {
                        deltaAngle[deltaFlag] = mue;
                        deltaPhiAngle[deltaFlag] = phie;
                        vectorCopy(position,&deltaPosition[deltaFlag]);
                        vectorCopy(tempAngle,&deltaInitAngle[deltaFlag]);
                        deltaEnergy[deltaFlag] = ee - M_ELECTRON*C_LIGHT*C_LIGHT/E_CHARGE;
                        //printf("Compton Scatter of %lf keV and cos angle of %lf\n",deltaEnergy[deltaFlag]/1e3,deltaAngle[deltaFlag]);
                        deltaFlag++;
                        i--;
                    }
 
                }
            }

        }
        // if outside volume or energy is 0 then get rid of
        if ((position.z > 0) || (position.z < -sensorthickness) ||
            (position.x > ((pixelsx/2 + activeBuffer) * pixsize)) ||
            (position.x < ((-pixelsx/2 - activeBuffer) * pixsize)) ||
            (position.y > ((pixelsy/2 + activeBuffer) * pixsize)) ||
            (position.y < ((-pixelsy/2 - activeBuffer) * pixsize))) {

            //if (counter > 0) printf("exited silicon %e %ld\n",position.z,counter);
            break;
        }
        if (energy <= mass) {
            //printf("lost all energy %e %e %ld\n",energy,mass,counter);
            break;
        }

        }


    }

    counterAdd(&localLog, &state.counterLog);
    counterAdd(&localLog, &state.globalLog);

}
