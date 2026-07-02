///
/// @package phosim
/// @file photonoptimization.cpp
/// @brief photon optimization routines (part of Image class)
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

int Image::dynamicTransmissionOptimization(int thread, int k, int *lastSurface, int *preGhost, int waveSurfaceIndex, int straylightcurrent, Photon *aph, int sourceType) {

    int newSurf;
    int status;
    int ambiguous1;
    int ambiguous2;
    int ambiguous;
    int ghost1;
    int ghost2;
    int run;
    int ccounter;
    int ghost;
    int absoluteCounter = 0;
    int temp = 0;

    // new event
    //redoghlp:;
    absoluteCounter++;
    if (absoluteCounter > 10000) {
        return(1);
    }
    *preGhost = 0;
    ghost1 = 0;
    ghost2 = 0;
    ambiguous1 = 0;
    ambiguous2 = 0;
    run = 0;
    ccounter = 0;

    // restart ambiguous cases
 redoamb:;
    aph->counter = -1;
    ambiguous = 0;
    newSurf = -1;
    aph->direction = 1;
    aph->preGhostFlag = 0;
    ghost = 0;
    int sc = 0;
    int ambsc = 0;
    for (int i = -1; i < natmospherefile; i++) {
        int iTemp = (i+1)*6 + waveSurfaceIndex;
        if (sourceType!=0) {
            if (atmospheremode == 2) {

                if (i >= 0) {
                    if (aph->counter >= (MAX_BOUNCE - 1)) {
                        ghost++;
                        ghost++;
                        goto maxbounce;
                    }

                    if ((sourceType==4) && (sc==0)) {
                        status = transmissionPreCheck(iTemp + 1, aph, run);
                        if (status == 1) {
                            status = transmissionPreCheck(iTemp + 4, aph, run);
                            if (status == 1) {
                                //
                            } else {
                                status = transmissionPreCheck(iTemp + 3, aph, run);
                                if (status == 1) {
                                    //
                                } else {
                                    sc += 1;
                                }
                            }
                        }
                    }

                    if ((sourceType==7) && (sc==0)) {
                        status = transmissionPreCheck(iTemp + 1, aph, run);
                        if (status == 1) {
                            status = transmissionPreCheck(iTemp + 6, aph, run);
                            if (status == 1) {
                                //
                            } else {
                                status = transmissionPreCheck(iTemp + 5, aph, run);
                                if (status == 1) {
                                    //
                                } else {
                                    sc += 1;
                                }
                            }
                        }
                    }

                    if (((sourceType!=4) && (sourceType!=7)) || (sc>0)) {
                        status = transmissionPreCheck(iTemp + 1, aph, run);
                        //can't optimize this because small chance of small angle scatter
                        //                if (status == 1) {
                        //if (k > 0) goto redoghlp;
                        //return(1);
                        //}
                        if (sourceType==8 && status==1) sc++;
                        if (sourceType==8 && status==2) ambsc=1;
                        if (status == 2) {
                            // assume pass
                        }
                    }
                }

                if ((sourceType==4) && (sc == 0)) {
                    status = transmissionPreCheck(iTemp + 2, aph, run);
                    if (status == 1) {
                        status = transmissionPreCheck(iTemp + 4, aph, run);
                        if (status == 1) {
                            //
                        } else {
                            status = transmissionPreCheck(iTemp + 3, aph, run);
                            if (status == 1) {
                                //
                            } else {
                                sc += 1;
                            }
                        }
                    }
                }

                if ((sourceType==6) && (sc == 0)) {
                    status = transmissionPreCheck(iTemp + 2, aph, run);
                    if (status == 1) {
                        status = transmissionPreCheck(iTemp + 6, aph, run);
                        if (status == 1) {
                            //
                        } else {
                            status = transmissionPreCheck(iTemp + 5, aph, run);
                            if (status == 1) {
                                //
                            } else {
                                sc += 1;
                            }
                        }
                    }
                }

                if (((sourceType!=4) && (sourceType!=7)) || (sc > 0)) {
                    if ((sourceType == 0) || (sourceType==1) ||
                        (sourceType == 2) || (sourceType==5) || (sourceType==8)) {
                        status = transmissionPreCheckSuppression(iTemp + 2, aph, run);
                    } else {
                        status = transmissionPreCheck(iTemp + 2, aph, run);
                    }
                    if (sourceType == 8) {
                        if (status == 1) sc++;
                        if (status == 2) ambsc=1;
                    } else {
                        if (status == 1) {
                            //                                                        if (k > 0) goto redoghlp;
                            return(1);
                        }
                    }
                    if (status == 2) {
                        // assume pass
                    }
                }
            }
        }
    }

    if ((sourceType==3) || (sourceType==4) || (sourceType==6) || (sourceType==7) ) {
        if (sc != 1) {
            //                        if (k > 0) goto redoghlp;
            return(1);
        }
    }
    if (sourceType==8) {
        if ((sc != 1) && (ambsc==0)) {
            //                        if (k > 0) goto redoghlp;
            return(1);
        }
    }

    temp = (natmospherefile+1)*6 + 2 + waveSurfaceIndex;
    while (1) {
        if (aph->direction == 1) {
            newSurf++;
        } else {
            newSurf--;
        }
        if (newSurf == -1) {
            if (ambiguous>=1) goto maxbounce;
            //                        if (k > 0) goto redoghlp;
            return(1);
        }

        if (aph->direction == 1) {
            *lastSurface = newSurf;
        }
        if (aph->counter >= (MAX_BOUNCE - 1)) goto maxbounce;
        status = transmissionPreCheck(temp + 2*newSurf, aph, run);
        if ((status==2) && (straylightcurrent==1)) {
            ambiguous++;
        }
        if (status == 1 || ((status==2) && (run==1))) {
            if (straylightcurrent == 1) {
                if (aph->counter >= (MAX_BOUNCE - 1)) {
                    ghost++;
                    ghost++;
                    goto maxbounce;
                }
                status = transmissionPreCheckNoRoll(temp + 2*newSurf + 1, aph);
                if (status == 1) {
                    if (ambiguous>=1) goto maxbounce;
                    //                                        if (k > 0) goto redoghlp;
                    return(1);
                } else if (status == 0) {
                    aph->direction = -aph->direction;
                    aph->preGhostFlag = aph->preGhostFlag + 1;
                    ghost++;
                } else if (status == 2) {
                    //assume pass
                    aph->direction = -aph->direction;
                    aph->preGhostFlag = aph->preGhostFlag + 1;
                    ghost++;
                }
            } else {
                if (ambiguous>=1) goto maxbounce;
                //                                if (k > 0) goto redoghlp;
                return(1);
            }
        }
        if (aph->direction == 1 && newSurf == nsurf - 1) {

            if (aph->counter >= (MAX_BOUNCE - 1)) {
                ghost++;
                ghost++;
                goto maxbounce;
            }
            status = transmissionPreCheck(temp + 2*nsurf, aph, run);
            if (status == 1) {
                if (ambiguous>=1) goto maxbounce;
                //                                if (k > 0) goto redoghlp;
                return(1);
            } else if (status == 0) {
                goto maxbounce;
            } else if (status == 2) {
                // assume pass
                goto maxbounce;
            }
        }

        if (aph->direction == -1 && newSurf == 0) {
            if (run==1) { //FIX!!!!!!
                //*preGhost=0;
                //return(0);
                goto maxbounce;
            }
            //                        if (k > 0) goto redoghlp;
            if (ambiguous>=1) goto maxbounce;
            return(1);
        }
    }
 maxbounce:;


    if (straylightcurrent==1) {
        if (run==0) ambiguous1=ambiguous; else ambiguous2=ambiguous;
        if (run==0) ghost1=ghost; else ghost2=ghost;
        if (ambiguous1 == 0) {
            //proceed
        } else {
            if (run == 0) {
                ccounter = aph->counter;
                run++;
                goto redoamb;
            }
        }
        *preGhost = ghost1;
        if (ghost2 > *preGhost) *preGhost=ghost2;
        if (ambiguous1 > *preGhost) *preGhost=ambiguous1;
        if (ambiguous2 > *preGhost) *preGhost=ambiguous2;
    }
    aph->maxcounter = aph->counter;
    if (ccounter > aph->maxcounter) aph->maxcounter=ccounter;

    return(0);

}




int Image::dynamicTransmissionOptimizationZero(int *lastSurface, int *preGhost, int waveSurfaceIndex, int straylightcurrent, Photon *aph, int sourceType) {

    int newSurf;
    int status;

    *preGhost = 0;

    aph->counter = -1;
    newSurf = -1;
    aph->direction = 1;
    aph->preGhostFlag = 0;
    int sc = 0;
    int ambsc = 0;
    for (int i = -1; i < natmospherefile; i++) {
        int iTemp = (i+1)*6 + waveSurfaceIndex;
        if (sourceType != 0) {
            if (atmospheremode == 2) {

                if (i >= 0) {
                    if ((sourceType==4) && (sc==0)) {
                        status = transmissionPreCheckZero(iTemp + 1, aph);
                        if (status == 1) {
                            status = transmissionPreCheckZero(iTemp + 4, aph);
                            if (status != 1) {
                                status = transmissionPreCheckZero(iTemp + 3, aph);
                                if (status != 1) {
                                    sc += 1;
                                }
                            }
                        }
                    }
                    if ((sourceType==7) && (sc==0)) {
                        status = transmissionPreCheckZero(iTemp + 1, aph);
                        if (status == 1) {
                            status = transmissionPreCheckZero(iTemp + 6, aph);
                            if (status != 1) {
                                status = transmissionPreCheckZero(iTemp + 5, aph);
                                if (status != 1) {
                                    sc += 1;
                                }
                            }
                        }
                    }
                    if (((sourceType!=4) && (sourceType!=7)) || (sc>0)) {
                        status = transmissionPreCheckZero(iTemp + 1, aph);
                        if (sourceType==8) {
                            if (status==1) sc++;
                            if (status==2) ambsc=1;
                        }
                            //can't optimize this because small chance of small angle scatter
                        // if (status == 1) {
                        //     return(1);
                        // }
                    }
                }

                if ((sourceType==4) && (sc == 0)) {
                    status = transmissionPreCheckZero(iTemp + 2, aph);
                    if (status == 1) {
                        status = transmissionPreCheckZero(iTemp + 4, aph);
                        if (status != 1) {
                            status = transmissionPreCheckZero(iTemp + 3,  aph);
                            if (status != 1) {
                                sc += 1;
                            }
                        }
                    }
                }

                if ((sourceType==7) && (sc == 0)) {
                    status = transmissionPreCheckZero(iTemp + 2, aph);
                    if (status == 1) {
                        status = transmissionPreCheckZero(iTemp + 6, aph);
                        if (status != 1) {
                            status = transmissionPreCheckZero(iTemp + 5, aph);
                            if (status != 1) {
                                sc += 1;
                            }
                        }
                    }
                }

                if (((sourceType!=4) && (sourceType!=7)) || (sc > 0)) {
                    if ((sourceType == 0) || (sourceType==1) ||
                        (sourceType == 2) || (sourceType==5) || (sourceType==8)) {
                        status = transmissionPreCheckSuppression(iTemp + 2, aph, 0);
                    } else {
                        status = transmissionPreCheckZero(iTemp + 2, aph);
                    }
                    if (sourceType==8) {
                        if (status == 1) {
                            sc++;
                        }
                        if (status == 2) {
                            ambsc=1;
                        }
                    } else {
                        if (status == 1) {
                            return(1);
                        }
                    }
                }
            }
        }
    }
    if ((sourceType==3) || (sourceType==4) || (sourceType==6) || (sourceType==7)) {
        if (sc != 1) {
            return(1);
        }
    }
    if (sourceType==8) {
        // printf("%d %d\n",sc,ambsc);
        if ((sc != 1) && (ambsc==0)) {
            return(1);
        }
    }

    int temp = (natmospherefile+1)*6 + 2 + waveSurfaceIndex;
    while (1) {
        newSurf++;
        *lastSurface = newSurf;
        status = transmissionPreCheckZero(temp + 2*newSurf, aph);
        if (status == 1) {
            return(1);
        }
        if (newSurf == nsurf - 1) {
            status = transmissionPreCheckZero(temp + 2*nsurf, aph);
            if (status == 1) {
                return(1);
            } else {
                break;
            }
        }

    }

    aph->maxcounter = aph->counter;
    return(0);

}
