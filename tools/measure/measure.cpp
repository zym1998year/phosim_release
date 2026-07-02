/*Created by JRP. Modified and updated by AD*/


#include <fitsio.h>
#include <fitsio2.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.1415926535897932846264


// function to swap elements
void swap(float *a, float *b) {
  float t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
int partition(float array[], int low, int high) {
  
  // select the rightmost element as pivot
  float pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swap(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSort(float array[], int low, int high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partition(array, low, high);
    
    // recursive call on the left of pivot
    quickSort(array, low, pi - 1);
    
    // recursive call on the right of pivot
    quickSort(array, pi + 1, high);
  }
}

// function to print array elements
void printArray(int array[], int size) {
  for (int i = 0; i < size; ++i) {
    printf("%d  ", array[i]);
  }
  printf("\n");
}


void calculate(char fits_file[100], double guess_flux, double guess_mux, double guess_muy, double guess_alphax, double guess_alphay, double guess_alphaxy, int counter_target, int distorted, int fixbackground, long *naxes, float *array){
	
    
    double *data;
    float *data2;
	
    
    //Adjust cut size. If array >100x100 reduce it to 100x100
	//axes[0] is dimension along x 
    int cut=1000;
	if(naxes[0] <30 or naxes[1]<30)
	{
		printf("Error:  Image too small \n");
        exit(1);
	}
    if(naxes[0] <= naxes[1])
    {
		cut = naxes[0]/2;
    }
    else
    {
		cut = naxes[1]/2;
    }

    
    // if (cut > 50)
    // {
    //     	cut = 50;
    // }

    //Resize the image if necessary
    if ((naxes[0] > (2*cut)) || (naxes[1] > (2*cut))) {
               
	       int maxx = naxes[0]/2;
	       int maxy = naxes[1]/2;
           int lowx=maxx-cut;
           if (lowx < 0) lowx=0;
           int lowy=maxy-cut;
           if (lowy < 0) lowy=0;
           int highx = maxx+cut;
           if (highx > naxes[0]-1) highx=naxes[0]-1;
           int highy = maxy+cut;
           if (highy > naxes[1]-1) highy=naxes[1]-1;
           data2=static_cast<float*>(calloc((highx-lowx)*(highy-lowy),sizeof(float)));
       		//Copy cut data to new array called data2
           for (int i = lowx; i < highx; i++) {
               for (int j = lowy; j< highy; j++) {
                   data2[(j-lowy)*(highx-lowx)+(i-lowx)]=array[j*naxes[0]+i];
               }
           }
	      //Copy cut data data2 to array. Update naxes
           naxes[0]=highx-lowx;
           naxes[1]=highy-lowy;
           for (int i = 0; i < naxes[0]; i++) {
               for (int j=0; j < naxes[1]; j++) {
                   array[j*naxes[0]+i]=data2[j*naxes[0]+i];
               }
           }
    }

    //Find total counts
    double tot3 = 0.0;
    for (long i=0;i<naxes[0];i++) {
        for (long j=0;j<naxes[1];j++) {
            tot3+=array[j*naxes[0]+i];
        }
    }
    //printf("Total counts in the image:  %lf\n",tot3);

    //Copy data from array to data
    double sizex=static_cast<double>(naxes[0]);
    double sizey=static_cast<double>(naxes[1]);
    data=static_cast<double*>(calloc(naxes[0]*naxes[1],sizeof(double)));
    for (int i = 0; i < naxes[0]; i++) {
        for (int j=0; j < naxes[1]; j++) {
            data[j*naxes[0]+i]=static_cast<double>(array[j*naxes[0]+i]);
        }
    }
    //Sort array. Used to find median
	
    long nn=naxes[0]*naxes[1];
	quickSort(array, 0, nn - 1);
	
	
	long n0;
	double med;
	if(nn %2 != 0)
	{
		n0 = (nn+1) / 2 - 1;
		med = array[n0];
	}
	else
	{
		n0 = (nn+1) / 2 - 1;
		med = (array[n0] +array[n0+1])/2.0;
	}
    printf("Median of image is %lf counts\n", med);
    
    long count=0;  //Count is basically a variable which denotes number of non zero elements in array.
    //We use only those 
    if (med == 0) {
        for (int i=0;i<sizex;i++) {
            for (int j=0;j<sizey;j++) {
               
                    count++;
                
            }
        }
    } else {
        count=naxes[0]*naxes[1];
    }
    double *x, *y, *w, *weight;
    long xx, yy;
    x=static_cast<double*>(calloc(count,sizeof(double)));
    y=static_cast<double*>(calloc(count,sizeof(double)));
    w=static_cast<double*>(calloc(count,sizeof(double)));
    weight=static_cast<double*>(calloc(count,sizeof(double)));
    count=0;
    for (int i=0;i<naxes[0];i++) {
        for (int j=0;j<naxes[1];j++) {
            
		    x[count]=i-(sizex/2.0-0.5);
		    y[count]=j-(sizey/2.0-0.5);
		    w[count]=data[j*naxes[0]+i];
		    count++;
            
        }
    }

    double total = 0.0;
    for (long i=0;i<count;i++) total+=w[i];
    //FILE *fp;
    double a=0.0;
    double b=0.0;
    double mux=guess_mux;
    double muy=guess_muy;
    double alphax=guess_alphax;
    double alphay=guess_alphay;
    double alphaxy=guess_alphaxy;
    double background=med;
    double t1, t2, t3, t4, t5, t6, t7, w_mod;
    double sigmaxx, sigmayy, sigmaxy, area, signal, e1, e2, prevsigmaxx, prevsigmayy, delsigmaxx, delsigmayy, size, q;
    int converged = 1;
    double ellip, pa;
    ellip=pa=0.0;
	delsigmaxx = delsigmayy = prevsigmaxx = prevsigmayy = 9999.0;
//	double guess_area = PI*sqrt(pow(guess_alphax,2) * pow(guess_alphay,2) - pow(guess_alphaxy,2)) ;
        double guess_area = PI*(0.5*pow(guess_alphax,2) + 0.5*pow(guess_alphay,2));
	long iteration = 0;
    
	while ( ((abs(delsigmaxx) > 1e-3) || (abs(delsigmayy) > 1e-3)) && (iteration<counter_target)){
		
		//printf("%lf and %lf \n", alphax, alphay);
		//Catching bad things i.e nans before happening
		if(isnan(sigmaxx) == 1 or isnan(sigmayy) == 1 or isnan(sigmaxy) == 1 or isnan(signal) == 1 or isinf(sigmaxx) == 1 or isinf(sigmayy) == 1 or isinf(sigmaxy) == 1 or isinf(signal) == 1)
		{
			converged =0; 
			break;
		}
		if(abs(alphax) >cut  or abs(alphay) > cut or abs(alphaxy) > cut)
		{
			converged =0; 
			break;
		}
        t1=0.0;
        t2=0.0;
        t3=0.0;
		t4=0.0;
        t5=0.0;
        t6=0.0;
        t7=0.0;

        //If gaussian is too sharply peaked restrict weight alphas
		if(alphax < 0.9)
	    	alphax = 0.9;
		if(alphay < 0.9)
	    	alphay = 0.9;
		// Making sure the sqrt in weight isnt giving nans if root of negative number
		while(pow((alphaxy/alphax/alphay),2.0) >=1.0 )
		{
			if(alphaxy > 0)
			alphaxy = alphaxy - 0.1;
			if(alphaxy < 0)
			alphaxy = alphaxy + 0.1;
			if(abs(alphaxy) > 5000)
			{
				converged = 0;
				break;
			}
		}
        for (long i = 0; i < count; i++) {
		        weight[i]=(exp(-(pow((x[i]-mux),2.0)/alphax/alphax-
		                         2.0*alphaxy/alphax/alphax/alphay/alphay*(x[i]-mux)*(y[i]-muy)+
		                         pow((y[i]-muy),2.0)/alphay/alphay)/2.0/(1-pow((alphaxy/alphax/alphay),2.0))))/
		            (2*PI*alphax*alphay*sqrt(1-pow((alphaxy/alphax/alphay),2.0)));
			
			if(distorted == 0)
			{
				w_mod = w[i]-background;
				//if(i== 0 or i == 1)
				//	printf("%lf and %lf and %lf \n", w_mod, w[i], background);
			}
			if(distorted == 1)
			{
				
				q = abs(w[i]-background)/sqrt(w[i]);
				if((w[i]-background) < 0)
				{
					w_mod = 1.41421*sqrt(w[i]) * (0.477/ exp(0.551*q - 0.06*q*q + 0.003*q*q*q));
				}
				else
				{
					w_mod = w[i]-background+ sqrt(w[i]) * (0.477/ exp(0.928*q + 0.247*q*q + 0.04*q*q*q))*1.41421;
				}
				//if(i== 0 or i == 1)
				//	printf("%lf and %lf and %lf \n", w_mod, w[i], background);
				
			}

	        t1 += (x[i]*y[i]*(w_mod)*weight[i]);
	        t2 += ((w_mod)*weight[i]);
	        t3 += (x[i]*(w_mod)*weight[i]);
	        t4 += (y[i]*(w_mod)*weight[i]);
	        t5 += (x[i]*x[i]*(w_mod)*weight[i]);
	        t6 += (y[i]*y[i]*(w_mod)*weight[i]);
	        t7 += (weight[i]*weight[i]);	
			
        }
        //printf("%lf and %lf and %lf and %lf and %lf and %lf and %lf \n" ,t1,t2,t3,t4,t5,t6, t7);
        sigmaxy=(t1/t2-t3*t4/t2/t2);
        sigmaxx=(t5/t2-t3*t3/t2/t2);
        sigmayy=(t6/t2-t4*t4/t2/t2);
	
		if(sigmaxx< 0 or sigmayy<0)
		{
			converged =0; 
			break;

		}
		    mux=(t3/t2);
		    muy=(t4/t2);
		    signal=t2/t7;
		    if (med != 0){
		        background=(total-signal)/sizex/sizey;
                        if (fixbackground!=0) {
                            background=med;
                            signal = total - background*count;
                        }
		        //total = signal + sizex*sizey*background;
		}
                    //        area = PI* sqrt(2*sigmaxx*2*sigmayy - 2*sigmaxy) ;
        area = PI* (sigmaxx+sigmayy) ;
        e1=(sigmaxx-sigmayy)/(sigmaxx+sigmayy);
        e2=(2.0*sigmaxy)/(sigmaxx+sigmayy);
        ellip=sqrt(e1*e1+e2*e2);
        pa=0.5*atan2(e2,e1);
		size = sqrt(sigmaxx+sigmayy);

	
	   
		delsigmaxx = prevsigmaxx - sigmaxx;
		delsigmayy = prevsigmayy - sigmayy;
		
		prevsigmaxx = sigmaxx;
		prevsigmayy = sigmayy;
	
	    alphax=sigmaxx*2.0;
        if (alphax < 0.9) alphax=0.9;
        alphax=sqrt(alphax);
        alphay=sigmayy*2.0;
        if (alphay < 0.9) alphay=0.9;
        alphay=sqrt(alphay);
        alphaxy=sigmaxy*2.0;
		iteration ++;
	

    }
    if (converged == 0) {
        printf("Error: Did not converge.\n");
    } 
	else {
		//This is to fix flux for super small objects(a few pixels)
		if(alphax<0.75 or alphay<0.75)
		{
			t2=0;
			t7=0;
			for (long i = 0; i < count; i++) {
		            weight[i]=(exp(-(pow((x[i]-mux),2.0)/alphax/alphax-
		                         2.0*alphaxy/alphax/alphax/alphay/alphay*(x[i]-mux)*(y[i]-muy)+
		                         pow((y[i]-muy),2.0)/alphay/alphay)/2.0/(1-pow((alphaxy/alphax/alphay),2.0))))/
		            (2*PI*alphax*alphay*sqrt(1-pow((alphaxy/alphax/alphay),2.0)));
			t2 += ((w[i]-background)*weight[i]);
			t7 += (weight[i]*weight[i]);
			}
			signal=t2/t7;
		}
		
		//Apply distortion corrections for 1 iteration
		if(distorted == 1 and counter_target == 1)
		{
			FILE *myFile = NULL;

			myFile = fopen("./calib_final.dat","r");

			if (myFile == NULL)
			{
				printf("\n Unable to read lookup table1 \n");
				exit(1);
			}

			if (myFile != NULL)
			{
				//printf("\n\nFile read succesfully a\n");
			}

			char buffer[80];
	 		float lut[218][3];
		    int row = 0;
		    int column = 0;

		    while (fgets(buffer,80, myFile)) {
		        column = 0; 
		        // Splitting the data
		        char* value = strtok(buffer, ", ");
		        while (value) {
		            
		            // Column 2
		            if (column == 1 or column == 2 or column==3) {
		                lut[row][column-1] = atof(value);	
		            }
	 
		            column++;
					value = strtok(NULL, ",");
		        }
	 			row++;
		        
		    }
	 
		    // Close the file
		    fclose(myFile);
			double sqrtba = sqrt(med)* guess_area;
			double fbysqrtba = log10(guess_flux/sqrtba);
			double sigma_factor, flux_bias;
			int index;
			if(fbysqrtba <-3.275)
			{
				index = 0;
				flux_bias = lut[index][1];
				sigma_factor = lut[index][2];
			}
			if(fbysqrtba > 2.1)
			{
				index = 217;
				flux_bias = lut[index][1];
				sigma_factor = lut[index][2];
			}
			if(fbysqrtba >-3.275 and fbysqrtba < 2.1)
			{
				index = int((fbysqrtba+3.3)/0.025);
				flux_bias = lut[index][1] + ((lut[index+1][1]-lut[index][1])/(lut[index+1][0]-lut[index][0])) * (fbysqrtba - lut[index][0]);
				sigma_factor = lut[index][2] + ((lut[index+1][2]-lut[index][2])/(lut[index+1][0]-lut[index][0])) * (fbysqrtba - lut[index][0]);
			}
			//printf("%lf and %lf and %lf", signal, med, flux_bias);
			signal = signal - flux_bias*sqrtba;
			sigmaxx = sigmaxx - sigma_factor*guess_alphax*guess_alphax;
			sigmayy = sigmayy - sigma_factor*guess_alphay*guess_alphay;
			sigmaxy = sigmaxy - sigma_factor*guess_alphaxy;
			size = sqrt(sigmaxx + sigmayy);
		
			
		}
		//Apply distortion corrections for convergence
		if(distorted == 1 and counter_target>1)
		{
			FILE *myFile = NULL;
			myFile = fopen("./calib_final_inf.dat","r");
			if (myFile == NULL)
			{
				printf("\n Unable to read lookup table2 \n");
				exit(1);
			}
			if (myFile != NULL)
			{
				
				//printf("\n\nFile read succesfully b\n");
			}
			char buffer[1024];
	 		float lut[35][3];
		    int row = 0;
		    int column = 0;
	 
		    while (fgets(buffer,
		                 1024, myFile)) {
		        column = 0;
		        
	 
		        // Splitting the data
		        char* value = strtok(buffer, ", ");
	 
		        while (value) {
		            
		            // Column 2
		            if (column == 1 or column == 2 or column==3) {
		                lut[row][column-1] = atof(value);
		            }
	 
		            column++;
					value = strtok(NULL, ",");
		        }
	 			row++;
				
		        
		    }
	 
		    // Close the file
		    fclose(myFile);
			double sqrtba = sqrt(med)* guess_area;
			double fbysqrtba = log10(guess_flux/sqrtba);
			double sigma_factor, flux_bias;
			int index;
			if(fbysqrtba <1.225)
			{
				index = 0;
				flux_bias = lut[index][1];
				sigma_factor = lut[index][2];
			}
			if(fbysqrtba > 2.825)
			{
				index = 34;
				flux_bias = lut[index][1];
				sigma_factor = lut[index][2];
			}
			if(fbysqrtba >1.225 and fbysqrtba < 2.825)
			{
				index = int((fbysqrtba-1.175)/0.05);
				flux_bias = lut[index][1] + ((lut[index+1][1]-lut[index][1])/(lut[index+1][0]-lut[index][0])) * (fbysqrtba - lut[index][0]);
				sigma_factor = lut[index][2] + ((lut[index+1][2]-lut[index][2])/(lut[index+1][0]-lut[index][0])) * (fbysqrtba - lut[index][0]);
			}
			//printf("%lf and %lf and %lf", signal, med, flux_bias);
			signal = signal - flux_bias*sqrtba;
			sigmaxx = sigmaxx - sigma_factor*guess_alphax*guess_alphax;
			sigmayy = sigmayy - sigma_factor*guess_alphay*guess_alphay;
			sigmaxy = sigmaxy - sigma_factor*guess_alphaxy;
			size = sqrt(sigmaxx + sigmayy);


		}	
		double backerr = (background*area*8)/signal/signal;
		printf("\n");
		printf("PhotometricFlux/Counts               %12lf +/- %12lf\n",signal,sqrt(signal+background*area*8));
		printf("AstrometricPositionX/Pixels          %12lf +/- %12lf\n",mux,sqrt(area/PI/signal + backerr*area/PI) );
		printf("AstrometricPositionY/Pixels          %12lf +/- %12lf\n",muy,sqrt(area/PI/signal + backerr*area/PI));
		printf("PsfSize/Pixels                       %12lf +/- %12lf\n",size, sqrt(0.5)*sqrt(area/PI/signal + backerr*area/PI));
		printf("PsfEllipticityComponent1             %12lf +/- %12lf\n",e1,sqrt(2.0)*sqrt(1/signal + backerr) );
		printf("PsfEllipticityComponent2             %12lf +/- %12lf\n",e2,sqrt(2.0)*sqrt(1/signal + backerr) );
		printf("BackgroundCountRate/(Counts/Pixels)  %12lf +/- %12lf\n",background,sqrt(background*count)/count);
		    
		/*A better estimation of elliptcity errors are of the following form 
		a = 8*sigmaxx**2 *sigmayy**2 /(sigmaxx + sigmayy)**4
		b = 16*sigmaxx**2 *sigmayy**2 /(sigmaxx + sigmayy)**4
		error in e1 = sqrt((a/N + b*backerr) * sqrt(1-e2) )
		error in e2 = sqrt((a/N + b*backerr) * ((1-e2)/sqrt(1-e1)) )
		*/
	   	//fprintf(fp, "PhotometricFlux/Counts               %12lf  +/- %12lf\n", signal,sqrt(signal)*backerr);
	   	//fprintf(fp, "AstrometricPositionX/Pixels          %12lf  +/- %12lf\n", mux,sqrt(2.0)*rms/sqrt(signal)*backerr);
		//fprintf(fp, "AstrometricPositionY/Pixels          %12lf  +/- %12lf\n",muy,sqrt(2.0)*rms/sqrt(signal)*backerr);
		//fprintf(fp, "PsfSize/Pixels                       %12lf  +/- %12lf\n",rms,rms/sqrt(signal)*backerr);
		//fprintf(fp, "PsfEllipticityComponent1             %12lf  +/- %12lf\n",e1,sqrt(2.0)/sqrt(signal)*backerr);
		//fprintf(fp, "PsfEllipticityComponent2    ,         %12lf  +/- %12lf\n",e2,sqrt(2.0)/sqrt(signal)*backerr);
		//fprintf(fp, "BackgroundCountRate/(Counts/Pixels)  %12lf  +/- %12lf\n",background,sqrt(background*count)/count);
	   	//fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", alphax, alphay, alphaxy, signal, pa, e1, e2, sigmaxx, sigmayy,sigmaxy, mux, muy, rms, background);
		//fprintf(fp, "%s \n", argv[1]);
    }
    //fclose(fp);
}

int main(int argc, char **argv) {
	double guess_flux, guess_alphax, guess_alphay, guess_alphaxy, guess_mux, guess_muy;
	int counter_target, distorted, fixbackground;
	//Below are the guess values of flux, centroid and sigmas
	//Guess flux is ONLY used for distorted measurement.

    if (argc > 2) {
        guess_flux = atof(argv[2]);
        guess_mux = atof(argv[3]);
        guess_muy = atof(argv[4]);
        guess_alphax = atof(argv[5]); 
        guess_alphay = atof(argv[6]);
        guess_alphaxy = atof(argv[7]);
        counter_target = atoi(argv[8]); //Change this to desired number of iteration
        distorted = atoi(argv[9]); //Change this 1 to force switch ON distortion
        fixbackground = atoi(argv[10]);
    } else {
        guess_flux = 9000.0;
        guess_mux = 0;
        guess_muy = 0;
        guess_alphax = 3.0; 
        guess_alphay = 3.0;
        guess_alphaxy = 0.0;
        counter_target = 100; //Change this to desired number of iteration
        distorted = 0; //Change this 1 to force switch ON distortion
        fixbackground = 0;
    }
	if(distorted != 0 && distorted != 1)
	{
		printf("Invalid Distortion number. Input either 0 or 1 \n");
		return 1;
	}
	fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    char keyname[4096];
    char comment[4096];
    char value[4096];
    float *array;
	// Reading the fits file into array Array is the image
    int status=0;
	if (fits_open_file(&faptr,argv[1],READONLY,&status)) {
		printf("Error:  Could not open %s\n",argv[1]);
		exit(1);
		}
		
	fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
	array=static_cast<float*>(calloc(naxes[0]*naxes[1],sizeof(float)));
	fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,array,&anynull,&status);
	fits_close_file(faptr,&status);
        calculate(argv[1], guess_flux, guess_mux, guess_muy, guess_alphax, guess_alphay, guess_alphaxy,counter_target, distorted, fixbackground, naxes, array );
    
		
	
}
