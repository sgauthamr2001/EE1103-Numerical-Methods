/*************************
Purpose:Real time data analysis of fluoroscence data
Authors:Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Date:28/10/2019
Group   : 2H
Inputs: None
Outputs:Plots of filterd data,locations anf FWHM of peaks,Mean of FWHM's for different number of cells.
Compilation:cc file_analysis.c 
**************************/

#include<stdio.h>
#include<stdlib.h>

#define w 10     //number of elements in the circular buffer
#define N 10000 //number of points being taken at a time
float Array_sort(float *array , int n)//function to find the median of array of data
{ 
	    // declare some local variables
	    int i=0 , j=0 ;
	    float temp=0;

	    for(i=0 ; i<n ; i++)
	    {
		for(j=0 ; j<n-1 ; j++)
		{
		    if(array[j]>array[j+1])
		    {
		        temp        = array[j];
		        array[j]    = array[j+1];
		        array[j+1]  = temp;
		    }
		}
	    }

	    float median=0;
	    if(n%2 == 0)
		median = (array[(n-1)/2] + array[n/2])/2.0;
	    else
		median = array[n/2];
	    
	    return median;  //function returns median
}  
void plot()//function to plot in gnuplot
{
        FILE *pipe = popen("gnuplot -persist", "w");
 	fprintf(pipe, "set terminal pdf\n");
	fprintf(pipe, "set output 'plot.pdf'\n");
	fprintf(pipe, "set title 'Plot of forward and side scatter data after filtering'\n");
        fprintf(pipe, "set xlabel 'X'\n"); 
        fprintf(pipe, "set ylabel 'Y'\n");
        fprintf(pipe, "plot 'FLdata2.txt' u 1:2 w l title 'forward scatter data','FLdata2.txt' u ($1):($3-2500) title 'side scatter data' w l\n");
        fprintf(pipe, "set title 'Plot of side scatter data before and after filtering'\n");
        fprintf(pipe, "set xlabel 'X'\n"); 
        fprintf(pipe, "set ylabel 'Y'\n");
        fprintf(pipe, "set xrange [0:10000]\n");
        fprintf(pipe, "plot 'FLdata.txt' u 0:2 w l title 'before filtering','FLdata2.txt' u ($1):($2) title 'after filtering' w l\n");
        pclose(pipe);
}
int main()
{
	FILE *fp=fopen("FLdata.txt","r"); //file pointer to the given data file
	FILE *fp1=fopen("FLdata1.txt","w+");//file pointer to the file FLdata1.txt
	FILE *fp2=fopen("FLdata2.txt","w+");//file pointer to the file FLdata2.txt

	float arr1[w][2]={0};//array to store values forward scatter and side scatter values in a circular buffer
	int arr2[w]={0};     //array to store values of locations of forward scatter and side scatter
	int i,k,j,h,c;           //defining local variables
	int peak[N/10];      //array to store values of maximum of forward scatter data
	float fwhm[N/10]={0};//array to store values of full width half maximum 
	int peakf[N/10];     //array to store values of maximum of side scatter data
	float xa1,xa2;       //xa1 is moving average value of forward scatter data,xa2 is moving average value of side scatter data
	float ya1,ya2;       //ya1,ya2 are base average values of forward scatter and side scstter data respectively
	float fs,ss;         //fs,ss are the values of forward scatter data,side scatter data read in a particular line in a file
	int lm1,lm2;         //lm1 and lm2 correspond to peak locations of forward scatter data and side scatter data in a loop
	int a=-1;            //a is a loop element 
	int time=0;          //time is a loop variable 
	float maxf,maxs;     //maxf and maxs are varibles to store maximum values of forward and side scatter in a circular buffer 
	int peakl=0,peakc=0; //peakl,peakc are variables used to verify threshhold conditions on locations of maxima's
	int cellcount=0;     //counts the number of cells
	int peakcount=0;     //counts the number of peaks
	float tmin,tmax;     //variables used to find full width half maximum
	float maxsa[N/10];   //array to store side scatter data at maximum locations
	float y[N/10]={0};   //array to store side scatter data after reducing noise
	int x[N/10]={0};     //array to store locations of side scatter data after reducing noise
	float diff[N/10];    //array to store time arrival between fluorosence peaks
	float diffc[N/10];   //array to store time arrival between cells
	float mean=0;        //variables used to find mean,median of the data;
	float meanhm=0;
	float median;
	float sumc=0;


	for(i=0;i<N/w;i++)   //Initialising the circular buffer
       	{             
                
                if(a<(w-1)) //ensures that the the value of a lies between 0 and w-1 ,where w is length of buffer to create a circular buffer
	            a++;
 	        else 
		    a=0;     
                     
                xa1=0; //setting initial values of each of four variables defined above to 0
		xa2=0;
		ya1=0;
		ya2=0;
          
		for(k=0; k<w; k++) //using moving average filter to filter the data

		{

			fscanf(fp,"%f %f",&fs,&ss);
                        xa1+=fs;
			xa2+=ss;
			time++;

		}
                
                
		arr2[a]=time;     //storing filtered data in respective arrays
		arr1[a][0]=xa1/w;
		arr1[a][1]=xa2/w;
                j=a-1;
		lm2=0;            //setting initial values of each varible to 0
		lm1=0;
                maxf=0;
		maxs=0;
		for(int b=0; b<w; b++)//for loop to find maximum value in a buffer
		{
                        if(j<(w-1)) //creating a circular buffer to find maximum values
			   j++;
		        else 
		           j=0;
			ya1 += (arr1[j][0]); //storing values to base average of a filter
                        ya2 += (arr1[j][1]);
                        
                        if(maxf<arr1[j][0]) //finding maxima of forward scatter data
			{
				 maxf=arr1[j][0];
                                 lm1=arr2[j];
			}
                              
			if(maxs<arr1[j][1] && b==w/2)//finding maxima of side scatter data
			{
                                maxs=arr1[j][1];
                                lm2=arr2[j];
                     	}
                         
                   
        
                }

		if((maxs-(ya2/w))>10 && (lm2-peakc)>30 )//checking if found maxima of side scatter data are greater than threshold and distance between them is not to less where 10 is threshold value;
		{
  
                         peakc=lm2;
                         peakcount++; //calculting the total count of peaks
                         fprintf(fp1,"%d %f\n",peakc,maxs);//writng the location and value of peaks in a file

			if((abs(lm2-lm1)<5)) //checking if there is a cell at found peak using seperation between maxima of forward and side scatter data
                        {
				peakl=(lm2+lm1)/2; //calulating the location of the cell
                                peakf[cellcount]=peakl;//array to store locations of cells
                              	cellcount++;
			}

		}
                fprintf(fp2,"%d %f %f\n",arr2[a],arr1[a][1],arr1[a][0]);//writing filtered data into a file

          
         }
      

         fclose(fp1);//closing both the files
         fclose(fp2);
         plot();
	 FILE *fp3=fopen("FLdata1.txt","r");//opening files to read the data
	 FILE *fp4=fopen("FLdata2.txt","r");

	 for(int i=0;i<N;i++) 
	 {   
         	if(i%10==0)         //since every 10th point is filtered in the given data
         	        fscanf(fp4,"%d %f %*f",&(x[i/10]),&y[i/10]);//reads the data file and feeds the corresponding filtered side scattered data 
         }	
         for(a=0;a<peakcount;a++)
         {  
        	 fscanf(fp3,"%d %f",&(peak[a]),&maxsa[a]);//reading the location and value of peaks from the data
	 }
   
 	 fclose(fp3);
 	 fclose(fp4);//closing the data files

         for(a=0;a<peakcount;a++) //finding full width half maxima at the peak
         {
  		 for(i=0;i<N;i++)
 		 {

    		       if(i%10==0)
		       {
      
     			     if(x[i/10]==peak[a])//asserting the location of a peak
                	     {
                		   for(h=-3;h<0;h++)//comparing in the neighbourhood of maxima to find full width half maxima
                   		   {
                                        if(y[(i/10)+h]<=(maxsa[a]/2)<=y[(i/10)+1+h])//linear interpolation to find the locations of half maximag
                                           tmin=x[(i/10)+h] + ((x[(i/10)+1+h]-x[(i/10)+h])/(y[(i/10)+1+h]-y[(i/10)+h]))*((maxsa[a]/2)-y[(i/10)+h]);
                                   }

                                   for(h=0;h<3;h++)
                                   {
                                        if(y[(i/10)+h]<=(maxsa[a]/2)<=y[(i/10)+1+h])
                                           tmax=x[(i/10)+h] + ((x[(i/10)+1+h]-x[(i/10)+h])/(y[(i/10)+1+h]-y[(i/10)+h]))*((maxsa[a]/2)-y[(i/10)+h]);
                                   }
                       
                             }
         
                       }
        
                 
		     if(0<(tmax-tmin) && (tmax-tmin)<50) //bypassing the values thathave a large error
		   		fwhm[a]=tmax-tmin;
		     else
		    		fwhm[a]=(fwhm[a-1]*(a-1)/a);
        
                }

	}

	FILE *fp5=fopen("FLdata1.txt","w+");//opening files to store the required data
	FILE *fp6=fopen("FLdata2.txt","w+");
 	for(a=0;a<peakcount;a++)
	{  
      		 fprintf(fp5,"%d %f\n",peak[a],fwhm[a]);//storing the location of peak and coress ponding full width half maxima associated with it
	}


	meanhm=fwhm[0];    //finding mean fullwidth half maxima for different number of cells
	for(a=1;a<peakcount;a++)
	{
 		 diff[a-1]=peak[a]-peak[a-1];
  		 meanhm = (meanhm*(a) + fwhm[a])/(a+1);
  		 mean= (mean*(a-1)+diff[a-1])/(a);    //finding mean of arrival time of peaks
  
 		 fprintf(fp6,"The mean of width of peaks for %d peaks is %f\n",a,meanhm);
 
	}
 
	for(a=1;a<cellcount;a++)   //finding sum arrival time  between cells
	{
  		diffc[a-1]=peakf[a]-peakf[a-1];
  		sumc += diffc[a-1];
	}
	   printf("The Final Report is (For number of points = %d ,(It can be changed in program)):\n",N);
           printf("\nThe number of peaks is %d\n",peakcount);//priting several values
	   printf("The number of cell is %d\n",cellcount);
           printf("The mean of peak arrival times is %f\n",mean);
	   printf("The mean of cell arrival times is %f\n",sumc/(cellcount-1));
	   median=Array_sort(diffc,cellcount);//function returns the value of median of arrival time of cells
	   printf("The median of cell arrival times is %f\n",median); 
           printf("The mean of width of peaks for %d peaks is %f\n",peakcount,meanhm);
           printf("\nThe locations of peaks and full width half maximas are stored in FLdata1.txt\n");
           printf("The mean fullwidth half maxima for different number of peaks is stored in FLdata2.txt\n");
	  fclose(fp5);//closing the files
	  fclose(fp6);
 
}

		
