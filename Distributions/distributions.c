/*****************
Purpose:Task 4
Authors:Gagan Deep,Santosh,Gautham
Date:28/08/2019
Inputs:None
Outputs:Different Distributions
**************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/******* mean1 is a function to calculate mean for values of normal distribution*******/
 
 
double mean1(double* values, int n)
  {
      int i;
      int m = n + n % 2;
      double s = 0;
      for ( i = 0; i < m; i+=2)
          s += values[i]+values[i+1];
      return s/m;
  }

/******* mean2 is a function to calculate mean for values of Rayleigh and Lorentzian distribution*******/

double mean2(double* values, int n)
  {
      int i;
      double s = 0;
      for ( i = 0; i < n; i+=1)
         s += values[i];
      return s/n;
  }
  
/******* stddev1 is a function to calculate Standard deviation for values of normal distribution*******/ 

double stddev1(double* values, int n)
  {
      int i;
      int m = n + n % 2; 
      double average = mean1(values,m);
      double s = 0;
      for ( i = 0; i < m; i+=2 )
          s += ((values[i] - average) * (values[i] - average))+((values[i+1] - average) * (values[i+1] - average));
      return sqrt(s / (m - 1));
  }

/******* stddev2 is a function to calculate Standard deviation for values of Rayleigh and Lorentzian distribution*******/ 

double stddev2(double* values, int n)
  {
      int i;
      double average = mean2(values,n);
      double s = 0;
      for ( i = 0; i < n; i++ )
          s += (values[i] - average) * (values[i] - average);
      return sqrt(s / (n - 1));
  } 
/*****************Rayleigh random number generator***************************/

double* generaterrm(int n,double *high,double *low)
    {
      int i;
      double* values = (double*)calloc(n,sizeof(double));
      double average, deviation;
      FILE*fptr=fopen("b.txt","w+");  
      if ( values )
      {
       	  for ( i = 0; i < n; i +=1 )
          {
                double x;
                x = 1.0* rand() / (double)RAND_MAX;
                values[i] =2*sqrt(-2*log(x));

           
                if(values[i]<=*low)
                   *low=values[i];
           
                if(values[i]>=*high)
                   *high=values[i];
                                                     
                       
                fprintf(fptr,"%f\n",values[i]);  
                 
          }
                fclose(fptr);   
    }
    return values;
   }
/*****************Lorentzian random number generator***************************/
double* generatelrm(int n,double *high,double *low)
  {
     int i;
     double* values = (double*)calloc(n,sizeof(double));
     double average, deviation;
     FILE*fptr=fopen("b.txt","w+");  
     if ( values )
     {
         for ( i = 0; i < n; i +=1 )
          {
                double x;
                x = 1.0* rand() / (double)RAND_MAX;
                values[i] = tan(M_PI*(x)-(M_PI/2));
             
                if(values[i]<=*low)
                   *low=values[i];
           
                if(values[i]>=*high)
                   *high=values[i];                
                                    
                fprintf(fptr,"%f\n",values[i]);  
                 
           }
                fclose(fptr);   
     }
         return values;
   }
 

/*****************Normal random numbers generator - Marsaglia algorithm********************/
 
double* generatenrm(int n,double *high,double *low)
 {
    int i;
    int m = n + n % 2;
    double* values = (double*)calloc(m,sizeof(double));
    double average, deviation;
    FILE*fptr=fopen("b.txt","w+");  
    if ( values )
    {
        for ( i = 0; i < m; i += 2 )
         {
            double x,y,rsq,f;
            do
            {
                x = 2.0 * rand() / (double)RAND_MAX - 1.0;
                y = 2.0 * rand() / (double)RAND_MAX - 1.0;
                rsq = x * x + y * y;
            }
            while( rsq >= 1. || rsq == 0. );
               f = sqrt( -2.0 * log(rsq) / rsq );
               values[i]   = x * f;
               values[i+1] = y * f;


            if(values[i]<*low)
               *low=values[i];
            if(values[i+1]<*low)
               *low=values[i+1];
            if(values[i]>*high)
               *high=values[i];
            if(values[i+1]>*high)
               *high=values[i];
           
                                  
            fprintf(fptr,"%f %f\n",values[i],values[i+1]);  
                 
         }
                fclose(fptr);   
    }
    return values;
 }
 
/***********printHistogram1 prints histogram for values of normal Distribution********/ 
void printHistogram1(double* values,int n,double high,double low)
 {
     const int width = 50;    
     int max = 0;
 
    
     int m = n + n % 2;
    
 
     int i,j,k;
     int nbins = ceil(1+3.322*log10((n+n%2)));
     double delta =(high-low)/nbins;
     int* bins = (int*)calloc(nbins,sizeof(int));
     if ( bins != NULL )
     {
         for ( i = 0; i < m; i++ )
         {
             int j = (int)( (values[i] - low) / delta );
             if ( 0 <= j  &&  j <=nbins )
                 bins[j]++;
         }
 
         for ( j = 0; j <=nbins; j++ )
             if ( max < bins[j] )
                 max = bins[j];
 
         for ( j = 0; j <=nbins; j++ )
         {
             printf("(%5.2f, %5.2f) |", low + j * delta, low + (j + 1) * delta );
             k = (int)( (double)width * (double)bins[j] / (double)max );
             while(k-- > 0) putchar('#');
               printf("  %-.1f%%", bins[j] * 100.0 / (double)m);
               putchar('\n');
         }
             free(bins);
     }
 }


/***********printHistogram2 prints histogram for values of Rayleigh and Lorentzian Distribution********/
 void printHistogram2(double* values,int n,double high,double low)
  {
       const int width = 50;    
       int max = 0;
    
    
       int i,j,k;
       int nbins = ceil(1+3.322*log10(n));
       double delta =(high-low)/nbins;
       int* bins = (int*)calloc(nbins,sizeof(int));
       if ( bins != NULL )
       {
           for ( i = 0; i < n; i++ )
           {
               int j = (int)( (values[i] - low) / delta );
               if ( 0 <= j  &&  j <=nbins )
                  bins[j]++;
           }
  
           for ( j = 0; j <=nbins; j++ )
               if ( max <= bins[j] )
                  max = bins[j];
 
           for ( j = 0; j <=nbins; j++ )
           {
               printf("(%5.2f, %5.2f) |", low + j * delta, low + (j + 1) * delta );
               k = (int)( (double)width * (double)bins[j] / (double)max );
               while(k-- > 0) putchar('#');
                  printf("  %-.1f%%", bins[j] * 100.0 / (double)n);
                  putchar('\n');
           }
 
               free(bins);
       }
  }
/***********Histogram plots values of Distributions in gnuplot********/
 void Histogram(char distribution[])
  {	
	FILE *pipe = popen("gnuplot -persist", "w");
 	fprintf(pipe, "set terminal pdf\n");
	fprintf(pipe, "set output 'data.pdf'\n");
	fprintf(pipe, "set title 'Plot Of Downsampled Data'\n");
        fprintf(pipe, "set xrange [0:502]\n");
        fprintf(pipe, "set xrange [0:0.1]\n");
        fprintf(pipe, "set xlabel 'X--->'\n"); 
        fprintf(pipe, "set ylabel 'Y--->'\n");
        fprintf(pipe, "plot 'downs.txt' u (hist($1, width)):(1.0) smooth freq w boxes lc rgb 'blue' notitle \n");
        pclose(pipe);
     
  }
int main(int argc ,char *argv[])
 {
    
  if(argc==1)
     printf("Please the type of distribution you want after ./a.out, Normal or Rayleigh or Lorentzian\n");   
  else
  if(argc>1)
  {
      double* seq;
      srand((unsigned int)time(NULL));
    
      int x;
      printf("Enter the number of values to be generted(in general large value):\n");
      scanf("%d",&x);
      double high=__DBL_MIN__,low=__DBL_MAX__;
    
   
      if(strcmp(argv[1],"Normal")==0)
      {
        
         (seq = generatenrm(x,&high,&low))!= NULL;
         printf("mean = %g, stddev =%g\n\n", mean1(seq,x), stddev1(seq,x));
         printHistogram1(seq,x,high,low);
         Histogram("Normal");
         free(seq);
                
        return EXIT_SUCCESS;
      }
  else
      if(strcmp(argv[1],"Rayleigh")==0)
      {
        
         (seq = generaterrm(x,&high,&low))!=NULL;
         printf("mean = %g, stddev = %g\n\n", mean2(seq,x), stddev2(seq,x));
         printHistogram2(seq,x,high,low);
         Histogram("Rayleigh");
         free(seq);
        
         return EXIT_SUCCESS;
      }

  else  
      if(strcmp(argv[1],"Lorentzian")==0)
      {
        
         (seq = generatelrm(x,&high,&low))!=NULL;
         printf("mean = %g, stddev = %g\n\n", mean2(seq,x), stddev2(seq,x));
         printHistogram2(seq,x,high,low);
         Histogram("Lorentzian");
         free(seq);
        
        
         return EXIT_SUCCESS;
    
      }
else printf("please check the entered name\n");
   }

    return EXIT_FAILURE;
 }
