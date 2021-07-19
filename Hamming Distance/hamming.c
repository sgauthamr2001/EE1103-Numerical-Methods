/*************************
Purpose:Comparing hamming distance
Authors:Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Date:12/09/2019
Group   : 2H
Inputs:  Number of bits transmitted,number of bits received,number of bits revealed,number of bits flipped
Outputs: The offset generated after comparision which gives minimum hamming distance and initial offset
Compilation:cc hamming.c 
Command line:./a.out N P Q (Preferably ./a.out 100000 10000 2500)
**************************/
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>

int * rndgen(int l,int u,int a[],int p)//The following function generates p distinct random numbers between u and l by using knuth algorithm
 {

	int in, im;

	im = 0;

	for (in = 0; in < (u-l+1) && im < p; ++in) 
	 {
	 	 int rn = (u-l+1) - in;
	 	 int rm = p - im;
	 	 if (rand() % rn < rm)    
	    	 	  a[im++] = in + l; 
	 }	
	assert(im == p);
	return a;
 }

int main(int argc,char** argv)
 {
         if(argc<4) //makes sure that the user enters the arguments after ./a.out
          {
	        printf("Please enter number of bits to be transmitted, number of bits to be received, number of bits to be flipped followed by ./a.out. Example: ./a.out 100000 10000 2500\n");
          }

         else
	  {
		srand(time(0));//seeding random number generator            
       		int i,N,M,P,Q,j;//Several variables
		N=atoi(argv[1]);//Taking number of bits transmitted from second argument of command line
       		P=atoi(argv[2]);//Taking number of bits received from third argument of command line
        	Q=atoi(argv[3]);//Taking number of bits to be flipped form fourth argument of command line
		int *ptr;//pointer to store N random bits
        	int *ptrn;//pointer to store P bits
 		int a[P];//a is an array which stores locations with respect to random offset M
		int c[P];//c is an array which stores locations with respect to 0
		int b[Q];//b stores Q locations from P  locations
		int t[P];//t stores the the loactions of bits that are revealed
		ptr=calloc(N,sizeof(int));
		int hamdist[N];//array to stor hamming distance values after each iteration.
	        int x;//x is the number of bits to be revealed
		printf("Please Enter the number of points to be revealed to compare the hamming distance\n");
		scanf("%d",&x);
   		printf("Please wait ....\n");
		for(i=0;i<N;i++)//generating N random bits
		 {
			*(ptr+i)=rand()%2;
                
     		 }
       			 M=rand()%((N/10)+1)+1;//Generating random offset between 1 and N/10
        
 		ptrn=calloc(P,sizeof(int));

                rndgen(M+1,N,a,P);//calls the fuction to generate a values
   
   		a[0]=M;
	
        
       		for(i=0;i<P;i++) //generating c values
		 {
                 	c[i]=a[i]-M;
	 	 }
	 	int cmax=c[0];//finding maximum value of c
      	        for(i=0; i<P; i++)  		
		 {
			
                	  if(c[i]>cmax)
		 	   {
				cmax=c[i];
			
       	            	   }
	
		 }
		for(i=0;i<P;i++)//storing P values into pointer ptrn
		 {
                 
			*(ptrn+i)=*(ptr+M+c[i]);
	        
        	 }

	        rndgen(0,P-1,b,Q); //Generating b values 

		for(j=0;j<Q;j++)//Flipping Q values in ptrn
		 {
                 
			*(ptrn+b[j])=1-*(ptrn+b[j]);
	          
        	 }
  		FILE *ft=fopen("bit.txt","w+"); //opening file

		rndgen(0,P,t,x);//generating t values
	
		for(int h=0; h<=N-cmax; h++) //comparing hamming distance for various offsets		
		{ 
			hamdist[h]=0;
		

			for(i=0; i<x;i++)
			 {
		
				if(*(ptrn + t[i])^*(ptr+h+c[t[i]])==1)
					hamdist[h]++;
			 }
			fprintf(ft,"%d %d\n",h,hamdist[h]);
 	
		}
       		int hammin;
        	int min=hamdist[0]; //finding minimum offset
		for(int h=0; h<=N-cmax; h++) 		
		 {
		
	
                 	 if(hamdist[h]<min)
		 	  {
				min=hamdist[h];
				hammin=h;
                          }
	
	         } 

                FILE *pipe = popen("gnuplot -persist", "w");              
			fprintf(pipe, "set terminal pdf\n");
			fprintf(pipe, "set output 'plot.pdf'\n");
			fprintf(pipe, "set title 'Plot of Hamming distance Vs Offsets'\n");
			fprintf(pipe, "set xlabel 'Offset'\n");
                        fprintf(pipe, "set ylabel 'Hamming distance'\n");
			fprintf(pipe, "plot 'bit.txt' u 1:2 w l title 'Plot'\n");
                
       		printf("The final report is :\n");
     	 	printf("1.Number of bits transmitted: %d\n",N);
		printf("2.Number of bits received: %d\n",P);
		printf("3.Number of bits flipped: %d\n",Q);
		printf("4.Number of bits revealed: %d\n",x);
		printf("5.The random offset generated at start is: %d\n",M);
		printf("6.The offset obtained after final comparision: %d\n",hammin);
                printf("7.The coressponding hamming distance is: %d\n",min);
		return 0;
		free(ptrn);
		free(ptr);
		fclose(ft);
		pclose(pipe);
         }
 }
