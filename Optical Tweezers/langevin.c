/*************************
Purpose:Solving the Langevin equation and ploting the graps of rms value of position and temperature
Authors:Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Date:31/09/2019
Group   : 2H
Inputs:  Number of steps,Starting Temperature,Ending Temperature,Size of particle in nanometer
Outputs: Plot of position vs time,Plot of rms value of position with temperature
Compilation:cc langevin.c -lm
Command line:./a.out N Tstart Tend D (Preferably ./a.out 10000 370 270 500)
**************************/
#include<stdio.h> 
#include<math.h>
#include<stdlib.h>
#define kb 0.01380 //Bolztmann constant in pN.nm(T)^-1
#define n 890   //The coefficient of viscosity of water in (pN/(nm)^2)*s

//function to return the value of db/dt in Langevin equation
//F is the thermal force
//b is the position
//g is gamma or drag coefficient
//k is kappa
//a is the value of db/dt

float f(float F,float b,float g,float k,float a) 
{
	 a=(F- k*b)/g;
	 return a;
}
void plot(int c,int d,int e)//function to plot in gnuplot
{
        
}
int main(int argc ,char *argv[])//Taking input from the command line
{
         if(argc!=5) //Verifying that user enters te required four arguments after ./a.out
         {
		 printf("Please enter the following in the respective order followed by ./a.out\n"); 
		 printf("1.Number of steps N\n");
		 printf("2.Starting value of temperature Tstart\n");
		 printf("3.Ending value of temperature Tstop\n") ;
		 printf("4.Size of the particle in nanometres\n");
		 printf("Example: ./a.out 10000 270 370 500\n");
         }
         else
         if(argc==5) //finding the values of b after the user enters four arguments after ./a.out 
         {
                 int k;       //variable used in loop
                 int N;       //Number of steps
                 float Tstart;//Starting value of temperature
                 float Tstop; //Ending value of temperature
                 float D;     //Size of the particle in nanometer
 
		 N=atoi(argv[1]);     //Taking N from the command line and converting to an integer  
                 Tstart=atof(argv[2]);//Taking Tstart from the command line and converting to a float value    
		 Tstop=atof(argv[3]); //Taking Tstop from the command line and converting to a float value
                 D=atof(argv[4]);     //Taking D from the command line and converting to a float value   
  
	
		 float dt=0.01;       //dt is stepsize and taking it as 0.001
		 float F,W;            //F is the thermal force and is generated randomly between -W to W
		 float kappa=0.1;      //kappa is a constant used in Langevein equation
		 float g;              //g is gamma and is a constant in Langevein equation
		 float b[N],t[N];//arrays to store b vales,t(time) values,u(slope) values
		 float k1b,k2b,k3b,k4b,k5b,k6b;//These are variables used in RK45 method
		 float sumsquare,rms; //sumsquare is the sum of squares and rms is the rms value of distance in nanometer


		 FILE *ft1=fopen("data.txt","w"); //Opening file to store rms values and temperature
		 FILE *ft2=fopen("data1.txt","w");//Opening file to store b and t values for one particular temperature
		 
                 g=3*(M_PI)*n*D; //Equating gamma to 3*PI*n*D
		 b[0]=0;//Taking initial value of b as 1
		 t[0]=0;//taking initial value of time as 0
         
                 int T;//Loop variable that corresponds to the value of temperature in each loop
                 int Tmean;//Calculating the mean value of temperature
                 Tmean=(int)((int)Tstart+(int)Tstop)/2;
                 for(T=(int)Tstart;T<(int)Tstop;T++)//Rounding off both Tstart and Tstop to integer values
                 {
                 	W=sqrt((2*kb*T*g));//Taking W a 2*Kb*T*(Gamma)
	                sumsquare=0;                         
	                for(k=1;k<N;k++)//using rk45 method to solve the equations
		        {
                                
               			 F=(rand()%(2*(int)W+1)) - (int)W;//Generating F randomly between -W to W;
		                 k1b=f(F,(b[k-1]),g,kappa,k1b)*dt;
				 k2b=f(F,(b[k-1]+(k1b/4)),g,kappa,k2b)*dt;
				 k3b=f(F,(b[k-1]+(3*k1b/32)+(9*k2b/32)),g,kappa,k3b)*dt;	
				 k4b=f(F,(b[k-1]+(1932*k1b/2197)-(7200*k2b/2197)+(7296*k3b/2197)),g,kappa,k4b)*dt;
				 k5b=f(F,(b[k-1]+(439*k1b/216)-(8*k2b)+(3680*k3b/513)-(845*k4b/4104)),g,kappa,k5b)*dt;
				 k6b=f(F,(b[k-1]-(8*k1b/27)+(2*k2b)-(3544*k3b/2565)+(1859*k4b/4104)-(11*k5b/40)),g,kappa,k6b)*dt;

				 b[k]= b[k-1] + ((16*k1b/135)+(6656*k3b/12825) + (28561*k4b/56430)-(9*k5b/50)+(2*k6b/55));
				 t[k]= k*dt;
				
                                 if(T==Tmean) //Writing values of b and t for mean Temperature
                                 {
                			 fprintf(ft2,"%f %f\n",t[k],b[k]); 
                                 }   
	                }
			for(k=0;k<N;k++)
		        {
		         	sumsquare+=(b[k]*b[k]);//Finding the sum of squares of b values
		  
		        }
	 
		        rms=sqrt(sumsquare/N);//Generating the rms value of b
		        fprintf(ft1,"%f %d\n",rms,T);//Writing rms values and temperatures in a file
	        }
                fclose(ft1);
		fclose(ft2);
                FILE *pipe = popen("gnuplot -persist", "w");//Calling gnuplot from the c program
 		fprintf(pipe, "set terminal pdf\n");
		fprintf(pipe, "set output 'plot.pdf'\n");
		fprintf(pipe, "set title 'Plot of rms values (in nm) vs Temperatue'\n");
        	fprintf(pipe, "set xlabel 'Temparature'\n"); 
      	  	fprintf(pipe, "set ylabel 'rms values'\n");
                fprintf(pipe, "plot 'data.txt' u 2:1 w l title 'Plot of rms data'\n");
        	fprintf(pipe, "set title 'Plot of b values vs Time'\n");
     		fprintf(pipe, "set xlabel 'Time'\n"); 
        	fprintf(pipe, "set ylabel 'b (in nm) values'\n");
        	fprintf(pipe, "plot 'data1.txt' u 1:2 w l title 'N= %d,D= %d T= %d'\n",N,(int)D,Tmean);
        	pclose(pipe);
                printf("The plots are saved as plot.pdf\n");

         }
 
}
