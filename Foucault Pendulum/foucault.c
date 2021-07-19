/*************************
Purpose:To solve the equations of foucault pendulum
Authors:Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Date:23/10/2019
Group   : 2H
Inputs:constants of the equations of foucault pendulum,intial coordinates and velocities,stepsize and number of points to be generated, the method to solve differential equations.
Outputs:Plots of foucault pendulum
Compilation:cc foucault.c -lm
**************************/

#include<stdio.h>
#include<math.h>
int euler(int n,double a,double b,double c,double d,double e,double f,double step)//euler method to solve differential equations
{
        int k;   
        double x[n],y[n],t[n],r[n]; //arrays to store x,y,r,t values
        double u[n];//arrays to store dx/dt values
        double v[n];//arrays to store dy/dt values
   
	x[0]=c;
	y[0]=d; 
	u[0]=e;
        v[0]=f;
	
	
	//Using euler's method to find solution
	for(k=1;k<n;k++)
	 {     
              x[k]=(u[k-1]*step) + x[k-1];
              y[k]=(v[k-1]*step )+ y[k-1];
	      u[k]=((a*v[k-1]-b*x[k-1])*step) + u[k-1];
	      v[k]=(((-a)*u[k-1]-b*y[k-1])*step) + v[k-1];
	      
         }


        for(k=0;k<n;k++) 
	{
	      r[k]=sqrt((x[k]*x[k]) + (y[k]*y[k]));
              t[k]=k*step;
        }
            


 	FILE *fp = fopen("data.txt", "w");
	//writing into a file
	for(k=0;k<n;k++)
	 {
		fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",x[k],y[k],u[k],v[k],r[k],t[k]);
         }
        fclose(fp);
        
}
int heuns(int n,double a,double b,double c,double d,double e,double f,double step)//huens method to solve differential equations
{
        int k;
	float x0[n],y0[n],x[n],y[n],r[n],t[n];//arrays to store different values
	float u[n];
	float v[n];
	
	x0[0]=c;
	y0[0]=d;
        x[0]=c;
	y[0]=d;	
	u[0]=e;
	v[0]=f;

        //Using heun's method to find solution
	for(k=1;k<n;k++)
	{
		x0[k]=(u[k-1]*step)+x0[k-1];
		u[k]=((a*v[k-1]-b*x0[k-1])*step)+u[k-1];
		x[k]=x0[k-1]+(step/2)*(u[k-1]+u[k]);
		y0[k]=(v[k-1]*step)+y0[k-1];
		v[k]=(((-a)*u[k-1]-b*y0[k-1])*step)+v[k-1];
		y[k]=y0[k-1]+(step/2)*(v[k-1]+v[k]);
	}	
	 for(k=0;k<n;k++)
	{
	      r[k]=sqrt((x[k]*x[k]) + (y[k]*y[k]));
              t[k]=k*step;
        }
            


 	FILE *fp = fopen("data.txt", "w");
	//writing into a file
	for(k=0;k<n;k++)
	 {
		fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",x[k],y[k],u[k],v[k],r[k],t[k]);
         }
	 fclose(fp);
}
int rk(int n,double a,double b,double c,double d,double e,double f,double dt)//runge-kutta method to solve differential equations
{

	int k;

	double x[n],y[n],u[n],v[n],r[n],t[n];//arrays to store different values
	double k1x,k2x,k3x,k4x;//coefficients in rk method
	double k1y,k2y,k3y,k4y;
	double k1u,k2u,k3u,k4u;
	double k1v,k2v,k3v,k4v;
	

	x[0]=c;//initial values 
	y[0]=d;
	u[0]=e;
	v[0]=f;


        //Using runge-kutta method to find solution
	for(k=1;k<n;k++)
	 {

		k1x=u[k-1]*dt;
		k1y=v[k-1]*dt;
		k1u=(a*v[k-1]-b*x[k-1])*dt;
		k1v=(-a*u[k-1]-b*y[k-1])*dt;

		k2x=(u[k-1]+(k1u/2))*dt;
		k2y=(v[k-1]+(k1v/2))*dt;
		k2u=(a*(v[k-1]+(k1v/2))-b*(x[k-1]+(k1x/2)))*dt;
		k2v=(-a*(u[k-1]+(k1u/2))-b*(y[k-1]+(k1y/2)))*dt;

		k3x=(u[k-1]+(k2u/2))*dt;
		k3y=(v[k-1]+(k2v/2))*dt;
		k3u=(a*(v[k-1]+(k2v/2))-b*(x[k-1]+(k2x/2)))*dt;
		k3v=(-a*(u[k-1]+(k2u/2))-b*(y[k-1]+(k2y/2)))*dt;

		k4x=(u[k-1]+ k3u)*dt;
		k4y=(v[k-1]+ k3v)*dt;
		k4u=(a*(v[k-1]+ k3v)-b*(x[k-1]+k3x))*dt;
		k4v=(-a*(u[k-1]+ k3u)-b*(y[k-1]+ k3y))*dt;

		x[k]= x[k-1] + (k1x+ 2*k2x +2*k3x + k4x)/6;
		y[k]= y[k-1] + (k1y+ 2*k2y +2*k3y + k4y)/6;
		u[k]= u[k-1] + (k1u+ 2*k2u +2*k3u + k4u)/6;
		v[k]= v[k-1] + (k1v+ 2*k2v +2*k3v + k4v)/6;
	 }
 	
	 for(k=0;k<n;k++)
	 {
	        r[k]=sqrt((x[k]*x[k]) + (y[k]*y[k]));
                t[k]=k*dt;
         }
            


 	FILE *fp = fopen("data.txt", "w");
	//writing into a file
	for(k=0;k<n;k++)
	 {
		fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",x[k],y[k],u[k],v[k],r[k],t[k]);
         }
         fclose(fp);
}

int rk45(int n,double a,double b,double c,double d,double e,double f,double dt)//rk45 method to solve differential equations
{

	int k;

	double x[n],y[n],u[n],v[n],r[n],t[n];//arrays to store different values
	double k1x,k2x,k3x,k4x,k5x,k6x;
	double k1y,k2y,k3y,k4y,k5y,k6y;
	double k1u,k2u,k3u,k4u,k5u,k6u;
	double k1v,k2v,k3v,k4v,k5v,k6v;


	x[0]=c;//initial values
	y[0]=d;
	u[0]=e;
	v[0]=f;

	//Using rk45 method to find solution
	for(k=1;k<n;k++)
	 {

		k1x=u[k-1]*dt;
		k1y=v[k-1]*dt;
		k1u=(a*v[k-1]-b*x[k-1])*dt;
		k1v=(-a*u[k-1]-b*y[k-1])*dt;

		k2x=(u[k-1]+(k1u/4))*dt;
		k2y=(v[k-1]+(k1v/4))*dt;
		k2u=(a*(v[k-1]+(k1v/4))-b*(x[k-1]+(k1x/4)))*dt;
		k2v=(-a*(u[k-1]+(k1u/4))-b*(y[k-1]+(k1y/4)))*dt;

		k3x=(u[k-1]+(3*k1u/32)+(9*k2u/32))*dt;
		k3y=(v[k-1]+(3*k1v/32)+(9*k2v/32))*dt;
		k3u=(a*(v[k-1]+(3*k1v/32)+(9*k2v/32))-b*(x[k-1]+(3*k1x/32)+(9*k2x/32)))*dt;
		k3v=(-a*(u[k-1]+(3*k1u/32)+(9*k2u/32))-b*(y[k-1]+(3*k1y/32)+(9*k2y/32)))*dt;

		k4x=(u[k-1]+(1932*k1u/2197)-(7200*k2u/2197)+(7296*k3u/2197))*dt;
		k4y=(v[k-1]+(1932*k1v/2197)-(7200*k2v/2197)+(7296*k3v/2197))*dt;
		k4u=(a*(v[k-1]+(1932*k1v/2197)-(7200*k2v/2197)+(7296*k3v/2197))-b*(x[k-1]+(1932*k1x/2197)-(7200*k2x/2197)+(7296*k3x/2197)))*dt;
		k4v=(-a*(u[k-1]+(1932*k1u/2197)-(7200*k2u/2197)+(7296*k3u/2197))-b*(y[k-1]+(1932*k1y/2197)-(7200*k2y/2197)+(7296*k3y/2197)))*dt;

		k5x=(u[k-1]+(439*k1u/216)-(8*k2u)+(3680*k3u/513)-(845*k4u/4104))*dt;
		k5y=(v[k-1]+(439*k1v/216)-(8*k2v)+(3680*k3v/513)-(845*k4v/4104))*dt;
		k5u=(a*(v[k-1]+(439*k1v/216)-(8*k2v)+(3680*k3v/513)-(845*k4v/4104))-b*(x[k-1]+(439*k1x/216)-(8*k2x)+(3680*k3x/513)-(845*k4x/4104)))*dt;
		k5v=(-a*(u[k-1]+(439*k1u/216)-(8*k2u)+(3680*k3u/513)-(845*k4u/4104))-b*(y[k-1]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4104)))*dt;

		k6x=(u[k-1]-(8*k1u/27)+(2*k2u)-(3544*k3u/2565)+(1859*k4u/4104)-(11*k5u/40))*dt;
		k6y=(v[k-1]-(8*k1v/27)+(2*k2v)-(3544*k3v/2565)+(1859*k4v/4104)-(11*k5v/40))*dt;
		k6u=(a*(v[k-1]+(439*k1v/216)-(8*k2v)+(3680*k3v/513)-(845*k4v/4101))-b*(x[k-1]-(8*k1x/27)+(2*k2x)-(3544*k3x/2565)+(1859*k4x/4104)-(11*k5x/40)))*dt;
		k6v=(-a*(u[k-1]-(8*k1u/27)+(2*k2u)-(3544*k3u/2565)+(1859*k4u/4104)-(11*k5u/40))-b*(y[k-1]+(439*k1y/216)-(8*k2y)+(3680*k3y/513)-(845*k4y/4101)))*dt;

		x[k]= x[k-1] + ((16*k1x/135)+(6656*k3x/12825) + (28561*k4x/56430)-(9*k5x/50)+(2*k6x/55));
		y[k]= y[k-1] + ((16*k1y/135)+(6656*k3y/12825) + (28561*k4y/56430)-(9*k5y/50)+(2*k6y/55));
		u[k]= u[k-1] + ((16*k1u/135)+(6656*k3u/12825) + (28561*k4u/56430)-(9*k5u/50)+(2*k6u/55));
		v[k]= v[k-1] + ((16*k1v/135)+(6656*k3v/12825) + (28561*k4v/56430)-(9*k5v/50)+(2*k6v/55));
	    }
	 
	   for(k=0;k<n;k++)
            {
		      r[k]=sqrt((x[k]*x[k]) + (y[k]*y[k]));
		      t[k]=k*dt;
	    }
		    


	   FILE *fp = fopen("data.txt", "w");
	   //writing into a file
	   for(k=0;k<n;k++)
	    {
	              fprintf(fp,"%lf %lf %lf %lf %lf %lf\n",x[k],y[k],u[k],v[k],r[k],t[k]);
            }
	    fclose(fp);
}
void plot(int n,double dt,char method[])//function to plot in gnuplot
 {
        FILE *pipe = popen("gnuplot -persist", "w");
 	fprintf(pipe, "set terminal pdf\n");
	fprintf(pipe, "set output 'plot.pdf'\n");
	fprintf(pipe, "set title 'Plot of x vs y for %s method'\n",method);
        fprintf(pipe, "set xlabel 'X'\n"); 
        fprintf(pipe, "set ylabel 'Y'\n");
        fprintf(pipe, "plot 'data.txt' w l title ' number of points = %d,stepsize = %.3lf \n",n,dt);
        pclose(pipe);
 }
int main()
{
	int num,n;
	double a,b,c,d,e,f,step;//inputs to functions
	
	printf("Enter the number of points to be generated (preferably 1,00,000)\n");
	scanf("%d",&n);
	printf("Enter a,b in the equation : d2y/dt2 = a*dx/dt - b*y preferably (a=0.256,b=1.452)\n");
	scanf("%lf %lf",&a,&b);
	printf("Enter the step size (preferably 0.001-0.01)\n");
	scanf("%lf",&step);
	printf("Enter initial x,y values and initial velocity in x direction and y direction respectively(preferably 0.67 0.67 0 0)\n");
	scanf("%lf %lf %lf %lf",&c,&d,&e,&f);
	printf("\n>Enter 1 if you want to use euler's method\n");
	printf(">Enter 2 if you want to use heun's method\n");
	printf(">Enter 3 if you want to use runge kutta method \n");
	printf(">Enter 4 if you want to use rk45 method\n");
	scanf("%d",&num);
	    switch (num)//switch is used here to use a particular method to solve the equations
	    {
		case 1:
			euler(n,a,b,c,d,e,f,step);
			plot(n,step,"euler");
                        printf("The graph is saved as plot.pdf\n");
		break;
		case 2:
			heuns(n,a,b,c,d,e,f,step);
			plot(n,step,"heuns");
			printf("The graph is saved as plot.pdf\n");
		break;
		case 3:
			rk(n,a,b,c,d,e,f,step);
			plot(n,step,"runge-kutta");
			printf("The graph is saved as plot.pdf\n");
		break;
		case 4:
			rk45(n,a,b,c,d,e,f,step);
			plot(n,step,"rk45");
			printf("The graph is saved as plot.pdf\n");
		break;
		default:
		        printf("Please enter a number valid number\n");
		break;
	    }
}

//gnuplot command for first graph in report
// plot 'data.txt'   u ($1):($2):($3/1000):($4/1000):4 with vectors head size 0.01,20,60 filled lc palette 






