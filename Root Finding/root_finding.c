/*************************
Purpose:To generate plots of logistic map and predator prey model
Authors:Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Date:26/09/2019
Group   : 2H
Inputs:Values of number of iterations, step size ,several variables
Outputs:Plots of logistic map and predator prey model
Compilation:cc assignment5.c -lm
**************************/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main()

{
        //variables of predator prey model
        int num;  //num is a variable to switch value
        int k;   
        int n;    //Number of iterations of predator-prey model
  
        printf("Please enter number of iterations for predator prey model preferably 10000\n");
        scanf("%d",&n);
        

	float a,b,c,d,e,f,step; //a,b,c,d,e,f are constants in differential equations,intial values of predator and prey,step size of time 
	float x[n],y[n],t[n];//arrays to store x values(prey),y values(predator),time values 
        float vx[n];   //arrays to store dx values
        float vy[n];   //arrays to store dy values
        
        //variables of logistic map                  
        double rs;   //start value of parameter r
	double re;   //end value of parameter r
	double de;   //delta r
	double rc;   //parameter r
	double stx = 0.5; //other variables
	double xc;
	unsigned int g, h;
	unsigned int isdone;
	unsigned int match;
	unsigned int maxrep = 100;
	printf("please enter the step size  for logistic map preferably 0.001\n");
        scanf("%lf",&de);

        printf(">enter 1 if you want Plots of Lotka-Volterra model\n>enter 2 if you want Logistic map\n");
        scanf("%d",&num);
 switch(num)
  {
     case 1:
 

	printf("dx/dt=ax-bxy ; dy/dt=cxy -dy ; x corresponds to prey , y coressponds to predator\n");
	printf("Enter the values of a,b,c,d,Initial value of x,Intial value of y with spaces in the above order;Preferably 1 2 3 4 0.2 0.3 0.01\n");

	scanf("%f %f %f %f %f %f %f",&a,&b,&c,&d,&e,&f,&step);


	x[0]=e;
	y[0]=f; 
	t[0]=0;
	vx[0]=0;
	vy[0]=0;
	//Using euler's method to find solution
	for(k=1;k<n;k++)
	 {
		x[k] = (a*x[k-1] - (b*x[k-1]*y[k-1]))*step+ x[k-1]; 
		y[k] = ((d*x[k-1]*y[k-1]) - c*y[k-1])*step + y[k-1];
                t[k]+=k*step;
     		vx[k] = (a*x[k-1] - (b*x[k-1]*y[k-1]))*step;
		vy[k] = ((d*x[k-1]*y[k-1]) - c*y[k-1])*step;
         }
 	FILE *fp = fopen("data.txt", "w");
	//writing into a file
	for(k=0;k<n;k++)
	 {
		fprintf(fp,"%f %f %f %f %f\n",t[k],x[k],y[k],vx[k],vy[k]);
         }
        fclose(fp);
 	//plotting in gnuplot
        FILE *pipe = popen("gnuplot -persist", "w");
        fprintf(pipe, "set terminal pdf\n");
	fprintf(pipe, "set output 'plots.pdf'\n");
        fprintf(pipe, "set title 'Plot of variation of predator and prey with time'\n");
        fprintf(pipe,"set xlabel 'Time'\n");
        fprintf(pipe,"set ylabel 'Number of species'\n");
        fprintf(pipe, "plot 'data.txt'  u 1:2 w l title 'Prey vs time' ,'data.txt' u 1:3 w l  title 'Predator vs Time' \n");
        fprintf(pipe, "set title 'Plot of  predator vs prey'\n");
        fprintf(pipe,"set xlabel 'Prey'\n");
        fprintf(pipe,"set ylabel 'Predator'\n");
        fprintf(pipe, "plot 'data.txt'   u ($2):($3):($4):($5):5 with vectors head size 0.1,20,60 filled lc palette title 'Predator vs Prey' ");
        
        pclose(pipe);
    break ;

    case 2:
       
	printf("please enter the Starting value of r,ending value of r with spaces for logistic map preferably 1 4 \n");
      	scanf("%lf %lf",&rs,&re);
	
        FILE *ofp_dat;
	
	
	if(rs > re)
         {
		rs = re;
	 }
	if(re > 4.0)
         {
		re = 4.0;
	 }

	double *states = calloc(maxrep,sizeof(double));
	
	rc = rs;
	
	//Creating blank file
	ofp_dat = fopen("res.dat","w");
	fclose(ofp_dat);
		
	//Start of main simulation
	while(rc <= re)
        {
	 	
		//steady-state 
		xc = stx;
		for(g = 0; g < 10000; g++)
                 {
			xc = rc*xc*(1.0-xc);
		 }
		
		// checking oscillations
		isdone = 0;
		g = 0;
		while(isdone != 1)
                 {
			xc = rc*xc*(1.0-xc);

			//Checking if oscillatory point reached
			if(g == 0)
                         {
				*(states + g) = xc;
			 } 
                        else if(g > 0)
                               {
				  match = 0;
				  for(h = 0; h < g; h++)
                                   {
					if(*(states + h) == xc)
                                         {
						match = h;
						isdone = 1;
						break;
					 }
			           }
			          if(isdone != 1)
                                   {
					*(states + g) = xc;
				   }
			       }
			g++;
			if(g == maxrep)
                         {
				isdone = 1;
			 }
		 }
		
		//Write out oscillatory points to file
		g--;
		ofp_dat = fopen("res.dat","a+");
		for(h = match; h < g; h++)
                 {
			fprintf(ofp_dat,"%0.10lf %0.10lf\n",rc,*(states + h));
		 }
		fclose(ofp_dat);
		
		rc += de;
	 }
	 free(states);

         //Plotting in gnuplot
         FILE *pipe1 = popen("gnuplot -persist", "w");
         fprintf(pipe1, "set terminal pdf\n");
	 fprintf(pipe1, "set output 'plots1.pdf'\n");
         fprintf(pipe1, "set title  'Logistic map'\n");
         fprintf(pipe1, "set xlabel 'Parameter r'\n");
         fprintf(pipe1, "set ylabel 'x values'\n");
         fprintf(pipe1,"plot 'res.dat' u 1:2 pt 7 ps 0.017\n");
        
         pclose(pipe1);
	 return 0;
      
     break ;
  }
 
}
