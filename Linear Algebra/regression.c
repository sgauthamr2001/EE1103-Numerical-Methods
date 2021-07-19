/*******************************************************************************************
Purpose :To obtain a function using linear regression
Authors : Muhammed Hisham(EE19B098),Sai Gautham Ravipati(EE19B053)
Group   : 2H
Date    :19/09/2019
Inputs  :None
Outputs :The value of surface area obtained at Height 187cm and Weight 78kg
Compilation: cc regression.c -lm
*******************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define e M_E

void LUdecomposer(int m, float A[m][m],float B[m], float X[m])//Function to perform LU decomposition
 {
       float L[m][m], U[m][m];
       int i,k,j;
       float Y[m];
       for(j=0; j<m; j++)           //Generation of both lower and upper diagonal matrices         
       {
           for(i=0; i<m; i++)
           {
              if(i<=j)
              {
                   U[i][j]=A[i][j];
                   for(k=0; k<i-1; k++)
                        U[i][j]-=L[i][k]*U[k][j];
                   if(i==j)
                        L[i][j]=1;
                   else
                        L[i][j]=0;
              }
            else
              {
                   L[i][j]=A[i][j];
                   for(k=0; k<=j-1; k++)
                        L[i][j]-=L[i][k]*U[k][j];
                        L[i][j]/=U[j][j];
                        U[i][j]=0;
              }
           }
       }
       for(i=0; i<m; i++)              //Generation of intermediate variable Y
       {
            Y[i]=B[i];
            for(j=0; j<i; j++)
            {
                 Y[i]-=L[i][j]*Y[j];
            }
       }
       for(i=m-1; i>=0; i--)            //Gives the values of variables involved
       {
            X[i]= Y[i];
            for(j=i+1; j<m; j++)
            {
                 X[i]-=U[i][j]*X[j];
            }
            X[i]/=U[i][i];
       }
    
   
 }
int main()
 {
      int h[9]; //array to store values of heights
      int w[9]; //array to store values of weights
      float a[9];//array to store values of area
      /*********************
       sosqh:Variable of sum of squares of h
       sosqw:Variable of sum of squares of w
       sh:Variable of sum of h
       sw:Variable of sum of w
       ps:Variable of sum of product of h and w
       sfa:Variable of sum of a values
       spaw:Variable of sum of product of w and a
       spah:Variable of sum of product of h and a     
       ******************/
 
       int sosqh=0,sosqw=0,sh=0,sw=0,ps=0;
       float sfa=0.00,spaw=0.00,spah=0.00;
       int i;
       float A[3][3],B[3];//Arrays that store values for LU decomposition
       float X[3];        //Values obtained from LU decomposition
       float func;        //Variable that gives value of surface area at H=187 cm and W =87 Kg

       FILE *file = fopen("data.txt","r"); //data.txt is a file with Heights in column 1,Weights in column 2,Areas in Column 3
          
       if(file==NULL)
       {
           printf("\n error!!!");
         return -1;
       }
       for(i=0;i<9;i++)
       {
           fscanf(file,"%d %d %f",&h[i],&w[i],&a[i])!=EOF;   //reading the file to obtain values of h,w,a    
       }
       fclose(file);
    
    
       for(i=0;i<9;i++)                 //evaluvating each variable
       {
	    sh+=h[i];
            sw+=w[i];
	    sosqh+=h[i]*h[i];
	    sosqw+=w[i]*w[i];
	    ps+=h[i]*w[i];
	    spah+=a[i]*h[i];
	    spaw+=a[i]*w[i];
	    sfa+=a[i];
       }
       //assigning each value of matrix to perform LU decomposition
       A[0][0]=9,A[0][1]=sh,A[0][2]=sw,A[1][0]=sh,A[1][1]=sosqh,A[1][2]=ps,A[2][0]=sw,A[2][1]=ps,A[2][2]=sosqw;
       B[0]=sfa,B[1]=spah,B[2]=spaw;  
       LUdecomposer(3,A,B,X);//performing LU decomposition
       func=X[0]+X[1]*187+X[2]*78;//gives value of surface area at H=187 cm and W =87 Kg

       //prints area as a function of height and weight
       printf("The function is A=%lf+%lf(H)+%lf(W)\n",X[0],X[1],X[2]);
       //prints value of surface area at H=187 cm and W =87 Kg
       printf("The value of surface area obtained at Height 187cm and Weight 78kg is %lf\n",func);
       return 0;    
 }
