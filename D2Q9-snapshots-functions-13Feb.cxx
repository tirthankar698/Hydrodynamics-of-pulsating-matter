#include<iostream>
#include<math.h>
#include<cmath>
#include<fstream>
#include<stdio.h>
#include <cstdio>
#include<stdlib.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include "pulsating.h"

//Defining the Laplacian

double Lap(double** array, int rows, int cols, int ii, int jj, double dxx) {

    // Laplacian with periodic boundary conditions
    return (4*(array[(ii+1)% rows][jj]+array[(ii+ rows-1)% rows][jj] + array[ii][(jj+1)% cols]+array[ii][(jj+cols-1)%cols])-20*array[ii][jj] +(1)*(array[(ii+1)%rows][(jj+1)%cols] + array[(ii-1+rows)%rows][(jj+1)%cols] + array[(ii+1)%rows][(jj-1+cols)%cols] + array[(ii-1+rows)%rows][(jj-1+cols)%cols]))/(6*dxx);
}

//Defining the divergence

double Div(double** arrayx, double** arrayy, int rows, int cols, int ii, int jj, double dxx) {

 // Divergence with periodic boundary conditions
   return (4*(arrayx[(ii+1)%rows][jj] - arrayx[(ii-1+rows)%rows][jj] + arrayy[ii][(jj+1)%cols] - arrayy[ii][(jj+cols-1)%cols]) +(arrayx[(ii+1)%rows][(jj+1)%cols] + arrayy[(ii+1)%rows][(jj+1)%cols] + arrayx[(ii+1)%rows][(jj+cols-1)%cols] + arrayy[(ii+rows-1)%rows][(jj+1)%cols] - arrayx[(ii-1+rows)%rows][(jj+1)%cols] - arrayx[(ii+rows-1)%rows][(jj+cols-1)%cols] - arrayy[(ii+rows-1)%rows][(jj+cols-1)%cols] - arrayy[(ii+1)%rows][(jj+cols-1)%cols]))/(12*dxx);
}


//main function starts here

int main()
{
long idum;
long sd1=-9781010,sd2=-9876010,sd3=-9013564,sd4=-9964431,sd5=-9057876,sd6=-9166654,sd7=-9277345,sd8=-9383314,sd9=-9493635,sd10=-9667324,sd11=-9554321,sd12=-9887654;

double Lx=100;  // extent in x direction
double Ly=100;  // extent in y direction

double dx=0.25;
double dx2=dx*dx;

int L= int (Lx/dx);  // extent in x direction
int L1=int (Ly/dx);  //extent in y direction

//int midx=0.5*L-1;
// int midy=0.5*L-1;
//int P=70;
//long int umax=10;
//long int nmax=50;


double** rho= new double*[L];
double** rhonew = new double*[L];
double** psi= new double*[L];
double** psinew= new double*[L];
double** sinpsi= new double*[L];
double** cospsi= new double*[L];
double** laprho=new double*[L];
//double** gradxpsi=new double*[L];
//double** gradypsi=new double*[L];
double** lappsi= new double*[L];
double** lapcospsi= new double*[L];
double** nlpsi=new double*[L];
double** etaxrho= new double*[L] ;
double** etayrho= new double*[L];
double** om= new double*[L];
double** divetarho= new double*[L];

 for (int i = 0; i < L; i++) {
        rho[i] = new double[L1];
        rhonew[i] = new double[L1];
        psi[i] = new double[L1];
         psinew[i] = new double[L1];
        sinpsi[i] = new double[L1];
        cospsi[i] = new double[L1];
        laprho[i] = new double[L1];
       // gradxpsi[i] = new double[L1];
       // gradypsi[i] = new double[L1];
       lapcospsi[i] = new double[L1];
        lappsi[i] = new double[L1];
        nlpsi[i] = new double[L1];
        etaxrho[i]= new double[L1];
        etayrho[i]= new double[L1];
        om[i]= new double[L1];
        divetarho[i]= new double[L1];
    }


double tPhi, Phi;
double* tc= new double[21000];

double pi=4*atan(1);

//double avgrho[L][L1]={};
//double avgpsi[L][L1]={};

double rho0=1;
double rhostar=1.04;
double kappa=1.0;
double mu=1.0;
double ep=0.6;
double lambda=1.0;
double a=1.0;
double gamma=1.0;
double Dr;
Dr=mu*kappa;
double beta;
double betap;
double alphap;
alphap=lambda/((gamma*rho0)*(Dr));
beta=(mu*rho0*lambda)/(gamma);
double Drho=0.01;
double Dpsi=Drho;
double omega=1.0/1.0;
betap=beta/omega;
double omegalarge=omega;
double rholess=rhostar;
double trelax=1e2;
double g=0.0;
double lc;
lc=sqrt(Dr/omega);
double ps0;
//double deno=mu*a*ep*rho0*rho0;
//double jj=-2*omega/(deno);
//double psi0=asin(jj);
//double u=0.001;

double dt=0.001;

//double dt=0.1*dx2/alphap;
if(alphap>1)
{dt=0.1*dx2/alphap;}
else if(alphap<1)
{dt=0.1*dx2;}

int hmax=10;

int count1=0;
double count=count1;
int cc=0;
double sumrho=0;
double sumpsi=0;
double A=(L*L1);

double minDensity=0;
double maxDensity=5;
double minPhase=-pi;
double maxPhase=pi;

double f1=0.0;
double f2=0.0;

count = count+1;
tc[0]=0;
for(int i=1;i<=hmax;i++)
{
double ii=i;
if(i==1){
tc[1]=tc[0]+dt;}
else{tc[i]=tc[i-1]+1;}
}

// Convert the parameter values to strings
    std::stringstream ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11, ss12, ss13, ss14, ss15, ss16, ss17, ss18, ss19, ss20, ss21, ss22 ;
    ss1 << Lx;
    ss2 << Ly;
    ss3 << rhostar;
    ss4 << mu;
    ss5 << alphap;
    ss6 << lambda;
    ss7 << kappa;
    ss8 << omega;
    ss9 << Drho;
    ss10 << Dpsi;
    ss11 << ep;
    ss12 << trelax;
    ss13 << dt;
    ss14 << Dr;
    ss15 << beta;
    ss16 << rho0;
    ss17 << betap;
    ss18 << rholess;
    ss19 << dx;
    ss20 << omegalarge;
     ss21 << f1;
      ss22 << f2;
    
    std::string param1Str = ss1.str();
    std::string param2Str = ss2.str();
    std::string param3Str = ss3.str();
    std::string param4Str = ss4.str();
    std::string param5Str = ss5.str();
    std::string param6Str = ss6.str();
    std::string param7Str = ss7.str();
    std::string param8Str = ss8.str();
    std::string param9Str = ss9.str();
    std::string param10Str = ss10.str();
    std::string param11Str = ss11.str();
    std::string param12Str = ss12.str();
    std::string param13Str = ss13.str();
    std::string param14Str = ss14.str();
    std::string param15Str = ss15.str();
    std::string param16Str = ss16.str();
    std::string param17Str = ss17.str();
    std::string param18Str = ss18.str();
    std::string param19Str = ss19.str();
    std::string param20Str = ss20.str();
     std::string param21Str = ss21.str();
      std::string param22Str = ss22.str();
    
    
    //create the directory to store files
    std::string dirName = "func-D2Q9-normal-hom-den-Lx_" + param1Str + "-Ly_" + param2Str + "-rhostar_" + param3Str + "-mu_" + param4Str + "-alphap_" + param5Str + "-lambda_" + param6Str + "-kappa_" + param7Str + "-Drho_" + param9Str + "-Dpsi_" + param10Str + "-ep_" + param11Str + "-dt_" + param13Str + "-Dr_" + param14Str + "-beta_" + param15Str + "-rho0_" + param16Str + "-dx_" + param19Str ;
    mkdir(dirName.c_str(), 0777);



//ofstream f2;
//f2.open("time-series-op-Lx_256-Ly_64-non-conserved-mu_1-a_1-lambda_1-kappa_1-omega_0.2-ep_1-t1e3-inhom-dt0.001.txt");

//initialization


/*for(int i=280;i<=300;i++)
 {
 for(int j=280;j<=300;j++)
  {
   om[i][j]=omegalarge;
  }  
 }*/

ps0=2.0*(pi)*ran2(&sd4);
for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
  psi[i][j]= ps0;
// psi[i][j]=2.10866+ 0.1*ran2(&sd5);
// psi[i][j]=0.0*ran2(&sd7);
  rho[i][j]=rhostar;
  etaxrho[i][j]=gaussian();
 etayrho[i][j]=gaussian();
  }
 }
 
 /*for(int i=140;i<=160;i++)
 {
 for(int j=140;j<=160;j++)
  {
   rho[i][j]=rholess;
  }  
 }*/


// defining non-linear functions
for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
 sinpsi[i][j]=sin(psi[i][j]);
 cospsi[i][j]=cos(psi[i][j]);
 nlpsi[i][j]= ((rho0*(1+ep*cospsi[i][j]))-rho[i][j])*sinpsi[i][j];
  }
  }
   
 
 // computing finite difference laplacian and divergence
 
 for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
  
  laprho[i][j]=Lap(reinterpret_cast<double**>(rho), L, L1, i, j, dx2);
   lappsi[i][j]=Lap(reinterpret_cast<double**>(psi), L, L1, i, j, dx2);
   lapcospsi[i][j]=Lap(reinterpret_cast<double**>(cospsi), L, L1, i, j, dx2);
   
   divetarho[i][j] = Div(reinterpret_cast<double**>(etaxrho), reinterpret_cast<double**>(etayrho), L, L1, i, j, dx);
   
  }
 }


//start dynamics

double t=0;
for(int h=1; h<=hmax; h++)
{
//ss17 << h;
//std::string param17Str = ss17.str();

 // Concatenate the parameter values to create the file name
    std::string fileName = dirName + "/Lx_" + param1Str + "-Ly_" + param2Str + "-rhostar_" + param3Str + "-mu_" + param4Str + "-a_" + param5Str + "-lambda_" + param6Str + "-kappa_" + param7Str + "-omega_" + param8Str + "-Drho_" + param9Str + "-Dpsi_" + param10Str + "-ep_" + param11Str + "-t_" + to_string(h) + "-dt_" + param13Str + "-Dr_" + param14Str + "-beta_" + param15Str + "-rho0_" + param16Str + ".txt";
    
while (t < tc[h])
 {
 
 // update steps
 for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  { 
  rhonew[i][j] = rho[i][j] + alphap*laprho[i][j]*dt - (alphap*ep)*dt*lapcospsi[i][j] + (sqrt(2*Drho*beta*dt/(dx2))/(Dr*rho0))*divetarho[i][j];
 psinew[i][j] = psi[i][j] +  0.0*dt + 1.0*lappsi[i][j]*dt + (ep)*dt*nlpsi[i][j] + sqrt(2*Dpsi*dt/(Dr*(dx2)))*gaussian();
  }
 }
 
 //updating non-linear terms
 for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
    rho[i][j]=rhonew[i][j];
    psi[i][j]=psinew[i][j]; 
   etaxrho[i][j]=gaussian();
     etayrho[i][j]=gaussian();
    sinpsi[i][j]=sin(psinew[i][j]); 
    cospsi[i][j]=cos(psinew[i][j]);
   nlpsi[i][j]= ((rho0*(1+(ep*cospsi[i][j])))-rhonew[i][j])*sinpsi[i][j];
  }
  }
  
  // updating finite difference laplacian and divergence
 for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  { 
   laprho[i][j]=Lap(reinterpret_cast<double**>(rho), L, L1, i, j, dx2);
   lappsi[i][j]=Lap(reinterpret_cast<double**>(psi), L, L1, i, j, dx2);
   lapcospsi[i][j]=Lap(reinterpret_cast<double**>(cospsi), L, L1, i, j, dx2);
   
   divetarho[i][j] = Div(reinterpret_cast<double**>(etaxrho), reinterpret_cast<double**>(etayrho), L, L1, i, j, dx);
   
   }
   }
   
   
   
   
//f2<<t<<"  "<<"  "<<sumrho/A<<"  "<<sumpsi/A<<endl;
   
   
 t=t+dt;
 }  // end of dynamics
 
 

 //print output file
std::ofstream f3(fileName);
for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
 // double f=psi[i][j]-psi0;
//  double p2=2*pi;
 // double kk=fmod(f,p2);
f3<<i*dx<<"   "<<j*dx<<"  "<<rho[i][j]<<"  "<<fmod(psi[i][j], 2*M_PI)<<endl;
  }
 }
//f3<<"#"<<count<<"  "<<sumrho<<"  "<<endl;
f3.close();

std::string fileOutden = dirName + "/hom-snapshot_density_L" + param1Str + "_omega_" + param8Str + "_time_" + to_string(h) + ".png";
std::string fileOutph = dirName + "/hom-snapshot_phase_L" + param1Str + "_omega_" + param8Str + "_time_" + to_string(h) + ".png";
std::string fileOutomc = dirName + "/hom-snapshot_omegac_L" + param1Str + "_omega_" + param8Str + "_time_" + to_string(h) + ".png";

// Open a pipe to Gnuplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    
    // Send Gnuplot commands to the pipe
    std::string plotCommand1 = "set origin 0.0,0.0 \n";
    plotCommand1 +="set key samplen 1 \n";
    plotCommand1 +="set size 1.0,1.0 \n";
    plotCommand1 +="set cbrange [0.1:2] \n";
    plotCommand1 += "unset border\n";
    //plotCommand1 += "se palette model HSV defined ( 0 0 1 1, 1 1 1 1 )\n";
    plotCommand1 += "set palette defined (0 0 0 0, 1 0 0 1, 3 0 1 0, 4 1 0 0, 6 1 1 1)\n";
    plotCommand1 +="set colorbox noborder vertical \n";
    plotCommand1 +="unse xtics \n";
    plotCommand1 +="unse ytics \n";
   // plotCommand1 +="set cbtics offset 0,2.8 \n";
    plotCommand1 +="se cbtics (0.5,1,1.5,2) font ',16' \n";
     plotCommand1 +="se cblabel '{/Symbol r}' rotate by 0 offset 0,1 font ',20' \n";
    
  //  plotCommand1 += "set palette defined ( 0 0 0 0, 1 1 1 1 ) \n";
  plotCommand1 += "plot \"" + fileName + "\" using 1:2:3 w p palette pt 7 ps 0.6 ti 'density; {/Symbol w}=" + std::to_string(omega) + "; t=" + std::to_string(t) +"'"; 
  //plotCommand2 = "set term pslatex \n";
  std::string  plotCommand2 = "\n set cbrange [0:2*pi] \n";
    plotCommand2 += "se palette model HSV defined ( 0 0 1 1, 1 1 1 1 )\n";
    plotCommand2 +="set colorbox noborder vertical \n";
    plotCommand2 +="unse xtics \n";
    plotCommand2 +="unse ytics \n";
  //  plotCommand1 +="set cbtics offset 0,2.2 \n";
   plotCommand2 +="se cbtics ('0' 0, '{/Symbol p}' pi, '2{/Symbol p}' 2*pi) font ',20' \n";
   
   plotCommand2 +="se cblabel '{/Symbol y}' offset 0,5 rotate by 0 font ',20' \n";
   plotCommand2  += "plot \"" + fileName + "\" using 1:2:4 w p palette pt 7 ps 0.6 ti 'phase; {/Symbol w}=" + std::to_string(omega) + "; t=" + std::to_string(t) +"'";
  //  plotCommand2 += "\n set cbrange [3.5:5.5] \n";
  //   std::string plotCommand3 = "plot \"" + fileName + "\" using 1:2:7 w p palette pt 7 ti '{/Symbol W}_c; {/Symbol W}=" + std::to_string(omega) + "; t=" + std::to_string(t) +"'";
    
    
    fprintf(gnuplotPipe, "set term png\n");
   // fprintf(gnuplotPipe, "set output 'plot-density.png'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", fileOutden.c_str());
    fprintf(gnuplotPipe, "set xlabel ''\n");
    fprintf(gnuplotPipe, "set ylabel ''\n");
   // fprintf(gnuplotPipe, "plot \"%s\" using 1:2:3 w p palette pt 7 ti 'density; {/Symbol W}=0.2'\n", fileName.c_str());
     fprintf(gnuplotPipe, "%s\n", plotCommand1.c_str());
    
   // fprintf(gnuplotPipe, "set output 'plot-phase.png'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", fileOutph.c_str());
    fprintf(gnuplotPipe, "set xlabel ''\n");
    fprintf(gnuplotPipe, "set ylabel ''\n");
    // fprintf(gnuplotPipe, "plot \"%s\" using 1:2:(atan2($4,$5)) w p palette pt 7 ti 'phase; {/Symbol W}=0.2'\n", fileName.c_str());
     fprintf(gnuplotPipe, "%s\n", plotCommand2.c_str());
    // Close the pipe to Gnuplot
    
  //   fprintf(gnuplotPipe, "set output '%s'\n", fileOutomc.c_str());
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    // fprintf(gnuplotPipe, "plot \"%s\" using 1:2:(atan2($4,$5)) w p palette pt 7 ti 'phase; {/Symbol W}=0.2'\n", fileName.c_str());
  //   fprintf(gnuplotPipe, "%s\n", plotCommand3.c_str());
    
    
    
    pclose(gnuplotPipe);
    
    // Delete the txt file at the end
    if (remove(fileName.c_str()) != 0) {
        std::cerr << "Error deleting file: " << fileName << std::endl;
    } /*else {
        std::cout << "File successfully deleted: " << fileName << std::endl;
    }*/


} // end of hmax loop
//f2.close();
 
 /*for(int i=0;i<L;i++)
 {
 for(int j=0;j<L1;j++)
  {
 sumpsi+=psi[i][j];
  }
 }*/

for (int i = 0; i < L; i++) {
        delete[] rho[i];
        delete[] rhonew[i];
        delete[] psi[i];
        delete[] psinew[i];
        delete[] sinpsi[i];
        delete[] cospsi[i];
        delete[] laprho[i];
      //  delete[] gradxpsi[i];
       // delete[] gradypsi[i];
       delete[] lapcospsi[i];
        delete[] lappsi[i];
        delete[] nlpsi[i];    
        delete[] etaxrho[i];
        delete[] etayrho[i];   
        delete[] om[i];
        delete[] divetarho[i];
    }


delete[] rho;
delete[] rhonew;
delete[] psi;
delete[] psinew;
delete[] sinpsi;
delete[] cospsi;
delete[] laprho;
//delete[] gradxpsi;
//delete[] gradypsi;
delete[] lapcospsi;
delete[] lappsi;
delete[] nlpsi;
delete[] etaxrho;
delete[] etayrho;
delete[] om;
delete[] divetarho;

delete[] tc;

return 0;
};
