/* 
*  model.h
*    Project: c2B
*    Author: Maxime Charlebois, Simon Verret
*    License: MIT License
*/

#pragma once

#include "basicMatrix.h"

typedef struct Model {
  // things that define the model (most can be read in para.dat)
  public:

    int verbose=1;

    //tk parameters:
    double MU=0.0;
    double ETA=0.1;
    double OMEGA=0.0;    
    double t=1.0;
    double tp=0.0;
    double tpp=0.0;
    double DELTA=0.0;

    int nbOfTerms=7;

    //cuba parameters:
    double EPSREL=0.0001;
    double EPSABS=0.0001;
    int MAXEVAL=1000000;
    int MINEVAL=15000;
    int VERBOSE=0;
    
    //density parameters:
    double smallScale=0.2;
    double largeScale=10.0;
    double w_domain=1.0;
    double pole=3.0;
    double beta=50.0;

    //
    int periodization = 0; //0=green, 1=cumulant, 2=compact tiling, 3=exact

    // loop info
    int nOmega=200;
    double omegaMin=-4.0;
    double omegaMax=4.0;

    //matrices
    BasicMatrix tc;
    BasicMatrix tc2;  // used only for the exact lattice and the compact tiling
    BasicMatrix dtk;  
    BasicMatrix dtk2; // used only for the exact lattice
    BasicMatrix dtk3; // used only for the compact tiling
    BasicMatrix green;
    BasicMatrix cumul;
    BasicMatrix sigma;
    complex<double> G_per;
    complex<double> M_per;

    Model():
      tc(8), tc2(8), dtk(8), dtk2(8), dtk3(8), green(8), cumul(8), sigma(8)
    {
      if (not exists("para.dat")) {printf("ERROR: couldn't find file 'para.dat'\n\n"); exit(1);}
      printf("reading parameters from para.dat\n\n") ;
      ifstream file;
      file.open("para.dat");

      //model parameters:
      readNumber(file,"MU",MU); 
      readNumber(file,"ETA",ETA);
      readNumber(file,"OMEGA",OMEGA);
      readNumber(file,"t",t);
      readNumber(file,"tp",tp);
      readNumber(file,"tpp",tpp);
      readNumber(file,"DELTA",DELTA);
      readNumber(file,"periodization",periodization); //0=green, 1=cumulant, 2=compact tiling, 3=exact

      //cuba parameters:
      readNumber(file,"EPSREL",EPSREL);
      readNumber(file,"EPSABS",EPSABS);
      readNumber(file,"MAXEVAL",MAXEVAL); 
      readNumber(file,"MINEVAL",MINEVAL);
      readNumber(file,"VERBOSE",VERBOSE);

      //dos parameters:
      readNumber(file,"nOmega",nOmega);
      readNumber(file,"omegaMin",omegaMin);
      readNumber(file,"omegaMax",omegaMax);
    }

    void calculate_dtk(const double kx, const double ky)
      //Hk = tk matrix
    {

      complex<double> ex(cos(-kx*2.+ky*2.), sin(-kx*2.+ky*2.));
      complex<double> emx = conj(ex);
      complex<double> ey(cos(-kx*2.-ky*2.), sin(-kx*2.-ky*2.));
      complex<double> emy = conj(ey);  
   
      //assignation of the left-half 8 by 8 tc:
      tc(0,0)=  0.;  tc(0,1)= -tp;  tc(0,2)= -t;   tc(0,3)= -tp;   
      tc(1,0)= -tp;  tc(1,1)=  0.;  tc(1,2)= -t;   tc(1,3)= -tpp;  
      tc(2,0)= -t;   tc(2,1)= -t;   tc(2,2)=  0.;  tc(2,3)= -t;    
      tc(3,0)= -tp;  tc(3,1)= -tpp; tc(3,2)= -t;   tc(3,3)=  0.;   
      tc(4,0)=  0.;  tc(4,1)= -t;   tc(4,2)= -tp;  tc(4,3)=  0.;   
      tc(5,0)= -tpp; tc(5,1)= -tp;  tc(5,2)= -t;   tc(5,3)= -tp;   
      tc(6,0)=  0.;  tc(6,1)=  0.;  tc(6,2)= -tp;  tc(6,3)= -t;   
      tc(7,0)=  0.;  tc(7,1)=  0.;  tc(7,2)= -tpp; tc(7,3)=  0.;   
      
      //right-half :   
      tc(0,4)=  0.;  tc(0,5)= -tpp; tc(0,6)=  0.;  tc(0,7)=  0.;  
      tc(1,4)= -t;   tc(1,5)= -tp;  tc(1,6)=  0.;  tc(1,7)=  0.;  
      tc(2,4)= -tp;  tc(2,5)= -t;   tc(2,6)= -tp;  tc(2,7)= -tpp;  
      tc(3,4)=  0.;  tc(3,5)= -tp;  tc(3,6)= -t;   tc(3,7)=  0.; 
      tc(4,4)=  0.;  tc(4,5)= -t;   tc(4,6)= -tpp; tc(4,7)= -tp;  
      tc(5,4)= -t;   tc(5,5)=  0.;  tc(5,6)= -t;   tc(5,7)= -t; 
      tc(6,4)= -tpp; tc(6,5)= -t;   tc(6,6)=  0.;  tc(6,7)= -tp;  
      tc(7,4)= -tp;  tc(7,5)= -t;   tc(7,6)= -tp;  tc(7,7)=  0.;  
      
     
      
      
      
      //assignation of the left-half 8 by 8 dtk: 
      dtk(0,0)=  0.;                      dtk(0,1)= -tp*emx;                    dtk(0,2)=  0.;                      dtk(0,3)= -tp*ey;       
      dtk(1,0)= -tp*ex;                   dtk(1,1)=  0.;                        dtk(1,2)=  0.;                      dtk(1,3)= -tpp*(ey + ex + ey*ex);     
      dtk(2,0)=  0.;                      dtk(2,1)=  0.;                        dtk(2,2)=  0.;                      dtk(2,3)=  0.;         
      dtk(3,0)= -tp*emy;                  dtk(3,1)= -tpp*(emx + emy + emx*emy); dtk(3,2)=  0.;                      dtk(3,3)=  0.;         
      dtk(4,0)= -t*ex;                    dtk(4,1)=  0.;                        dtk(4,2)= -tp*ex;                   dtk(4,3)= -t*ex;      
      dtk(5,0)= -tpp*(ex + emy + ex*emy); dtk(5,1)= -tp*emy;                    dtk(5,2)=  0.;                      dtk(5,3)= -tp*ex;         
      dtk(6,0)= -t*emy;                   dtk(6,1)= -t*emy;                     dtk(6,2)= -tp*emy;                  dtk(6,3)=  0.;         
      dtk(7,0)= -t*ex*emy;                dtk(7,1)= -t*emy;                     dtk(7,2)= -tpp*(ex + emy + ex*emy); dtk(7,3)= -t*ex;         
            
      //Right-half :
      dtk(0,4)= -t*emx;                     dtk(0,5)= -tpp*(emx + ey + emx*ey); dtk(0,6)= -t*ey;                  dtk(0,7)= -t*emx*ey;     
      dtk(1,4)=  0.;                        dtk(1,5)= -tp*ey;                   dtk(1,6)= -t*ey;                  dtk(1,7)= -t*ey;     
      dtk(2,4)= -tp*emx;                    dtk(2,5)=  0.;                      dtk(2,6)= -tp*ey;                 dtk(2,7)= -tpp*(ey + emx + emx*ey);     
      dtk(3,4)= -t*emx;                     dtk(3,5)= -tp*emx;                  dtk(3,6)=  0.;                    dtk(3,7)= -t*emx; 
      dtk(4,4)=  0.;                        dtk(4,5)=  0.;                      dtk(4,6)= -tpp*(ex + ey + ex*ey); dtk(4,7)= -tp*ey;  
      dtk(5,4)=  0.;                        dtk(5,5)=  0.;                      dtk(5,6)=  0.;                    dtk(5,7)=  0.;     
      dtk(6,4)= -tpp*(emx + emy + emx*emy); dtk(6,5)=  0.;                      dtk(6,6)=  0.;                    dtk(6,7)= -tp*emx;     
      dtk(7,4)= -tp*emy;                    dtk(7,5)=  0.;                      dtk(7,6)= -tp*ex;                 dtk(7,7)=  0.;     
            
      
      
      
      
      
      //used for exact lattice and compact tiling (left-half) :
      tc2(0,0)=  0.;  tc2(0,1)= -tp;  tc2(0,2)=  t;   tc2(0,3)= -tp;   
      tc2(1,0)= -tp;  tc2(1,1)=  0.;  tc2(1,2)=  t;   tc2(1,3)= -tpp;  
      tc2(2,0)=  t;   tc2(2,1)=  t;   tc2(2,2)=  0.;  tc2(2,3)=  t;    
      tc2(3,0)= -tp;  tc2(3,1)= -tpp; tc2(3,2)=  t;   tc2(3,3)=  0.;   
      tc2(4,0)=  0.;  tc2(4,1)=  t;   tc2(4,2)= -tp;  tc2(4,3)=  0.;   
      tc2(5,0)= -tpp; tc2(5,1)= -tp;  tc2(5,2)=  t;   tc2(5,3)= -tp;   
      tc2(6,0)=  0.;  tc2(6,1)=  0.;  tc2(6,2)= -tp;  tc2(6,3)=  t;   
      tc2(7,0)=  0.;  tc2(7,1)=  0.;  tc2(7,2)= -tpp; tc2(7,3)=  0.;   
      
      //right-half :   
      tc2(0,4)=  0.;  tc2(0,5)= -tpp; tc2(0,6)=  0.;  tc2(0,7)=  0.;  
      tc2(1,4)=  t;   tc2(1,5)= -tp;  tc2(1,6)=  0.;  tc2(1,7)=  0.;  
      tc2(2,4)= -tp;  tc2(2,5)=  t;   tc2(2,6)= -tp;  tc2(2,7)= -tpp;  
      tc2(3,4)=  0.;  tc2(3,5)= -tp;  tc2(3,6)=  t;   tc2(3,7)=  0.; 
      tc2(4,4)=  0.;  tc2(4,5)=  t;   tc2(4,6)= -tpp; tc2(4,7)= -tp;  
      tc2(5,4)=  t;   tc2(5,5)=  0.;  tc2(5,6)=  t;   tc2(5,7)=  t; 
      tc2(6,4)= -tpp; tc2(6,5)=  t;   tc2(6,6)=  0.;  tc2(6,7)= -tp;  
      tc2(7,4)= -tp;  tc2(7,5)=  t;   tc2(7,6)= -tp;  tc2(7,7)=  0.;  





      //used for exact lattice (Left-half) :
      dtk2(0,0)=  0.;                      dtk2(0,1)= -tp*emx;                    dtk2(0,2)=  0.;                      dtk2(0,3)= -tp*ey;       
      dtk2(1,0)= -tp*ex;                   dtk2(1,1)=  0.;                        dtk2(1,2)=  0.;                      dtk2(1,3)= -tpp*(ey + ex + ey*ex);     
      dtk2(2,0)=  0.;                      dtk2(2,1)=  0.;                        dtk2(2,2)=  0.;                      dtk2(2,3)=  0.;         
      dtk2(3,0)= -tp*emy;                  dtk2(3,1)= -tpp*(emx + emy + emx*emy); dtk2(3,2)=  0.;                      dtk2(3,3)=  0.;         
      dtk2(4,0)=  t*ex;                    dtk2(4,1)=  0.;                        dtk2(4,2)= -tp*ex;                   dtk2(4,3)=  t*ex;      
      dtk2(5,0)= -tpp*(ex + emy + ex*emy); dtk2(5,1)= -tp*emy;                    dtk2(5,2)=  0.;                      dtk2(5,3)= -tp*ex;         
      dtk2(6,0)=  t*emy;                   dtk2(6,1)=  t*emy;                     dtk2(6,2)= -tp*emy;                  dtk2(6,3)=  0.;         
      dtk2(7,0)=  t*ex*emy;                dtk2(7,1)=  t*emy;                     dtk2(7,2)= -tpp*(ex + emy + ex*emy); dtk2(7,3)=  t*ex;         
            
      //Right-half :
      dtk2(0,4)=  t*emx;                     dtk2(0,5)= -tpp*(emx + ey + emx*ey); dtk2(0,6)=  t*ey;                  dtk2(0,7)=  t*emx*ey;     
      dtk2(1,4)=  0.;                        dtk2(1,5)= -tp*ey;                   dtk2(1,6)=  t*ey;                  dtk2(1,7)=  t*ey;     
      dtk2(2,4)= -tp*emx;                    dtk2(2,5)=  0.;                      dtk2(2,6)= -tp*ey;                 dtk2(2,7)= -tpp*(ey + emx + emx*ey);     
      dtk2(3,4)=  t*emx;                     dtk2(3,5)= -tp*emx;                  dtk2(3,6)=  0.;                    dtk2(3,7)=  t*emx; 
      dtk2(4,4)=  0.;                        dtk2(4,5)=  0.;                      dtk2(4,6)= -tpp*(ex + ey + ex*ey); dtk2(4,7)= -tp*ey;  
      dtk2(5,4)=  0.;                        dtk2(5,5)=  0.;                      dtk2(5,6)=  0.;                    dtk2(5,7)=  0.;     
      dtk2(6,4)= -tpp*(emx + emy + emx*emy); dtk2(6,5)=  0.;                      dtk2(6,6)=  0.;                    dtk2(6,7)= -tp*emx;     
      dtk2(7,4)= -tp*emy;                    dtk2(7,5)=  0.;                      dtk2(7,6)= -tp*ex;                 dtk2(7,7)=  0.;     






      //used for compact tiling:
      dtk3(0,0)=  0.;                      dtk3(0,1)= -tp*emx;                    dtk3(0,2)=  0.;                      dtk3(0,3)= -tp*ey;       
      dtk3(1,0)= -tp*ex;                   dtk3(1,1)=  0.;                        dtk3(1,2)=  0.;                      dtk3(1,3)= -tpp*(ey + ex + ey*ex);     
      dtk3(2,0)=  0.;                      dtk3(2,1)=  0.;                        dtk3(2,2)=  0.;                      dtk3(2,3)=  0.;         
      dtk3(3,0)= -tp*emy;                  dtk3(3,1)= -tpp*(emx + emy + emx*emy); dtk3(3,2)=  0.;                      dtk3(3,3)=  0.;         
      dtk3(4,0)=  t*ex;                    dtk3(4,1)=  0.;                        dtk3(4,2)= -tp*ex;                   dtk3(4,3)=  t*ex;      
      dtk3(5,0)= -tpp*(ex + emy + ex*emy); dtk3(5,1)= -tp*emy;                    dtk3(5,2)=  0.;                      dtk3(5,3)= -tp*ex;         
      dtk3(6,0)=  t*emy;                   dtk3(6,1)=  t*emy;                     dtk3(6,2)= -tp*emy;                  dtk3(6,3)=  0.;         
      dtk3(7,0)=  t*ex*emy;                dtk3(7,1)=  t*emy;                     dtk3(7,2)= -tpp*(ex + emy + ex*emy); dtk3(7,3)=  t*ex;         
            
      //Right-half :
      dtk3(0,4)=  t*emx;                     dtk3(0,5)= -tpp*(emx + ey + emx*ey); dtk3(0,6)=  t*ey;                  dtk3(0,7)=  t*emx*ey;     
      dtk3(1,4)=  0.;                        dtk3(1,5)= -tp*ey;                   dtk3(1,6)=  t*ey;                  dtk3(1,7)=  t*ey;     
      dtk3(2,4)= -tp*emx;                    dtk3(2,5)=  0.;                      dtk3(2,6)= -tp*ey;                 dtk3(2,7)= -tpp*(ey + emx + emx*ey);     
      dtk3(3,4)=  t*emx;                     dtk3(3,5)= -tp*emx;                  dtk3(3,6)=  0.;                    dtk3(3,7)=  t*emx; 
      dtk3(4,4)=  0.;                        dtk3(4,5)=  0.;                      dtk3(4,6)= -tpp*(ex + ey + ex*ey); dtk3(4,7)= -tp*ey;  
      dtk3(5,4)=  0.;                        dtk3(5,5)=  0.;                      dtk3(5,6)=  0.;                    dtk3(5,7)=  0.;     
      dtk3(6,4)= -tpp*(emx + emy + emx*emy); dtk3(6,5)=  0.;                      dtk3(6,6)=  0.;                    dtk3(6,7)= -tp*emx;     
      dtk3(7,4)= -tp*emy;                    dtk3(7,5)=  0.;                      dtk3(7,6)= -tp*ex;                 dtk3(7,7)=  0.;     
    };

    void calculate_sigma(const complex<double> z)
      //S = sigma matrix
    {
      // the self-energy matrix is S = 1/(z+mu-tc2)
      for(int i=0;i<sigma.dim;i++)
        for(int j=0;j<sigma.dim;j++)
        {
          sigma(i,j) = -tc2(i,j);
          if (i==j) sigma(i,j) += z + MU;  
          if (periodization==2) sigma(i,j) += -dtk3(i,j);
          if (periodization==3) sigma(i,j) += -dtk2(i,j);
        }

      sigma.invert();

      for(int i=0;i<sigma.dim;i++)
        for(int j=0;j<sigma.dim;j++)
        {
          sigma(i,j) = DELTA*DELTA*sigma(i,j);
        }
    }

    void calculate_cumulant(const complex<double> z)
      //M = cumul matrix
    {
      // the cumulant matrix is M = 1/(z+mu-Sigma)
      for(int i=0;i<cumul.dim;i++)
        for(int j=0;j<cumul.dim;j++)
        {
          cumul(i,j) = -sigma(i,j);
          if(i==j) cumul(i,j) += z + MU;
        }
      cumul.invert();
    }


    void calculate_Gk(const complex<double> z)
      //Gk = Green matrix
    {   
      // the Green matrix is Gk = 1/(z-Hk)
      for(int i=0;i<green.dim;i++)
        for(int j=0;j<green.dim;j++)
        {
          green(i,j) = -tc(i,j) -dtk(i,j) -sigma(i,j);
          if(i==j) green(i,j) += z + MU;
        }
      green.invert();
    }


    void calculate_Gperiodized(const double px, const double py)
    {
      double Ry[8] = {0.,1.,1.,1.,2.,2.,2.,3.};
      double Rx[8] = {0.,-1.,0.,1.,-1.,0.,1.,0.};
      complex<double> z(OMEGA,ETA);
      
      calculate_dtk(px, py);
      calculate_sigma(z);

      if ((periodization==1) || (periodization==4)){
        calculate_cumulant(z);
        M_per = 0;
        for (int ii=0; ii<8; ++ii) {
          for (int jj=0; jj<8; ++jj) {
            double arg = ((Rx[jj]-Rx[ii])*px + (Ry[jj]-Ry[ii])*py);
            complex<double> phase(cos(arg), sin(arg));
            M_per += 0.125 * cumul(ii,jj) * phase; 
          }
        }
        double epsilon_k = -2*t*(cos(px)+cos(py)) -4*tp*cos(px)*cos(py) -2*tpp*(cos(2*px)+cos(2*py));
        if(periodization==1) {G_per = 1./((1./M_per) - epsilon_k);}
        else {G_per = (z + MU - (1./M_per) );}      //calculate self-energy for M_per
      }
      else {
        calculate_Gk(z);
        G_per = 0;
        for (int ii=0; ii<8; ++ii) {
          for (int jj=0; jj<8; ++jj) {
            double arg = ((Rx[jj]-Rx[ii])*px + (Ry[jj]-Ry[ii])*py);
            complex<double> phase(cos(arg), sin(arg));
            G_per += 0.125 * green(ii,jj) * phase; 
          }
        }
      }
    }
} Model;
