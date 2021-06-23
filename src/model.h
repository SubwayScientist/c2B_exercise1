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
      tc(16), tc2(16), dtk(16), dtk2(16), dtk3(16), green(16), cumul(16), sigma(16)
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

      complex<double> ex(cos(-kx*4.), sin(-kx*4.));
      complex<double> emx = conj(ex);
      complex<double> ey(cos(-ky*4.), sin(-ky*4.));
      complex<double> emy = conj(ey);  
   
      //assignation of the left-half 16 by 16 tc:
      tc(0,0)=  0.;  tc(0,1)= -t;    tc(0,2)= -tpp;  tc(0,3)=  0.;   tc(0,4)= -t;    tc(0,5)= -tp;   tc(0,6)=  0.;   tc(0,7)=  0.;
      tc(1,0)= -t;   tc(1,1)=  0.;   tc(1,2)= -t;    tc(1,3)= -tpp;  tc(1,4)= -tp;   tc(1,5)=  -t;   tc(1,6)= -tp;   tc(1,7)=  0.;
      tc(2,0)= -tpp; tc(2,1)= -t;    tc(2,2)=  0.;   tc(2,3)= -t;    tc(2,4)=  0.;   tc(2,5)= -tp;   tc(2,6)=  -t;   tc(2,7)= -tp;
      tc(3,0)=  0.;  tc(3,1)= -tpp;  tc(3,2)= -t;    tc(3,3)=  0.;   tc(3,4)=  0.;   tc(3,5)=  0.;   tc(3,6)= -tp;   tc(3,7)= -t;
      tc(4,0)=  -t;  tc(4,1)= -tp;   tc(4,2)=  0.;   tc(4,3)=  0.;   tc(4,4)=  0.;   tc(4,5)= -t;    tc(4,6)= -tpp;  tc(4,7)=  0.;
      tc(5,0)= -tp;  tc(5,1)=  -t;   tc(5,2)= -tp;   tc(5,3)=  0.;   tc(5,4)= -t;    tc(5,5)=  0.;   tc(5,6)= -t;    tc(5,7)= -tpp;
      tc(6,0)=  0.;  tc(6,1)= -tp;   tc(6,2)= -t;    tc(6,3)= -tp;   tc(6,4)= -tpp;  tc(6,5)= -t;    tc(6,6)=  0.;   tc(6,7)= -t;
      tc(7,0)=  0.;  tc(7,1)=  0.;   tc(7,2)= -tp;   tc(7,3)=  -t;   tc(7,4)=  0.;   tc(7,5)= -tpp;  tc(7,6)= -t;    tc(7,7)=  0.;
      tc(8,0)= -tpp; tc(8,1)=  0.;   tc(8,2)=  0.;   tc(8,3)=  0.;   tc(8,4)=  -t;   tc(8,5)= -tp;   tc(8,6)=  0.;   tc(8,7)=  0.;
      tc(9,0)=  0.;  tc(9,1)=  -tpp; tc(9,2)=  0.;   tc(9,3)=  0.;   tc(9,4)= -tp;   tc(9,5)=  -t;   tc(9,6)= -tp;   tc(9,7)=  0.;
      tc(10,0)= 0.;  tc(10,1)= 0.;   tc(10,2)= -tpp; tc(10,3)= 0.;   tc(10,4)= 0.;   tc(10,5)= -tp;  tc(10,6)=  -t;  tc(10,7)= -tp;
      tc(11,0)= 0.;  tc(11,1)= 0.;   tc(11,2)= 0.;   tc(11,3)= -tpp; tc(11,4)= 0.;   tc(11,5)= 0.;   tc(11,6)= -tp;  tc(11,7)=  -t;
      tc(12,0)= 0.;  tc(12,1)= 0.;   tc(12,2)= 0.;   tc(12,3)= 0.;   tc(12,4)= -tpp; tc(12,5)= 0.;   tc(12,6)= 0.;   tc(12,7)= 0.;
      tc(13,0)= 0.;  tc(13,1)= 0.;   tc(13,2)= 0.;   tc(13,3)= 0.;   tc(13,4)= 0.;   tc(13,5)= -tpp; tc(13,6)= 0.;   tc(13,7)= 0.;
      tc(14,0)= 0.;  tc(14,1)= 0.;   tc(14,2)= 0.;   tc(14,3)= 0.;   tc(14,4)= 0.;   tc(14,5)= 0.;   tc(14,6)= -tpp; tc(14,7)= 0.;
      tc(15,0)= 0.;  tc(15,1)= 0.;   tc(15,2)= 0.;   tc(15,3)= 0.;   tc(15,4)= 0.;   tc(15,5)= 0.;   tc(15,6)= 0.;   tc(15,7)= -tpp;
      //right-half of the 16 by 16 tc:    
      tc(0,8)= -tpp;  tc(0,9)=  0.;   tc(0,10)=  0.;  tc(0,11)=  0.;  tc(0,12)=  0.;   tc(0,13)=  0.;   tc(0,14)=  0.;   tc(0,15)=  0.;
      tc(1,8)=  0.;   tc(1,9)= -tpp;  tc(1,10)=  0.;  tc(1,11)=  0.;  tc(1,12)=  0.;   tc(1,13)=  0.;   tc(1,14)=  0.;   tc(1,15)=  0.;
      tc(2,8)=  0.;   tc(2,9)=  0.;   tc(2,10)= -tpp; tc(2,11)=  0.;  tc(2,12)=  0.;   tc(2,13)=  0.;   tc(2,14)=  0.;   tc(2,15)=  0.;
      tc(3,8)=  0.;   tc(3,9)=  0.;   tc(3,10)=  0.;  tc(3,11)= -tpp; tc(3,12)=  0.;   tc(3,13)=  0.;   tc(3,14)=  0.;   tc(3,15)=  0.;
      tc(4,8)= -t;    tc(4,9)= -tp;   tc(4,10)=  0.;  tc(4,11)=  0.;  tc(4,12)= -tpp;  tc(4,13)=  0.;   tc(4,14)=  0.;   tc(4,15)=  0.;
      tc(5,8)= -tp;   tc(5,9)=  -t;   tc(5,10)= -tp;  tc(5,11)=  0.;  tc(5,12)=  0.;   tc(5,13)= -tpp;  tc(5,14)=  0.;   tc(5,15)=  0.;
      tc(6,8)=  0.;   tc(6,9)= -tp;   tc(6,10)= -t;   tc(6,11)= -tp;  tc(6,12)=  0.;   tc(6,13)=  0.;   tc(6,14)= -tpp;  tc(6,15)=  0.;
      tc(7,8)=  0.;   tc(7,9)=  0.;   tc(7,10)= -tp;  tc(7,11)=  -t;  tc(7,12)=  0.;   tc(7,13)=  0.;   tc(7,14)=  0.;   tc(7,15)= -tpp;
      tc(8,8)=  0.;   tc(8,9)= -t;    tc(8,10)= -tpp; tc(8,11)=  0.;  tc(8,12)= -t;    tc(8,13)= -tp;   tc(8,14)=  0.;   tc(8,15)=  0.;
      tc(9,8)= -t;    tc(9,9)=  0.;   tc(9,10)= -t;   tc(9,11)= -tpp; tc(9,12)= -tp;   tc(9,13)=  -t;   tc(9,14)= -tp;   tc(9,15)=  0.;
      tc(10,8)= -tpp; tc(10,9)= -t;   tc(10,10)= 0.;  tc(10,11)= -t;  tc(10,12)= 0.;   tc(10,13)= -tp;  tc(10,14)= -t;   tc(10,15)= -tp;
      tc(11,8)= 0.;   tc(11,9)= -tpp; tc(11,10)= -t;  tc(11,11)=  0.; tc(11,12)= 0.;   tc(11,13)= 0.;   tc(11,14)= -tp;  tc(11,15)=  -t;
      tc(12,8)= -t;   tc(12,9)= -tp;  tc(12,10)= 0.;  tc(12,11)=  0.; tc(12,12)= 0.;   tc(12,13)= -t;   tc(12,14)= -tpp; tc(12,15)= 0.;
      tc(13,8)= -tp;  tc(13,9)=  -t;  tc(13,10)= -tp; tc(13,11)=  0.; tc(13,12)= -t;   tc(13,13)= 0.;   tc(13,14)= -t;   tc(13,15)= -tpp;
      tc(14,8)= 0.;   tc(14,9)= -tp;  tc(14,10)= -t;  tc(14,11)= -tp; tc(14,12)= -tpp; tc(14,13)= -t;   tc(14,14)=  0.;  tc(14,15)= -t;
      tc(15,8)= 0.;   tc(15,9)= 0.;   tc(15,10)= -tp; tc(15,11)=  -t; tc(15,12)= 0.;   tc(15,13)= -tpp; tc(15,14)= -t;   tc(15,15)= 0.;
     
      
      
      
      //assignation of the left-half 4 by 4 dtk: 
      dtk(0,0)=  0.;          dtk(0,1)=  0.;      dtk(0,2)= -tpp*ex;   dtk(0,3)= -t*ex;       dtk(0,4)=  0.;       dtk(0,5)=  0.;       dtk(0,6)=  0.;       dtk(0,7)= -tp*ex;
      dtk(1,0)=  0.;          dtk(1,1)=  0.;      dtk(1,2)=  0.;       dtk(1,3)= -tpp*ex;     dtk(1,4)=  0.;       dtk(1,5)=  0.;       dtk(1,6)=  0.;       dtk(1,7)=  0.;
      dtk(2,0)= -tpp*emx;     dtk(2,1)=  0.;      dtk(2,2)=  0.;       dtk(2,3)=  0.;         dtk(2,4)=  0.;       dtk(2,5)=  0.;       dtk(2,6)=  0.;       dtk(2,7)=  0.;
      dtk(3,0)= -t*emx;       dtk(3,1)= -tpp*emx; dtk(3,2)=  0.;       dtk(3,3)=  0.;         dtk(3,4)= -tp*emx;   dtk(3,5)=  0.;       dtk(3,6)=  0.;       dtk(3,7)=  0.;
      dtk(4,0)=  0.;          dtk(4,1)=  0.;      dtk(4,2)=  0.;       dtk(4,3)= -tp*ex;      dtk(4,4)=  0.;       dtk(4,5)=  0.;       dtk(4,6)= -tpp*ex;   dtk(4,7)= -t*ex;
      dtk(5,0)=  0.;          dtk(5,1)=  0.;      dtk(5,2)=  0.;       dtk(5,3)=  0.;         dtk(5,4)=  0.;       dtk(5,5)=  0.;       dtk(5,6)=  0.;       dtk(5,7)= -tpp*ex;
      dtk(6,0)=  0.;          dtk(6,1)=  0.;      dtk(6,2)=  0.;       dtk(6,3)=  0.;         dtk(6,4)= -tpp*emx;  dtk(6,5)=  0.;       dtk(6,6)=  0.;       dtk(6,7)=  0.;
      dtk(7,0)= -tp*emx;      dtk(7,1)=  0.;      dtk(7,2)=  0.;       dtk(7,3)=  0.;         dtk(7,4)= -t*emx;    dtk(7,5)= -tpp*emx;  dtk(7,6)=  0.;       dtk(7,7)=  0.;
      dtk(8,0)= -tpp*emy;     dtk(8,1)=  0.;      dtk(8,2)=  0.;       dtk(8,3)=  0.;         dtk(8,4)=  0.;       dtk(8,5)=  0.;       dtk(8,6)=  0.;       dtk(8,7)= -tp*ex;
      dtk(9,0)=  0.;          dtk(9,1)= -tpp*emy; dtk(9,2)=  0.;       dtk(9,3)=  0.;         dtk(9,4)=  0.;       dtk(9,5)=  0.;       dtk(9,6)=  0.;       dtk(9,7)=  0.;
      dtk(10,0)= 0.;          dtk(10,1)= 0.;      dtk(10,2)= -tpp*emy; dtk(10,3)= 0.;         dtk(10,4)= 0.;       dtk(10,5)= 0.;       dtk(10,6)= 0.;       dtk(10,7)= 0.;
      dtk(11,0)= 0.;          dtk(11,1)= 0.;      dtk(11,2)= 0.;       dtk(11,3)= -tpp*emy;   dtk(11,4)= -tp*emx;  dtk(11,5)= 0.;       dtk(11,6)= 0.;       dtk(11,7)= 0.;
      dtk(12,0)= -t*emy;      dtk(12,1)= -tp*emy; dtk(12,2)= 0.;       dtk(12,3)= -tp*ex*emy; dtk(12,4)= -tpp*emy; dtk(12,5)= 0.;       dtk(12,6)= 0.;       dtk(12,7)= 0.;
      dtk(13,0)= -tp*emy;     dtk(13,1)= -t*emy;  dtk(13,2)= -tp*emy;  dtk(13,3)= 0.;         dtk(13,4)= 0.;       dtk(13,5)= -tpp*emy; dtk(13,6)= 0.;       dtk(13,7)= 0.;
      dtk(14,0)= 0.;          dtk(14,1)= -tp*emy; dtk(14,2)=  -t*emy;  dtk(14,3)= -tp*emy;    dtk(14,4)= 0.;       dtk(14,5)= 0.;       dtk(14,6)= -tpp*emy; dtk(14,7)= 0.;
      dtk(15,0)= -tp*emy*emx; dtk(15,1)= 0.;      dtk(15,2)= -tp*emy;  dtk(15,3)= -t*emy;     dtk(15,4)= 0.;       dtk(15,5)= 0.;       dtk(15,6)= 0.;       dtk(15,7)= -tpp*emy;
      //Right-half of the 16 by 16 dtk
      dtk(0,8)= -tpp*ey;   dtk(0,9)=  0.;       dtk(0,10)=  0.;     dtk(0,11)=  0.;     dtk(0,12)= -t*ey;      dtk(0,13)= -tp*ey;    dtk(0,14)=  0.;      dtk(0,15)= -tp*ex*ey;
      dtk(1,8)=  0.;       dtk(1,9)= -tpp*ey;   dtk(1,10)=  0.;     dtk(1,11)=  0.;     dtk(1,12)= -tp*ey;     dtk(1,13)= -t*ey;     dtk(1,14)= -tp*ey;   dtk(1,15)=  0.;
      dtk(2,8)=  0.;       dtk(2,9)=  0.;       dtk(2,10)= -tpp*ey; dtk(2,11)=  0.;     dtk(2,12)=  0.;        dtk(2,13)= -tp*ey;    dtk(2,14)= -t*ey;    dtk(2,15)= -tp*ey;
      dtk(3,8)=  0.;       dtk(3,9)=  0.;       dtk(3,10)=  0.;     dtk(3,11)= -tpp*ey; dtk(3,12)= -tp*emx*ey; dtk(3,13)=  0.;       dtk(3,14)= -tp*ey;   dtk(3,15)=  -t*ey;
      dtk(4,8)=  0.;       dtk(4,9)=  0.;       dtk(4,10)=  0.;     dtk(4,11)= -tp*ex;  dtk(4,12)= -tpp*ey;    dtk(4,13)=  0.;       dtk(4,14)=  0.;      dtk(4,15)=  0.;
      dtk(5,8)=  0.;       dtk(5,9)=  0.;       dtk(5,10)=  0.;     dtk(5,11)=  0.;     dtk(5,12)=  0.;        dtk(5,13)= -tpp*ey;   dtk(5,14)=  0.;      dtk(5,15)=  0.;
      dtk(6,8)=  0.;       dtk(6,9)=  0.;       dtk(6,10)=  0.;     dtk(6,11)=  0.;     dtk(6,12)=  0.;        dtk(6,13)=  0.;       dtk(6,14)= -tpp*ey;  dtk(6,15)=  0.;
      dtk(7,8)= -tp*emx;   dtk(7,9)=  0.;       dtk(7,10)=  0.;     dtk(7,11)=  0.;     dtk(7,12)=  0.;        dtk(7,13)=  0.;       dtk(7,14)=  0.;      dtk(7,15)= -tpp*ey;
      dtk(8,8)=  0.;       dtk(8,9)=  0.;       dtk(8,10)= -tpp*ex; dtk(8,11)= -t*ex;   dtk(8,12)=  0.;        dtk(8,13)=  0.;       dtk(8,14)=  0.;      dtk(8,15)= -tp*ex;
      dtk(9,8)=  0.;       dtk(9,9)=  0.;       dtk(9,10)=  0.;     dtk(9,11)= -tpp*ex; dtk(9,12)=  0.;        dtk(9,13)=  0.;       dtk(9,14)=  0.;      dtk(9,15)=  0.;
      dtk(10,8)= -tpp*emx; dtk(10,9)= 0.;       dtk(10,10)= 0.;     dtk(10,11)= 0.;     dtk(10,12)= 0.;        dtk(10,13)= 0.;       dtk(10,14)= 0.;      dtk(10,15)= 0.;
      dtk(11,8)= -t*emx;   dtk(11,9)= -tpp*emx; dtk(11,10)= 0.;     dtk(11,11)= 0.;     dtk(11,12)= -tp*emx;   dtk(11,13)= 0.;       dtk(11,14)= 0.;      dtk(11,15)= 0.;
      dtk(12,8)= 0.;       dtk(12,9)= 0.;       dtk(12,10)= 0.;     dtk(12,11)= -tp*ex; dtk(12,12)= 0.;        dtk(12,13)= 0.;       dtk(12,14)= -tpp*ex; dtk(12,15)= -t*ex;
      dtk(13,8)= 0.;       dtk(13,9)= 0.;       dtk(13,10)= 0.;     dtk(13,11)= 0.;     dtk(13,12)= 0.;        dtk(13,13)= 0.;       dtk(13,14)= 0.;      dtk(13,15)= -tpp*ex;
      dtk(14,8)= 0.;       dtk(14,9)= 0.;       dtk(14,10)= 0.;     dtk(14,11)= 0.;     dtk(14,12)= -tpp*emx;  dtk(14,13)= 0.;       dtk(14,14)= 0.;      dtk(14,15)= 0.;
      dtk(15,8)= -tp*emx;  dtk(15,9)= 0.;       dtk(15,10)= 0.;     dtk(15,11)= 0.;     dtk(15,12)= -t*emx;    dtk(15,13)= -tpp*emx; dtk(15,14)= 0.;      dtk(15,15)= 0.;
      
      
      
      
      
      
      //used for exact lattice and compact tiling (left-half) :
      tc2(0,0)=  0.;  tc2(0,1)=  t;    tc2(0,2)= -tpp;  tc2(0,3)=  0.;   tc2(0,4)=  t;    tc2(0,5)= -tp;   tc2(0,6)=  0.;   tc2(0,7)=  0.;
      tc2(1,0)=  t;   tc2(1,1)=  0.;   tc2(1,2)=  t;    tc2(1,3)= -tpp;  tc2(1,4)= -tp;   tc2(1,5)=   t;   tc2(1,6)= -tp;   tc2(1,7)=  0.;
      tc2(2,0)= -tpp; tc2(2,1)=  t;    tc2(2,2)=  0.;   tc2(2,3)=  t;    tc2(2,4)=  0.;   tc2(2,5)= -tp;   tc2(2,6)=   t;   tc2(2,7)= -tp;
      tc2(3,0)=  0.;  tc2(3,1)= -tpp;  tc2(3,2)=  t;    tc2(3,3)=  0.;   tc2(3,4)=  0.;   tc2(3,5)=  0.;   tc2(3,6)= -tp;   tc2(3,7)=  t;
      tc2(4,0)=  t;   tc2(4,1)= -tp;   tc2(4,2)=  0.;   tc2(4,3)=  0.;   tc2(4,4)=  0.;   tc2(4,5)=  t;    tc2(4,6)= -tpp;  tc2(4,7)=  0.;
      tc2(5,0)= -tp;  tc2(5,1)=   t;   tc2(5,2)= -tp;   tc2(5,3)=  0.;   tc2(5,4)=  t;    tc2(5,5)=  0.;   tc2(5,6)=  t;    tc2(5,7)= -tpp;
      tc2(6,0)=  0.;  tc2(6,1)= -tp;   tc2(6,2)=  t;    tc2(6,3)= -tp;   tc2(6,4)= -tpp;  tc2(6,5)=  t;    tc2(6,6)=  0.;   tc2(6,7)=  t;
      tc2(7,0)=  0.;  tc2(7,1)=  0.;   tc2(7,2)= -tp;   tc2(7,3)=   t;   tc2(7,4)=  0.;   tc2(7,5)= -tpp;  tc2(7,6)=  t;    tc2(7,7)=  0.;
      tc2(8,0)= -tpp; tc2(8,1)=  0.;   tc2(8,2)=  0.;   tc2(8,3)=  0.;   tc2(8,4)=   t;   tc2(8,5)= -tp;   tc2(8,6)=  0.;   tc2(8,7)=  0.;
      tc2(9,0)=  0.;  tc2(9,1)=  -tpp; tc2(9,2)=  0.;   tc2(9,3)=  0.;   tc2(9,4)= -tp;   tc2(9,5)=   t;   tc2(9,6)= -tp;   tc2(9,7)=  0.;
      tc2(10,0)= 0.;  tc2(10,1)= 0.;   tc2(10,2)= -tpp; tc2(10,3)= 0.;   tc2(10,4)= 0.;   tc2(10,5)= -tp;  tc2(10,6)=   t;  tc2(10,7)= -tp;
      tc2(11,0)= 0.;  tc2(11,1)= 0.;   tc2(11,2)= 0.;   tc2(11,3)= -tpp; tc2(11,4)= 0.;   tc2(11,5)= 0.;   tc2(11,6)= -tp;  tc2(11,7)=   t;
      tc2(12,0)= 0.;  tc2(12,1)= 0.;   tc2(12,2)= 0.;   tc2(12,3)= 0.;   tc2(12,4)= -tpp; tc2(12,5)= 0.;   tc2(12,6)= 0.;   tc2(12,7)= 0.;
      tc2(13,0)= 0.;  tc2(13,1)= 0.;   tc2(13,2)= 0.;   tc2(13,3)= 0.;   tc2(13,4)= 0.;   tc2(13,5)= -tpp; tc2(13,6)= 0.;   tc2(13,7)= 0.;
      tc2(14,0)= 0.;  tc2(14,1)= 0.;   tc2(14,2)= 0.;   tc2(14,3)= 0.;   tc2(14,4)= 0.;   tc2(14,5)= 0.;   tc2(14,6)= -tpp; tc2(14,7)= 0.;
      tc2(15,0)= 0.;  tc2(15,1)= 0.;   tc2(15,2)= 0.;   tc2(15,3)= 0.;   tc2(15,4)= 0.;   tc2(15,5)= 0.;   tc2(15,6)= 0.;   tc2(15,7)= -tpp;
      //right-half :   
      tc2(0,8)= -tpp;  tc2(0,9)=  0.;   tc2(0,10)=  0.;  tc2(0,11)=  0.;  tc2(0,12)=  0.;   tc2(0,13)=  0.;   tc2(0,14)=  0.;   tc2(0,15)=  0.;
      tc2(1,8)=  0.;   tc2(1,9)= -tpp;  tc2(1,10)=  0.;  tc2(1,11)=  0.;  tc2(1,12)=  0.;   tc2(1,13)=  0.;   tc2(1,14)=  0.;   tc2(1,15)=  0.;
      tc2(2,8)=  0.;   tc2(2,9)=  0.;   tc2(2,10)= -tpp; tc2(2,11)=  0.;  tc2(2,12)=  0.;   tc2(2,13)=  0.;   tc2(2,14)=  0.;   tc2(2,15)=  0.;
      tc2(3,8)=  0.;   tc2(3,9)=  0.;   tc2(3,10)=  0.;  tc2(3,11)= -tpp; tc2(3,12)=  0.;   tc2(3,13)=  0.;   tc2(3,14)=  0.;   tc2(3,15)=  0.;
      tc2(4,8)=  t;    tc2(4,9)= -tp;   tc2(4,10)=  0.;  tc2(4,11)=  0.;  tc2(4,12)= -tpp;  tc2(4,13)=  0.;   tc2(4,14)=  0.;   tc2(4,15)=  0.;
      tc2(5,8)= -tp;   tc2(5,9)=   t;   tc2(5,10)= -tp;  tc2(5,11)=  0.;  tc2(5,12)=  0.;   tc2(5,13)= -tpp;  tc2(5,14)=  0.;   tc2(5,15)=  0.;
      tc2(6,8)=  0.;   tc2(6,9)= -tp;   tc2(6,10)=  t;   tc2(6,11)= -tp;  tc2(6,12)=  0.;   tc2(6,13)=  0.;   tc2(6,14)= -tpp;  tc2(6,15)=  0.;
      tc2(7,8)=  0.;   tc2(7,9)=  0.;   tc2(7,10)= -tp;  tc2(7,11)=   t;  tc2(7,12)=  0.;   tc2(7,13)=  0.;   tc2(7,14)=  0.;   tc2(7,15)= -tpp;
      tc2(8,8)=  0.;   tc2(8,9)=  t;    tc2(8,10)= -tpp; tc2(8,11)=  0.;  tc2(8,12)=  t;    tc2(8,13)= -tp;   tc2(8,14)=  0.;   tc2(8,15)=  0.;
      tc2(9,8)=  t;    tc2(9,9)=  0.;   tc2(9,10)=  t;   tc2(9,11)= -tpp; tc2(9,12)= -tp;   tc2(9,13)=   t;   tc2(9,14)= -tp;   tc2(9,15)=  0.;
      tc2(10,8)= -tpp; tc2(10,9)=  t;   tc2(10,10)= 0.;  tc2(10,11)=  t;  tc2(10,12)= 0.;   tc2(10,13)= -tp;  tc2(10,14)=  t;   tc2(10,15)= -tp;
      tc2(11,8)= 0.;   tc2(11,9)= -tpp; tc2(11,10)=  t;  tc2(11,11)=  0.; tc2(11,12)= 0.;   tc2(11,13)= 0.;   tc2(11,14)= -tp;  tc2(11,15)=   t;
      tc2(12,8)=  t;   tc2(12,9)= -tp;  tc2(12,10)= 0.;  tc2(12,11)=  0.; tc2(12,12)= 0.;   tc2(12,13)=  t;   tc2(12,14)= -tpp; tc2(12,15)= 0.;
      tc2(13,8)= -tp;  tc2(13,9)=   t;  tc2(13,10)= -tp; tc2(13,11)=  0.; tc2(13,12)=  t;   tc2(13,13)= 0.;   tc2(13,14)=  t;   tc2(13,15)= -tpp;
      tc2(14,8)= 0.;   tc2(14,9)= -tp;  tc2(14,10)=  t;  tc2(14,11)= -tp; tc2(14,12)= -tpp; tc2(14,13)=  t;   tc2(14,14)=  0.;  tc2(14,15)=  t;
      tc2(15,8)= 0.;   tc2(15,9)= 0.;   tc2(15,10)= -tp; tc2(15,11)=   t; tc2(15,12)= 0.;   tc2(15,13)= -tpp; tc2(15,14)=  t;   tc2(15,15)= 0.;





      //used for exact lattice (Left-half) :
      dtk2(0,0)=  0.;          dtk2(0,1)=  0.;      dtk2(0,2)= -tpp*ex;   dtk2(0,3)=  t*ex;       dtk2(0,4)=  0.;       dtk2(0,5)=  0.;       dtk2(0,6)=  0.;       dtk2(0,7)= -tp*ex;
      dtk2(1,0)=  0.;          dtk2(1,1)=  0.;      dtk2(1,2)=  0.;       dtk2(1,3)= -tpp*ex;     dtk2(1,4)=  0.;       dtk2(1,5)=  0.;       dtk2(1,6)=  0.;       dtk2(1,7)=  0.;
      dtk2(2,0)= -tpp*emx;     dtk2(2,1)=  0.;      dtk2(2,2)=  0.;       dtk2(2,3)=  0.;         dtk2(2,4)=  0.;       dtk2(2,5)=  0.;       dtk2(2,6)=  0.;       dtk2(2,7)=  0.;
      dtk2(3,0)=  t*emx;       dtk2(3,1)= -tpp*emx; dtk2(3,2)=  0.;       dtk2(3,3)=  0.;         dtk2(3,4)= -tp*emx;   dtk2(3,5)=  0.;       dtk2(3,6)=  0.;       dtk2(3,7)=  0.;
      dtk2(4,0)=  0.;          dtk2(4,1)=  0.;      dtk2(4,2)=  0.;       dtk2(4,3)= -tp*ex;      dtk2(4,4)=  0.;       dtk2(4,5)=  0.;       dtk2(4,6)= -tpp*ex;   dtk2(4,7)=  t*ex;
      dtk2(5,0)=  0.;          dtk2(5,1)=  0.;      dtk2(5,2)=  0.;       dtk2(5,3)=  0.;         dtk2(5,4)=  0.;       dtk2(5,5)=  0.;       dtk2(5,6)=  0.;       dtk2(5,7)= -tpp*ex;
      dtk2(6,0)=  0.;          dtk2(6,1)=  0.;      dtk2(6,2)=  0.;       dtk2(6,3)=  0.;         dtk2(6,4)= -tpp*emx;  dtk2(6,5)=  0.;       dtk2(6,6)=  0.;       dtk2(6,7)=  0.;
      dtk2(7,0)= -tp*emx;      dtk2(7,1)=  0.;      dtk2(7,2)=  0.;       dtk2(7,3)=  0.;         dtk2(7,4)=  t*emx;    dtk2(7,5)= -tpp*emx;  dtk2(7,6)=  0.;       dtk2(7,7)=  0.;
      dtk2(8,0)= -tpp*emy;     dtk2(8,1)=  0.;      dtk2(8,2)=  0.;       dtk2(8,3)=  0.;         dtk2(8,4)=  0.;       dtk2(8,5)=  0.;       dtk2(8,6)=  0.;       dtk2(8,7)= -tp*ex;
      dtk2(9,0)=  0.;          dtk2(9,1)= -tpp*emy; dtk2(9,2)=  0.;       dtk2(9,3)=  0.;         dtk2(9,4)=  0.;       dtk2(9,5)=  0.;       dtk2(9,6)=  0.;       dtk2(9,7)=  0.;
      dtk2(10,0)= 0.;          dtk2(10,1)= 0.;      dtk2(10,2)= -tpp*emy; dtk2(10,3)= 0.;         dtk2(10,4)= 0.;       dtk2(10,5)= 0.;       dtk2(10,6)= 0.;       dtk2(10,7)= 0.;
      dtk2(11,0)= 0.;          dtk2(11,1)= 0.;      dtk2(11,2)= 0.;       dtk2(11,3)= -tpp*emy;   dtk2(11,4)= -tp*emx;  dtk2(11,5)= 0.;       dtk2(11,6)= 0.;       dtk2(11,7)= 0.;
      dtk2(12,0)=  t*emy;      dtk2(12,1)= -tp*emy; dtk2(12,2)= 0.;       dtk2(12,3)= -tp*ex*emy; dtk2(12,4)= -tpp*emy; dtk2(12,5)= 0.;       dtk2(12,6)= 0.;       dtk2(12,7)= 0.;
      dtk2(13,0)= -tp*emy;     dtk2(13,1)=  t*emy;  dtk2(13,2)= -tp*emy;  dtk2(13,3)= 0.;         dtk2(13,4)= 0.;       dtk2(13,5)= -tpp*emy; dtk2(13,6)= 0.;       dtk2(13,7)= 0.;
      dtk2(14,0)= 0.;          dtk2(14,1)= -tp*emy; dtk2(14,2)=   t*emy;  dtk2(14,3)= -tp*emy;    dtk2(14,4)= 0.;       dtk2(14,5)= 0.;       dtk2(14,6)= -tpp*emy; dtk2(14,7)= 0.;
      dtk2(15,0)= -tp*emy*emx; dtk2(15,1)= 0.;      dtk2(15,2)= -tp*emy;  dtk2(15,3)=  t*emy;     dtk2(15,4)= 0.;       dtk2(15,5)= 0.;       dtk2(15,6)= 0.;       dtk2(15,7)= -tpp*emy;
      //Right-half :
      dtk2(0,8)= -tpp*ey;   dtk2(0,9)=  0.;       dtk2(0,10)=  0.;     dtk2(0,11)=  0.;     dtk2(0,12)=  t*ey;      dtk2(0,13)= -tp*ey;    dtk2(0,14)=  0.;      dtk2(0,15)= -tp*ex*ey;
      dtk2(1,8)=  0.;       dtk2(1,9)= -tpp*ey;   dtk2(1,10)=  0.;     dtk2(1,11)=  0.;     dtk2(1,12)= -tp*ey;     dtk2(1,13)=  t*ey;     dtk2(1,14)= -tp*ey;   dtk2(1,15)=  0.;
      dtk2(2,8)=  0.;       dtk2(2,9)=  0.;       dtk2(2,10)= -tpp*ey; dtk2(2,11)=  0.;     dtk2(2,12)=  0.;        dtk2(2,13)= -tp*ey;    dtk2(2,14)=  t*ey;    dtk2(2,15)= -tp*ey;
      dtk2(3,8)=  0.;       dtk2(3,9)=  0.;       dtk2(3,10)=  0.;     dtk2(3,11)= -tpp*ey; dtk2(3,12)= -tp*emx*ey; dtk2(3,13)=  0.;       dtk2(3,14)= -tp*ey;   dtk2(3,15)=   t*ey;
      dtk2(4,8)=  0.;       dtk2(4,9)=  0.;       dtk2(4,10)=  0.;     dtk2(4,11)= -tp*ex;  dtk2(4,12)= -tpp*ey;    dtk2(4,13)=  0.;       dtk2(4,14)=  0.;      dtk2(4,15)=  0.;
      dtk2(5,8)=  0.;       dtk2(5,9)=  0.;       dtk2(5,10)=  0.;     dtk2(5,11)=  0.;     dtk2(5,12)=  0.;        dtk2(5,13)= -tpp*ey;   dtk2(5,14)=  0.;      dtk2(5,15)=  0.;
      dtk2(6,8)=  0.;       dtk2(6,9)=  0.;       dtk2(6,10)=  0.;     dtk2(6,11)=  0.;     dtk2(6,12)=  0.;        dtk2(6,13)=  0.;       dtk2(6,14)= -tpp*ey;  dtk2(6,15)=  0.;
      dtk2(7,8)= -tp*emx;   dtk2(7,9)=  0.;       dtk2(7,10)=  0.;     dtk2(7,11)=  0.;     dtk2(7,12)=  0.;        dtk2(7,13)=  0.;       dtk2(7,14)=  0.;      dtk2(7,15)= -tpp*ey;
      dtk2(8,8)=  0.;       dtk2(8,9)=  0.;       dtk2(8,10)= -tpp*ex; dtk2(8,11)=  t*ex;   dtk2(8,12)=  0.;        dtk2(8,13)=  0.;       dtk2(8,14)=  0.;      dtk2(8,15)= -tp*ex;
      dtk2(9,8)=  0.;       dtk2(9,9)=  0.;       dtk2(9,10)=  0.;     dtk2(9,11)= -tpp*ex; dtk2(9,12)=  0.;        dtk2(9,13)=  0.;       dtk2(9,14)=  0.;      dtk2(9,15)=  0.;
      dtk2(10,8)= -tpp*emx; dtk2(10,9)= 0.;       dtk2(10,10)= 0.;     dtk2(10,11)= 0.;     dtk2(10,12)= 0.;        dtk2(10,13)= 0.;       dtk2(10,14)= 0.;      dtk2(10,15)= 0.;
      dtk2(11,8)=  t*emx;   dtk2(11,9)= -tpp*emx; dtk2(11,10)= 0.;     dtk2(11,11)= 0.;     dtk2(11,12)= -tp*emx;   dtk2(11,13)= 0.;       dtk2(11,14)= 0.;      dtk2(11,15)= 0.;
      dtk2(12,8)= 0.;       dtk2(12,9)= 0.;       dtk2(12,10)= 0.;     dtk2(12,11)= -tp*ex; dtk2(12,12)= 0.;        dtk2(12,13)= 0.;       dtk2(12,14)= -tpp*ex; dtk2(12,15)=  t*ex;
      dtk2(13,8)= 0.;       dtk2(13,9)= 0.;       dtk2(13,10)= 0.;     dtk2(13,11)= 0.;     dtk2(13,12)= 0.;        dtk2(13,13)= 0.;       dtk2(13,14)= 0.;      dtk2(13,15)= -tpp*ex;
      dtk2(14,8)= 0.;       dtk2(14,9)= 0.;       dtk2(14,10)= 0.;     dtk2(14,11)= 0.;     dtk2(14,12)= -tpp*emx;  dtk2(14,13)= 0.;       dtk2(14,14)= 0.;      dtk2(14,15)= 0.;
      dtk2(15,8)= -tp*emx;  dtk2(15,9)= 0.;       dtk2(15,10)= 0.;     dtk2(15,11)= 0.;     dtk2(15,12)=  t*emx;    dtk2(15,13)= -tpp*emx; dtk2(15,14)= 0.;      dtk2(15,15)= 0.;






      //used for compact tiling:
      dtk3(0,0)=  0.;          dtk3(0,1)=  0.;      dtk3(0,2)= -tpp*ex;   dtk3(0,3)=  t*ex;       dtk3(0,4)=  0.;       dtk3(0,5)=  0.;       dtk3(0,6)=  0.;       dtk3(0,7)= -tp*ex;
      dtk3(1,0)=  0.;          dtk3(1,1)=  0.;      dtk3(1,2)=  0.;       dtk3(1,3)= -tpp*ex;     dtk3(1,4)=  0.;       dtk3(1,5)=  0.;       dtk3(1,6)=  0.;       dtk3(1,7)=  0.;
      dtk3(2,0)= -tpp*emx;     dtk3(2,1)=  0.;      dtk3(2,2)=  0.;       dtk3(2,3)=  0.;         dtk3(2,4)=  0.;       dtk3(2,5)=  0.;       dtk3(2,6)=  0.;       dtk3(2,7)=  0.;
      dtk3(3,0)=  t*emx;       dtk3(3,1)= -tpp*emx; dtk3(3,2)=  0.;       dtk3(3,3)=  0.;         dtk3(3,4)= -tp*emx;   dtk3(3,5)=  0.;       dtk3(3,6)=  0.;       dtk3(3,7)=  0.;
      dtk3(4,0)=  0.;          dtk3(4,1)=  0.;      dtk3(4,2)=  0.;       dtk3(4,3)= -tp*ex;      dtk3(4,4)=  0.;       dtk3(4,5)=  0.;       dtk3(4,6)= -tpp*ex;   dtk3(4,7)=  t*ex;
      dtk3(5,0)=  0.;          dtk3(5,1)=  0.;      dtk3(5,2)=  0.;       dtk3(5,3)=  0.;         dtk3(5,4)=  0.;       dtk3(5,5)=  0.;       dtk3(5,6)=  0.;       dtk3(5,7)= -tpp*ex;
      dtk3(6,0)=  0.;          dtk3(6,1)=  0.;      dtk3(6,2)=  0.;       dtk3(6,3)=  0.;         dtk3(6,4)= -tpp*emx;  dtk3(6,5)=  0.;       dtk3(6,6)=  0.;       dtk3(6,7)=  0.;
      dtk3(7,0)= -tp*emx;      dtk3(7,1)=  0.;      dtk3(7,2)=  0.;       dtk3(7,3)=  0.;         dtk3(7,4)=  t*emx;    dtk3(7,5)= -tpp*emx;  dtk3(7,6)=  0.;       dtk3(7,7)=  0.;
      dtk3(8,0)= -tpp*emy;     dtk3(8,1)=  0.;      dtk3(8,2)=  0.;       dtk3(8,3)=  0.;         dtk3(8,4)=  0.;       dtk3(8,5)=  0.;       dtk3(8,6)=  0.;       dtk3(8,7)= -tp*ex;
      dtk3(9,0)=  0.;          dtk3(9,1)= -tpp*emy; dtk3(9,2)=  0.;       dtk3(9,3)=  0.;         dtk3(9,4)=  0.;       dtk3(9,5)=  0.;       dtk3(9,6)=  0.;       dtk3(9,7)=  0.;
      dtk3(10,0)= 0.;          dtk3(10,1)= 0.;      dtk3(10,2)= -tpp*emy; dtk3(10,3)= 0.;         dtk3(10,4)= 0.;       dtk3(10,5)= 0.;       dtk3(10,6)= 0.;       dtk3(10,7)= 0.;
      dtk3(11,0)= 0.;          dtk3(11,1)= 0.;      dtk3(11,2)= 0.;       dtk3(11,3)= -tpp*emy;   dtk3(11,4)= -tp*emx;  dtk3(11,5)= 0.;       dtk3(11,6)= 0.;       dtk3(11,7)= 0.;
      dtk3(12,0)=  t*emy;      dtk3(12,1)= -tp*emy; dtk3(12,2)= 0.;       dtk3(12,3)= -tp*ex*emy; dtk3(12,4)= -tpp*emy; dtk3(12,5)= 0.;       dtk3(12,6)= 0.;       dtk3(12,7)= 0.;
      dtk3(13,0)= -tp*emy;     dtk3(13,1)=  t*emy;  dtk3(13,2)= -tp*emy;  dtk3(13,3)= 0.;         dtk3(13,4)= 0.;       dtk3(13,5)= -tpp*emy; dtk3(13,6)= 0.;       dtk3(13,7)= 0.;
      dtk3(14,0)= 0.;          dtk3(14,1)= -tp*emy; dtk3(14,2)=   t*emy;  dtk3(14,3)= -tp*emy;    dtk3(14,4)= 0.;       dtk3(14,5)= 0.;       dtk3(14,6)= -tpp*emy; dtk3(14,7)= 0.;
      dtk3(15,0)= -tp*emy*emx; dtk3(15,1)= 0.;      dtk3(15,2)= -tp*emy;  dtk3(15,3)=  t*emy;     dtk3(15,4)= 0.;       dtk3(15,5)= 0.;       dtk3(15,6)= 0.;       dtk3(15,7)= -tpp*emy;
      //Right-half 
      dtk3(0,8)= -tpp*ey;   dtk3(0,9)=  0.;       dtk3(0,10)=  0.;     dtk3(0,11)=  0.;     dtk3(0,12)=  t*ey;      dtk3(0,13)= -tp*ey;    dtk3(0,14)=  0.;      dtk3(0,15)= -tp*ex*ey;
      dtk3(1,8)=  0.;       dtk3(1,9)= -tpp*ey;   dtk3(1,10)=  0.;     dtk3(1,11)=  0.;     dtk3(1,12)= -tp*ey;     dtk3(1,13)=  t*ey;     dtk3(1,14)= -tp*ey;   dtk3(1,15)=  0.;
      dtk3(2,8)=  0.;       dtk3(2,9)=  0.;       dtk3(2,10)= -tpp*ey; dtk3(2,11)=  0.;     dtk3(2,12)=  0.;        dtk3(2,13)= -tp*ey;    dtk3(2,14)=  t*ey;    dtk3(2,15)= -tp*ey;
      dtk3(3,8)=  0.;       dtk3(3,9)=  0.;       dtk3(3,10)=  0.;     dtk3(3,11)= -tpp*ey; dtk3(3,12)= -tp*emx*ey; dtk3(3,13)=  0.;       dtk3(3,14)= -tp*ey;   dtk3(3,15)=   t*ey;
      dtk3(4,8)=  0.;       dtk3(4,9)=  0.;       dtk3(4,10)=  0.;     dtk3(4,11)= -tp*ex;  dtk3(4,12)= -tpp*ey;    dtk3(4,13)=  0.;       dtk3(4,14)=  0.;      dtk3(4,15)=  0.;
      dtk3(5,8)=  0.;       dtk3(5,9)=  0.;       dtk3(5,10)=  0.;     dtk3(5,11)=  0.;     dtk3(5,12)=  0.;        dtk3(5,13)= -tpp*ey;   dtk3(5,14)=  0.;      dtk3(5,15)=  0.;
      dtk3(6,8)=  0.;       dtk3(6,9)=  0.;       dtk3(6,10)=  0.;     dtk3(6,11)=  0.;     dtk3(6,12)=  0.;        dtk3(6,13)=  0.;       dtk3(6,14)= -tpp*ey;  dtk3(6,15)=  0.;
      dtk3(7,8)= -tp*emx;   dtk3(7,9)=  0.;       dtk3(7,10)=  0.;     dtk3(7,11)=  0.;     dtk3(7,12)=  0.;        dtk3(7,13)=  0.;       dtk3(7,14)=  0.;      dtk3(7,15)= -tpp*ey;
      dtk3(8,8)=  0.;       dtk3(8,9)=  0.;       dtk3(8,10)= -tpp*ex; dtk3(8,11)=  t*ex;   dtk3(8,12)=  0.;        dtk3(8,13)=  0.;       dtk3(8,14)=  0.;      dtk3(8,15)= -tp*ex;
      dtk3(9,8)=  0.;       dtk3(9,9)=  0.;       dtk3(9,10)=  0.;     dtk3(9,11)= -tpp*ex; dtk3(9,12)=  0.;        dtk3(9,13)=  0.;       dtk3(9,14)=  0.;      dtk3(9,15)=  0.;
      dtk3(10,8)= -tpp*emx; dtk3(10,9)= 0.;       dtk3(10,10)= 0.;     dtk3(10,11)= 0.;     dtk3(10,12)= 0.;        dtk3(10,13)= 0.;       dtk3(10,14)= 0.;      dtk3(10,15)= 0.;
      dtk3(11,8)=  t*emx;   dtk3(11,9)= -tpp*emx; dtk3(11,10)= 0.;     dtk3(11,11)= 0.;     dtk3(11,12)= -tp*emx;   dtk3(11,13)= 0.;       dtk3(11,14)= 0.;      dtk3(11,15)= 0.;
      dtk3(12,8)= 0.;       dtk3(12,9)= 0.;       dtk3(12,10)= 0.;     dtk3(12,11)= -tp*ex; dtk3(12,12)= 0.;        dtk3(12,13)= 0.;       dtk3(12,14)= -tpp*ex; dtk3(12,15)=  t*ex;
      dtk3(13,8)= 0.;       dtk3(13,9)= 0.;       dtk3(13,10)= 0.;     dtk3(13,11)= 0.;     dtk3(13,12)= 0.;        dtk3(13,13)= 0.;       dtk3(13,14)= 0.;      dtk3(13,15)= -tpp*ex;
      dtk3(14,8)= 0.;       dtk3(14,9)= 0.;       dtk3(14,10)= 0.;     dtk3(14,11)= 0.;     dtk3(14,12)= -tpp*emx;  dtk3(14,13)= 0.;       dtk3(14,14)= 0.;      dtk3(14,15)= 0.;
      dtk3(15,8)= -tp*emx;  dtk3(15,9)= 0.;       dtk3(15,10)= 0.;     dtk3(15,11)= 0.;     dtk3(15,12)=  t*emx;    dtk3(15,13)= -tpp*emx; dtk3(15,14)= 0.;      dtk3(15,15)= 0.;
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
      double Ry[16] = {0.,0.,0.,0.,1.,1.,1.,1.,2.,2.,2.,2.,3.,3.,3.,3.};
      double Rx[16] = {0.,1.,2.,3.,0.,1.,2.,3.,0.,1.,2.,3.,0.,1.,2.,3.};
      complex<double> z(OMEGA,ETA);
      
      calculate_dtk(px, py);
      calculate_sigma(z);

      if (periodization==1){
        calculate_cumulant(z);
        M_per = 0;
        for (int ii=0; ii<16; ++ii) {
          for (int jj=0; jj<16; ++jj) {
            double arg = ((Rx[jj]-Rx[ii])*px + (Ry[jj]-Ry[ii])*py);
            complex<double> phase(cos(arg), sin(arg));
            M_per += 0.0625 * cumul(ii,jj) * phase; 
          }
        }
        double epsilon_k = -2*t*(cos(px)+cos(py)) -4*tp*cos(px)*cos(py) -2*tpp*(cos(2*px)+cos(2*py));
        G_per = 1./((1./M_per) - epsilon_k);
      }
      else {
        calculate_Gk(z);
        G_per = 0;
        for (int ii=0; ii<16; ++ii) {
          for (int jj=0; jj<16; ++jj) {
            double arg = ((Rx[jj]-Rx[ii])*px + (Ry[jj]-Ry[ii])*py);
            complex<double> phase(cos(arg), sin(arg));
            G_per += 0.0625 * green(ii,jj) * phase; 
          }
        }
      }
    }
} Model;


