#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

#include "../include/LSODA.h"

#include "../include/asa047.hpp"

#include "../boost/math/interpolators/pchip.hpp"

#define G 6.6732e-8
#define C 2.9979e+10
#define MSUN 1.989e+33
#define PI 3.1415926535
#define MB 1.66e-24
#define N 45001  // grid size

// #define __USE_MINGW_ANSI_STDIO  (__MINGW_FEATURES__ & __MINGW_ANSI_STDIO__)


char out_name[80];
char gr_name[80];

double Length = G*MSUN/pow(C,2),
Time = Length/C,
Density = MSUN/pow(Length,3),
M0,
DI;
std::vector<double> M,
		M2,
		R,
		R2,
		q,
		I,
		I2;              


void load_out( char out_file[])
{
 std::size_t i;                    /* counter */

int n_tab;   
      
 FILE *f_out;              /* pointer to eos_file */
double Me,Re,qe,je,Ie; 

    /* OPEN FILE TO READ */
    char out[100];
    strcpy(out,"outr/");  /* add path */
    strcat(out,out_file);
	
    if((f_out=fopen(out,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",out_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_out,"%d\n",&n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(int i=1;i<=n_tab;i++) {  
      /*  input line for _four_column equation of state file */
       fscanf(f_out,"%lf %lf %lf %lf %lf\n",&Me,&Re,&Ie,&je,&qe) ; 
       M.push_back(Me);
       R.push_back(Re);
       I.push_back(Ie);
       q.push_back(qe);

    }
}

void load_gr( char gr_file[])
{
 std::size_t i;                    /* counter */

int ntab;
      
 FILE *f_gr;              /* pointer to eos_file */
double Mg,Rg,Ig,Jg; 

    /* OPEN FILE TO READ */
    char gr[100];
    strcpy(gr,"outr/");  /* add path */
    strcat(gr,gr_file);
	
    if((f_gr=fopen(gr,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",gr_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_gr,"%d\n",&ntab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(int i=1;i<=ntab;i++) {  
      /*  input line for _four_column equation of state file */
       fscanf(f_gr,"%lf %lf %lf %lf \n",&Mg,&Rg,&Ig,&Jg) ; 
       M2.push_back(Mg);
       R2.push_back(Rg);
       I2.push_back(Ig);
    }
}



double mr(double rr)
{
 	auto spline = boost::math::interpolators::pchip<decltype(R)>(R,M);
 	return spline(rr);
}
double im(double mm)
{
 	auto spline = boost::math::interpolators::pchip<decltype(M)>(M,I);
 	return spline(mm);
}

double mr2(double rr2)
{
 	auto spline = boost::math::interpolators::pchip<decltype(R2)>(R2,M2);
 	return spline(rr2);
}
double im2(double mm2)
{
 	auto spline = boost::math::interpolators::pchip<decltype(M2)>(M2,I2);
 	return spline(mm2);
}
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

int main(int argc, char const *argv[])
{
   for(int i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
        switch(argv[i][1]){
          
          case 'f':
              sscanf(argv[i+1],"%s",out_name);
              break;
          case 's':
              sscanf(argv[i+1],"%s",gr_name);
              break;
          case 'p':
              sscanf(argv[i+1],"%lf",&M0);
              break;
  		
  		}
  	}
  	
  load_out(out_name);
  
    load_gr(gr_name);

 DI=(im(M0)-im2(M0))/im2(M0);
 if(DI<0){
    printf("t");
     DI=-DI;
  }
 

    printf("\n");
    printf("  %s                  \n",out_name);    
    printf("  %s                  \n",gr_name);   
    printf("  I_m= %6.5e    \n",im(M0));
    printf("DI= %6.5e  \n", DI);
    cout<<fixed<<"test: "<<DI<<endl;

//for(int i=3;i<=70;i++){
//    printf("%6.5e  %6.5e  \n",im(M[i]),M[i]);
//    }
   
  return 0;
}
