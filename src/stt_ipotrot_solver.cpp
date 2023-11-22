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


char eos_name[80];

double central_density,
alpha_phi = -1/pow(3,0.5),
Length = G*MSUN/pow(C,2),
Time = Length/C,
Density = MSUN/pow(Length,3),
scalar_charge,
Mass,                                // gravitational mass
Radius,
W_kepler,
W_last,
ratio,
J,
I;                            



std::vector<double> log_e_tab,               /* rho points in tabulated EOS */
	log_p_tab,               /* p points in tabulated EOS */
	log_h_tab,              /* h points in EOS file */
	log_n0_tab,         /* number density in EOS file */  
	log_rho_tab;


double coupling;

int n_tab,                           /* number of tabulated EOS points */
    n_nearest=1,                     /* nearest grid point, used in interp. */ 
    print_option;                    /* select print out */ 

void load_eos( char eos_file[])
{
 std::size_t i;                    /* counter */

 double p,                 /* pressure */
        rho,               /* density */
        h,                 /* enthalpy */
        n0;                /* number density */    
      
 FILE *f_eos;              /* pointer to eos_file */
 

    /* OPEN FILE TO READ */
    char eos[100];
    strcpy(eos,"EoS/");  /* add path */
    strcat(eos,eos_file);
	
    if((f_eos=fopen(eos,"r")) == NULL ) {    
       printf("cannot open file:  %s\n",eos_file); 
       exit(0);
    }

 
    /* READ NUMBER OF TABULATED POINTS */

    fscanf(f_eos,"%d\n",&n_tab);


    /* READ EOS, H, N0 AND MAKE THEM DIMENSIONLESS (EXCEPT N0) */
 
    for(int i=1;i<=n_tab;i++) {  
      /*  input line for _four_column equation of state file */
       fscanf(f_eos,"%lf %lf %lf %lf\n",&rho,&p,&h,&n0) ; 

/*  input line for _five_ column equation of state file         
    fscanf(f_eos,"%lf %lf %lf %lf %*lf\n",&rho,&p,&h,&n0) ; */

       log_e_tab.push_back(log10(rho/Density));       /* multiply by C^2 to get */ 
       log_p_tab.push_back(log10(p/(Density*pow(C,2))));             /* energy density. */
       log_h_tab.push_back(log10(h/(C*C)));        
       log_n0_tab.push_back(log10(n0));
       log_rho_tab.push_back(log10(MB*n0*pow(Length,3)/MSUN));             /* STILL IN CGS ! */
    }
}

/*******************************************************************/
double e_at_p(double pp)
{
 	auto spline = boost::math::interpolators::pchip<decltype(log_p_tab)>(log_p_tab,log_e_tab);
	 return pow(10.0,spline(log10(pp)));
}
/*******************************************************************/
double e_at_rho(double rhorho)
{
	auto spline = boost::math::interpolators::pchip<decltype(log_rho_tab)>(log_rho_tab,log_e_tab);
	return pow(10.0,spline(log10(rhorho)));
}
/*******************************************************************/
double p_at_e(double ee)
{
	auto spline = boost::math::interpolators::pchip<decltype(log_e_tab)>(log_e_tab,log_p_tab);
	return pow(10.0,spline(log10(ee)));
}
/*******************************************************************/
double p_at_rho(double rhorho)
{
	auto spline = boost::math::interpolators::pchip<decltype(log_rho_tab)>(log_rho_tab,log_p_tab);
	return pow(10.0,spline(log10(rhorho)));
}
/*******************************************************************/
double p_at_h(double hh)
{ 
	auto spline = boost::math::interpolators::pchip<decltype(log_h_tab)>(log_h_tab,log_p_tab);
	return pow(10.0,spline(log10(hh)));
}
/*******************************************************************/
double h_at_p(double pp)
{ 
	auto spline = boost::math::interpolators::pchip<decltype(log_p_tab)>(log_p_tab,log_h_tab);
	return pow(10.0,spline(log10(pp)));
}
/*******************************************************************/
double n0_at_e(double ee)
{ 
	auto spline = boost::math::interpolators::pchip<decltype(log_e_tab)>(log_e_tab,log_n0_tab);
	return pow(10.0,spline(log10(ee)));
}
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
double A_phi(double phi,double alpha){
  return exp(alpha_phi*phi);
}
/*******************************************************************/

double V_phi(double phi, double alpha){
   return (1/(4*alpha))*pow(1 - exp(2*alpha_phi*phi),2);
  
}
/*******************************************************************/
double dV_phi(double phi,double alpha){
  return (-alpha_phi/alpha)*(exp(alpha_phi*2*phi) - exp(alpha_phi*4*phi));
  
}
/*******************************************************************/
/*******************************************************************/
/*******************************************************************/
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);

  return linspaced;
}
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
static void system0_inside(double ru, double* u, double* uprime,double alpha)
{
    (void)ru;
    
    double A = A_phi(u[4],alpha);
    double V = V_phi(u[4],alpha);
    double dV = dV_phi(u[4],alpha);

    double eps = e_at_p(u[3]);


    uprime[0] = 4*PI*pow(A,4)*pow(ru,2.0)*eps + 0.5*ru*(ru-2*u[0])*pow(u[2],2) + 0.25*V*pow(ru,2);
    uprime[1] = 2.0*(u[0] + 4.0*PI*pow(A,4)*pow(ru,3.0)*u[3] - 0.25*V*pow(ru,3.0))/( ru*(ru-2.0*u[0])) + ru*pow(u[2],2);
    uprime[2] = 4*PI*pow(A,4)*ru*(alpha_phi*(eps-3.0*u[3]) + ru*(eps-u[3])*u[2] )/(ru-2.0*u[0]) - 2*u[2]*(1 - u[0]/ru)/(ru - 2.0*u[0]) + (0.5*V*u[2]*pow(ru,2.0) + 0.25*dV*ru)/(ru - 2.0*u[0]);
    uprime[3] = -(eps + u[3])*((u[0] + 4.0*PI*pow(A,4)*pow(ru,3.0)*u[3] - 0.25*V*pow(ru,3.0))/(ru*(ru-2.0*u[0])) + 0.5*ru*pow(u[2],2) + alpha_phi*u[2]);
    uprime[4] = u[2];
            
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

static void system0_outside(double ru, double* u, double* uprime,double alpha)
{
    (void)ru;
    

    double A = A_phi(u[4],alpha);
    double V = V_phi(u[4],alpha);
    double dV = dV_phi(u[4],alpha);

    
    uprime[0] = 0.5*ru*(ru-2*u[0])*pow(u[2],2) + 0.25*V*pow(ru,2);
    uprime[1] = 2.0*(u[0] - 0.25*V*pow(ru,3.0))/(ru*(ru-2.0*u[0])) + ru*pow(u[2],2);   
    uprime[2] = -2.0*u[2]*(1 - u[0]/ru)/(ru - 2.0*u[0]) + (0.5*V*u[2]*pow(ru,2.0) + 0.25*dV*ru)/(ru - 2*u[0]);
    uprime[3] = 0;
    uprime[4] = u[2];  
            
}

/**************************************************************************************************************************************/


static void system_inside(double r, double* y, double* yprime,double alpha)
{
    (void)r;
    
    double A = A_phi(y[5],alpha);
    double V = V_phi(y[5],alpha);
    double dV = dV_phi(y[5],alpha);

    double eps = e_at_p(y[3]);


    yprime[0] = 4*PI*pow(A,4)*pow(r,2.0)*eps + 0.5*r*(r-2.0*y[0])*pow(y[2],2) + 0.25*V*pow(r,2);   // y[0]=m
   
    yprime[1] = 2.0*(y[0] + 4.0*PI*pow(A,4)*pow(r,3.0)*y[3] - 0.25*V*pow(r,3.0))/( r*(r-2.0*y[0])) + r*pow(y[2],2);  //y[1]=nu

    yprime[2] =  4*PI*pow(A,4)*r*(alpha_phi*(eps-3.0*y[3]) + r*(eps-y[3])*y[2] )/(r-2.0*y[0]) - 2*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]) + (0.5*V*y[2]*pow(r,2.0) + 0.25*dV*r)/(r - 2.0*y[0]);                                                                   //y[2]=dphi

    yprime[3] = -(eps + y[3])*((y[0] + 4.0*PI*pow(A,4)*pow(r,3.0)*y[3] - 0.25*V*pow(r,3.0))/(r*(r-2.0*y[0])) + 0.5*r*pow(y[2],2) + alpha_phi*y[2]);                                                                                                     //y[3]=P

    yprime[4] = (4.0*PI*pow(A,4)*(eps+y[3])*(4.0*y[6]+r*y[4])*pow(r,2) + (r-2.0*y[0])*(-4.0 + pow(r,2)*pow(y[2],2))*y[4])/(r*(r-2.0*y[0]));//y[4]=dw
    
    yprime[5] = y[2];   //y[5]=phi
       
    yprime[6] = y[4];  //y[6]=w
                
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

static void system_outside(double r, double* y, double* yprime,double alpha)
{
    (void)r;
    

    double A = A_phi(y[5],alpha);
    double V = V_phi(y[5],alpha);
    double dV = dV_phi(y[5],alpha);

    
    yprime[0] = 0.5*r*(r-2*y[0])*pow(y[2],2) + 0.25*V*pow(r,2);

    yprime[1] = 2.0*(y[0] - 0.25*V*pow(r,3.0))/(r*(r-2.0*y[0])) + r*pow(y[2],2);
    
    yprime[2] = -2.0*y[2]*(1 - y[0]/r)/(r - 2.0*y[0]) + (0.5*V*y[2]*pow(r,2.0) + 0.25*dV*r)/(r - 2*y[0]);

    yprime[3] = 0;

    yprime[4] = (r-2.0*y[0])*(-4.0 + pow(r,2)*pow(y[2],2))*y[4]/(r*(r-2.0*y[0]));


    yprime[5] = y[2];
    
    yprime[6] = y[4];

}
/**************************************************************************************************************************************/

std::vector<double> gradient(std::vector<double> input, double h){
    if (input.size() <= 1) return input;
    std::vector<double> res;
    for(int j=0; j<input.size(); j++) {
        int j_left = j - 1;
        int j_right = j + 1;
        if (j_left < 0) {
            j_left = 0;
            j_right = 1;
        }
        if (j_right >= input.size()){
            j_right = input.size() - 1;
            j_left = j_right - 1;
        }
        // gradient value at position j
        double dist_grad = (input[j_right] - input[j_left]) / (2.0*h);
        res.push_back(dist_grad);
    }
    return res;
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

double function_to_minimize(double args[2])
{
  double nu_c_obj,
  phi_c_obj,
  eps_c_obj,
  p_c_obj,
  r_max_obj;

  nu_c_obj = args[0];
  phi_c_obj = args[1];

  
  
  eps_c_obj = central_density*pow(10,15)/Density;
  p_c_obj = p_at_e(eps_c_obj);

	
  if (coupling >= 1.0 && coupling < 20.0)
  {
    r_max_obj = 20.0;
  }

  else if (coupling>pow(10,3))
  {
    r_max_obj = 500.0;
  }
  
  else 
  {
    r_max_obj = coupling/2;
  }


  if (r_max_obj < pow(3.0/(4.0*PI*central_density),1.0/3.0)){

    r_max_obj = pow(3.0/(4.0*PI*central_density),1.0/3.0);

  }



  int neq0 = 5;
  double t_obj, tout_obj;
  int istate_obj = 1;

	std::vector<double> r_vec_obj = linspace(0.0,r_max_obj,N);

	double dr_obj = r_vec_obj[1]-r_vec_obj[0];

	std::vector<double> y_obj = {0.0,nu_c_obj,0.0,p_at_e(eps_c_obj),phi_c_obj};

	t_obj = dr_obj;
	tout_obj = t_obj+dr_obj;

  LSODA lsoda;

  std::vector<double> yout_obj,
    m_obj,
    nu_obj,
    dfi_obj,
    pr_obj,
    fi_obj;

  m_obj.push_back(y_obj[0]);
  nu_obj.push_back(y_obj[1]);
  dfi_obj.push_back(y_obj[2]);
  pr_obj.push_back(y_obj[3]);
  fi_obj.push_back(y_obj[4]);

  int idx_obj=2;

  while(y_obj[3]>=1e-13 && t_obj<r_vec_obj[N-2]){

    lsoda.lsoda_update(system0_inside, neq0, y_obj, yout_obj, &t_obj, tout_obj, &istate_obj, coupling,1e-06,1e-20);
      m_obj.push_back(yout_obj[1]);
      nu_obj.push_back(yout_obj[2]);
      dfi_obj.push_back(yout_obj[3]);
      pr_obj.push_back(yout_obj[4]);
      fi_obj.push_back(yout_obj[5]);

      if(istate_obj==-3 || isnan(nu_obj[idx_obj-2]) || isnan(fi_obj[idx_obj-2])){
        break;
      }
      
      tout_obj += dr_obj;

      y_obj[0] = yout_obj[1];
      y_obj[1] = yout_obj[2];
      y_obj[2] = yout_obj[3];
      y_obj[3] = yout_obj[4];
      y_obj[4] = yout_obj[5];

      idx_obj += 1;

  }

  int idxlast_obj = idx_obj-1;

  while(t_obj<r_vec_obj[N-2]){

    lsoda.lsoda_update(system0_outside, neq0, y_obj, yout_obj, &t_obj, tout_obj, &istate_obj, coupling,1e-06,1e-20);
        m_obj.push_back(yout_obj[1]);
        nu_obj.push_back(yout_obj[2]);
        dfi_obj.push_back(yout_obj[3]);
        pr_obj.push_back(yout_obj[4]);
        fi_obj.push_back(yout_obj[5]);

      if(istate_obj==-3 || isnan(nu_obj[idx_obj-2]) || isnan(fi_obj[idx_obj-2])){
        break;
      }
   
      tout_obj += dr_obj;

      y_obj[0] = yout_obj[1];
      y_obj[1] = yout_obj[2];
      y_obj[2] = yout_obj[3];
      y_obj[3] = yout_obj[4];
      y_obj[4] = yout_obj[5];

      idx_obj += 1;    	
  }

return pow(nu_obj[idx_obj-3],2)+pow(fi_obj[idx_obj-3],2);
}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

void find_min(double (&min_vals)[2])
{
  int number_of_variables = 2;

  double *starting_points;
  starting_points = new double[number_of_variables];

  double *xmin;
  double ynewlo;
  double *step;

  xmin = new double[number_of_variables];
  step = new double[number_of_variables];

  int icount;
  int konvge,kcount;
  int numres;
  int ifault;

  double reqmin = 1.0E-08;

  starting_points[0] = -1.0;
  starting_points[1] = 0.1;

  step[0] = 1.0;
  step[1] = 0.1;

  konvge = 10;
  kcount = 500;
  
  nelmin (function_to_minimize, 2, starting_points, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault);

 if (abs(ynewlo)>1e-06){

    double *starting_points2;

    starting_points2 = new double[number_of_variables];

    starting_points2[0] = -1.0;
    starting_points2[1] = -0.1;

    double *xmin2;
    double ynewlo2;
    double *step2;

    double reqmin2 = 1.0E-08;

    xmin2 = new double[number_of_variables];
    step2 = new double[number_of_variables];

    int icount2;
    int konvge2,kcount2;
    int numres2;
    int ifault2;

    step2[0] = 1.0;
    step2[1] = 0.1;

    konvge2 = 10;
    kcount2 = 500;

    nelmin (function_to_minimize, 2, starting_points2, xmin2, &ynewlo2, reqmin2, step2,
    konvge2, kcount2, &icount2, &numres2, &ifault2);

    std::cout << "\n";
    std::cout << "\n";

    std::cout << "  Estimate of minimizing value X*:\n";
    std::cout << "\n";
    for (std::size_t i = 0; i < number_of_variables; i++ )
    {
      if (i==0)
      {
        cout << "   "<<"metric_phi_central = " << setw(7) << xmin2[i] << "\n";
      }
      else{
        cout << "   "<<"scalar_central = " << setw(7) << xmin2[i] << "\n";
      }
    }

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "  F(X*) = " << ynewlo2 << "\n";

    std::cout << "\n";
    std::cout << "  Number of Nelder-Mead iterations = " << icount2 << "\n";
    std::cout << "  Number of Nelder-Mead restarts =   " << numres2 << "\n";
    std::cout << "\n";
    std::cout << "\n";

    min_vals[0] = xmin2[0];
    min_vals[1] = xmin2[1];
    printf("-------------------------------------------------------------------------------\n");  

    delete [] starting_points2;
    delete [] step2;
    delete [] xmin2;

  }

  else{

    std::cout << "\n";
    std::cout << "\n";

    std::cout << "  Estimate of minimizing value X*:\n";
    std::cout << "\n";
    for (std::size_t i = 0; i < number_of_variables; i++ )
    {
      if (i==0)
      {
        std::cout << "   "<<"metric_phi_central = " << setw(7) << xmin[i] << "\n";
      }
      else{
        std::cout << "   "<<"scalar_central = " << setw(7) << xmin[i] << "\n";
      }
    }

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "  F(X*) = " << ynewlo << "\n";

    std::cout << "\n";
    std::cout << "  Number of Nelder-Mead iterations = " << icount << "\n";
    std::cout << "  Number of Nelder-Mead restarts =   " << numres << "\n";
    std::cout << "\n";
    std::cout << "\n";

    min_vals[0] = xmin[0];
    min_vals[1] = xmin[1];
    printf("-------------------------------------------------------------------------------\n");  

  }

  delete [] starting_points;
  delete [] step;
  delete [] xmin;

}

/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/
/**************************************************************************************************************************************/

int main(int argc, char *argv[])
{

  for(int i=1;i<argc;i++) 
      if(argv[i][0]=='-'){
        switch(argv[i][1]){
          
          case 'f':
              sscanf(argv[i+1],"%s",eos_name);
              break;
          case 'c':
              sscanf(argv[i+1],"%lf",&coupling);
              break;
          case 'e':
              sscanf(argv[i+1],"%lf",&central_density);
              break;
          case 'p':
              sscanf(argv[i+1],"%i",&print_option);
              break;
          }
      }

  load_eos(eos_name);

  double min_val[2];

  find_min(min_val);

  double nu_c,
          phi_c,
          eps_c,
          p_c,
          r_max;

  nu_c = min_val[0];
  phi_c = min_val[1];

  eps_c = central_density*pow(10,15)/Density;
  p_c = p_at_e(eps_c); 


  if (coupling >= 1.0 && coupling < 20.0)
  {
    r_max = 20.0;
  }

  else if (coupling>pow(10,3))
  {
    r_max = 500.0;
  }

  else 
  {
    r_max = coupling/2;
  }


  if (r_max < pow(3.0/(4.0*PI*central_density),1.0/3.0)){

    r_max = pow(3.0/(4.0*PI*central_density),1.0/3.0);

  }

  int neq = 7;
  int itol, itask, istate, iopt, jt, iout;
  double atol, rtol, t, tout;

  istate = 1;

  std::vector<double> r_vec = linspace(0.0,r_max,N);
  double dr = r_vec[1]-r_vec[0];

  std::vector<double> y = {0.0,nu_c,0.0,p_at_e(eps_c),0.0,phi_c,1.0};

  t = dr;
  tout = t+dr;

  LSODA lsoda;
  
  std::vector<double> yout,
                 res_m,
                 res_nu,
                 res_dphi,
                 res_pr,
                 res_dw,
                 res_fi,
                 res_w;

  res_m.push_back(y[0]);
  res_nu.push_back(y[1]);
  res_dphi.push_back(y[2]);
  res_pr.push_back(y[3]);
  res_dw.push_back(y[4]);
  res_fi.push_back(y[5]);
  res_w.push_back(y[6]);

  int idx=2;

  while(y[3]>=1e-13 && t<r_vec[N-2]){

    lsoda.lsoda_update(system_inside, 7, y, yout, &t, tout, &istate, coupling,1e-06,1e-20);
      res_m.push_back(yout[1]);
      res_nu.push_back(yout[2]);
      res_dphi.push_back(yout[3]);
      res_pr.push_back(yout[4]);
      res_dw.push_back(yout[5]);
      res_fi.push_back(yout[6]);
      res_w.push_back(yout[7]);
         
      tout += dr;

      y[0] = yout[1];
      y[1] = yout[2];
      y[2] = yout[3];
      y[3] = yout[4];
      y[4] = yout[5];
      y[5] = yout[6];
      y[6] = yout[7];

      idx += 1;
  }

  int idxlast = idx-1;

  while(t<r_vec[N-2]){

    lsoda.lsoda_update(system_outside, 7, y, yout, &t, tout, &istate, coupling,1e-06,1e-20);
      res_m.push_back(yout[1]);
      res_nu.push_back(yout[2]);
      res_dphi.push_back(yout[3]);
      res_pr.push_back(yout[4]);
      res_dw.push_back(yout[5]);
      res_fi.push_back(yout[6]);
      res_w.push_back(yout[7]);
            
      tout += dr;

      y[0] = yout[1];
      y[1] = yout[2];
      y[2] = yout[3];
      y[3] = yout[4];
      y[4] = yout[5];
      y[5] = yout[6];
      y[6] = yout[7];

      idx += 1;     

  }

  std::vector<double> dnudr;
  std::vector<double> dphidr;

  dnudr = gradient(res_nu,dr);
  dphidr = gradient(res_fi,dr);

  Mass = dnudr[N-3]*pow(r_vec[N-3],2)/(2*(1+r_vec[N-3]*dnudr[N-3]));
  
  Radius = r_vec[idxlast]*A_phi(res_fi[idxlast],coupling)*Length/pow(10,5);
  
  W_kepler = sqrt(Mass/pow(Radius,3.0));
  
  W_last = res_w[N-5];
  
  ratio = W_kepler/W_last;
  
  J = (ratio*(pow(r_vec[N-5],4.0)*res_dw[N-5]))/6.0;
  
  I = J/W_kepler;
  
 /* printf("%5.4e  %5.4e  %5.4e  %5.4e  %5.4e   \n", W_kepler, W_last, W_last*ratio, ratio, J, I);*/
        
  //Mass = dphdr[N-3]*pow(r_vec[N-3],2)/(2*(1+r_vec[N-3]*dphdr[N-3]));
  
  scalar_charge = -dphidr[N-3]*pow(r_vec[N-3],2);

  switch(print_option){

    case 0:

    printf("\n");
    printf("  %s                  \n",eos_name);
    printf("  %i         N             \n",N);
    printf("  %6.5e  e_c           (10^15 gr/cm^3)\n",central_density);
    printf("\n");
    printf("  %2.2e     alpha             \n",coupling);
    printf("\n");
    printf("  %6.5e  M             (M_sun)\n",Mass);
    printf("  %6.5e  R             (km)\n",Radius);
    printf("  %6.5e  I             \n",I);
    printf("  %6.5e  W_kepler             \n",W_kepler);
    printf("  %6.5e  dw             \n",res_w[N-5]);
    printf("  %6.5e  scalar charge          \n",scalar_charge);
    break;

    case 1:

    printf("\n");
    printf("  %s                  \n",eos_name);
    printf("  %i         N             \n",N);
    printf("  %6.5e  e_c           (10^15 gr/cm^3)\n",central_density);
    printf("\n");
    printf("  %2.2e     alpha             \n",coupling);
    printf("\n");
    printf("%5.4e ratio \n",ratio);
    printf("  %6.5e  M             (M_sun)\n",Mass);
    printf("  %6.5e  R             (km)\n",Radius);
    printf("  %6.5e  scalar charge          \n",scalar_charge);
    printf("\n");


    printf("-------------------------------------------------------------------------------\n");   
    printf("r        mu         nu        phi        epsilon        presure	w	dw	pow(r_vec[w],4)*dw\n");
    printf("-------------------------------------------------------------------------------\n");   


    for(std::size_t w =0;w<=N;w++){

        printf("%5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e %5.4e\n",
          r_vec[w]*A_phi(res_fi[w],coupling)*Length/pow(10,5), res_m[w], res_nu[w], res_fi[w], e_at_p(res_pr[w]), res_pr[w], ratio*res_w[w],res_dw[w],pow(r_vec[w],4)*res_dw[w]*ratio);
        }

    break;
  }

  return 0;
}
