# include "spectral.h"

// *****************************************************************************
double const Fluids::d2pc(double s)
{
    double aux = -1.0/theta - 2.0;
    
    if(s<0.0)s=0.0;
    else if(s>0.99)s=0.99;
    return ( pd/theta/theta*(1.0+theta)*pow((1.0-s),aux));
};
// *****************************************************************************
double const Fluids::dkrn(double s)
{
    double aux = (2.0+theta)/theta;
    
    if(s<0.0)s=0.0;
    else if(s>1.0)s=1.0;
    
    double sum = 2.0*s*(1.0-pow((1.0-s),aux));
    sum += s*s*aux*pow((1.0-s),(aux-1.0));
    return (sum);
};
// *****************************************************************************
double const Fluids::dkrw(double s)
{
    double aux = (2.0+3.0*theta)/theta;
    
    if(s<0.0)s=0.0;
    else if(s>1.0)s=1.0;
    
    return (-aux*pow((1.0-s),(aux-1.0)));
};
// *****************************************************************************
double const Fluids::dpc(double s)
{
    double aux = -1.0/theta - 1.0;
    
    if(s<0.0)s=0.0;
    else if(s>0.99)s=0.99;
    
    return ((pd/theta) * pow((1.0-s),aux));
};
// *****************************************************************************
double const Fluids::Krn(double s)
{
    double aux=(2.0+theta)/theta;
    
    if(s<0.0)s=0.0;
    else if(s>1.0)s=1.0;
    
    return (s*s*(1.0-pow((1.0-s),aux)));
};
// *****************************************************************************
double const Fluids::Krw(double s)
{ 
  double aux=(2.0+3.0*theta)/theta;
  
  if(s<0.0)s=0.0;
  else if(s>1.0)s=1.0;
  
  return (pow((1.0-s),aux));
};
// *****************************************************************************
double const Fluids::pressao_capilar(double s)
{
  double aux = -1.0/theta;
  if(s<0.0)s=0.0;
  else if(s>0.99)s=0.99;

  return (pd*pow((1.0-s),aux));
};
// ****************************************************************************
// FUNCOES GLOBAIS
// ****************************************************************************
extern  double p_in, p_out, sn_in, sn_ini, pw_ini;

double funcao_pdir(double x,double y,const int type)
{
    //extern double p_in,p_out;

    if(type==-1) return(p_in);
    else if(type==1)return(p_out);
    else return(0);
};
// ***************************************************************************
double funcao_sdir(double x,double y)
{
    //extern double sn_in;
    return (sn_in);
};
// ***************************************************************************
double sn_inicial(double x, double y,double z)
{
    //extern double sn_ini;
    return (sn_ini);
};
// ***************************************************************************
double pw_inicial(double x, double y,double z)
{
    //extern double pw_ini;
    return (pw_ini);
};
