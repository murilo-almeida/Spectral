#include "spectral.h"

// ****************************************************************************
//
//  Modified Principal Functions of the first type (Section 3.2.3 Karnadiakis)
//
// ****************************************************************************
double Psia(int I, int i, double z)
{
  double value;
  double y;
  
  if(i==I+1) value = 1.0; //Pontos colapsados
  else
    if(i==0) value = (1.0 - z) / 2.0;
    else
      if (i==I) value = (1.0 + z) / 2.0;
      else {
	Jacobi_P( i-1, 1.0, 1.0, z, y);
	value= (1.0-z)*(1.0+z)/4.0*y;
      }
  return value;
}

// ************************************************************************
//  Derivative of the
//  Modified Principal Functions of the first type
//
// ************************************************************************
double DPsia(int I, int i, double z)
{
  double value;
  double y,dy;
  if(i==I+1) value = 0.0;// Ponto colapsado
  else 
    if(i==0) value = -0.5;
    else
      if (i==I) value = 0.5;
      else {
	Jacobi_P( i-1, 1.0, 1.0, z, y, dy);
	value= (1.0-z)*(1.0+z)/4.0*dy-z/2.0*y;
      }
  return value;
}

// ************************************************************************
//
//  Modified Principal Functions of the second type
//
// ************************************************************************
double Psib(int I, int J, int i, int j, double z)
{
  double value, y;
  if(i==I+1) value = Psia(J,j,z); // Ponto colapsado
  else
    if(i==0 || i == I) value = Psia(J,j,z);
    else
      if(j==0) value = pow((1.0-z)/2.0,(double)i+1.0);
      else {
	Jacobi_P(j-1,2.0*(double)i+1.0, 1.0,z,y);
	value=pow( (1.0-z)/2.0 , (double)i+1.0)* (1.0+z)/2.0*y;
      }
  return value;
}

// ************************************************************************
//  Derivative of the
//  Modified Principal Functions of the second type
//
// ************************************************************************
double DPsib(int I, int J, int i, int j, double z)
{
  double value, y, dy;
  if(i==I+1) value=DPsia(J,j,z);// Ponto colapsado
  else
    if(i==0 || i == I) value = DPsia(J,j,z);
    else
      if(j==0) value = -((double)i+1.0) * pow((1.0-z)/2.0,(double)i)/2.0;
      else {
	Jacobi_P(j-1,2.0*(double)i+1.0, 1.0,z,y,dy);
	value= y * ( pow((1.0-z)/2.0,(double)i+1.0)*0.5 +  
		     (1.0+z)/2.0*((double)i+1.0)* pow((1.0-z)/2.0,(double)i) )
	  + pow((1.0-z)/2.0, (double)i+1.0) * (1.0+z)/2.0 * dy;
      }
  return value;
}

// ************************************************************************
//
//  Modified Principal Functions of the third type
//
// ************************************************************************
double Psic(int I, int J, int K, int i, int j, int k, double z)
{
  double value, y;
  if(i==I+1) value = Psia(K,k,z);
  else 
    if(i==0 || i == I) value = Psib(J,K,j,k,z);
    else
      if(j==0 || j == J) value= Psib(I,K,i,k,z);
      else 
	if (k == 0)  value = pow((1.0-z)/2.0,(double)(i+j)+1.0);
	else {
	  Jacobi_P(k-1, 2.0*(double)(i+j)+1.0, 1.0,z,y);
	  value=pow((1.0-z)/2.0, (double)(i+j)+1.0) * (1.0+z)/2.0 * y;
	}
  return value;
}
