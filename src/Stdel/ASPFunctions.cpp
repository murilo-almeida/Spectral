# include "spectral.h"

// ****************************************************************************

void jacobi_poly ( int n, double alpha, double beta, double x, double cx[] )

// ****************************************************************************
//
//  Purpose:
//
//    JACOBI_POLY evaluates the Jacobi polynomials at X.
//
//  Differential equation:
//
//    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
//
//  Recursion:
//
//    P(0,ALPHA,BETA,X) = 1,
//
//    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
//
//    P(N,ALPHA,BETA,X)  =
//      (
//        (2*N+ALPHA+BETA-1)
//        * ((ALPHA**2-BETA**2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X)
//        * P(N-1,ALPHA,BETA,X)
//        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
//      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
//
//  Restrictions:
//
//    -1 < ALPHA
//    -1 < BETA
//
//  Norm:
//
//    Integral ( -1 <= X <= 1 ) ( 1 - X )**ALPHA * ( 1 + X )**BETA
//      * P(N,ALPHA,BETA,X)**2 dX
//    = 2**(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
//      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
//
//  Special values:
//
//    P(N,ALPHA,BETA)(1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.  Note
//    that polynomials 0 through N will be computed.
//
//    Input, double ALPHA, one of the parameters defining the Jacobi
//    polynomials, ALPHA must be greater than -1.
//
//    Input, double BETA, the second parameter defining the Jacobi
//    polynomials, BETA must be greater than -1.
//
//    Input, double X, the point at which the polynomials are to be evaluated.
//
//    Output, double CX[N+1], the values of the first N+1 Jacobi
//    polynomials at the point X.
//
{
  double c1;
  double c2;
  double c3;
  double c4;
  int i;
  
  if ( alpha <= -1.0 )
    {
      cout << "\n";
      cout << "JACOBI_POLY - Fatal error!\n";
      cout << "  Illegal input value of ALPHA = " << alpha << "\n";
      cout << "  But ALPHA must be greater than -1.\n";
      exit ( 1 );
    }
  
  if ( beta <= -1.0 )
    {
      cout << "\n";
      cout << "JACOBI_POLY - Fatal error!\n";
      cout << "  Illegal input value of BETA = " << beta << "\n";
      cout << "  But BETA must be greater than -1.\n";
      exit ( 1 );
    }
  
  if ( n < 0 )
    {
      return;
    }
  
  cx[0] = 1.0;
  
  if ( n == 0 )
    {
      return;
    }
  
  cx[1] = ( 1.0 + 0.5 * ( alpha + beta ) ) * x 
    + 0.5 * ( alpha - beta );
  
  for ( i = 2; i <= n; i++ )
    {
      c1 = 2.0 * ( double ) ( i ) * ( ( double ) ( i ) + alpha + beta ) 
	* ( ( double ) ( 2 * i - 2 ) + alpha + beta );
      
      c2 = ( ( double ) ( 2 * i - 1 ) + alpha + beta ) 
	* ( ( double ) ( 2 * i ) + alpha + beta ) 
	* ( ( double ) ( 2 * i - 2 ) + alpha + beta );
      
      c3 = ( ( double ) ( 2 * i - 1 ) + alpha + beta ) 
	* ( alpha + beta ) * ( alpha - beta );
      
      c4 = - ( double ) ( 2 ) * ( ( double ) ( i - 1 ) + alpha ) 
	* ( ( double ) ( i - 1 ) + beta )  
	* ( ( double ) ( 2 * i ) + alpha + beta );
      
      cx[i] = ( ( c3 + c2 * x ) * cx[i-1] + c4 * cx[i-2] ) / c1;
      
    }
  
  return;
}
// ******************************************************************
void Jacobi_P (int n, double a, double b, double x, 
	       double &y, double &dy, double &d2y)
// ******************************************************************
{
  // check parameters
  if (n < 0 || a <= -1.0 || b <= -1.0 || fabs(x) > 1.0)
    printf("%s: %s n=%d a=%lf b=%lf x=%lf\n", "Jacobi_P",   "bad argument 1\n",n,a,b,x);
  
  double c0, c1, c2, c3, c4, ab, di, ym, yp, dym, dyp, d2ym, d2yp;
  
  y = 1.0;
  dy = 0.0;
  d2y = 0.0;
  if (n == 0) return;
  
  ab = a + b;
  y = (ab + 2.0) * 0.5 * x + (a - b) * 0.5;
  dy = (ab + 2.0) * 0.5;
  d2y = 0.0;
  if (n == 1) return;
  
  yp = 1.0;
  dyp = 0.0;
  d2yp = 0.0;
  for (int i = 2; i <= n; ++i) {
    di = (double) i;
    c0 = di * 2.0 + ab;
    c1 = di * 2.0 * (di + ab) * (c0 - 2.0);
    c2 = (c0 - 1.0) * (c0 - 2.0) * c0;
    c3 = (c0 - 1.0) * (a - b) * ab;
    c4 = (di + a - 1.0) * 2.0 * c0 * (di + b - 1.0);
    ym = y;
    y = ((c2 * x + c3) * y - c4 * yp) / c1;
    yp = ym;
    dym = dy;
    dy = ((c2 * x + c3) * dy - c4 * dyp + c2 * yp) / c1;
    dyp = dym;
    d2ym = d2y;
    d2y = ((c2 * x + c3) * d2y - c4 * d2yp + c2 * 2.0 * dyp) / c1;
    d2yp = d2ym;
  }
}
// *****************************************************************
// Polymorphic form
// *****************************************************************
void Jacobi_P (int n, double a, double b, double x, 
              double &y)
// *****************************************************************
{
  // check parameters
  if (n < 0 || a <= -1.0 || b <= -1.0 || fabs(x) > 1.0)
    printf("%s: %s", "Jacobi_P",   "bad argument 2");
  
  double c0, c1, c2, c3, c4, ab, di, ym, yp;
  
  y = 1.0;
  // dy = 0.0;
  //d2y = 0.0;
  if (n == 0) return;
  
  ab = a + b;
  y = (ab + 2.0) * 0.5 * x + (a - b) * 0.5;
  //dy = (ab + 2.0) * 0.5;
  //d2y = 0.0;
  if (n == 1) return;
  
  yp = 1.0;
  //dyp = 0.0;
  //d2yp = 0.0;
  for (int i = 2; i <= n; ++i) {
    di = (double) i;
    c0 = di * 2.0 + ab;
    c1 = di * 2.0 * (di + ab) * (c0 - 2.0);
    c2 = (c0 - 1.0) * (c0 - 2.0) * c0;
    c3 = (c0 - 1.0) * (a - b) * ab;
    c4 = (di + a - 1.0) * 2.0 * c0 * (di + b - 1.0);
    ym = y;
    y = ((c2 * x + c3) * y - c4 * yp) / c1;
    yp = ym;
    //dym = dy;
    //dy = ((c2 * x + c3) * dy - c4 * dyp + c2 * yp) / c1;
    //dyp = dym;
    //d2ym = d2y;
    //d2y = ((c2 * x + c3) * d2y - c4 * d2yp + c2 * 2.0 * dyp) / c1;
    //d2yp = d2ym;
  }
}
// *****************************************************************
// Polymorphic form
// *****************************************************************
void Jacobi_P (int n, double a, double b, double x, 
              double &y,double &dy)
// *****************************************************************
{
  // check parameters
  if (n < 0 || a <= -1.0 || b <= -1.0 || fabs(x) > 1.0)
    printf("%s: %s", "Jacobi_P",   "bad argument 3");
  
  double c0, c1, c2, c3, c4, ab, di, ym, yp, dym, dyp;
  
  y = 1.0;
  dy = 0.0;
  //d2y = 0.0;
  if (n == 0) return;
  
  ab = a + b;
  y = (ab + 2.0) * 0.5 * x + (a - b) * 0.5;
  dy = (ab + 2.0) * 0.5;
  //d2y = 0.0;
  if (n == 1) return;
  
  yp = 1.0;
  dyp = 0.0;
  //d2yp = 0.0;
  for (int i = 2; i <= n; ++i) {
    di = (double) i;
    c0 = di * 2.0 + ab;
    c1 = di * 2.0 * (di + ab) * (c0 - 2.0);
    c2 = (c0 - 1.0) * (c0 - 2.0) * c0;
    c3 = (c0 - 1.0) * (a - b) * ab;
    c4 = (di + a - 1.0) * 2.0 * c0 * (di + b - 1.0);
    ym = y;
    y = ((c2 * x + c3) * y - c4 * yp) / c1;
    yp = ym;
    dym = dy;
    dy = ((c2 * x + c3) * dy - c4 * dyp + c2 * yp) / c1;
    dyp = dym;
    //d2ym = d2y;
    //d2y = ((c2 * x + c3) * d2y - c4 * d2yp + c2 * 2.0 * dyp) / c1;
    //d2yp = d2ym;
  }
}

// **********************************************************************
void Jacobi_roots(int m, double alpha, double beta, double x[])
// **********************************************************************
//
// Purpose:
//
//         Calculate the m roots of the Jacobi polynomial
//         P(m,alpha,beta)
//
// *******************************************************************
{
  double r,y,dy;
  double dm= m;
  double s, delta;
  const double pi     = 3.1415926535897932384626433832795;
  for(int k=0; k < m; k++){
    r=-cos((2.0*(double) k + 1.0)*pi/2.0/dm);
    if(k>0) r=(r+x[k-1])/2.0;
    do{
      s=0;
      for(int i=0;i<k;i++) s+=1.0/(r-x[i]);
      Jacobi_P(m, alpha, beta,r,y,dy);
      delta=-y/(dy-y*s);
      r+=delta;
    }while(delta>1.0e-10);
    x[k]=r;
  }
}
// ****************************************************************************

double my_gamma ( double x )
  
// ****************************************************************************
//
//  Purpose:
//
//    GAMMA returns the value of the Gamma function at X.
//
//  Definition:
//
//    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
//
//  Recursion:
//
//    GAMMA(X+1) = X*GAMMA(X)
//
//  Restrictions:
//
//    0 < X ( a software restriction).
//
//  Special values:
//
//    GAMMA(0.5) = sqrt(PI)
//
//    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the point at which the Gamma function is desired.
//
//    Output, double GAMMA, the Gamma function of X.
//
{
  return ( exp ( gamma_log ( x ) ) );
}
// ****************************************************************************

double gamma_log ( double x )

// ****************************************************************************
//
//  Purpose:
//
//    GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    Computation is based on an algorithm outlined in references 1 and 2.
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
//    approximation for 12 < X is from reference 3, while approximations
//    for X < 12.0 are similar to those in reference 1, but are unpublished.
//    The accuracy achieved depends on the arithmetic system, the compiler,
//    intrinsic functions, and proper selection of the machine-dependent
//    constants.
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    W. J. Cody and L. Stoltz
//    Argonne National Laboratory
//
//  Reference:
//
//    # 1)
//    W. J. Cody and K. E. Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, 1967, pages 198-203.
//
//    # 2)
//    K. E. Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    # 3)
//    Hart, Et. Al.,
//    Computer Approximations,
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable floating point number.
//
// ****************************************************************************
//
//  Explanation of machine-dependent constants
//
//  BETA   - radix for the floating-point representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
//
{
  double c[7] = {
    -1.910444077728E-03, 
    8.4171387781295E-04, 
    -5.952379913043012E-04, 
    7.93650793500350248E-04, 
    -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
  //
  //  Return immediately if the argument is out of range.
  //
  if ( x <= 0.0 || xbig < x )
    {
      return xbig ;//d_huge ( );
    }
  
  if ( x <= d_epsilon ( ) )
    {
      res = - log ( x );
    }
  else if ( x <= 1.5 )
    {
      if ( x < pnt68 )
	{
	  corr = - log ( x );
	  xm1 = x;
	}
      else
	{
	  corr = 0.0;
	  xm1 = ( x - 0.5 ) - 0.5;
	}
      
      if ( x <= 0.5 || pnt68 <= x )
	{
	  xden = 1.0;
	  xnum = 0.0;
	  
	  for ( i = 0; i < 8; i++ )
	    {
	      xnum = xnum * xm1 + p1[i];
	      xden = xden * xm1 + q1[i];
	    }
	  
	  res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
	}
      else
	{
	  xm2 = ( x - 0.5 ) - 0.5;
	  xden = 1.0;
	  xnum = 0.0;
	  for ( i = 0; i < 8; i++ )
	    {
	      xnum = xnum * xm2 + p2[i];
	      xden = xden * xm2 + q2[i];
	    }
	  
	  res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );
	  
	}
    }
  else if ( x <= 4.0 )
    {
      xm2 = x - 2.0;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
	{
	  xnum = xnum * xm2 + p2[i];
	  xden = xden * xm2 + q2[i];
	}
      
      res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
    }
  else if ( x <= 12.0 )
    {
      xm4 = x - 4.0;
      xden = - 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
	{
	  xnum = xnum * xm4 + p4[i];
	  xden = xden * xm4 + q4[i];
	}
      
      res = d4 + xm4 * ( xnum / xden );
    }
  else
    {
      res = 0.0;
      
      if ( x <= frtbig )
	{
	  
	  res = c[6];
	  xsq = x * x;
	  
	  for ( i = 0; i < 6; i++ )
	    {
	      res = res / xsq + c[i];
	    }
	  
	}
      
      res = res / x;
      corr = log ( x );
      res = res + sqrtpi - 0.5 * corr;
      res = res + x * ( corr - 1.0 );
      
    }
  
  return res;
}

// ***********************************************************;
//
//
//   Funcao C(Q,alpha,beta,xi)
//
// ***********************************************************;
double C(int Q, double alpha, double beta, double xi)
{
  double res;
  double aux=1.0;
  double y;
  
  for(int i=1 ; i<(Q+2); i++) aux*= (double) i;
  aux*=(Q+1);
  Jacobi_P (Q+1, alpha,beta, xi, y);
  res= pow(2.0, alpha+beta+1.0) * my_gamma(alpha+Q+2.0) 
    * my_gamma(beta+Q+2.0)/aux/
    my_gamma(alpha+beta+Q+3.0)/y/y;
  
  return res;
}

// ****************************************************************************
//
// Calcula os pontos e os pesos para o metodo de 
//    integracao de Gauss-Lobatto-Jacobi
//    e os pesos da derivada por colocacao
//
// ****************************************************************************
void Gauss_Lobatto_Jacobi_parameters(int Q, double alpha, double beta, 
				     double x[], double w[], 
				     double D[MAXQ][MAXQ])
{
  int Qm1 = Q-1;
  int Qm2 = Q-2;
  int i;
  double dp1[MAXQ], dp2[MAXQ];
 // double sgn=1.0;
  double y, dy, d2y, aux, aaux;
  
  // Calcula os pontos de quadratura
  if(Qm2>0) {
    Jacobi_roots(Qm2,alpha+1.0, beta+1.0, x);
    for ( i=Qm2; i>0; i--) x[i]=x[i-1];
  }
  x[0]   = -1.0;
  x[Qm1] =  1.0;
  
  // Calcula os pesos
  w[0]   =  (beta+1.0) * C(Qm2,alpha,beta,x[0]);
  
  if(Qm2>0)
    for( i=1; i<Qm1; i++)
      w[i]= C (Qm2, alpha, beta,x[i]);
  
  w[Qm1] = (alpha+1.0) * C(Qm2,alpha,beta,x[Qm1]);
  
  // Collocation Differentiation
  for(i=0; i<Q;i++){
    aux=x[i];
    aaux=1.0-(aux*aux);
    Jacobi_P(Qm2, alpha+1.0, beta+1.0, x[i], y, dy,d2y);
    dp1[i]=aaux*dy-2.0*aux*y;
    dp2[i]=aaux*d2y-4.0*aux*dy-2.0*y;
  }
  
  for(i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
 	D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}

// ****************************************************************************
//
// Calcula os pontos e os pesos para o metodo de 
//    integracao de Gauss-Radau-Jacobi
//    e os pesos da derivada por colocacao
//
// ****************************************************************************
void Gauss_Radau_Jacobi_parameters(int Q, double alpha, double beta, 
				   double x[], double w[], 
				   double D[MAXQ][MAXQ])
{
  int Qm1 = Q-1;
  int i;
  double dp1[MAXQ], dp2[MAXQ];
  double y, dy, d2y, aux, aaux;
  
  // Calcula os pontos de quadratura
  if(Qm1>0){
    Jacobi_roots(Qm1,alpha, beta+1.0, x);
    for ( i=Qm1; i>0; i--) x[i]=x[i-1];
  }
  x[0]   = -1.0;
  
  // Calcula os pesos
  w[0]   = (beta+1.0) * B(Qm1,alpha,beta,x[0]);
  for( i=1; i<Q; i++)
    w[i]= B (Qm1, alpha, beta,x[i]);
  
  //Collocation differentiation
  for(i=0; i<Q;i++){
    aux=x[i];
    aaux=1.0+aux;
    Jacobi_P(Qm1, alpha, beta+1.0, aux, y, dy,d2y);
    dp1[i]=aaux*dy+y;
    dp2[i]=aaux*d2y+2.0*dy;
  }
  
  for(i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
        D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}
// ****************************************************************************
//
// Calcula os pontos e os pesos para o metodo de 
//    integracao de Gauss-Jacobi
//    e os pesos da derivada por colocacao
//
// ****************************************************************************
void Gauss_Jacobi_parameters(int Q, double alpha, double beta, 
				   double x[], double w[], 
				   double D[MAXQ][MAXQ])
{
  int i;
  double dp1[MAXQ], dp2[MAXQ];
  double y, dy, d2y, aux, aaux;
  double H, Qfat;

  Qfat=1.0;
  for (i=2;i<=Q;i++) Qfat *= (double)i;

  // Calcula os pontos de quadratura
  Jacobi_roots(Q,alpha, beta, x);

  H = pow(2.0, alpha+beta+1.0) 
    * my_gamma(alpha+Q+1.0) 
    * my_gamma(beta+Q+1.0) 
    / Qfat 
    / my_gamma(alpha+beta+Q+1.0);
  
  // Calcula os pesos
  //Collocation differentiation
  for(i=0; i<Q;i++){
    aux=x[i];
    aaux=1.0-aux*aux;
    Jacobi_P(Q, alpha, beta, aux, y, dy,d2y);
    w[i] = H /aaux / dy / dy;
    dp1[i]=dy;
    dp2[i]=d2y;
  }
  
  for(i=0; i<Q; i++)
    for(int j=0; j<Q; j++){
      if(i==j) D[i][j] = dp2[i]/2.0/dp1[i];
      else
        D[i][j] = dp1[i]/dp1[j]/(x[i]-x[j]);
    }
}
// ***********************************************************;
//
//
//   Funcao B(Q,alpha,beta,xi)
//
// ***********************************************************;
double B(int q, double alpha, double beta, double xi)
{
  double res;
  double aux=1.0;
  double y;
  
  for(int i=1 ; i<(q+1); i++) aux*= (double) i;
  Jacobi_P (q, alpha,beta, xi, y);
  res= pow(2.0, alpha+beta) * my_gamma(alpha+(double)q+1.0) 
    * my_gamma(beta+(double)q+1.0)
    * (1.0-xi)/aux/(beta+(double)q+1.0)/my_gamma(alpha+beta+q+2.0)/y/y;
  return res;
}
// ****************************************************************************

double d_epsilon ( void )

// ****************************************************************************
//
//  Purpose:
//
//    D_EPSILON returns the round off unit for floating arithmetic.
//
//  Discussion:
//
//    D_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double D_EPSILON, the floating point round-off unit.
//
{
  double r;
  
  r = 1.0;
  
  while ( 1.0 < ( double ) ( 1.0 + r )  )
    {
      r = r / 2.0;
    }
  
  return ( 2.0 * r );
}
