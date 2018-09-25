#ifndef _Fluids_headers
#define _Fluids_headers

class Fluids

{
 public:

  void set_theta(double x) {theta=x;};
  void set_pd(double x){pd=x;};
  void set_mun(double x){mun=x;};
  void set_muw(double x){muw=x;};
  void set_rhon(double x){rhon=x;};
  void set_rhow(double x){rhow=x;};
  
  double const show_theta(){return theta;};
  double const show_pd(){return pd;};
  double const show_mun(){return mun;};
  double const show_muw(){return muw;};
  double const show_rhon(){return rhon;};
  double const show_rhow(){return rhow;};

  double const Krw(double s);
  double const Krn(double s);
  double const dkrw(double s);
  double const dkrn(double s);
  double const pressao_capilar(double s);
  double const dpc(double s);
  double const d2pc(double s);


 private:
  
  double theta,pd,mun,muw,rhon,rhow;
};
#endif
