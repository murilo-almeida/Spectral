/*****************************************************************************/
/*****************************************************************************/
# include "spectral.h"
// ****************************************************************************

// ****************************************************************************

void somavv(int q,const double * ptx,const double * pty,double * ptr)
{
  for(int i=0;i<q;i++){
    ptr[i]=ptx[i]+pty[i];
  }
  /*
  int count=0;
  while((q--)>0){
    (*ptr++)=(*ptx++) + (*pty++);
    //    printf("soma[%d] = %4g\n",count++,*(ptr-1));
  }
  */
};

// ****************************************************************************
void soma3v(int q,const double *x1,const double *x2,const double *x3,
	    double * ptr)
{
//  int count=0;
  while((q--)>0){
    (*ptr++)=(*x1++) + (*x2++) + (*x3++);
    // printf("soma3v[%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void soma4v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,double * ptr)
{
 // int count=0;
  while((q--)>0){
    (*ptr++)=(*x1++) + (*x2++) + (*x3++) + (*x4++);
    // printf("soma4v[%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void soma5v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,const double *x5,
	    double * ptr)
{
 // int count=0;
  while((q--)>0){
    (*ptr++)=(*x1++) + (*x2++) + (*x3++) + (*x4++) + (*x5++);
    // printf("soma5v[%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void subtvv(int q,const double * ptx,const double * pty,double * ptr)
{
//  int count=0;
  while((q--)>0){
    (*ptr++)=(*ptx++) - (*pty++);
    // printf("subtvv[%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void prodvv(int q,const double * ptx,const double * pty,double * ptr)
{
//  int count=0;
  while((q--)>0){
    (*ptr++)=(*ptx++) * (*pty++);
    // printf("prodvv [%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void prod3v(int q,const double * ptx,const double * pty, const double * ptz,
	    double * ptr)
{
 // int count=0;
  while((q--)>0){
    (*ptr++)=(*ptx++) * (*pty++) * (*ptz++);
    // printf("prod3v [%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void prod4v(int q,const double *x1,const double *x2, const double *x3,
	    const double *x4,double * ptr)
{
 // int count=0;
  while((q--)>0){
    (*ptr++)=(*x1++) * (*x2++) * (*x3++) * (*x4++);
    // printf("prod4v [%d] = %4g\n",count++,*(ptr-1));
    
  }
};
// ****************************************************************************
void prod5v(int q,const double *x1,const double *x2,const double *x3,
	    const double *x4,const double *x5,
	    double * ptr)
{
 // int count=0;
  while((q--)>0){
    (*ptr++)=(*x1++) * (*x2++) * (*x3++) * (*x4++) *(*x5++);
    // printf("prod5v [%d] = %4g\n",count++,*(ptr-1));
    
  }
};

// ****************************************************************************
void prodev(int q, const double fator, const double * ptx, double *ptr)
{ 
  for(int i=0; i < q; i++){
    ptr[i]= fator * ptx[i];
  }
  
  //while((q--)>0){
  // (*ptr++)=fator * (*ptx++);
    // printf("prodev [%d] = %4g\n",count++,*(ptr-1));  
  // }
};
// ****************************************************************************
void prodgg(int q,const double *ptx0, const double *ptx1,const double *pty0,
	    const double *pty1, 
	    double *ptr)
{
  double s0[q],s1[q]; 
  prodvv(q,ptx0,pty0,s0);
  prodvv(q,ptx1,pty1,s1);
  somavv(q,s0,s1,ptr);
};
// ****************************************************************************
void prodgn(int q,const double *ptx0, const double *ptx1,const double *n, 
	    double *ptr)
{
  double s0[q],s1[q]; 
  prodev(q,n[0],ptx0,s0);
  prodev(q,n[1],ptx1,s1);
  somavv(q,s0,s1,ptr);
};
// ****************************************************************************
void prodgn(const int ndim, const int qmax, double **grad, const double *n, 
	    double *ptr)
{
  for(int qtemp=0;qtemp<qmax;qtemp++) {
    ptr[qtemp]=0.0;
    for(int ndir=0;ndir<ndim;ndir++){
      ptr[qtemp] += grad[ndir][qtemp]*n[ndir];
    }
  }
};
// ****************************************************************************
void Kgradn(const int ndim,const int qmax,const double *K,double **grad, 
	    const double *n,double *ptr)
{
  for(int qtemp=0;qtemp<qmax;qtemp++) {
    ptr[qtemp]=0.0;
  }

  for(int nd =0; nd < ndim; nd++){
    double aux = K[nd]*n[nd];
    for(int qtemp=0;qtemp<qmax;qtemp++) {
      //  printf("AQUI ndir = %d q = %d grad = %e\n",ndir,qtemp,grad[ndir][qtemp]);
      ptr[qtemp] += aux*grad[nd][qtemp];
    }
  }
};
// ****************************************************************************
void Kgradn(int q,const double *K,const double *ptx0, const double *ptx1, 
	    const double *n,double *ptr)
{
  /*
  for(int qtemp=0;qtemp<q;qtemp++) {
    ptr[qtemp]=0.0;
    printf("AQUI  q = %d grad[0] = %e grad[1] = %e\n",qtemp,ptx0[qtemp],ptx1[qtemp]);
  }
  */
  double s0[q],s1[q]; 
  double aux= K[0]*n[0];
  prodev(q,aux,ptx0,s0);
  aux= K[1]*n[1];
  prodev(q,aux,ptx1,s1);
  somavv(q,s0,s1,ptr);
};
// ****************************************************************************
void Kgradgrad(int q,const double *K,const double *ptx0, const double *ptx1, 
	       const double *pty0, 
	       const double *pty1,double *ptr)
{
  double s0[q],s1[q]; 
  prodvv(q,ptx0,pty0,s0);
  prodvv(q,ptx1,pty1,s1);
  prodev(q,K[0],s0,s0);
  prodev(q,K[1],s1,s1);
  somavv(q,s0,s1,ptr);
};
// ****************************************************************************
void integral(int q,const double *ptw, const double *ptx, double &r)
{
  r=0.0;
  while((q--)>0){ 
    r += ((*ptw++) * (*ptx++));
  }
  // printf("integral = %4g\n",r);  
};
// ****************************************************************************
void soma_vetor(int q,const double *ptx, double & r)
{
  r=0.0;
  while((q--)>0){ 
    r+=(*ptx++);
  }
  // printf("soma_vetor = %4g\n",r);  
};
