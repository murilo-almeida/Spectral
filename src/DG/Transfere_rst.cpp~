# include "spectral.h"
// #include "PhElem.cpp"
/******************************************************************************/
/*   Transfere dados do arquivo para restart                                 */
/******************************************************************************/
void DG_Prob::DG_Transfere_rst(int * b_in,double * b_do,double * buff)
{
  int count;
  MPI::COMM_WORLD.Bcast(b_in,3,MPI::INT,0);
  passo= b_in[0];
  tnc  = b_in[1];
  count= b_in[2]; // tamanho do buff
  MPI::COMM_WORLD.Bcast(b_do,10,MPI::DOUBLE,0);
  Dt         =b_do[0];
  t          =b_do[1];
  injrate_w0 =b_do[2];
  injrate_n0 =b_do[3];
  Iw         =b_do[4];
  In         =b_do[5];
  prodrate_w0=b_do[6];
  prodrate_n0=b_do[7];
  Qw         =b_do[8];
  Qn         =b_do[9];
  MPI::COMM_WORLD.Bcast(buff,count,MPI::DOUBLE,0);
  count=0;
  for(int i=0;i<NELEM;i++)el[i].restart_element(buff,count);
}
