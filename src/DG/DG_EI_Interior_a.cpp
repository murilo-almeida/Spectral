#include "spectral.h"
#include "DG_Prob.h"

// ****************************************************************************
void DG_Prob::DG_EI_Epetra_Interior(const EDGE border,
                                    Teuchos::RCP<Epetra_FECrsMatrix> A,
                                    Teuchos::RCP<Epetra_FEVector> RHS)
{
	MyElem e0,e1;
	e0=el[border.elemento[0]];
	e1=el[border.elemento[1]];
	const int ntot =  e0.show_ptr_stdel(sat)->nn_val() +
		e0.show_ptr_stdel(pres)->nn_val() +
		e1.show_ptr_stdel(sat)->nn_val() +
		e1.show_ptr_stdel(pres)->nn_val();
	
  double  mx  [ntot*ntot];
  double  B   [ntot];
  int     indx[ntot];

	const int ntot1 = DG_EI_Interior(border,mx,B,indx);
	
	if(ntot1 > ntot){
		cout<<"DG_EI_Epetra_Interior: ntot1 > ntot\n";
		exit(1);
	}
 
	A->InsertGlobalValues(ntot1,&indx[0],&mx[0],Epetra_FECrsMatrix::ROW_MAJOR);
  RHS->SumIntoGlobalValues(ntot1,&indx[0],&B[0]);
 
};
// alterado em 30/11/2013
