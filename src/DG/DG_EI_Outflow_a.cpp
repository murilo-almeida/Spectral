#include "spectral.h"
#include "DG_Prob.h"
// ****************************************************************************
void DG_Prob::DG_EI_Epetra_Outflow(const EDGE border,
																	 Teuchos::RCP<Epetra_FECrsMatrix> A,
																	 Teuchos::RCP<Epetra_FEVector> RHS)
{
	MyElem e0;
	e0=el[border.elemento[0]];
	const int ntot =  e0.show_ptr_stdel(sat)->nn_val() +
	e0.show_ptr_stdel(pres)->nn_val();
	
  double  mx  [ntot*ntot];
  double  B   [ntot];
  int     indx[ntot];

	const int ntot1 = DG_EI_Outflow(border,mx,B,indx);
	
	if(ntot1 > ntot){
		cout<<"DG_EI_Epetra_Outflow: ntot1 > ntot\n";
		exit(1);
	}
  
  A->InsertGlobalValues(ntot1,&indx[0],&mx[0],Epetra_FECrsMatrix::ROW_MAJOR);
  RHS->SumIntoGlobalValues(ntot1,&indx[0],&B[0]);

};

// corrigido em 24/05/2008
// alterado em 21/09/2011
// alterado em 21/10/2011
// alterado em 13/09/2012
// alterado em 30/11/2013
