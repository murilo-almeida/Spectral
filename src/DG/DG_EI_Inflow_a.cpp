#include "spectral.h"
#include "DG_Prob.h"
/*****************************************************************************/
/*     USO com Trilinos  */
/*****************************************************************************/
// ****************************************************************************
void DG_Prob::DG_EI_Epetra_Inflow(const EDGE border,
                                  Teuchos::RCP<Epetra_FECrsMatrix> A,
                                  Teuchos::RCP<Epetra_FEVector> RHS)
{
	MyElem e0;
	e0=el[border.elemento[0]];
	const int ntot =  e0.show_ptr_stdel(sat)->nn_val() +
	e0.show_ptr_stdel(pres)->nn_val();
	
  double* mx  = new double [ntot*ntot];
  double*  B  = new double [ntot];
  int*     indx = new int[ntot];
	
  // Chama DG_EI_Inflow para calcular mx e B
	const int ntot1 = DG_EI_Inflow(border,mx,B,indx);
	
	if(ntot1 > ntot){
		cout<<"DG_EI_Epetra_Inflow: ntot1 > ntot\n";
		exit(1);
	}
	
  A->SumIntoGlobalValues(ntot1,&indx[0],&mx[0],Epetra_FECrsMatrix::ROW_MAJOR);
  RHS->SumIntoGlobalValues(ntot1,&indx[0],&B[0]);
  
  delete [] mx; mx = nullptr;
  delete [] B; B = nullptr;
  delete [] indx; indx = nullptr;
};

// alterado em 21/09/2011
// alterado em 21/10/2011
// alterado em 23/10/2011
// alterado em 13/09/2012
// alterado em 30/11/2013
// *****************************************************
void DG_Prob::DG_EI_Epetra_Inflow_map(const EDGE border,
                                  Teuchos::RCP<Epetra_FECrsGraph> A,
                                  Teuchos::RCP<Epetra_FEVector> RHS)
{
    MyElem e0;
    e0=el[border.elemento[0]];
    const int ntot =  e0.show_ptr_stdel(sat)->nn_val() +
    e0.show_ptr_stdel(pres)->nn_val();
    
    double* mx  = new double [ntot*ntot];
    double*  B  = new double [ntot];
    int*     indx = new int[ntot];
    
    // Chama DG_EI_Inflow para calcular mx e B
    const int ntot1 = DG_EI_Inflow(border,mx,B,indx);
    
    if(ntot1 > ntot){
        cout<<"DG_EI_Epetra_Inflow: ntot1 > ntot\n";
        exit(1);
    }
    
    A->InsertGlobalIndices(ntot1,&indx[0],ntot1,&indx[0]);
    //RHS->SumIntoGlobalValues(ntot1,&indx[0],&B[0]);
    
    delete [] mx; mx = nullptr;
    delete [] B; B = nullptr;
    delete [] indx; indx = nullptr;
};
