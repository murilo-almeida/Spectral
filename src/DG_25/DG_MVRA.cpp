#include "spectral.h"
#include <AztecOO_ConditionNumber.h>
#include "DG_Prob.h"
// *******************************************************************
// *******************************************************************
void DG_Prob::Ge_MVRA0(const double Dt,/* Epetra_Map Map,*/
                       double & valor,
                       double_t & norm_delta_X
                       )
{
  Comm->Barrier();

  A->PutScalar(0.0);
  RHS->PutScalar(0.0);
  DG_MatrizVetor_Epetra(Dt,A,RHS);
	
  double valor_temp;
  RHS->Norm2( &valor_temp );valor_temp/= RHS->GlobalLength();
  //RHS->NormInf( &valor_temp );
  valor = valor_temp;
  //cout << "Ge_MVRA0: valor_temp = "<< valor_temp << endl; 
  if(std::isnan(valor_temp)) cout << "Ge_MVRA0  Float was Not a Number: valor_temp " << valor_temp << endl;
  ResolverComTrilinos(TrilinosSolver,*StandardMap,A,RHS,norm_delta_X);
}
// *******************************************************************
// *******************************************************************
void DG_Prob::Ge_MVRA(const double Dt,/*Epetra_Map Map,*/
                      double & valor,
                      int & token,
                      double_t & norm_delta_X
                      )
{
  Comm->Barrier();
  // *********************************************************************************************
  // Precisamos criar a matrix A e o vetor RHS aqui porque em DG_MatrizVetor_Epetra
  // A e RHS vao ser preeenchidos
  // e terminadas suas construcoes com RHS->GlobalAssemble(Add) e A->FillComplete()
  // *********************************************************************************************
  //Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp (new Epetra_FECrsMatrix (Copy, *AA));
  //Teuchos::RCP<Epetra_FEVector>  RHS = Teuchos::rcp(new Epetra_FEVector(*StandardMap,1));

  A->PutScalar(0.0);
  RHS->PutScalar(0.0);
  DG_MatrizVetor_Epetra(Dt,A,RHS);
     
  double valor_temp;
  RHS->NormInf( &valor_temp );// ou RHS->Norm2(&valor_temp); valor_temp/= RHS->GlobalLength();
  if(std::isnan(valor_temp)) cout << "Ge_MVRA  Float was Not a Number: valor_temp " << valor_temp << endl;
  //cout << "Ge_MVRA: valor_temp = "<< valor_temp << endl;  
  if ( token == 0 ) {
    valor = valor_temp; 
    if ( valor_temp <= precisao ) { // aceita o incremento
      token = 2;
    }
    else {
      ResolverComTrilinos(TrilinosSolver,*StandardMap,A,RHS,norm_delta_X);
      token = 1; // passa para a iteracao
    }
  }
  else {
    // if(token != 0) {
    
		if ( valor_temp <= precisao || norm_delta_X <= precisao) { // aceita o incremento
				token = 2;
				valor = valor_temp;
    }
		else if(valor_temp > valor) { // restaura
            token = 3;
        }
    else { // continua a simulacao
      valor = valor_temp;
      ResolverComTrilinos(TrilinosSolver,*StandardMap,A,RHS,norm_delta_X);
    }
  }
}

/*****************************************************************************/
// DG_MatrizVetor para usar com Trilinos
/*****************************************************************************/
void DG_Prob::DG_MatrizVetor_Epetra(const double Dt,
                                    Teuchos::RCP<Epetra_FECrsMatrix> A,
                                    Teuchos::RCP<Epetra_FEVector> RHS)
{
  // *******************************************************
  //int p = myid;
 
	// elementos da particao
	for (int k=0;k<Particao[myid].nele;++k) {
		int i=Particao[myid].ele[k];
		// printf(" p =%d nele %d i = %d\n",p,Particao[myid].nele,i);

		el[i].VolumeIntegrals(Dt,fluids,A,RHS,gbtrsn,gbtrpw);
			
		// printf("Passou \n");
		/* Usar se precisar atualizar so o termo de DT nas integrais de volume
		 da equacao de sn
		 el[i].VolumeIntegralsT(Dt,Ti,Tj,Tx,count,B);
		*/
		
		//printf("terminou VolumeIntegralsI %2d/%2d  count= %2d\n",i+1,Particao[myid].nele,count);
	}

	// Elementos ghosts
	for (int k=0;k<Particao[myid].ngho;++k) {
		int i=Particao[myid].gho[k];
		el[i].VolumeTracos(Dt,fluids,gbtrsn,gbtrpw);
	}
	
	//printf("Integrais de bordas\n");
	//DG_EI(A,RHS);
	
	EDGE border1;
  
	for (int k=0;k<Particao[myid].nbor;++k) {
		// bordas da particao
		int i=Particao[myid].bor[k];
		//printf(" EdgeIntegrals %3d/%d\n",i,Particao[myid].nbor);
		border1=border[i];
		int t=border1.tipo;
		//printf("tipo %d\n",t);
		switch(t) {
			case -1:
        DG_EI_Epetra_Inflow(border1,A,RHS);
				break;
				
			case 1:
        DG_EI_Epetra_Outflow(border1,A,RHS);
				break;
	
			case 2:
        DG_EI_Epetra_Interior(border1,A,RHS);
				break;
		}
	}

	//MPI::COMM_WORLD.Barrier();
  Comm->Barrier();
  // -----------------------------------------------------------
  // Montagem da Matriz e do Vetor deve ser feito uma unica vez
  // -----------------------------------------------------------
  A->GlobalAssemble(true,Add);
  RHS->GlobalAssemble(Add);
  A->FillComplete();
	//printf("Criou Matriz e Vetor\n");
};
// ************************************************
// ************************************************
void DG_Prob::DG_FECrsGraph_generate(Teuchos::RCP<Epetra_FECrsGraph> A)
{
    // *******************************************************
    //int p = myid;
    Teuchos::RCP<Epetra_FEVector>  RHS = Teuchos::rcp(new Epetra_FEVector(*StandardMap,1));
    // elementos da particao
    for (int k=0;k<Particao[myid].nele;++k) {
        int i=Particao[myid].ele[k];
        // printf(" p =%d nele %d i = %d\n",p,Particao[myid].nele,i);
        
        el[i].VolumeIntegrals_map(A,RHS);
        
        // printf("Passou \n");
        /* Usar se precisar atualizar so o termo de DT nas integrais de volume
         da equacao de sn
         el[i].VolumeIntegralsT(Dt,Ti,Tj,Tx,count,B);
         */
        
        //printf("terminou VolumeIntegralsI %2d/%2d  count= %2d\n",i+1,Particao[myid].nele,count);
    }
    
    // Elementos ghosts
    /* for (int k=0;k<Particao[myid].ngho;++k) {
        int i=Particao[myid].gho[k];
        el[i].VolumeTracos(Dt,fluids,gbtrsn,gbtrpw);
    } */
    
    //printf("Integrais de bordas\n");
    //DG_EI(A,RHS);
    
    EDGE border1;
    
    for (int k=0;k<Particao[myid].nbor;++k) {
        // bordas da particao
        int i=Particao[myid].bor[k];
        //printf(" EdgeIntegrals %3d/%d\n",i,Particao[myid].nbor);
        border1=border[i];
        int t=border1.tipo;
        //printf("tipo %d\n",t);
        switch(t) {
            case -1:
                DG_EI_Epetra_Inflow_map(border1,A,RHS);
                break;
                
            case 1:
                DG_EI_Epetra_Outflow_map(border1,A,RHS);
                break;
                
            case 2:
                DG_EI_Epetra_Interior_map(border1,A,RHS);
                break;
        }
    }
    
    //MPI::COMM_WORLD.Barrier();
    Comm->Barrier();
    // -----------------------------------------------------------
    // Montagem da Matriz e do Vetor deve ser feito uma unica vez
    // -----------------------------------------------------------
    A->GlobalAssemble(true);
    //RHS->GlobalAssemble(Add);
    //A->FillComplete();
    //printf("Criou Matriz e Vetor\n");
};

// ************************************************
// ************************************************
void DG_Prob::DG_conditionNumber(Epetra_Map Map)
{
  char  arquivo_nome [80];
  MPI::COMM_WORLD.Barrier();
  Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,Map,0));//&NNz[0]);
  Teuchos::RCP<Epetra_FEVector> RHS = Teuchos::rcp(new Epetra_FEVector(Map,1));
  AztecOOConditionNumber COND;
  printf("Entre com o nome do arquivo contendo valores de sigma e sigma1\n");
  scanf("%s",arquivo_nome);
  strcat(arquivo_nome,"\0");
  FILE * fin0 = fopen(arquivo_nome,"rb");
  printf("Entre com o nome do arquivo saida\n");
  scanf("%s",arquivo_nome);
  strcat(arquivo_nome,"\0");
  FILE * fou0 = fopen(arquivo_nome,"wb");
  int nvalues;
  fscanf(fin0,"%d",&nvalues);
  int ind[nvalues];
  double s[nvalues][2];
  for(int i = 0; i<nvalues; i++) {
    fscanf(fin0,"%d %lf %lf",&ind[i],&s[i][0],&s[i][1]);
  }
  fclose(fin0);
  Dt=2.0;
  for(int k=0;k<15;k++){
    Dt/=2.0;
    for(int j=0;j<nvalues;j++){
      sigma=s[j][0];
      sigma1=s[j][1];
      DG_MatrizVetor_Epetra(Dt,A,RHS);
      COND.initialize(*(A.ptr()));
      COND.computeConditionNumber(1000,1.0E-8);
      double condnum = COND.getConditionNumber();
      fprintf(fou0,"%2d %14.5e %14.5e %14.5e %14.5e\n",ind[j],sigma,sigma1,Dt,condnum);
      printf ("%2d (%14.5e %14.5e) %14.5e %14.5e\n",ind[j],sigma,sigma1,Dt,condnum);
    }
  }
  fclose(fou0);
};
