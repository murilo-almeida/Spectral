// ********************************************************
//---------------------------------------------------------
#include "spectral.h"

//# include "Problem_Interface.HXX"
//# include "My_Nox_Problem.HXX"

//#include "TesteHeader.h"

//---------------------------------------------------------
int main(int argc, char* argv[])
//---------------------------------------------------------
{ // Seq: 00  main
  int rank,size;
  
#ifdef HAVE_MPI
  MPI::Init(argc, argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
#else
  rank=0;
  size=1;
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  
  if(argc!=2) {
      if(rank==0){
        printf("Numero de argumentos de chamada diferente de 2!\n");
        printf("Forneca o nome do arquivo de dados!\n");
      }
#ifdef HAVE_MPI
    MPI::Finalize();
#endif
      exit(0);
    }
  
  //cout<< "ALO"<< endl;
  if(rank==0){
    printf("Numero de processadores %d\n",size);
  }
  
  // ***********************************
  // Tipo de Problema a ser resolvido  *
  // ***********************************
  
  // unica linha a ser editada para escolher o problema
  typedef DG_Prob My_Prob;
  //      ^^^^^^^
  // **********************************************
  
  My_Prob  * p = new My_Prob (Comm); // simulador
  
  // **********************************************
 
  if (p)
	{
    // passe o  (myid=) rank  para My_Prob
    
    p->set_id(rank,size);
   
    // ****************************************************
    // Ler arquivos de entrada e montar estrutura de dados
    // ****************************************************
    // preamble deve ser especificado em cada problema
    // ****************************************************
    
    p->My_Prob::preamble(argv[1]);
    
    // ******** Inicio do Driver ************************************
    // Driver deve ser especificado em cada problema
    // Driver contem as especificacoes dos algoritmos a serem usados
    // Deve ser usado para mudar os parâmetros como o tipo de solver
    
    p->Driver();
    
    // Criar o nox solver my_solver
    // Teuchos::RCP<My_Nox_Problem<DG_Prob>>  nox_problem =
    // Teuchos::rcp(new My_Nox_Problem<DG_Prob> (*p,Comm));
    
    // ********************************************
    // Terminando a simulacao;
    // Limpar a memoria usada;
    // ********************************************
    delete p; p=nullptr;       // delete simulator
    
    cout << "\nFinal da simulacao: myid = " << rank <<"\nLimpou memoria\n";
	}
	else {
		printf("No simulator created");
	}

#ifdef HAVE_MPI
  MPI::Finalize();
#endif
  
  return 0;
}
