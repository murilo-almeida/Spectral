# include "spectral.h"
#include "DG_Prob.h"

// ****************************************************************************
// ****************************************************************************
void DG_Prob::Driver(char * str)
{
  // Seq: 02  DG_Prob::Driver

  // ****************************************************
  // Teste de numeracao de mesh com elementos continous
/*
   int count = 0;
   gbnmap_continuous(count);
   if(myid==0) DG_eco();
   FILE * fsaida;
   fsaida=fopen("saida_proj_C0","wb");
   projetar_C0(fsaida,funcao,0);
   fclose(fsaida);
   exit(0);
*/
  // fim do teste
  // ****************************************************

	//tnc=0;

  // *****************************************************
  // Condicoes iniciais
  // *****************************************************

cout << "Entrando em initial_conditions\n";
  DG_initial_conditions();
    cout << "Passou DG_initial_conditions();\n";
    NumGlobalElements = NumD;
    StandardMap = Teuchos::rcp(new Epetra_Map(NumGlobalElements, 0, *Comm));

  if(myid==0) {
    DG_eco(); // Escreve arquivo de eco dos dados
      printf("Driver apos eco: NumD = %d\n",NumD);
  }


  //int ncount[NumD];// so para testar: ha repeticoes quando ncount[i] > NumNz[i]

  // ******************************************************************
  // ********************
  // Cerne do problema  *
  // ********************
  // Especifica o algoritmo de resolucao de sistema linear
  TrilinosSolver="Amesos";
  //TrilinosSolver="AztecOO";
  cout << "\nUsando TrilinosSolver = "<< TrilinosSolver << "\n";

  //"Amesos_Mumps; // nao existe
  //"Amesos_Superludist";  // nao funciona
  //"Amesos_Superlu; // nao funciona

  //"Amesos_Umfpack" ;
  //"Amesos_Klu"; // lento
  //"Amesos_Lapack";

  AmesosSolverType="Amesos_Umfpack";
  if(TrilinosSolver=="Amesos")
    cout << "Usando AmesosSolverType = "<< AmesosSolverType << "\n\n";
  // Faz iteracao temporal
  // Escolha o algoritmo a ser usado

  DG_Iterate();

  // ********************
  // ******************************************************************

#ifdef HAVE_MPI
 // MPI::COMM_WORLD.Barrier();
  Comm->Barrier();
#endif

  if(myid==0) {
    cout << "Final de DG_driver\n";
    printf("FINAL: clock acumulado =%u\n",(unsigned) clock());
    printf("Terminou com sucesso!!!\n");
  }//if(myid==0)
};
