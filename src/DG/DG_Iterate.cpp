# include "spectral.h"
#include "DG_Prob.h"

//# include "My_Nox_Problem_Interface.H"
//# include "My_Nox_Problem.H"

// ***************************************************************************
// *******************
// Iteracao temporal *
// *******************
// Minha implementacao do Metodo de Newton-Raphson em paralelo
// ***************************************************************************
void DG_Prob::DG_Iterate()
{
  
  Epetra_Map Map(NumD,0,*Comm);
 
  tnc=0;
  // *****************************************
  // Arquivos de saida
  FILE *fout,*fout3;
  if(myid==0) {
    fout=fopen(arq_sai,"ab");
    
    if (restart_flag==1){
      fout3=fopen(arq_flu,"ab");
    }
    else{
      passo=0;
      fout3=fopen(arq_flu,"wb");
    } 
    
  }
  printf("Iterate entrada\n"); 
   //***************************************************
  // Construcao alternativa do Mapa  
  // Programa distribui as variaveis automaticamente
  cout << "DG_Iterate: NumD = "<< NumD << endl;

  printf("Iterate ponto 1\n"); 
  
  double tcount=0.0;
  double_t norm_delta_X=0.0;
  double valor, valor0, valor1;
  int iter,cut;
  int nprint=0;
  unsigned tmatriz,tresolver;//t1,t2,
  
  if(myid == 0) {
    printf("\nIniciou iteracao temporal DG_Prob::DG_Iterate; tempo t=%e\n\n",t); 
    // cout << " tempo t= " << t << "\n\n";
  }// if(myid==0)   
  
  tmatriz=0.0;
  tresolver=0.0;
  
  // **********************************
  // Avanco temporal
  // **********************************
  // salvar Condicao inicial;
  // usave <- u0
  for (int i = 0; i < NELEM; ++i)
    el[i].Salvar_u0();  // usave <- u0

  // **********************************
  // token = 0 : Calcula o primeiro ponto para fazer a iteracao Gaussiana
  //       = 1 : Faz a Iteracao Gaussiana para avancar no tempo
  //       = 2 : Termina a Iteracao Gaussiana e avanca o tempo
  //       = 3 : Rejeita a Iteracao, reduz o Dt e recomeca a Iteracao
  //       = 4 : Termina Simulacao (atingiu tempo maximo)
	int token;
 
  if(t < t_final)
    token=0;
  else 
    token=4;
  
  cut=0;
  iter=0;
  while(token != 4) {
    
		switch (token) {
			
			case 0: // primeiro ponto
				
        Ge_MVRA0(Dt,Map,valor0,norm_delta_X);
				iter = 1;
        Ge_MVRA(Dt,Map,valor,token,norm_delta_X);
				valor1 = valor;
        if(cut==0 && myid==0) {
          cout<< "Inicio do Passo de tempo = "<< passo << endl;;
#ifdef saida_detalhada
          cout << "DG_Iterate: case 0: Dt = "<< Dt << endl;
          cout << "DG_Iterate: case 0: iteracao = "<< iter<< " valor (residual) = "<< valor << endl;
#endif
        }
        break;
			
			case 1: // Iteracao de Newton
        Ge_MVRA(Dt,Map,valor,token,norm_delta_X);
				++iter;
#ifdef saida_detalhada
				cout << "DG_Iterate: case 1: iteracao = "<< iter<< " valor (residual) = "<< valor << endl;
#endif
      break;
			
			case 2: // valor <= precisao
				// aceita o novo ponto; avanca o tempo
				//incrementa o t
				t+=Dt;
				tcount+=Dt;
				if( ( (tcount >= tprint) || (t >= t_final) ) && myid == 0 ) {
					// salva configuracao em fout
					fprintf(fout,"# t = %e\n",t);
					// imprimir dados
					for(int i = 0; i < NELEM; ++i)
						el[i].Atualizar_valores(fout);
					fprintf(fout,"\n\n");
					tcount=0.0;
					DG_Escrever_rst(++nprint);
					//printf("numero de impressao = %d\n",++nprint);
				} //if(myid==0)
      
				else { // (myid != 0)
					for(int i = 0; i < NELEM; ++i)
						el[i].Atualizar_valores();
					//atualizar sna e pwa com os novos valores calculados com u0
				}
				for (int i = 0; i < NELEM; ++i)
					el[i].Salvar_u0(); // usave <- u0
				DG_imprimir_taxas_de_producao(fout3,valor,valor0,valor1,iter);
        //std::cout << "cut = "<< cut << endl;
				cut=0;
				if( (tnc >= tncut) || (iter > 10) ) {
					tnc=0;
					Dt*=fatordt;// restaura o Dt padrao
					Dt = Dt > Dt0 ? Dt0 : Dt;
				} //if(tnc>=tncut)
      
				if( t >= t_final)
					token = 4; // termina iteracao
				else
					token = 0;// comeca nova iteracao no novo tempo.
				
			break;
			
			case 3: // Rejeita incremento
				// Recomecar iteracao com outro Dt
				tnc=0;
				++cut;
				Dt/=fatordt; // Ajuste o Dt
//#ifdef saida_detalhada
        if(myid==0) {
				//printf("Rejeitou o incremento cut=%3d ",cut);
        printf("Tentativa %3d ",cut);
				printf("Fez %3d iteracoes\n",iter);
        }
//#endif
        // restaurar os valores de u0 para os de antes da iteracao
				for (int i = 0; i < NELEM; ++i)
					el[i].Restaurar_u0(); // u0 <- usave
      
				token = 0;//reinicia o processo com Dt menor
			break; // sai da iteracao
		} //switch(token)
	} // while(token!=4)
  
  // ***************************************
  // encontrar os autovalores e autovetores
  //Eigenvectors(Dt,Map);
  
  // *****************************************************************
  // Atingiu tempo final. Escrever restart file e terminar.
  // *****************************************************************
  if(myid == 0) {
    FILE * frestart;
    frestart=fopen(arq_rst,"wb");
    fprintf(frestart,"%d %d %e %e\n",passo,tnc,Dt,t);
    fprintf(frestart,"%e %e %e %e ",injrate_w0,injrate_n0,Iw,In);
    fprintf(frestart,"%e %e %e %e\n",prodrate_w0,prodrate_n0,Qw,Qn);
    for(int i=0;i<NELEM;++i)
      el[i].escrever_restart(frestart);
    fclose(frestart);
    printf("\nFim da iteracao temporal\n");
    printf("\ntempo total para:\ncalcular matriz %u\nResolver sistema %u\n",tmatriz,tresolver);
    //cout << "Liberando memoria: Final de DG_Iterate\n";
		cout << "Final de DG_Prob::DG_Iterate\n";
  } //if(myid==0)

  // DG_conditionNumber(Map);
	
};


// ***************************************************************************
// *******************
// Iteracao temporal *
// *******************
// ***************************************************************************
