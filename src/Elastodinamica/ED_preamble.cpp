//
//  EG_preamble.cpp
//
//  Created by Murilo Pereira de Almeida on 10/29/16.
//  Copyright (c) 2016 Murilo Almeida. All rights reserved.
//
# include "spectral.h"
//#include "ED_Prob.h"

// ****************************************************************************
void ED_Prob::preamble(char * arq_entrada)
// ****************************************************************************
{
  char arq_geo[256]; // arquivo de geometria
  char arq_par[256]; // arquivo de parametros
  
  // Seq: 01  ED_Prob::ED_preamble
  
  if(myid==0) {
    
    printf("Numero de processos = %d\n",comm_size);
    printf("PONTO 1: Comeco. clock acumulado =%u\n",(unsigned) clock());
    
    // Alteracoes no preamble - parte trazida de Ler_arquivos
    cout << "\nArquivo de Entrada: "<<arq_entrada<< '\n';
    FILE *finput0;
    finput0=fopen(arq_entrada,"rb");
    // **********************
    // Arquivo de Entrada   *
    // **********************
    // Ler o nome do arquivo de dados da geometria grade
    fscanf(finput0,"%s",arq_geo);  strcat(arq_geo,"\0");  // Linha 1
    cout << "Arquivo de Geometria: "<<arq_geo<< '\n';
    
    // Ler o nome do arquivo de saida
    fscanf(finput0,"%s",arq_sai);  // Linha 2
    strcpy(arq_rst,arq_sai);
    strcpy(arq_flu,arq_sai);
    strcat(arq_sai,"\0");
    cout << "Arquivo de Saida: "<<arq_sai<< '\n';
    strcat(arq_rst,".rst\0");
    printf("Arquivo para re-inicio %s\n",arq_rst);
    strcat(arq_flu,".fluxes\0");
    cout << "Arquivo de Saida de fluxos: "<<arq_flu<< '\n';
    
    // Ler o nome do arquivo de eco
    fscanf(finput0,"%s",arq_eco);  strcat(arq_eco,"\0");  // Linha 3
    cout << "Arquivo de Eco: "<<arq_eco<< "\n\n";
    
    // Ler o numero de espacos interpolantes (N_FIELDS)
    //fscanf(finput0,"%d %*s",&N_FIELDS);
    fscanf(finput0,"%*s"); // Ler e descarta a linha contendo "NFIELDS"
    for(int i =0; i<N_FIELDS;++i) {
			// Ler os dados de cada Field
      fscanf(finput0,"%d %d %d",&Field[i].ordem, &Field[i].P, &Field[i].Q);
    }
    
    // Ler o numero de variaveis (N_var)
    fscanf(finput0,"%*s");// Ler string e descarta
    
    for(int i =0; i<N_VAR;++i) {
      //Ler os dados de cada variavel
      fscanf(finput0,"%d",&FieldOfVar[i]);
    }
    
    fscanf(finput0,"%s",arq_par);  // nome do arquivo de parametros
    strcat(arq_par,"\0");
    fclose(finput0);
  }
  // terminou myid ==0

#ifdef HAVE_MPI
  Comm->Barrier();
  if(comm_size > 1) {
    for(int ii=0;ii<N_FIELDS;++ii) {
      MPI::COMM_WORLD.Bcast(&Field[ii].ordem,1,MPI::INT,0);
      MPI::COMM_WORLD.Bcast(&Field[ii].P,1,MPI::INT,0);
      MPI::COMM_WORLD.Bcast(&Field[ii].Q,1,MPI::INT,0);
    }
    
    for(int ii=0;ii<N_VAR;++ii) {
      MPI::COMM_WORLD.Bcast(&FieldOfVar[ii],1,MPI::INT,0);
    }
  }
#endif
  
  // *************************************************************************
  // Ler os parametros especificos do problema
  // *************************************************************************
 
  Ler_Arquivo_Dados_ED(arq_par); // formato novo do arquivo de entrada
 
  // Ler o arquivo com os dados da malha no formato .meu  montar a estrutura 
  // de dados da geometria do problema
  Ler_e_Processar_malha(arq_geo);
 
 
  // *************************************************
  // Criar vetores globais para tracos de pw e sn    *
  // *************************************************
  int ntr=0;// total number of traces
  int count=0;
  for(int i=0;i<NELEM;++i) {
    ntr+=el[i].show_ptr_stdel(0)->nborder_val();
    el[i].inicia_gbtrbmap(count);
  }
  
  gbtrpw = new double [count];//trpw;
  gbtrsn = new double [count];//trsn;
  
  // ***************************************************
  // Iniciar os mapas nos elementos                    *
  // Preencher os mapas de valores internos e externos *
  // nos contornos dos elementos                       *
  // ***************************************************  
  for ( int i = 0; i < NBORDER; ++i ) {
    int t  = border[i].tipo;
    int e0 = border[i].elemento[0];
    int a0 = border[i].num_local[0];
    //   int s0 = border[i].sinal[0];
    int bind0 = el[e0].get_trace_border_map(a0);
    border[i].gbtrbind[0]=bind0;
    
    if(t==2) {// interior border; two adjacent elements 
      int e1 = border[i].elemento[1];
      int a1 = border[i].num_local[1];
      //   int s1 = border[i].sinal[1];
      int bind1 = el[e1].get_trace_border_map(a1);
      border[i].gbtrbind[1]=bind1;
      //   int s = s0 * s1;
      el[e0].set_stgbtrbmap(a0,bind0,bind1);
      el[e1].set_stgbtrbmap(a1,bind1,bind0);
    }
    else // boundary border; vmapM=vmapP 
      el[e0].set_stgbtrbmap(a0,bind0,bind0);
  }
}

// *****************************************************************************
void ED_Prob::MPI_Recebe_Dados_ED()
// *****************************************************************************
{
  MPI::COMM_WORLD.Bcast(&restart_flag,1,MPI::INT,0);
  MPI::COMM_WORLD.Bcast(&mun,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&muw,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&rhon,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&rhow,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&theta,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&pd,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(gravidade,3,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&perm,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&porosidade,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&sigma,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&sigma1,1,MPI::DOUBLE,0);//para termos com pc; penaliza descontinuidade em sn
  MPI::COMM_WORLD.Bcast(&beta,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&Dt0,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&fatordt,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&t_final,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&tprint,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&precisao,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&qdado,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&tncut,1,MPI::INT,0);
// termos de fontes
  MPI::COMM_WORLD.Bcast(&p_in,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&p_out,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&sn_in,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&sn_ini,1,MPI::DOUBLE,0);
  MPI::COMM_WORLD.Bcast(&pw_ini,1,MPI::DOUBLE,0);
} 


// ****************************************************************************
void ED_Prob::Ler_Arquivo_Dados_ED(char *arq_par)
// ****************************************************************************
//(char *str,double *coord,int *Ed,int *BC)
{
  if(myid==0)
  {
    FILE * finput0 = fopen(arq_par,"rb");
    // Dados da simulacao
    // Restart flag: 0=inicio; 1=continuar
    fscanf(finput0,"%d %*s",&restart_flag);
    fscanf(finput0,"%lf %*s",&mun);
    fscanf(finput0,"%lf %*s",&muw);
    fscanf(finput0,"%lf %*s",&rhon);
    fscanf(finput0,"%lf %*s",&rhow);
    fscanf(finput0,"%lf %*s",&theta);
    fscanf(finput0,"%lf %*s",&pd);
    fscanf(finput0,"%lf %*s",&gravidade[0]);
    fscanf(finput0,"%lf %*s",&gravidade[1]);
    fscanf(finput0,"%lf %*s",&gravidade[2]);
    fscanf(finput0,"%lf %*s",&perm);
    fscanf(finput0,"%lf %*s",&porosidade);
    fscanf(finput0,"%lf %*s",&sigma);
    fscanf(finput0,"%lf %*s",&sigma1);//para termos com pc; penaliza descontinuidade em sn
    fscanf(finput0,"%lf %*s",&beta);
    fscanf(finput0,"%lf %*s",&Dt0);
    fscanf(finput0,"%lf %*s",&fatordt);
    fscanf(finput0,"%lf %*s",&t_final);
    fscanf(finput0,"%lf %*s",&tprint);
    fscanf(finput0,"%lf %*s",&precisao);
    fscanf(finput0,"%lf %*s",&qdado);// termos de fontes
    fscanf(finput0,"%d %*s",&tncut);
    fscanf(finput0,"%lf %*s",&p_in);
    fscanf(finput0,"%lf %*s",&p_out);
    fscanf(finput0,"%lf %*s",&sn_in);
    fscanf(finput0,"%lf %*s",&sn_ini);
    fscanf(finput0,"%lf %*s",&pw_ini);
    fclose(finput0);
  }
#ifdef HAVE_MPI
  //MPI::COMM_WORLD.Barrier();
  Comm->Barrier();
  if(comm_size > 1) MPI_Recebe_Dados_ED();
#endif
 };

