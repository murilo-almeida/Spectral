//
//  DG_preamble.cpp
//  SDG
//
//  Created by Murilo Pereira de Almeida on 11/1/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//
# include "spectral.h"
# include "DG_Prob.h"

// Global variables
double p_in, p_out, sn_in, sn_ini, pw_ini;

// ****************************************************************************
void DG_Prob::preamble(char * arq_entrada)
// ****************************************************************************
{
  //FILE *finput_geo,*finput_par;
  char arq_geo[256]; // arquivo de geometria
  char arq_par[256]; // arquivo de parametros

  // Seq: 01  DG_Prob::preamble

  if(myid==0) {

    printf("Numero de processos = %d\n",comm_size);
    printf("PONTO 1: Comeco. clock acumulado =%u\n",(unsigned) clock());

    // Alteracoes no preamble - parte trazida de Ler_arquivos
    cout << "\nArquivo de Entrada: "<<arq_entrada<< '\n';
    FILE *finput0; // arquivo de entrada
    finput0=fopen(arq_entrada,"rb");
    // **********************
    // Arquivo de Entrada   *
    // **********************

    // Ler o nome do arquivo de dados da geometria -- grade
    fscanf(finput0,"%s",arq_geo); strcat(arq_geo,"\0");// Linha 1
    cout << "Arquivo de Geometria: "<<arq_geo<< '\n';

    // Ler o nome do arquivo de saida
    fscanf(finput0,"%s",arq_sai); strcat(arq_sai,"\0");// Linha 2
    strcpy(arq_rst,arq_sai);
    strcpy(arq_flu,arq_sai);
    cout << "Arquivo de Saida: "<<arq_sai<< '\n';
    strcat(arq_rst,".rst\0");
    printf("Arquivo para re-inicio %s\n",arq_rst);
    strcat(arq_flu,".fluxes\0");
    cout << "Arquivo de Saida de fluxos: "<<arq_flu<< '\n';

    // Ler o nome do arquivo de eco
    fscanf(finput0,"%s",arq_eco); strcat(arq_eco,"\0");  // Linha 3
    cout << "Arquivo de Eco: "<<arq_eco<< "\n\n";

    // Ler o numero de espacos interpolantes (NumFIELDS)
    //fscanf(finput0,"%d %*s",&NumFIELDS);
    fscanf(finput0,"%*s"); // Ler e descarta a linha contendo "NumFIELDS"; NumFIELDS especificado em GeProb<N_VAR,NumFIELDS>
    for(int i =0; i<NumFIELDS;++i) {
			// Ler os dados de cada Field
      fscanf(finput0,"%d %d %d",&Field[i].ordem, &Field[i].P, &Field[i].Q);
    }

    // Ler o numero de variaveis (NumVAR); N_VAR eh especificado em GeProb<N_VAR,NumFIELDS> e pode ser maior que NumFIELDS;
    fscanf(finput0,"%*s");// Ler string e descarta

    for(int i =0; i<NumVAR;++i) {
      //Ler os dados de cada variavel
      fscanf(finput0,"%d",&FieldOfVar[i]);
    }

    fscanf(finput0,"%s",arq_par); strcat(arq_par,"\0"); // nome do arquivo de parametros
    fclose(finput0); // fechou o arquivo de entrada
  } // terminou if(myid == 0)

#ifdef HAVE_MPI
  Comm->Barrier();
  if(comm_size > 1) {
    for(int ii=0;ii<NumFIELDS;++ii) {
      MPI::COMM_WORLD.Bcast(&Field[ii].ordem,1,MPI::INT,0);
      MPI::COMM_WORLD.Bcast(&Field[ii].P,1,MPI::INT,0);
      MPI::COMM_WORLD.Bcast(&Field[ii].Q,1,MPI::INT,0);
    }

    for(int ii=0;ii<NumVAR;++ii) {
      MPI::COMM_WORLD.Bcast(&FieldOfVar[ii],1,MPI::INT,0);
    }
  }
#endif

  // *************************************************************************
  // Ler os parametros especificos do problema
  // *************************************************************************

  Ler_Arquivo_Dados_DG(arq_par); // formato novo do arquivo de entrada

  Ler_e_Processar_malha(arq_geo);
    cout << "Passou Ler_e_Processar_malha(arq_geo); "<< std::endl;
  // ****************************
  // Impor as Condicoes_contorno
  // ****************************

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
  // Iniciar os mapas nos elementos;                   *
  // Preencher os mapas de valores internos e externos *
  // nos contornos dos elementos;                      *
  // Impor as condicoes de contorno                    *
  // ***************************************************
  for ( int i = 0; i < NBORDER; ++i ) {
    int t  = border[i].tipo;
    int e0 = border[i].elemento[0];
    int a0 = border[i].num_local[0];

    //   int s0 = border[i].sinal[0];
    int bind0 = el[e0].get_trace_border_map(a0);
    border[i].gbtrbind[0]=bind0;

      // *****************************************************
      // Tentativa de 02/10/2018
      int qlocal = el[e0].show_ptr_stdel(0)->qborder_val();
      double x[qlocal], y[qlocal], z[qlocal];
      const int nvert=el[e0].show_ptr_stdel(0)->nv_val();
     // cout << "nvert em Preamble "<< nvert << std::endl;
      int map[nvert];
      for (int i=0;i<nvert;++i){
          map[i]=el[e0].show_Vert_map(i);
      }
      el[e0].show_ptr_stdel(0)->face_GQCoord(V,map,a0,qlocal,x,y,z);

      if(t == -1 || t == 1){
          border[i].pdir = new double [qlocal];
          for(int q=0;q<qlocal;++q){
              border[i].pdir[q]=funcao_pdir(x[q],y[q],t);
          }
          if(t==-1){
              nin++;
              in_borders.push_back(i);
              border[i].sdir = new double[qlocal];
              for(int q=0;q<qlocal;++q){
                  border[i].sdir[q]=funcao_sdir(x[q],y[q]);
              }
          }
          else {
              out_borders.push_back(i);
              nout++;
          }
      }

      // Fim de tentativa de 02/10/2018
    // **************************************
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
    // Alterado em 2/10/2018
  //Processa_condicoes_contorno();
    cout << "Saindo de Preamble\n";
}


// ****************************************************************************
void DG_Prob::Ler_Arquivo_Dados_DG(char *arq_par)
// ****************************************************************************
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
  if(comm_size > 1) MPI_Recebe_Dados_DG();
  Comm->Barrier();
#endif
 };

// *****************************************************************************
void DG_Prob::MPI_Recebe_Dados_DG()
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
};
