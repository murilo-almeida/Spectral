# include "spectral.h"
#include "DG_Prob.h"
// Variaveis globais
extern  double p_in, p_out, sn_in, sn_ini, pw_ini;

void DG_Prob::DG_eco()
{
  printf("Entrando em DG_Prob::DG_eco\n");
  printf("Echo myid = %d\n", myid);
  if(myid==0) {
    FILE * fout1;
    fout1=fopen(arq_eco,"wb");
#ifdef _NEWMAT
    fprintf(fout1,"PROGRAMA USA NEWMAT \n\n");
#endif
    // ************************************************************************ 
    fprintf(fout1,"******************************************************");
    fprintf(fout1,"**************************\n");    
    fprintf(fout1,"Eco do arquivo dos dados\n");
    fprintf(fout1,"******************************************************");
    fprintf(fout1,"**************************\n");  
    fprintf(fout1,"\nNumFIELDS     = %d\n",NumFIELDS);
    fprintf(fout1,"P0           = %d\n",(Field[0].P));
		cout <<"P0           = "<<Field[0].P<<'\n';
		cout <<"P1           = "<<Field[1].P<<'\n';
    fprintf(fout1,"Q0           = %d\n",(Field[0].Q));
    fprintf(fout1,"P1           = %d\n",(Field[1].P));
    fprintf(fout1,"Q1           = %d\n",(Field[1].Q));
    
    fprintf(fout1,"mun          = %e\n",mun);
    fprintf(fout1,"muw          = %e\n",muw);
    fprintf(fout1,"rhon         = %e\n",rhon);
    fprintf(fout1,"rhow         = %e\n",rhow);
    fprintf(fout1,"theta        = %e\n",theta);
    fprintf(fout1,"pd           = %e\n",pd);
    fprintf(fout1,"gravidade[0] = %e\n",gravidade[0]);
    fprintf(fout1,"gravidade[1] = %e\n",gravidade[1]);
    fprintf(fout1,"gravidade[2] = %e\n",gravidade[2]);
    fprintf(fout1,"perm         = %e\n",perm);
    fprintf(fout1,"porosidade   = %e\n",porosidade);
    fprintf(fout1,"sigma        = %e\n",sigma);
    fprintf(fout1,"sigma1       = %e\n",sigma1);
    fprintf(fout1,"beta         = %e\n",beta);
    fprintf(fout1,"Dt0          = %e\n",Dt0);
    fprintf(fout1,"fatordt      = %e\n",fatordt);
    fprintf(fout1,"t_final      = %e\n",t_final);
    fprintf(fout1,"tprint       = %e\n",tprint);
    fprintf(fout1,"precisao     = %e\n",precisao);
    fprintf(fout1,"fonte        = %e\n",qdado);
    fprintf(fout1,"tncut        = %d\n",tncut);
    fprintf(fout1,"p_in         = %e\n",p_in);
    fprintf(fout1,"p_out        = %e\n",p_out);
    fprintf(fout1,"sn_in        = %e\n",sn_in);
    fprintf(fout1,"sn_ini       = %e\n",sn_ini);
    fprintf(fout1,"pw_ini       = %e\n",pw_ini);

    // ****************************************************************************************    
    fprintf(fout1,"\n******************************************************");
    fprintf(fout1,"**************************\n");
    fprintf(fout1,"Eco do arquivo de Geometria\n");
    fprintf(fout1,"******************************************************");
    fprintf(fout1,"**************************\n");
    fprintf(fout1,"Dimensao geometrica = %5d\n",dim);
    fprintf(fout1,"Numero de Vertices  = %5d\n",NUMNP);
    fprintf(fout1,"Numero de Elementos = %5d\n",NELEM);
    fprintf(fout1,"Numero de Arestas   = %5d\nNumero de Faces     = %5d\n",NL,NF);
    fprintf(fout1,"Numero de modos     = %5d\n",NumD); 
    fprintf(fout1,"\nEco dos Vertices: NUMNP= %4d\n",NUMNP);
    fprintf(fout1,"         i            x                   y                   z\n");
    for(int i=0;i<NUMNP;++i) {
      fprintf(fout1,"%10d   %17.10e   %17.10e   %17.10e\n",i,V[i].x,V[i].y,V[i].z);
    }
    
    fprintf(fout1,"\n******************************************************");
    fprintf(fout1,"**************************");
    fprintf(fout1,"\nEco dos dados processados\n");
    fprintf(fout1,"******************************************************");
    fprintf(fout1,"**************************\n");
    fprintf(fout1,"\nEco dos Elementos: NELEM= %4d\n",NELEM);
    fprintf(fout1,"  num  tipo");
    for(int i=0;i<4;i++)
      fprintf(fout1,"     V%1d",i);
    for(int i=0;i<4;i++)
      fprintf(fout1,"  sinal%1d",i);
    fprintf(fout1,"\n");
    for(int i=0;i<NELEM;i++) {
      int tipo = el[i].type_val();
      fprintf(fout1,"%5d    %2d",i,tipo);
      for(int j=0;j<el[i].show_ptr_stdel(0)->nv_val();j++)
	fprintf(fout1," %6d",el[i].show_Vert_map(j));
      if(tipo==2)fprintf(fout1,"        ");
      for(int j=0;j<el[i].show_ptr_stdel(0)->ne_val();j++)
	fprintf(fout1," %7d",el[i].show_sinal(j));
      fprintf(fout1,"\n");
    }
    
     // ************************************************************************
     // Eco das bordas
     fprintf(fout1,"\n******************************************************");
     fprintf(fout1,"**************************");
     fprintf(fout1,"\nDG_Prob::DG_eco\n");
     fprintf(fout1,"******************************************************");
     fprintf(fout1,"**************************\n");
     fprintf(fout1,"\nEco das bordas: NBORDER= %4d\n",NBORDER);
     fprintf(fout1,"   i tipo   Na   Nb comprimento");
     fprintf(fout1,"  ele0 lado sinal ele1 lado sinal");
     fprintf(fout1,"      n0      n1\n");
     
     for (int i=0; i < NBORDER; ++i) {
       int e=border[i].elemento[0];
       int a=border[i].num_local[0];
       int s=border[i].sinal[0];
       int t=border[i].tipo;
       double l=border[i].comprimento;
       fprintf(fout1,"%4d   %2d %4d %4d ",i,t,border[i].Na,border[i].Nb);
       // fprintf(fout1,"",t);
       fprintf(fout1,"%e ",l);
       fprintf(fout1,"%4d %4d %5d ",e,a,s);
       
       if(t!=2) {
 	fprintf(fout1,"  EXTERIOR t=%2d ",t);
 	//	if(t==-1) nin++;
 	//if(t==1) nout++;
 	//fprintf(fout1," ARESTA EXTERNA t=%2d",t);
       }
       else {
         e=border[i].elemento[1];
         a=border[i].num_local[1];
         s=border[i].sinal[1];
         fprintf(fout1,"%4d %4d %5d ",e,a,s);
       }
       fprintf(fout1,"%7.4lf %7.4lf",border[i].normal[0],border[i].normal[1]);  
       fprintf(fout1,"\n");
     }
 		/*
     for (int i=0;i<nin;i++) {
       printf("in_border %d %d\n",i,in_borders[i]);
     } 
     for(int i=0;i<nout;i++) {
       printf("out_border %d %d\n",i,out_borders[i]);
     }
     printf("nin = %d nout = %d\n",nin,nout);
     printf("\nMapeamento global dos nos: mode gbnmap sgn\n");
 		*/
    
    //Facet_eco(fout1);
    
     for(int i=0;i<NELEM;++i) {
 			fprintf(fout1,"\n");
       for(int k=0;k<NumVAR;++k) {
         int n=el[i].show_ptr_stdel(k)->nn_val();
         int b=el[i].show_ptr_stdel(k)->nb_val();
         int part=el[i].show_part_num();
         fprintf(fout1,"Elemento %d variavel %d nn = %3d nb = %3d part = %2d\n",i,k,n,b,part);
         fprintf(fout1,"mode  gbnmap  sgn\n");
         for(int j=0;j<n;++j) {
           int jj=el[i].map(k,j);
           fprintf(fout1,"%4d    %4d   %2d\n",j,jj,el[i].show_sgn(k,j));
           //	printf("AQUI i = %d j = %d jj = %d\n",i,j,jj);
         }
       }
     }
     //printf("\nMapeamento global das arestas: gbtrbmap tipo vmapM vmapP\n");
 		
     for(int i=0;i<NELEM;++i) {
       fprintf(fout1,"\nElemento %d\n",i);
       fprintf(fout1,"Mapa global das arestas:%4d\n",el[i].show_gbtrbmap());
       fprintf(fout1,"aresta border tipo vmapM vmapP\n");
       for(int k=0; k < el[i].show_ptr_stdel(0)->nborder_val(); k++) {
         int a = el[i].show_border_num(k);
         int t = border[a].tipo;
         int M = el[i].show_stgbtrbmapM(k);
         int P = el[i].show_stgbtrbmapP(k);
         fprintf(fout1,"%2d       %4d %4d %5d %5d\n",k,a,t,M,P);
       }
     }

    printf("Ao final de DG_Prob::DG_eco\nNumD = %d\n",NumD);
    fprintf(fout1,"Ao final de DG_Prob::DG_eco\nNumD = %d\n",NumD);
  
    fclose(fout1);
  }
};

