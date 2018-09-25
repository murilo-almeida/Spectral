//
//  DG_EI_Header.h
//  SDG
//
//  Created by Murilo Almeida on 11/30/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//

    #ifndef DG_EI_Header_h
    #define DG_EI_Header_h
    #include "DG_Prob.h"
//#define indice_interior(nnparc,iE,var,i) (nnparc[(2*iE + var)] + i)
//#define npos(ntot,Ti,Tj) ((ntot*Ti) + Tj)

// border.tipo = 2  => interior => Numero de elementos = 2
const int Num_elem = (border.tipo==2) ? 2 : 1;

double aux,aux1,aux2;
int m,iE,pos;

// elementos vizinhos ao border
MyElem e[2];
e[0]=el[border.elemento[0]];
e[1]=el[border.elemento[1]];

// Variaveis para dimensionar os arrays abaixo
const int ndim =e[0].show_ptr_stdel(0)->ndim_val();// Numero de dimensoes
const int qmax =e[0].show_ptr_stdel(0)->qborder_val();// Pontos de gauss nas arestas
const int n0   =e[0].show_ptr_stdel(0)->nn_val();// Num de modos para a variavel 0
const int nsat =e[0].show_ptr_stdel(sat)->nn_val();// Num de modos para saturacao
const int npres=e[0].show_ptr_stdel(pres)->nn_val();// Num de modos para pressao

// substituir por uma funcao que retorne so os w's
double length;
double Jac[qmax],w[qmax],Dtemp[MAXQ][MAXQ];
if(ndim==1){
    w[0]=1.0;
    length=1.0;
}
else {
    Gauss_Jacobi_parameters(qmax,0.0,0.0,Jac,w,Dtemp);
    length=border.modulo;
}
// ******************************************************
for(int q=0;q<qmax;q++)w[q]*=length/2.0;// salva em w = w*length/2.0
// *******************************************************
const double penalty =sigma/pow(length,beta);
const double penalty1=sigma1/pow(length,beta);// 07/08/08

// componentes da normal ao border
double n_e[ndim];
for(int i=0;i<ndim;i++)
    n_e[i]=border.normal[i]; // direcao i

// Valores das permeabilidades nos elementos vizinhos 0 e 1
double perm[2][3];
e[0].perm_val(perm[0]);
e[1].perm_val(perm[1]);

int a[2];
double * sn[2], * pw[2];
double lambdaw[2][qmax],lambdan[2][qmax],lambdat[2][qmax];
double d_lambdaw[2][qmax],d_lambdan[2][qmax],d_lambdat[2][qmax];
double pc[2][qmax],d_pc[2][qmax],d2_pc[2][qmax];
double mun=fluids.show_mun();
double muw=fluids.show_muw();
double K_g_sn_n[2][qmax],K_g_pw_n[2][qmax],K_g_pc_n[2][qmax];
double phi_r[qmax];
double phi_l[qmax];

int iEr, iEl;
double kn,kln,klw,klt;
double aaux[2];
// *************************************************************************
// Elemento vizinho iE
// *************************************************************************
for(iE=0;iE<Num_elem;iE++) {
    
    a[iE]=border.num_local[iE]; // aresta do elemento "iE" equivalente ao "border"
    // teste dos vetores globais
    // ************************************
    // Usa os tracos globais de sn e pw   *
    // ************************************
    int offset=border.gbtrbind[iE];
    sn[iE]=gbtrsn+offset;
    pw[iE]=gbtrpw+offset;
    
    for(m=0;m<qmax;m++) {
        // Calculo dos escalares
        double aux=sn[iE][m];
        lambdan[iE][m]=fluids.Krn(aux)/mun;
        lambdaw[iE][m]=fluids.Krw(aux)/muw;
        lambdat[iE][m]=lambdan[iE][m]+lambdaw[iE][m];
        d_lambdan[iE][m]=fluids.dkrn(aux)/mun;
        d_lambdaw[iE][m]=fluids.dkrw(aux)/muw;
        d_lambdat[iE][m]=d_lambdan[iE][m]+d_lambdaw[iE][m];
        pc[iE][m]=fluids.pressao_capilar(aux);
        d_pc[iE][m]=fluids.dpc(aux);
        d2_pc[iE][m]=fluids.d2pc(aux);
    }
    e[iE].Traco_Kgrad_sn(a[iE],Kgsn[iE]);
    e[iE].Traco_Kgrad_pw(a[iE],Kgpw[iE]);
    e[iE].Traco_Kgrad_pc(a[iE],Kgpc[iE]);
    
    for(int q=0;q<qmax;q++){
        for(int i=0;i<ndim;i++) {
            K_g_sn_n[iE][q]=0.0;
            K_g_pw_n[iE][q]=0.0;
            K_g_pc_n[iE][q]=0.0;
        }
        for(int i=0;i<ndim;i++) {
            K_g_sn_n[iE][q]+=(n_e[i]*Kgsn[iE][i][q]);
            K_g_pw_n[iE][q]+=(n_e[i]*Kgpw[iE][i][q]);
            K_g_pc_n[iE][q]+=(n_e[i]*Kgpc[iE][i][q]);
        }
    }
    
    // Iniciar K_g_phi_n
    
    // inicia as variaveis para sat (=0)
    for(int iM=0;iM<nsat;iM++) {
        e[iE].Traco_Kgrad_phi_n(a[iE],sat,iM,K_g_phi_n[iE][sat][iM]);
    }
    // inicia as variaveis para pres (=1)
    for(int iM=0;iM<npres;iM++) {
        e[iE].Traco_Kgrad_phi_n(a[iE],pres,iM,K_g_phi_n[iE][pres][iM]);
    }
    // fator dimensional que multiplica os sigmas (K/mu/L)
    aaux[iE]=0.0;
    for(int i=0;i<ndim;i++) {
        aux=perm[iE][i]*n_e[i];
        aaux[iE]+=aux*aux;
    }
}

if(Num_elem==1)
    aux=aaux[0];
else
    aux = aaux[0] < aaux[1] ? aaux[0] : aaux[1];

kn=sqrt(aux);
kln=1.0; //kn/mun;
klw=1.0; //kn/muw;
klt=1.0; //kln+klw;

    #endif
