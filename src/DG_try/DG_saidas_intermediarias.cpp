//
//  DG_saidas_intermediarias.cpp
//  SDG
//
//  Created by Murilo Almeida on 11/7/13.
//  Copyright (c) 2013 Murilo Almeida. All rights reserved.
//

#include "spectral.h"
#include "DG_Prob.h"

// ****************************************
// Imprimir as taxas de producao          *
// ****************************************
void DG_Prob::DG_imprimir_taxas_de_producao(FILE * fout3,
                                            double valor,
                                            double valor0,
                                            double valor1,
                                            int iter)
{
	double auxw,auxn,pauxw,pauxn,iauxw,iauxn;
    //cout << "nin = "<< nin << "nout = "<< nout<<std::endl;
	//if(myid==0) {
		injrate_w = 0.0;
		injrate_n = 0.0;
		for(int ii = 0;ii < nin; ++ii) {
			int in_b = in_borders[ii];
			if ( border[in_b].part_num != myid) {
				int eln = border[in_b].elemento[0];
				el[eln].VolumeTracos(Dt,fluids,gbtrsn,gbtrpw);
			}
			DG_EI_flux(border[in_b],auxw,auxn);
			injrate_w += auxw;
			injrate_n += auxn;
		}
		iauxw = 0.5*(injrate_w0 + injrate_w);
		iauxn = 0.5*(injrate_n0 + injrate_n);
		Iw += iauxw*Dt;
		In += iauxn*Dt;
		injrate_w0 = injrate_w;
		injrate_n0 = injrate_n;
	
		prodrate_w = 0.0;
		prodrate_n = 0.0;
		for(int ii = 0; ii < nout; ++ii) {
			int out_b = out_borders[ii];
			if ( border[out_b].part_num != myid) {
				int eln = border[out_b].elemento[0];
				el[eln].VolumeTracos(Dt,fluids,gbtrsn,gbtrpw);
			}
			DG_EI_flux(border[out_b],auxw,auxn);
			prodrate_w -= auxw;
			prodrate_n -= auxn;
		}
		pauxw = 0.5*(prodrate_w0 + prodrate_w);
		pauxn = 0.5*(prodrate_n0 + prodrate_n);
		Qw += pauxw*Dt;
		Qn += pauxn*Dt;
		prodrate_w0 = prodrate_w;
		prodrate_n0 = prodrate_n;
	
		printf("%5d ",passo);
		printf("v0=%e v1=%e ",valor0,valor1);
		printf("vf=%e ",valor);
		printf("ite=%3d ",iter);
		printf("tnc=%3d ",tnc);tnc++;
		printf("Dt=%e ",Dt);
		printf("t=%e\n",t);
	
		fprintf(fout3,"%5d ",passo++);
		fprintf(fout3,"%e ",t);
		fprintf(fout3,"%e ",iauxw);
		fprintf(fout3,"%e ",iauxn);
		fprintf(fout3,"%e ",Iw);
		fprintf(fout3,"%e ",In);
		fprintf(fout3,"%e ",pauxw);
		fprintf(fout3,"%e ",pauxn);
		fprintf(fout3,"%e ",Qw);
		fprintf(fout3,"%e\n",Qn);
	//} //if(myid==0)
}
