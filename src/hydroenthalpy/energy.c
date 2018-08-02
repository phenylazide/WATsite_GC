#include "energy.h"
//=================================================================================================
// CalculateEnthalpy
//=================================================================================================
void CalculateFreeEnergy(energy *engy)
{	

	water		*wtr;
	int			clustern;
	float		tmp, enthalpy, TdeltaS;
	float		deltaE, deltaS;
	wtr = ReadWaterIndexFile();
	
	ReadEnthalpyFiles(wtr, engy);
	
	ReadEntropyFile(engy);

	engy->free_energy = (float*)calloc(par.numclusters, sizeof(float));
	
	for(clustern = 0; clustern < par.numclusters; clustern++)
	{
		//engy->free_energy[clustern] = (engy->enthalpy[clustern] - BulkSolE) - T0*0.001*(engy->entropy[clustern] - BulkSolS);

		//use cutoff value from Freisner's paper
/*		tmp = engy->enthalpy[clustern] - par.Ecutoff;
		enthalpy = 0.0;
		if(tmp > 0)
		{
			enthalpy = par.Ereward;
		}

		tmp = engy->entropy[clustern] - par.Scutoff;

		TS = 0.0;
		if(tmp > 0)
		{
			TS = par.TSreward;
		}*/

		enthalpy = engy->enthalpy[clustern];	//the cutoff is from Freisner's paper
		TdeltaS = T0*engy->entropy[clustern]/1000.0;
		
		deltaE = enthalpy - par.E_bulk;
		//deltaS = -TS;					//S_bulk - S_prot
		engy->free_energy[clustern] = deltaE - TdeltaS;	//
		printf("%d %f\n", clustern, engy->free_energy[clustern]);
	}
	
}

//=================================================================================================
// MapontoSAMfiles
//=================================================================================================
// void MapontoSAMfiles(energy *engy)
// {
// 	grid		wtrgrd;
// 	
// 	//clusters = (hcluster*) calloc(par.numclusters, sizeof(hcluster));
// 	
// 	ReadWaterOccupancyGrid(engy, &wtrgrd);		//Read wateroccupancy grid file and map free energy on to grid
// 	printf("kkk1\n");
// 	ReadSAMGP(&wtrgrd);
// 	printf("kkk2\n");
// //	OutputNewgridFile(clusters, &grd);
// }

