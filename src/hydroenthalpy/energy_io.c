#include "energy_io.h"
//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadUsersParameters(int argc, char *argv[])
{
    char    filecomm[kString];
    int     j;

    // Read command line arguments
    // Standard name of command file
    // default: 
    sprintf(filecomm, "HydroEnthalpy.bcf");
    sprintf(par.outfolder, "Enthalpy_Output");
	printf("command file %s\n", filecomm);
    
    for(j = 1; j < argc; j+=2)
    {
        if(!strcmp(argv[j], "-c"))
        {
            strcpy(filecomm, argv[j+1]);
			printf("command file %s\n", filecomm);
        }
        else if(!strcmp(argv[j], "-H") || !strcmp(argv[j], "-h"))
        {
            printf("Usage: hydroenthalpy \n  -c       HydroEnthalpy.bcf \n  -o       Output folder \n");
            exit(0);
        }
        else if(!strcmp(argv[j], "-O") || !strcmp(argv[j], "-o"))
        {
        	strcpy(par.outfolder, argv[j+1]);
        }
    }
	ReadCommandFile(filecomm);

}

//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadCommandFile(char *filename)
{
    FILE    *fp;
    char    line[lString];
	
	fp = fopen(filename, "r");

    if (!fp)
    {
		printf("Command file %s is missing.\n", filename);
		exit(0);
    }
    else
    {
		// input files
		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$Input_Files:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%s", par.waterindex);
		printf("waterindex file: %s\n", par.waterindex);
        fgets(line, lString, fp);
		sscanf(line, "%s", par.wtrocupygrd);
		printf("wateroccupancygrid file: %s\n", par.wtrocupygrd);
        fgets(line, lString, fp);
		sscanf(line, "%s", par.entropy);
		printf("entropy file: %s\n", par.entropy);
        fgets(line, lString, fp);
		sscanf(line, "%s", par.enthalpy_dir);
		printf("enthalpy folder: %s\n", par.enthalpy_dir);
// 		fgets(line, lString, fp);
// 		sscanf(line, "%s", par.SAMGP0);
// 		printf("SAM file: %s\n", par.SAMGP0);
// 		fgets(line, lString, fp);
// 		sscanf(line, "%s", par.SAMGP1);
// 		fgets(line, lString, fp);
// 		sscanf(line, "%s", par.SAMGP2);

		rewind(fp);

		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$Input_Parameters:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%d", &(par.outputcluster));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.E_bulk));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.Ecutoff));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.Scutoff));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.Ereward));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.TSreward));

	}
	fclose(fp);
}
//=================================================================================================
// ReadUsersParameters
//=================================================================================================
water* ReadWaterIndexFile(void)
{
	FILE	*fi;
	char	line[lString];
	int		wtrn, numpresent, np;
	int		cluster, frame;
	water	*wtr;

	//printf("3 %p\n", wtr);
	fi = fopen(par.waterindex, "r");
	fgets(line, lString, fi);
	fgets(line, lString, fi);
	fgets(line, lString, fi);
	sscanf(line, "%d%d%d", &par.numwaters, &par.numclusters, &par.numframes);
	printf( "%d  %d  %d\n", par.numwaters, par.numclusters, par.numframes );
	
	wtr = (water*) calloc(par.numwaters, sizeof(water));
	//printf("4 %p\n", wtr);
	for(wtrn = 0; wtrn < par.numwaters; wtrn++)
	{
		fgets(line, lString, fi);
		sscanf(line, "%*d%*d%*d%d", &numpresent);
		
		wtr[wtrn].numpresent = numpresent;
		wtr[wtrn].clusterlist = (int*)calloc(numpresent, sizeof(int));
		wtr[wtrn].framelist = (int*)calloc(numpresent, sizeof(int));
		for(np = 0; np < numpresent; np++)
		{
			fgets(line, lString, fi);
			sscanf(line, "%d%d", &cluster, &frame);
/*			if (cluster == 33 || cluster == 36) 
			{
				printf("cluster %d %d\n", cluster, frame);
			}*/
			wtr[wtrn].clusterlist[np] = cluster;
			wtr[wtrn].framelist[np] = frame;
		}
	}
	fclose(fi);

/*	for(wtrn = 0; wtrn < par.numwaters; wtrn++)
	{
		printf("wtrnup %d\n", wtr[wtrn].numpresent);
		for(np = 0; np < wtr[wtrn].numpresent; np++)
		{
			printf("%d %d\n",wtr[wtrn].clusterlist[np], wtr[wtrn].framelist[np] ); 
		}
	}
*/
	return wtr;


}

//=================================================================================================
// ReadEnthalpyFiles
//=================================================================================================	
void ReadEnthalpyFiles(water *wtr, energy *engy)
{
	FILE	*fi, *fo;
	char	fname[kString], line[lString];
	int		wtrn, i, frmn, nump, clustern;
	float	*frame, *score;
	int		*clwtrn;
	float	tmps;
	float	tmpf;
	int		itmpf;

	fo = fopen("Enthalpy.txt", "w"); //YY
	frame = (float*)calloc(par.numframes, sizeof(float));
	score = (float*)calloc(par.numframes, sizeof(float));
	clwtrn = (int*)calloc(par.numclusters, sizeof(int));
	engy->enthalpy = (float*)calloc(par.numclusters, sizeof(float));
	engy->occupancy = (int*)calloc(par.numclusters, sizeof(int));
	for(i = 0; i < par.numclusters; i++)
	{
		engy->enthalpy[i] = 0.0;
		clwtrn[i] = 0;
	}
	
	for(wtrn = 0; wtrn < par.numwaters; wtrn++)
	{
		sprintf(fname, "%s/Energy_output_%d.xvg", par.enthalpy_dir, wtrn);
		//printf("%s\n", fname);
		fi = fopen(fname, "r");
		printf("fname %s\n", fname);
		
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(!strstr(line, "@") && !strstr(line, "#"))
				break;
		}
		sscanf(line, "%f%f", &frame[0], &score[0]);
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			sscanf(line, "%f%f", &tmpf, &tmps);
			itmpf = ((int)tmpf);
			//frmn = itmpf/10;
			//frmn = itmpf; //YY commented out
			frmn = itmpf-frame[0]; //YY
			//if(itmpf%10 == 0)
			//{
				frame[frmn] = tmpf;
				score[frmn] = tmps;
			//}
		}
/*		for(frmn = 1; frmn < par.numframes; frmn++)
		{
			fgets(line, lString, fi);
			sscanf(line, "%f%f", &tmpf, &tmps);
			if(tmpf%10 == 0)
			{
				frame[frmn] = tmpf;
				score[frmn] = tmps;
			}
		//printf("%d %f \n", (int)frame[frmn], score[frmn]);
		}
*/
//		printf("%d: %d\n", wtrn, wtr[wtrn].numpresent);
		// assign value to each cluster 
		for(nump = 0; nump < wtr[wtrn].numpresent; nump++)
		{
			clustern = wtr[wtrn].clusterlist[nump];
			frmn = wtr[wtrn].framelist[nump];
			
			engy->enthalpy[clustern] += score[frmn];
			fprintf(fo, "%d\t%d\t%f\n", clustern, frmn, score[frmn]); //YY
/*			if (clustern == 33 || clustern == 36) 
			{
				printf("clustern%d wtrn %d frmn %d score %f\n", clustern, wtrn, frmn, score[frmn]);
			}*/
			if (score[frmn] > -0.00001 && score[frmn] < 0.00001) 
			{
				printf("some of the water enthalpy values are incorrect! Check the g_energy run!\n");
				printf("cl%d fr%d score %f\n", clustern, frmn, score[frmn]);
				exit(0);
			}
				
			clwtrn[clustern]++;
		}
		fclose(fi);
	}
	printf("done reading energy files...\n");
/*	for(wtrn = 2000; wtrn < par.numwaters; wtrn++)
	{
		sprintf(fname, "%s/Energy_output_%d.xvg", par.enthalpy_dir, wtrn);
		printf("%s\n", fname);
		fi = fopen(fname, "r");
		//printf("fname %s\n", fname);
		
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(!strstr(line, "@") && !strstr(line, "#"))
				break;
		}
		sscanf(line, "%f%f", &frame[0], &score[0]);
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			sscanf(line, "%f%f", &tmpf, &tmps);
			itmpf = ((int)tmpf);
			frmn = itmpf/10;

			frame[frmn] = tmpf;
			score[frmn] = tmps;

		}
		// assign value to each cluster 
		for(nump = 0; nump < wtr[wtrn].numpresent; nump++)
		{
			clustern = wtr[wtrn].clusterlist[nump];
			frmn = wtr[wtrn].framelist[nump];
			
			engy->enthalpy[clustern] += score[frmn];
			//printf("cl%d fr%d score %f\n", clustern, frmn, score[frmn]);
			clwtrn[clustern]++;
		}
		fclose(fi);
	}
	*/
	for(clustern = 0; clustern < par.numclusters; clustern++)
	{
		engy->enthalpy[clustern] /= clwtrn[clustern]*4.184;		//4.184: joules->cal
		printf("enthalpy %f %d\n", engy->enthalpy[clustern], clwtrn[clustern]);
		engy->occupancy[clustern] = clwtrn[clustern];
	}
	free(frame);
	free(score);
	free(clwtrn);
	fclose(fo);
}

//=================================================================================================
// ReadEntropyFile
//=================================================================================================	
void ReadEntropyFile(energy *engy)
{
	FILE	*fi;
	char	line[lString];
	float	etop;	
	int		clustern;
	float	x, y, z, r, oc;

	engy->entropy = (float*) calloc(par.numclusters, sizeof(float));
	engy->coor = (float**) calloc(par.numclusters, sizeof(float*));
	engy->radius = (float*) calloc(par.numclusters, sizeof(float));
	engy->occupancy_probability = (float*)calloc(par.numclusters,sizeof(float));
	fi = fopen(par.entropy, "r");
	
    if (!fi)
    {
		printf("entropy file %s is missing.\n", par.entropy);
		exit(0);
    }
    else
    {
		// input files
		while(!feof(fi))
		{
	        fgets(line, lString, fi);
			if(!strstr(line, "@"))
				break;
		}
		
		for(clustern = 0; clustern < par.numclusters; clustern++)
		{
			engy->coor[clustern] = (float*) calloc(3, sizeof(float));
			fgets(line, lString, fi);
			sscanf(line, "%*d%f%*s%*f%*f%*f%*f%*f%*f%f%f%f%f%f", &etop, &x, &y, &z, &r, &oc);
			engy->entropy[clustern] = etop*1000/T0;		//mind the sign!!! read TS from entropy.egy...
			engy->coor[clustern][0] = x;
			engy->coor[clustern][1] = y;
			engy->coor[clustern][2] = z;
			engy->radius[clustern] = r;
			engy->occupancy_probability[clustern] = oc;
			printf("x y z r %f %f %f %f\n", x, y, z, r);
			//printf("%d %f\n", clustern, engy->entropy[clustern]);
		}
	}
	fclose(fi);
}
//=================================================================================================
// ReadUsersParameters
//=================================================================================================	
/*void ReadWaterOccupancyGrid(energy *engy, grid *wtrgrd)
{
	FILE	*fi;
	char	line[lString];
	int		numpoints, np, cln, nump, i, kk;
	float	x, y, z, occupancy, maxoccupancy;
	int		gx, gy, gz, nx, ny, nz;
	
	fi = fopen(par.wtrocupygrd, "r"); 
	
    if (!fi)
    {
		printf("wateroccupancy file %s is missing.\n", par.wtrocupygrd);
		exit(0);
    }
    else
    {
		// input files
	    fgets(line, lString, fi);
	    sscanf(line, "%f%f%f%d%d%d%d%f", &wtrgrd->ZeroPoint[0], &wtrgrd->ZeroPoint[1], &wtrgrd->ZeroPoint[2], &nx, &ny, &nz, &wtrgrd->numpoints, &wtrgrd->griddelta);
	    wtrgrd->numGP[0] = nx;
		wtrgrd->numGP[1] = ny;
		wtrgrd->numGP[2] = nz;
		wtrgrd->numGP1 = nz;
		wtrgrd->numGP2 = ny*nz;
		wtrgrd->numGP3 = nx*ny*nz;
		wtrgrd->clusternum = (int*)calloc(wtrgrd->numGP3, sizeof(int));
		wtrgrd->free_energy = (float*)calloc(wtrgrd->numGP3, sizeof(float));
		printf("%d %d %d %d %d %f\n", nx, ny, nz, wtrgrd->numpoints, wtrgrd->numGP3, wtrgrd->griddelta);
		for(i = 0; i < wtrgrd->numGP3; i++)
		{
			wtrgrd->clusternum[i] = -1;
			wtrgrd->free_energy[i] = 0.0;
		}
		
	    for(cln = 0; cln < par.outputcluster; cln++)
	    {
			fgets(line, lString, fi);
			sscanf(line, "%d%f", &numpoints, &maxoccupancy);

			for(np = 0; np < numpoints; np++)
			{
				//clusters[cln].hydropoint[np] = (float*)calloc(3, sizeof(float));
				fgets(line, lString, fi);
				sscanf(line, "%d%f%f%f%f", &cln, &x, &y, &z, &occupancy);
					
				// get point number
				gx = round((x - wtrgrd->ZeroPoint[0])/wtrgrd->griddelta);
				gy = round((y - wtrgrd->ZeroPoint[1])/wtrgrd->griddelta);
				gz = round((z - wtrgrd->ZeroPoint[2])/wtrgrd->griddelta);
				nump = gx*wtrgrd->numGP2 + gy*wtrgrd->numGP1 + gz;
				//printf("np %d %d\n", np, nump);
				wtrgrd->clusternum[nump] = cln;
				wtrgrd->free_energy[nump] = engy->free_energy[cln]*occupancy/maxoccupancy;				
			}
	    }
	}
	fclose(fi);
		int cc=0;
		for(i = 0; i < wtrgrd->numGP3; i++)
		{
			//printf("%7d %3d %10.3f\n", i, wtrgrd->clusternum[i], wtrgrd->free_energy[i]);
			if(wtrgrd->clusternum[i] == -1)	cc++;
		}
		printf("%d %d %d\n", wtrgrd->numGP3, wtrgrd->numpoints, cc);
}
*/
//=================================================================================================
// Read SAMGP files for grid points
//=================================================================================================
/*void ReadSAMGP(grid *grd)
{
	
	FILE	*fi[3], *fo;
	int		fn, i, j, k, *numG, nx, ny, nz, gx, gy, gz, type, gridp;
	int		searchgridp;	//search flag
	int		st;				//multiply griddelta in searching for nearest grid.
	float	xo, yo, zo, delta, xtmp[3];
	char	line[lString], property[kString], fname[kString];
	int		*numflag;
	//fo = fopen("OUT_SAM", "w");
	fi[0] = fopen(par.SAMGP0, "r");
	fi[1] = fopen(par.SAMGP1, "r");
	fi[2] = fopen(par.SAMGP2, "r");
	//fprintf(par.inputlog, "\nSAMfiles:\n");
	
	for(fn = 0; fn < 3; fn++)
	{
		if(!fi[fn])
		{
			printf("Cannot open SAMfile_%d\n", fn);
			exit(0);
		}
		sprintf(fname, "wtrgrd_%d", fn);
		fo = fopen(fname, "w"); 

		fgets(line, lString, fi[fn]);
		sscanf(line, "%f%f%f%d%d%d%f", &xo, &yo, &zo, &nx, &ny, &nz, &delta); 	
		
		j = 0;
		while(!feof(fi[fn]))
		{
			fgets(line, lString, fi[fn]);
			
			sscanf(line, "%*s%f%f%f", &xtmp[0], &xtmp[1], &xtmp[2]);	
			searchgridp = 1;
			st = 0;
			while(searchgridp)
			{
				numG = GetGridPoint(xtmp, grd, st);
				
				gridp = numG[0];
				printf("G_num out %d\n", gridp);
				for(k = 0; k < 8; k++)
				{
					if(grd->clusternum[numG[k]] != -1)
					{
						gridp = numG[k];
						searchgridp = 0;
						break;
					}
				}
				if(st > 3)		break;
				st++;
			}
			//if(st>3) printf("%d\n", st);
			fprintf(fo, "%7d %10.3f %10.3f %10.3f %10.3f %7d %7d\n", j, xtmp[0], xtmp[1], xtmp[2], grd->free_energy[gridp], grd->clusternum[gridp], gridp);
			j++;
 		}
		fclose(fo);
	}

	
	for(i = 0; i < 3; i++)
	{
		fclose(fi[i]);
	}	
}
*/
//=================================================================================================
// Read SAMGP files for grid points
//=================================================================================================
void OutputClusterInfo(energy *engy)
{
	FILE	*fo;
	int		clustern;
	fo = fopen("cluster.egy", "w");

	fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z radius occupancy\n");
	for(clustern = 0; clustern < par.numclusters; clustern++)
	{
		fprintf(fo, "%4d %10.3f %10.3f %10.3f %7.4f %7.4f %7.4f %7.2f %4d   %5.2f\n", clustern, -engy->entropy[clustern]*T0/1000, engy->enthalpy[clustern], engy->free_energy[clustern], engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], engy->radius[clustern], engy->occupancy[clustern], engy->occupancy_probability[clustern]);

		//printf("%d %f\n", clustern, engy->entropy[clustern]);
	}
	fclose(fo);

	fo = fopen("HydrationSites.pdb", "w");

	fprintf(fo, "@cluster# entropy(-TS_e) enthalpy total(E-TS) center_x center_y center_z radius occupancy\n");
	for(clustern = 0; clustern < par.numclusters; clustern++)
	{
		fprintf(fo, "ATOM    %3d  O   WAT A   1    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n", clustern+1, engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2], -engy->entropy[clustern]*T0/1000, engy->free_energy[clustern]);
	}
	fclose(fo);

	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", par.numclusters, 0);
	fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");	
	for(clustern = 0; clustern < par.numclusters; clustern++)
	{
		fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>       % 8.4f %5.2f\n",
                        clustern+1, "WAT",
                        engy->coor[clustern][0], engy->coor[clustern][1], engy->coor[clustern][2],
                        "O", engy->free_energy[clustern], -engy->entropy[clustern]*T0/1000);
	}
	fclose(fo);
}

//=================================================================================================
// GetGridPoint
//=================================================================================================
/*int* GetGridPoint(float *xtmp, grid *grd, int st)
{
	int		gx, gy, gz, closestGP, G_num[8], i;
	float	tx, ty, tz, alpha, beta, gama;

	alpha = modff( (xtmp[0] - grd->ZeroPoint[0])/grd->griddelta+Null, &tx);
	beta  = modff( (xtmp[1] - grd->ZeroPoint[1])/grd->griddelta+Null, &ty);
	gama = modff( (xtmp[2] - grd->ZeroPoint[2])/grd->griddelta+Null, &tz);
	gx = (int) (tx+Null-st);
	gy = (int) (ty+Null-st);
	gz = (int) (tz+Null-st);

	if (gx < 0 || gx >= grd->numGP[0]-1 || gy < 0 || gy >= grd->numGP[1]-1 || gz < 0 || gz >= grd->numGP[2]-1)
	{
		for(i = 0; i < 8; i++)
			G_num[i] = -1;
	}
	else
	{
		if (alpha < Null)
		{
			alpha = Null;
		}
		if (beta < Null)
		{
			beta = Null;
		}
		if (gama < Null)
		{
			gama = Null;
		}
		closestGP = (int)(2*alpha) + 2* (int)(2*beta) + 4* (int)(2*gama);
		switch(closestGP)
		{
			case 0:

				G_num[0] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 1:

				G_num[0] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[1] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 2:

				G_num[0] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 3:

				G_num[0] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 4:

				G_num[0] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 5:

				G_num[0] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 6:

				G_num[0] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[7] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));

				break;

			case 7:

				G_num[0] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[1] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + gz;
				G_num[2] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[3] = (gx+(2*st+1))*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + gz;
				G_num[4] = gx*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[5] = (gx+(2*st+1))*grd->numGP2 + gy*grd->numGP1 + (gz+(2*st+1));
				G_num[6] = gx*grd->numGP2 + (gy+(2*st+1))*grd->numGP1 + (gz+(2*st+1));
				G_num[7] = gx*grd->numGP2 + gy*grd->numGP1 + gz;


				break;

		}

	}
	printf("G_num in %d\n", G_num[0]);
	return G_num;
}
*/
//=================================================================================================
// OutputNewgridFile
//=================================================================================================

