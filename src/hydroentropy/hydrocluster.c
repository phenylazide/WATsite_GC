/***************************************************************************
 *	 This script 
 *	 1) reads ligand pdb file and defines the binding pocket			   *
 *	 (grid) extending from the minimum and maximum point of ligand;		   *
 *	 2) calculates water occupancy by going through the MD snapshot and	   *
 *	 mapping occupancy onto grid using gaussian distribution.			   *
 *	 3) clusters grip points by modified k-mean cluster algorithm		   *
 ***************************************************************************/
#define _FILE_OFFSET_BITS 64
#include "hydrocluster.h"
//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadLigand(char *infile, ligand *lig)
{
	FILE	*fi;
	int		nl = 0, nh = 0, i, j;	// nh: number of heavy atoms
	char	line[lString], atmn[8], coor[kString];
	float	x, y, z;
	
	fi = fopen(infile, "r");
	if(!fi)
	{
		printf("ligand file %s is missing.\n", infile);
		exit(0);
	}

	else
	{
		rewind(fi);
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(strstr(line, "ATOM") || strstr(line, "HETATM"))
			{
				nl++;
				GetChars( line, 75, 81, atmn ); 
				if(!strstr(atmn, "H"))
					nh++;
			}

		}	
		lig->atm = (atom *) calloc(nh, sizeof(atom));
		printf("number of heavy atoms in ligand: %d\n", nh);
		lig->numheavyatms = nh;

		rewind(fi);
		j = 0;
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(strstr(line, "ATOM") || strstr(line, "HETATM"))
			{
				GetChars( line, 75, 81, atmn ); 
				if(!strstr(atmn, "H"))
				{
					GetChars( line, 30, 37, coor ); 
					x = atof(coor);
					GetChars( line, 38, 45, coor ); 
					y = atof(coor);
					GetChars( line, 46, 54, coor ); 
					z = atof(coor);
					strcpy(lig->atm[j].name, atmn);
					lig->atm[j].x[0] = x;
					lig->atm[j].x[1] = y;
					lig->atm[j].x[2] = z;
					lig->atm[j].Num = j;
					j++;
					printf("atmn %s x y z: %f %f %f\n",atmn, x, y, z);	
				}
			}
		}
	}
	printf("Finish reading ligand\n");
	fclose(fi);
}

//=================================================================================================
// DefineBindingPocket
//=================================================================================================
void DefineBindingPocket(ligand *lig, grid *grd)
{
	int		i, j, k, l, m, num;
	float	longest, y[3];
	FILE	*fo;


	for(j = 0; j < 3; j++)
	{
		grd->xmin[j] = 99999.9;
		grd->xmax[j] = -99999.9;
	}

	// determine size of binding pocket
	for(k = 0; k < lig->numheavyatms; k++)
	{
		for(j = 0; j < 3; j++)
		{
			if(lig->atm[k].x[j] < grd->xmin[j])
			{
				grd->xmin[j] = lig->atm[k].x[j];	
			}
			if(lig->atm[k].x[j] > grd->xmax[j])
			{
				grd->xmax[j] = lig->atm[k].x[j];
			}
		}
	}
			
	// FarestDist is user-defined parameter (e.g. 10A)
	// griddelta is grid spacing (e.g. 0.25A)
	printf("min %f %f %f max %f %f %f\n", grd->xmin[0], grd->xmin[1], grd->xmin[2], grd->xmax[0], grd->xmax[1], grd->xmax[2]);
	for(j = 0; j < 3; j++)
	{
		grd->xmin[j] -= par.FarestDist + 2;
		grd->xmax[j] += par.FarestDist + 2;
		longest = grd->xmax[j] - grd->xmin[j];
		grd->GridSize[j] = par.griddelta * floor(longest/par.griddelta);
	}
	grd->griddelta = par.griddelta;
	printf("Farest %f longest %f\n", par.FarestDist, longest);
	// Determine zeropoint of Grid
	for(j = 0; j < 3; j++)
	{
		 grd->ZeroPoint[j] = grd->xmin[j];
	}

	// Number of grid points
	for(j = 0; j < 3; j++)
	{
		grd->NumGP_r[j] = (int) (grd->GridSize[j]/par.griddelta) + 1;
	}
	printf("%f %f %f %d %d %d\n", grd->GridSize[0]/par.griddelta, grd->GridSize[1]/par.griddelta, grd->GridSize[2]/par.griddelta, grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2]);
	grd->numGP1 = grd->NumGP_r[2];
	grd->numGP2 = grd->NumGP_r[1]*grd->NumGP_r[2];
	grd->numGP3 = grd->NumGP_r[0]*grd->NumGP_r[1]*grd->NumGP_r[2];

	// Coordinates of grid points
	grd->GP_r = (float **) calloc(grd->numGP3,sizeof(float *));
	for(k = 0; k < grd->NumGP_r[0]; k++) //YY NumGP_r[0] = the number of grid points on x axis
	{
		y[0] = k*grd->griddelta + grd->ZeroPoint[0];
		for(l = 0; l < grd->NumGP_r[1]; l++)
		{
			y[1] = l*grd->griddelta + grd->ZeroPoint[1];
			for(m = 0; m < grd->NumGP_r[2]; m++)
			{
				y[2] = m*grd->griddelta + grd->ZeroPoint[2];
				num = k*grd->numGP2 + l*grd->numGP1 + m; //YY num = index of the grid point

				grd->GP_r[num] = (float *) calloc(3,sizeof(float));

				grd->GP_r[num][0] = y[0];
				grd->GP_r[num][1] = y[1];
				grd->GP_r[num][2] = y[2]; //YY store the coordinates of grid points: grd->GP_r[num][0-2] x,y,x coordinates of the num_th grid point
			}
		}
	}
	
	printf("ZeroPoint %f %f %f %f %f %f GridSize %f %f %f\n", grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2], grd->xmax[0], grd->xmax[1], grd->xmax[2], grd->GridSize[0] , grd->GridSize[1], grd->GridSize[2]);

	//YY_8-29-14: GridSize (grd->GridSize[j]) = distance; grd->GridSize[j]/par.griddelta = number of grid point;  grd->numGP3 = 3D number of grid points. 

	printf("numGP %d %d %d %d\n", grd->numGP3, grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2]);
	printf("Binding Pocket Defined!\n");	

}


//=================================================================================================
// WaterOccupancy // read in NMRsnapshots, calculate wateroccupancy on grid: nmrpdb in vmd format
// modified for the water distribution: radius = 1A = 3*sigma (more steep distribution than 
// WaterOccupancy2)
//=================================================================================================
void WaterOccupancy3(char *infile, ligand *lig, grid *grd, fgrid *fgrd, water *wtr)
{
	FILE	*fi, *fo;
	char	line[lString], atmn[8], coor[kString];
	int		numfr = 0, numatm = 0, numWAT = 0, cc, dd, countAtm, fstWATatm; //YY 4.21.2015
	float	x, y, z, alpha, beta, gamma, tx, ty, tz, delta, dis, dis2, dx, dy, dz, xx, yy, zz, dd2;
	int		gx, gy, gz, num, gnum, del, atn, resn;
	int		i, j, k, k2, l2, m2;
	int		watinside = 0;
	int		G_num, closestGP;
	float	size_water = par.size_water;
	float	size_water2 = pow(size_water, 2);
	float	sigma = size_water/3.0;
	float	sigma2 = pow(sigma, 2);
	float	denom;
	int   grid_counter, new_counter;

	delta = grd->griddelta;		// grid spacing
	dis = size_water;			// cutoff radius around center of water molecule (3*sigma will account for 99.7% of the set)
	del = (int) (dis/delta) + 1;		// number of gridpoints around center of water molecule up to cutoff
	fi = fopen(infile, "r");
	printf("start reading protein snapshots...\n");
	if(!fi)
	{
		printf("protein file %s is missing.\n", infile);
		exit(0);
	}
	else
	{
		grd->WaterOccupancy = (float*)calloc(grd->numGP3, sizeof(float));
		for(i = 0; i < grd->numGP3; i++)
		{
			grd->WaterOccupancy[i] = 0;
		}
	
		// get number of atoms and waters inside of binding pocket
		rewind(fi);
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(strstr(line, "ATOM"))
			{
				numatm++;
			}
			if(strstr(line, "O   WAT") || strstr(line, "OW") )
			{
				numWAT++;
				GetChars( line, 30, 37, coor ); 
				x = atof(coor);
				GetChars( line, 38, 45, coor ); 
				y = atof(coor);
				GetChars( line, 46, 54, coor ); 
				z = atof(coor);

				if((grd->xmin[0] < x && x < grd->xmax[0]) && (grd->xmin[1] < y && y < grd->xmax[1]) && (grd->xmin[2] < z && z < grd->xmax[2]))
					watinside++;
			}
			if(strstr(line, "END"))
				break;
		}
		//wtr->totalresn = resn;
		wtr->numatms = numatm;
		wtr->numWATs = numWAT;
		printf("number of atoms:%d waters: %d\n", numatm, numWAT);
		printf("number of waters inside pocket in 1st frame: %d\n", watinside);
		// get number of frames
/*		rewind(fi);
		while(!feof(fi))
		{
			fgets(line, lString, fi);
			if(strstr(line, "MODEL"))
				numfr++;
		}*/
		wtr->numfr = par.numfr;
		numfr = par.numfr;
		rewind(fi);
		printf("number of frames %d\n",wtr->numfr);
		wtr->O = (atom**)calloc(numfr, sizeof(atom*));
		wtr->H1 = (atom**)calloc(numfr, sizeof(atom*));
		wtr->H2 = (atom**)calloc(numfr, sizeof(atom*));
		wtr->H3 = (atom**)calloc(numfr, sizeof(atom*));
		wtr->atn = (int**)calloc(numfr, sizeof(int*));
		wtr->resn = (int**)calloc(numfr, sizeof(int*));
		wtr->WATinside = (int*)calloc(numfr, sizeof(int));
		
		for(i = 0; i < numfr; i++)
		{
			wtr->O[i] = (atom*)calloc(watinside*2, sizeof(atom));
			wtr->H1[i] = (atom*)calloc(watinside*2, sizeof(atom));
			wtr->H2[i] = (atom*)calloc(watinside*2, sizeof(atom));
			wtr->H3[i] = (atom*)calloc(watinside*2, sizeof(atom));
			wtr->atn[i] = (int*)calloc(watinside*2, sizeof(int));
			wtr->resn[i] = (int*)calloc(watinside*2, sizeof(int));
		}
		
		// go over every frame of MD simulation, assign water occupancy to each grid point and store water coordinates
		for(i = 0; i < numfr; i++)
		{
			dd = 0;
			cc = 0;
			fstWATatm = 0; //YY 4.21.2015
			countAtm = 0; //YY 4.21.2015

			while(!feof(fi))
			{
				fgets(line, lString, fi);
				
				if(strstr(line, "O   WAT") || strstr(line, "OW"))
				{
					sscanf(line, "%*s%d%s",&atn, atmn); 
					if (countAtm == 0) //YY 4.21.2015
						fstWATatm = atn; //YY 4.21.2015
						//printf ("First Water Atom number is %d\n", fstWATatm); //YY 4.21.2015
					resn = 0;
					GetChars( line, 30, 37, coor ); 
					x = atof(coor);
					GetChars( line, 38, 45, coor ); 
					y = atof(coor);
					GetChars( line, 46, 54, coor ); 
					z = atof(coor);
					
							
					// calculate occupancy for waters within defined binding pocket
					alpha = modff( (x - grd->ZeroPoint[0])/delta+Null, &tx);
					beta  = modff( (y - grd->ZeroPoint[1])/delta+Null, &ty);
					gamma = modff( (z - grd->ZeroPoint[2])/delta+Null, &tz);
					//printf("%f %f %f\n", x, y, z);
					gx = (int) (tx+Null);
					gy = (int) (ty+Null);
					gz = (int) (tz+Null);


					//YY grd->NumGP_r[0] = number of grid points on x axis

					if(gx < 0 || gx >= grd->NumGP_r[0]-1 || gy < 0 || gy >= grd->NumGP_r[1]-1 || gz < 0 || gz >= grd->NumGP_r[2]-1) //YY water outside binding pocket
					{
						cc++;
						// skip the next two hydrogen
						fgets(line, lString, fi);
						fgets(line, lString, fi);
						if(par.TIP4P) fgets(line, lString, fi);//TIP4P water
					}


					else  //YY water inside binding pocket
					{
						for(k2 = gx - del; k2 <= gx + del; k2++)
						{
							for(l2 = gy - del; l2 <= gy + del; l2++)
							{
								for(m2 = gz - del; m2 <= gz + del; m2++)
								{
									if(k2 >= 0 && k2 < grd->NumGP_r[0] && l2 >= 0 && l2 < grd->NumGP_r[1] && m2 >= 0 && m2 < grd->NumGP_r[2])
									{
										num = k2*grd->numGP2 + l2*grd->numGP1 + m2;
										xx = grd->GP_r[num][0] - x;
										yy = grd->GP_r[num][1] - y;
										zz = grd->GP_r[num][2] - z;
										dd2 = xx*xx + yy*yy + zz*zz;
										if(dd2 <= size_water2)
										{
											grd->WaterOccupancy[num] += exp(-dd2/(2*sigma2));
										}
									}
								}
							}
						}
						// assign nearest grid point to the water
						if (alpha < Null)
						{
							alpha = Null;
						}
						if (beta < Null)
						{
							beta = Null;
						}
						if (gamma < Null)
						{
							gamma = Null;
						}
						closestGP = (int)(2*alpha) + 2* (int)(2*beta) + 4* (int)(2*gamma);
						switch(closestGP)
						{
							case 0:
								G_num = gx*grd->numGP2 + gy*grd->numGP1 + gz;
								break;				
							case 1:
								G_num = (gx+1)*grd->numGP2 + gy*grd->numGP1 + gz;
								break;			
							case 2:
								G_num = gx*grd->numGP2 + (gy+1)*grd->numGP1 + gz;
								break;				
							case 3:
								G_num = (gx+1)*grd->numGP2 + (gy+1)*grd->numGP1 + gz;
								break;
							case 4:
								G_num = gx*grd->numGP2 + gy*grd->numGP1 + (gz+1);
								break;			
							case 5:
								G_num = (gx+1)*grd->numGP2 + gy*grd->numGP1 + (gz+1);
								break;			
							case 6:
								G_num = gx*grd->numGP2 + (gy+1)*grd->numGP1 + (gz+1);
								break;
							case 7:
								G_num = (gx+1)*grd->numGP2 + (gy+1)*grd->numGP1 + (gz+1);
								break;
						}
						wtr->O[i][dd].nearestgrdp = G_num;
						// store Oxygen coordinates
						wtr->O[i][dd].x[0] = x;
						wtr->O[i][dd].x[1] = y;
						wtr->O[i][dd].x[2] = z;
						// store Hydrogen coordinates//TIP4 water

						fgets(line, lString, fi);
						//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
						GetChars( line, 30, 37, coor ); 
						x = atof(coor);
						GetChars( line, 38, 45, coor ); 
						y = atof(coor);
						GetChars( line, 46, 54, coor ); 
						z = atof(coor);

						wtr->H1[i][dd].x[0] = x;
						wtr->H1[i][dd].x[1] = y;
						wtr->H1[i][dd].x[2] = z;

						fgets(line, lString, fi);
						//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
						GetChars( line, 30, 37, coor ); 
						x = atof(coor);
						GetChars( line, 38, 45, coor ); 
						y = atof(coor);
						GetChars( line, 46, 54, coor ); 
						z = atof(coor);

						wtr->H2[i][dd].x[0] = x;
						wtr->H2[i][dd].x[1] = y;
						wtr->H2[i][dd].x[2] = z;

						if(par.TIP4P)
						{
							fgets(line, lString, fi);
							//sscanf(line, "%*s%*d%*s%*s%*d%f%f%f", &x, &y, &z);
							GetChars( line, 30, 37, coor ); 
							x = atof(coor);
							GetChars( line, 38, 45, coor ); 
							y = atof(coor);
							GetChars( line, 46, 54, coor ); 
							z = atof(coor);
							
							wtr->H3[i][dd].x[0] = x;
							wtr->H3[i][dd].x[1] = y;
							wtr->H3[i][dd].x[2] = z;
						}
						//wtr->atn[i][dd] = atn;
						wtr->atn[i][dd] = fstWATatm + countAtm; //YY 4.21.2015
						/*if (wtr->atn[i][dd] > 99999)//YY 4.21.2015
						{
							printf ("countAtm is %d\n", countAtm);
							printf ("fstWATatm is %d\n", fstWATatm);
							printf ("atom number now is %d = %d\n", wtr->atn[i][dd], fstWATatm + countAtm);							
						}*/
						wtr->resn[i][dd] = resn;
						dd++;
					}			
					if(par.TIP4P) 
					{
						countAtm += 4;
					}
					else 
					{
						countAtm += 3; //YY 4.21.2015
					}
				}
				
				if(strstr(line, "END")) 
					break;	
			}
			//printf("waters inside binding pocket: %d-%d\n", dd, cc);		
			wtr->WATinside[i] = dd; 
		}
		
		wtr->largest_wat = 0;
		for(i = 0; i < numfr; i++)
		{
			if(wtr->WATinside[i] > wtr->largest_wat)
				wtr->largest_wat = wtr->WATinside[i];
		}
		
		fo = fopen("WATinside.mol2", "w");
		for(i = 0; i < numfr; i++)
		{
			fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
			fprintf(fo, "% 5d % 5d	   0	 0	   0\n", wtr->WATinside[i]*3, wtr->WATinside[i]*2);
			fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
			
			for(j = 0; j < wtr->WATinside[i]; j++)
			{
				fprintf(fo, "% 7d %-4s	  % 10.4f% 10.4f% 10.4f %-8s  1 <0>		  % 8.4f\n",
						3*j+1, "WAT",
						wtr->O[i][j].x[0], wtr->O[i][j].x[1], wtr->O[i][j].x[2],
						"O", 0.0);
				fprintf(fo, "% 7d %-4s	  % 10.4f% 10.4f% 10.4f %-8s  1 <0>		  % 8.4f\n",
						3*j+2, "WAT",
						wtr->H1[i][j].x[0], wtr->H1[i][j].x[1], wtr->H1[i][j].x[2],
						"H", 0.0);
				fprintf(fo, "% 7d %-4s	  % 10.4f% 10.4f% 10.4f %-8s  1 <0>		  % 8.4f\n",
						3*j+3, "WAT",
						wtr->H2[i][j].x[0], wtr->H2[i][j].x[1], wtr->H2[i][j].x[2],
						"H", 0.0);
			}
			fprintf(fo, "@<TRIPOS>BOND\n");
			for(j = 0; j < wtr->WATinside[i]; j++)
			{
				fprintf(fo, "%6d %5d %5d %4s\n", 2*j+1, 3*j+1, 3*j+2, "1");
				fprintf(fo, "%6d %5d %5d %4s\n", 2*j+2, 3*j+1, 3*j+3, "1");
			}
		}
		fclose(fo);
		
		fo = fopen("WaterGrid.mol2", "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d	   0	 0	   0\n", grd->numGP3, 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");

		denom = 1.0/(sqrt_2*sqrtPi*sigma*wtr->numfr);
		printf("denom: %f\n", denom);
		for(num = 0; num < grd->numGP3; num++)
		{
			grd->WaterOccupancy[num] *= denom;		// normalize occupancy
			fprintf(fo, "% 7d %-4s	  % 10.4f% 10.4f% 10.4f %-8s  1 <0>		  % 8.4f\n",
						num, "WAT",
						grd->GP_r[num][0], grd->GP_r[num][1], grd->GP_r[num][2],
						"O", grd->WaterOccupancy[num]);
		}
		fclose(fo);

// filter grid density > par.griddensity
		grid_counter = 0;
		for(gx = 0; gx < grd->NumGP_r[0]; gx++)
		{
			for(gy = 0; gy < grd->NumGP_r[1]; gy++)
			{
				for(gz = 0; gz < grd->NumGP_r[2]; gz++)
				{
					gnum = gx*grd->numGP2 + gy*grd->numGP1 + gz;
					if (grd->WaterOccupancy[gnum] > par.griddensity) 
					{
						grid_counter ++;
					}
				}
			}
		}
		fgrd->gridnum = grid_counter;


		fgrd->GP_r = (float **) calloc(fgrd->gridnum, sizeof(float *));
		fgrd->WaterOccupancy = (float*)calloc(fgrd->gridnum, sizeof(float));
		fgrd->ndx_in_grid = (int*)calloc(fgrd->gridnum, sizeof(float));
		for (i = 0; i < fgrd->gridnum; i++)
		{
			fgrd->GP_r[i] = (float *)calloc (3, sizeof(float) );
			fgrd->WaterOccupancy[i] = 0;
			fgrd->ndx_in_grid[i] = 0;
		}

		new_counter = 0; 
		for(gx = 0; gx < grd->NumGP_r[0]; gx++)
		{
			for(gy = 0; gy < grd->NumGP_r[1]; gy++)
			{
				for(gz = 0; gz < grd->NumGP_r[2]; gz++)
				{
					gnum = gx*grd->numGP2 + gy*grd->numGP1 + gz;
					if (grd->WaterOccupancy[gnum] > par.griddensity) 
					{
						fgrd->GP_r[new_counter][0] = grd->GP_r[gnum][0];
						fgrd->GP_r[new_counter][1] = grd->GP_r[gnum][1];
						fgrd->GP_r[new_counter][2] = grd->GP_r[gnum][2];
						fgrd->WaterOccupancy[new_counter] = grd->WaterOccupancy[gnum];
						fgrd->ndx_in_grid[new_counter] = gnum;
						new_counter++;
					}
				}
			}
		}
		fo = fopen("FilterGrid.mol2", "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d     0   0     0\n", fgrd->gridnum, 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");

		for(num = 0; num < fgrd->gridnum; num++)
		{
			fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>     % 8.4f\n", num, "WAT", fgrd->GP_r[num][0], fgrd->GP_r[num][1], fgrd->GP_r[num][2], "O", fgrd->WaterOccupancy[num]);
		}
		fclose(fo);


	}
	
}




//=================================================================================================
// Cluster Hydrophilic grid using QT clustering
//=================================================================================================
void ClusterHydrogrid_QT(ligand *lig, grid *grd, water *wtr, hcluster *cluster_head)
{
	hcluster	*cluster_0;
	QTClustering(grd, cluster_head, wtr);
	cluster_0 = cluster_head->next;


	OutputWaterIndex(cluster_0, par.occupancy_cutoff, wtr);	
 	ClusterFilter_QT(cluster_0, par.occupancy_cutoff, wtr);
	OutputWaterOccupancyGrid(cluster_0, grd);

}
//=================================================================================================
// Cluster Hydrophilic grid
//=================================================================================================
void QTClustering(grid *grd, hcluster *cluster_head, water *wtr)
{
	hcluster	*cluster_0, *cluster, *prev_cluster;
    float	size_water = par.size_water;
    float	size_water2 = pow(size_water, 2);
	float	sigma = size_water/3.0;
	float	sigma2 = pow(sigma, 2);
	float	denom;
    float	delta = grd->griddelta;		// grid spacing
    int		del = (int) (size_water/delta) + 1;  // number of gridpoints around center of water molecule up to 
	int		grid_cluster_dim = pow(2*del, 3); 	//maximum nnumber of points could be included in each grid-centered cluster
	int		cluster_found = 1;
	int		gx, gy, gz, gnum;
	int		k2, l2, m2, num;
	float	x, y, z, xx, yy, zz, dd2;
	int		mi, i;
	int		largest_cluster_gnum;	//grip point num of the cluster with the largest occupancy
	int		num_largest_cluster_members;	//number of grid points belong to the largest cluster
	float	largest_cluster_occupancy;		//occupancy of the largest cluster
	float	largest_cluster_avg_occupancy;
	int		*largest_cluster_members;
	int		*flag_in_cluster;				//flag of whether the grid has been clustered
	float	*grid_cluster_occupancy;
	int		**grid_cluster_members;
	int		num_clusters;
	int		mm, j, k, p, f, n, numwaters;
	float	max_occu;
	float	occupancy_cutoff = 8.0;
	int		occupancy_cutoff2 = par.occupancy_cutoff;




	cluster_0 = NULL;
	flag_in_cluster = (int*)calloc(grd->numGP3, sizeof(int));
	grid_cluster_occupancy = (float*)calloc(grd->numGP3, sizeof(float));
	grid_cluster_members = (int**)calloc(grd->numGP3, sizeof(int*));
	for(i = 0; i < grd->numGP3; i++)
	{
		flag_in_cluster[i] = 0;
		grid_cluster_members[i] = (int*)calloc(grid_cluster_dim, sizeof(int));
	}
	largest_cluster_members = (int*)calloc(grid_cluster_dim, sizeof(int));

	num_clusters = 0;
    while(cluster_found)
    {
    	largest_cluster_occupancy = -999.0;
    	largest_cluster_avg_occupancy = -999.0;
		for(gx = 0; gx < grd->NumGP_r[0]; gx++)
		{
			for(gy = 0; gy < grd->NumGP_r[1]; gy++)
			{
				for(gz = 0; gz < grd->NumGP_r[2]; gz++)
				{
					gnum = gx*grd->numGP2 + gy*grd->numGP1 + gz;
					if(flag_in_cluster[gnum] == 0)
					{
						x = grd->GP_r[gnum][0];
						y = grd->GP_r[gnum][1];
						z = grd->GP_r[gnum][2];
						grid_cluster_occupancy[gnum] = 0.0;
						mi = 0;
						for(k2 = gx - del; k2 <= gx + del; k2++)
						{
							for(l2 = gy - del; l2 <= gy + del; l2++)
							{
								for(m2 = gz - del; m2 <= gz + del; m2++)
								{
									if(k2 >= 0 && k2 < grd->NumGP_r[0] && l2 >= 0 && l2 < grd->NumGP_r[1] && m2 >= 0 && m2 < grd->NumGP_r[2])
									{
										num = k2*grd->numGP2 + l2*grd->numGP1 + m2;
										if(flag_in_cluster[num] == 0)
										{
											xx = grd->GP_r[num][0] - x;
											yy = grd->GP_r[num][1] - y;
											zz = grd->GP_r[num][2] - z;
											dd2 = xx*xx + yy*yy + zz*zz;
											if(dd2 <= size_water2)
											{
												grid_cluster_occupancy[gnum] += grd->WaterOccupancy[num];
												grid_cluster_members[gnum][mi] = num;
												mi++;
											}
										}
									}
								}
							}
						}
						if(grid_cluster_occupancy[gnum] > largest_cluster_occupancy)
						//if(grid_cluster_occupancy[gnum]/mi > largest_cluster_avg_occupancy)
						{
							largest_cluster_gnum = gnum;
							largest_cluster_occupancy = grid_cluster_occupancy[gnum];
							largest_cluster_avg_occupancy = grid_cluster_occupancy[gnum]/mi;
							num_largest_cluster_members = mi;
							for(i = 0; i < num_largest_cluster_members; i++)
								largest_cluster_members[i] = grid_cluster_members[gnum][i]; 
						}
					}
				}
			}
		}
    	printf("largest_cluster_occupancy: %f %f %d\n", largest_cluster_occupancy, largest_cluster_avg_occupancy, num_largest_cluster_members);
    	//add into cluster list
		if(largest_cluster_occupancy > par.QT_occupancy_cutoff)
		{
			cluster = (hcluster *) calloc(1, sizeof(hcluster));
			if(cluster_0 == NULL)
			{
				cluster_0 = cluster;
			}
			else
			{
				prev_cluster->next = cluster;
			}
			cluster->numpoints = num_largest_cluster_members;
			cluster->center_occupancy = grd->WaterOccupancy[largest_cluster_gnum];
			cluster->overall_occupancy = largest_cluster_occupancy;
			cluster->center[0] = grd->GP_r[largest_cluster_gnum][0];
			cluster->center[1] = grd->GP_r[largest_cluster_gnum][1];
			cluster->center[2] = grd->GP_r[largest_cluster_gnum][2];

            cluster->hydropoint = (float **) calloc(num_largest_cluster_members, sizeof(float *));
            cluster->hydroscore = (float *) calloc(num_largest_cluster_members, sizeof(float));
            cluster->grdpn = (int *) calloc(num_largest_cluster_members, sizeof(int));
            cluster->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));

            for(j = 0; j < num_largest_cluster_members; j++)
            {
            	cluster->hydropoint[j] = (float *) calloc(3, sizeof(float));
			} 
			for(f = 0; f < wtr->numfr; f++)
			{
				cluster->wtrindex[f] = -1;
			}
			max_occu = -999.0;
            for(j = 0; j < num_largest_cluster_members; j++)
            {
				num = largest_cluster_members[j];
            	for(p = 0; p < 3; p++)
            	{
                	cluster->hydropoint[j][p] = grd->GP_r[num][p];				// atom numbers
                }
                cluster->hydroscore[j] = grd->WaterOccupancy[num];
                if(grd->WaterOccupancy[num] > max_occu)
                	max_occu = grd->WaterOccupancy[num];
                cluster->grdpn[j] = num;
                
                //go over all water molecules in each frame, assign to the cluster
                for(f=0; f<wtr->numfr; f++)
                {
                	//nearestwtr = -1;
                	for(n=0; n<wtr->WATinside[f]; n++)
                	{
                		if(cluster->grdpn[j] == wtr->O[f][n].nearestgrdp)
                		{	
                    		cluster->wtrindex[f] = n;	//water index
                    		break;		
                    	}	
                    }
                }
            }
            cluster->max_occupancy = max_occu;	
           
            mm = 0;
            for(j = 0; j < wtr->numfr; j++)
            {
            	if(cluster->wtrindex[j] >= 0)
            	{
            		mm++;
            	}
            }
            cluster->wtrnum = mm;
			numwaters = mm;
            cluster->average_occupancy = largest_cluster_occupancy/(num_largest_cluster_members);
            printf("max %f %f %f %d %d %d\n", max_occu, cluster->center_occupancy, cluster->average_occupancy, cluster->numpoints, num_largest_cluster_members, numwaters);
				
            prev_cluster = cluster;
            prev_cluster->next = NULL;
			
			//mark the grid point as clustered
			flag_in_cluster[largest_cluster_gnum] = 1;
			for(i = 0; i < num_largest_cluster_members; i++)
			{
				num = largest_cluster_members[i];
				flag_in_cluster[num] = 1;
			}
			num_clusters++;
		}
		else
		{
			cluster_found = 0;
		}
	}
	par.occupancy_cutoff = 200;

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0;		//give address of cluster_0 to the input pointer !!!this is very important!!!

    // free memory
    for(i = 0; i < grd->numGP3; i++)
    {
        free(grid_cluster_members[i]);
    }
    free(grid_cluster_members);
    free(flag_in_cluster);
    free(grid_cluster_occupancy);
    free(largest_cluster_members);

	//-----------------------OUTPUT cluster MOL2------------------------//
	FILE	*fo;
	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", num_clusters, 0);
	fprintf(fo, "SMALL min %f max %f\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", 0.0, 0.0);	
	int ii = 0;
	for(cluster = cluster_0, i=0; i < cluster_0->numclusters; cluster = cluster->next, i++)
	{
		fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
                        ii, "CLT",
                        cluster->center[0], cluster->center[1], cluster->center[2],
                        "O", i+1, cluster->average_occupancy, (float)cluster->wtrnum/par.numfr);
		ii++;

	}
	fclose(fo);



}

//=================================================================================================
// OutputWaterIndex: output waterindex.wtndx file for calc enthalpy
//=================================================================================================
void OutputWaterIndex(hcluster *cluster_0, int cutoff, water *wtr)
{
	int		wtrln, frm, i, atn, add, resn, wtn, ci, cl, j;
	FILE	*fo;
	int		*atomlist, *resnlist, *cltnum, **cltlist, **frmlist;
	hcluster	*cluster;
	int		*cflag;
	printf("wtr->numfr %d*%d = %d\n", wtr->numfr, wtr->WATinside[0], wtr->WATinside[0]*wtr->numfr);
	atomlist = (int*)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int));
	resnlist = (int*)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int));
	cltnum = (int*)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int));

	cltlist = (int **)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int *));
	frmlist = (int **)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int *));
	
//	for(i = 0; i < wtr->WATinside[0]*wtr->numfr; i++)
//	{
//		cltlist[i] = (int *)calloc(wtr->numfr, sizeof(int));
//		frmlist[i] = (int *)calloc(wtr->numfr, sizeof(int));
//	}

	cflag = (int *)calloc(cluster_0->numclusters, sizeof(int));
	int	clusternum = cluster_0->numclusters;
	for(ci = 0; ci < clusternum; ci++)
	{
		cflag[ci] = 1;
	}
/*	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		if(cluster->wtrnum < cutoff)
		{
			cflag[cl] = 0;
			//printf("cutoff!\n");
		}
	}*/
	wtrln = 0;
	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		
		if(cflag[cl])
		{
			printf("cluster %d\n", cl);
			for(frm = 0; frm < wtr->numfr; frm++)
			{
				if(cluster->wtrindex[frm] >= 0)
				{
					cltlist[wtrln] = (int *)calloc(wtr->numfr, sizeof(int));
					frmlist[wtrln] = (int *)calloc(wtr->numfr, sizeof(int));

					wtn = cluster->wtrindex[frm];
					atn = wtr->atn[frm][wtn];
					resn = wtr->resn[frm][wtn];
					add = 1;
					for(i = 0; i < wtrln; i++)
					{
						if(atn == atomlist[i])
						{
							add = 0;
							//if(cl == 33 || cl == 36) printf("clutnum %d %d %d %d\n", i, cltnum[i], cl, frm);
							cltlist[i][cltnum[i]] = cl;
							frmlist[i][cltnum[i]] = frm;
							cltnum[i]++;
	
							break;
						}
					}
					if(add)
					{
						atomlist[wtrln] = atn;
						resnlist[wtrln] = resn;
						//if(cl == 33 || cl == 36) printf("wtrln %d %d %d\n", wtrln, cl, frm);
						cltlist[wtrln][0] = cl;
						frmlist[wtrln][0] = frm;
						cltnum[wtrln] = 1;
						
						wtrln++;
					}
				}
			}	
 		}			
	}
	fo = fopen("waterindex.wtndx", "w");
	fprintf(fo, "@wtr#	atom#	resn#	#of_clusters\n@    	     	cluster#	frame#\n%6d %6d %6d\n", wtrln, cluster_0->numclusters, wtr->numfr);

	for(i=0; i < wtrln; i++)
	{
		fprintf(fo, "%5d	%6d	%5d	  %4d\n", i, atomlist[i], resnlist[i], cltnum[i]);	
		for(j = 0; j < cltnum[i]; j++)
		{
			fprintf(fo, "     	     %5d	%6d\n", cltlist[i][j], frmlist[i][j]);	
		}
	}	
	for(i = 0; i < wtr->WATinside[0]*wtr->numfr; i++)
	{
		free(frmlist[i]); 
		free(cltlist[i]);
	}
	
	free(cltlist);
	free(frmlist);
	free(atomlist);
	
	free(resnlist);
	free(cltnum);
	fclose(fo);

}


//=================================================================================================
// ClusterFilter1: when maximum value in the cluster < cutoff, remove cluster
//=================================================================================================
void ClusterFilter_QT(hcluster *cluster_0, int cutoff, water *wtr)
{
	float	maxscore;
	int		*cflag;
	int		i, j, k, ii, jj, clusternum, nump = 0, cl, ci;
	FILE	*fo, *fi;
	char	outfname[kString], wtrfolder[kString], line[lString];
	hcluster	*cluster;
	int		*wtrlist, wtrln, add, output, frm, wtn, atn, *wtrlistflag;

	wtrlist = (int*)calloc(wtr->WATinside[0]*wtr->numfr, sizeof(int));
	for(i = 0; i < wtr->WATinside[0]*wtr->numfr; i++)
		wtrlist[i] = -1;

	cflag = (int *)calloc(cluster_0->numclusters, sizeof(int));
	clusternum = cluster_0->numclusters;
	for(ci = 0; ci < clusternum; ci++)
	{
		cflag[ci] = 1;
	}
	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
/*		if(cluster->wtrnum < cutoff)
		{
			cflag[cl] = 0;
			//printf("cutoff!\n");
		}
		else
		{*/
			nump += cluster->numpoints+1;
//		}
	}

	wtrln = 0;
	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		if(cflag[cl])
		{

			//generate waterindex list for output waterindex files
			for(frm = 0; frm < wtr->numfr; frm++)
			{
				if(cluster->wtrindex[frm] >= 0)
				{
					wtn = cluster->wtrindex[frm];
					atn = wtr->atn[frm][wtn];
					add = 1;
					for(i = 0; i < wtrln; i++)
					{
						if(atn == wtrlist[i])
						{
							add = 0;
							break;
						}
					}
					if(add)
					{
						wtrlist[wtrln] = atn;
						//frmlist[wtrln] = frm;
						//resnlist[wtrln] = wtn;
						wtrln++;
					}
				}
			}			
		}
	}
	//fclose(fo);

	//generate folder for output waterindex
	sprintf(wtrfolder, "wtrindex_folder");
	mkdir(wtrfolder, S_IRWXU | S_IRWXG | S_IRWXO);
	chdir(wtrfolder);
	
	printf("wtrln %d\n", wtrln);

	sprintf(outfname, "wtrindex_0000.ndx");
	fo = fopen(outfname, "w");
	fi = fopen(par.indexfile, "r");
	if(!fi)
    {
		printf("index file %s is missing.\n", par.indexfile);
		exit(0);
    }
	wtrlistflag = (int*)calloc(wtr->numatms, sizeof(int));
	
	for(j = 0; j < wtr->numatms; j++)
		wtrlistflag[j] = 1;
	
	for(i=0; i < wtrln; i++)
	{
		if(i%200==0)
		{
			while(fgets(line, lString, fi) != NULL)
			{
				fprintf(fo, line);	
			}
			rewind(fi);
		}
		
		//printf("%d\n", wtrlist[i]);
		if(par.TIP4P) 
		{
			fprintf(fo, "[ Water_%04d ]\n", i);
			fprintf(fo, "%6d %6d %6d %6d\n", wtrlist[i], wtrlist[i]+1, wtrlist[i]+2, wtrlist[i]+3);
			wtrlistflag[wtrlist[i]-1] = 0;
			wtrlistflag[wtrlist[i]] = 0;
			wtrlistflag[wtrlist[i]+1] = 0;
			wtrlistflag[wtrlist[i]+2] = 0;
		}
		else
		{
			fprintf(fo, "[ Water_%04d ]\n", i);
			fprintf(fo, "%6d %6d %6d\n", wtrlist[i], wtrlist[i]+1, wtrlist[i]+2);			
			wtrlistflag[wtrlist[i]-1] = 0;
			wtrlistflag[wtrlist[i]] = 0;
			wtrlistflag[wtrlist[i]+1] = 0;
		}
		if((i+1)%200==0 || i == wtrln-1)
		{
			fprintf(fo, "[ Environment ]\n");
			output = 0;
			for(j = 0; j < wtr->numatms; j++)
			{
				if(wtrlistflag[j])
				{
					fprintf(fo, "%6d ", j+1);
					output++;
					if(output%15==0)
						fprintf(fo, "\n");

				}
			}
			fprintf(fo, "\n");
			
			fclose(fo);
			if(i != wtrln-1)	
			{
				sprintf(outfname, "wtrindex_%04d.ndx", (i+2)/200);
				fo = fopen(outfname, "w");
			}
			
			// reset flag for new file
			for(j = 0; j < wtr->numatms; j++)
				wtrlistflag[j] = 1;

		}
		
	}
	//fclose(fo);
	fclose(fi);
	chdir("../");
}


//=================================================================================================
// OutputWaterOccupancyGrid: output wateroccupancy.grd file
//=================================================================================================
void OutputWaterOccupancyGrid(hcluster *cluster_0, grid *grd)
{
	FILE	*fo;
	int		totalnumpoints = 0;	// total number of grid points in all the clusters
	int		cl, nump;
	hcluster	*cluster;
	
	fo = fopen("wateroccupancy.grd", "w");
	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		totalnumpoints += cluster->numpoints;
	}	

	printf("totalnumpoints = %d \n", totalnumpoints);

	fprintf(fo, "%10.3f %10.3f %10.3f %10d %10d %10d %10d %10.3f\n", grd->ZeroPoint[0], grd->ZeroPoint[1], grd->ZeroPoint[2], grd->NumGP_r[0], grd->NumGP_r[1], grd->NumGP_r[2], totalnumpoints, grd->griddelta);
	for(cluster = cluster_0, cl=0; cl < cluster_0->numclusters; cluster = cluster->next, cl++)
	{
		fprintf(fo, "%5d %10.3f\n", cluster->numpoints, cluster->max_occupancy);
		for(nump = 0; nump < cluster->numpoints; nump++)
		{
			fprintf(fo, "%3d %10.3f %10.3f %10.3f %10.3f\n", cl, cluster->hydropoint[nump][0], cluster->hydropoint[nump][1], cluster->hydropoint[nump][2], cluster->hydroscore[nump]);
		}
	}	
	fclose(fo);
}
