/***************************************************************************
 *   Copyright (C) 2009 by Bingjie Hu, Markus A. Lill                      *
 *   hub@purdue.edu                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************
 ***************************************************************************
 *   This program is for calculating the entropy score of waters inside    *
 *   of the binding pocket after MD simulation. The entropy is calculated  *
 *   by integrating waters' freedom in each rotation and translation 	   *
 *   dimension using PCA and simpson's integral
 *   																	   *
 *   Input 	->1. MD snapshots: for calculating water occupancy             *
 *   		->2. ligand pdb file: for defining binding pocket	           *
 *   											                           *
 *   Output	->1. MOL2 files containing information of clusters             *
 *          ->2. files contain principle component analysis and simpson's  *
 *               integral information                                      *
 *          ->3. waterindex files for rerun MD simulation to get enthalpy  *
 *               values                                                    *
 *          ->4. water occupancy grid file (wateroccupancy.grd)            *
 *          ->5. entropy value (entropy.egy): entropy*Temp(298.15K)        *
 ***************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include "main.h"
int main(int argc, char *argv[])
{
	ligand  	lig;
	grid		grd;
	fgrid		fgrd;
	water		wtr;
	hcluster	cluster_head;

	clock_t		t0, t1;
	float		ratio;
	par.occupancy_cutoff = 150;
	ratio = 1./CLOCKS_PER_SEC;

	ReadUsersParameters(argc, argv);
	
	ReadLigand(par.LigandFile, &lig);
	
	DefineBindingPocket(&lig, &grd);
	//BulkSolvent(&grd);
	mkdir(par.outfolder, S_IRWXU | S_IRWXG | S_IRWXO);
	chdir(par.outfolder);

	
	t0 = clock();
	WaterOccupancy3(par.ProtFile, &lig, &grd, &fgrd, &wtr);
	t1 = clock();
	printf("Time for calcwateroccupaney = %f\n", ratio*(long)t1-ratio*(long)t0);

	if(par.cluster_method == 1) // using DBSCAN clustering
	{
		printf("start DBSCAN clustering ...\n");
		t0 = clock();
		Cluster_DB(&fgrd, &lig, &grd, &wtr, &cluster_head);
		t1 = clock();
		printf("Time for DBSCAN clustering = %f\n", ratio*(long)t1-ratio*(long)t0);
	}

	if(par.cluster_method == 2) // using QT clustering
	{
		printf("start clustering hydrogrid...\n");
		t0 = clock();
		ClusterHydrogrid_QT(&lig, &grd, &wtr, &cluster_head);
		t1 = clock();
		printf("Time for clustering = %f\n", ratio*(long)t1-ratio*(long)t0);
	}

	printf("start calculating entropy...\n");
	CalcEntropy(&wtr, cluster_head.next);
	printf("Entropy calculation finished\nNow rerun Gromacs to get enthalpy values...\n");

	return EXIT_SUCCESS;
}
