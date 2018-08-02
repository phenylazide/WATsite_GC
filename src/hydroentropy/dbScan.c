#include "dbscan.h"
#include <stdio.h>  
#include <stdlib.h>  
#include <unistd.h>
#include <sys/stat.h>

void Cluster_DB(fgrid *fgrd, ligand *lig, grid *grd, water *wtr, hcluster *cluster_head)
{

	// &fgrd, &lig, &grd, &wtr, &cluster_head;

	dbScan( fgrd, grd, wtr, cluster_head );

	hcluster	*cluster_0;
	cluster_0 = cluster_head->next;

	par.occupancy_cutoff = 200;

	OutputWaterIndex(cluster_0, par.occupancy_cutoff, wtr);	
	ClusterFilter_QT(cluster_0, par.occupancy_cutoff, wtr);
	OutputWaterOccupancyGrid(cluster_0, grd);

}


//=================================================================================================
// DBSCAN clustering algorithm: output Cluster_DB.mol2 & HydrationSites.mol2 
//=================================================================================================

void dbScan(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head)
{

	printf("wtr->numfr %d*%d = %d\n", wtr->numfr, wtr->WATinside[0], wtr->WATinside[0]*wtr->numfr);

	int	DATASET_SIZE = fgrd->gridnum;

	int	FEATURES = 3; // dimensions of data

	float **data;
	int *clusters;
	int *clustered;
	int *visited;
	int *neigh_points;
	int i, j;

	data = (float**)calloc(DATASET_SIZE, sizeof(float*));		
	for(i = 0; i < DATASET_SIZE; i++)
	{
		data[i] = (float*)calloc(FEATURES, sizeof(float));
	}

	for(i = 0; i < DATASET_SIZE; i++)
	{
		for(j = 0; j < FEATURES; j++)
		{
			data[i][j] = fgrd->GP_r[i][j];
		}
	}

	int next_cluster = 1, num_npoints;
	int MIN_POINTS = par.dbscan_START;
	clusters = (int*)calloc(DATASET_SIZE, sizeof(int));
	clustered = (int*)calloc(DATASET_SIZE, sizeof(int));

	while( MIN_POINTS >= par.dbscan_END ) //CHANGE ME
	{
		printf("%d\n", MIN_POINTS);
		visited = (int*)calloc(DATASET_SIZE, sizeof(int));
		neigh_points = (int*)calloc(DATASET_SIZE*DATASET_SIZE, sizeof(int));	

		for(i = 0; i < DATASET_SIZE; i++)
		{
			if(!visited[i] & !clustered[i])
			{
				visited[i] = 1;
				num_npoints = regionQuery(0, i, DATASET_SIZE, FEATURES, data, neigh_points, MIN_POINTS);
				
				if(num_npoints > MIN_POINTS)
				{
					expandCluster(next_cluster, num_npoints, i,  DATASET_SIZE, FEATURES, data, neigh_points, clusters, visited, MIN_POINTS, clustered);
					next_cluster++;
				}
			}
		}
		printf("%d\n", next_cluster );
		MIN_POINTS -= 10;
	}	

	int		num_clusters;
	int		mm, qq, k, p, f, n, num, numwaters;
	float	max_occu;
	float	size_water = par.size_water;
	float	size_water2 = pow(size_water, 2);
	float	delta = grd->griddelta;		// grid spacing
	int		del = (int) (size_water/delta) + 1;  // number of gridpoints around center of water molecule up to 
	int		grid_cluster_dim = pow(2*del, 3); 	//maximum nnumber of points could be included in each grid-centered cluster
	int		*largest_cluster_members;
	int		*WAT;

	int counter;
	float x,y,z, occu;
	float **save_cluster = (float**)calloc(next_cluster, sizeof(float*));
	for(i = 0; i < next_cluster; i++)
	{
		save_cluster[i] = (float*)calloc(5, sizeof(float));
	}


	hcluster *cluster_0, *prev_cluster, *cluster_new;
	cluster_0 = NULL;
	num_clusters = 0;

	for (j=1; j<next_cluster; j++)
	{
//		printf("cluster index  j = %d\n", j);

		counter = 0;
		x = 0;
		y = 0;
		z = 0;
		occu = 0;
		largest_cluster_members = (int*)calloc(grid_cluster_dim, sizeof(int));

		for(i = 0; i < DATASET_SIZE; i++)
		{
			if (clusters[i] == j)
			{
				x += fgrd->GP_r[i][0];
				y += fgrd->GP_r[i][1];
				z += fgrd->GP_r[i][2];
				occu += fgrd->WaterOccupancy[i];
				//printf("ndx in grid   %d\n", fgrd->ndx_in_grid[i]);
				largest_cluster_members[counter] = fgrd->ndx_in_grid[i];
				counter ++;
			}
		}		
//		printf("number of gridpoints(%d) for counter j=%d\n", counter,j );	

//		if( occu > par.QT_occupancy_cutoff)
//		{
		cluster_new = (hcluster *) calloc(1, sizeof(hcluster));
		if(cluster_0 == NULL)
		{
			cluster_0 = cluster_new;
		}
		else
		{
			prev_cluster->next = cluster_new;
		}

		cluster_new->numpoints = counter;

		save_cluster[num_clusters][4] = counter;
		save_cluster[num_clusters][0] = x/counter;
		save_cluster[num_clusters][1] = y/counter;
		save_cluster[num_clusters][2] = z/counter;
		save_cluster[num_clusters][3] = occu;

		cluster_new->overall_occupancy = occu;
		cluster_new->average_occupancy = occu/counter;
		cluster_new->center[0] = x/counter;
		cluster_new->center[1] = y/counter;
		cluster_new->center[2] = z/counter;

		//printf("occu %f, ave_occu %f, x %f, y %f, z %f\n", cluster_new->overall_occupancy, cluster_new->average_occupancy, cluster_new->center[0], cluster_new->center[1], cluster_new->center[2]);

		//printf("wtr->numfr = %d\n", wtr->numfr);
		WAT = (int*)calloc(counter, sizeof(int));


		cluster_new->wtrindex = (int *) calloc(wtr->numfr, sizeof(int));
		cluster_new->hydropoint = (float **) calloc(counter, sizeof(float *));
		cluster_new->hydroscore = (float *) calloc(counter, sizeof(float));
		cluster_new->grdpn = (int *) calloc(counter, sizeof(int));

		for (i = 0; i < counter; i++)
		{
			cluster_new->hydropoint[i] = (float *) calloc(3, sizeof(float));
		}
		//printf("iii = %d\n", i);
		for(f = 0; f < wtr->numfr; f++)
		{
			cluster_new->wtrindex[f] = -1;
		}

		max_occu = -999.0;
		//printf("max_occu %f\n", max_occu);

		printf("NOW cluster = %d\n", num_clusters );

		for(i = 0; i < counter; i++)
		{
			//printf("gridpoint index = %d\n", i);
			num = largest_cluster_members[i];
			//printf("gridpoint number = %d\n", num);

			for(p = 0; p < 3; p++)
			{
				//printf("xyz  %f\n", grd->GP_r[num][p]);
				cluster_new->hydropoint[i][p] = grd->GP_r[num][p];				// atom numbers
			}
			cluster_new->hydroscore[i] = grd->WaterOccupancy[num];
			if(grd->WaterOccupancy[num] > max_occu)
			{
				max_occu = grd->WaterOccupancy[num];
				printf("%f\n", max_occu);
			}
			//printf("REACHED HERE\n");
			cluster_new->grdpn[i] = num;

			//go over all water molecules in each frame, assign to the cluster
//			mm = 0;
			for(f=0; f<wtr->numfr; f++)
			{
				//nearestwtr = -1;
				mm=0;

				for(n=0; n<wtr->WATinside[f]; n++)
				{
					if(cluster_new->grdpn[i] == wtr->O[f][n].nearestgrdp)
					{	
						cluster_new->wtrindex[f] = n;	//water index
						WAT[mm] = n;
						mm++;
						break;		
					}	
				}
/*				qq = mm;
				printf("qq  %d\n", qq);
				if (qq > 0)
				{
					printf("Cluster ndx %d   ", j-1);
					printf("Grid ndx %d   ", i);
					printf("frame %d   gnum = %d\n", f, num);

					//printf("number of water near this gridpoint %d\n", mm);
					//for (p=0; p<qq; p++) 
					//{
					//	printf("WAT[%d] index %d\n",p, WAT[p]);
					//}
						
				}   
*/				 
			}
		}
		cluster_new->max_occupancy = max_occu;	

		mm = 0;
		for(f = 0; f < wtr->numfr; f++)
		{
			if(cluster_new->wtrindex[f] >= 0)
			{
				mm++;
			}
		}
		cluster_new->wtrnum = mm;
		numwaters = mm;
		printf("max %f %f %d %d %d\n", cluster_new->overall_occupancy, cluster_new->average_occupancy, cluster_new->numpoints, counter, numwaters);
			
		prev_cluster = cluster_new;
		prev_cluster->next = NULL;

		num_clusters ++;
//		}

	}

	printf("num_clusters %d\n", num_clusters);

	cluster_0->numclusters = num_clusters;
	cluster_head->next = cluster_0;		//give address of cluster_0 to the input pointer !!!this is very important!!!
	// free memory
	free(largest_cluster_members);
	free(WAT);




	//-----------------------OUTPUT cluster MOL2------------------------//
	FILE  *fo;
	fo = fopen("HydrationSites.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0     0     0\n", num_clusters, 0);
	fprintf(fo, "SMALL min %f max %f\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n", 0.0, 0.0);	
	int ii = 0;
	for(cluster_new = cluster_0, i=0; i < cluster_0->numclusters; cluster_new = cluster_new->next, i++)
	{
		printf("% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
						ii, "CLT",
						cluster_new->center[0], cluster_new->center[1], cluster_new->center[2],
						"O", i+1, cluster_new->overall_occupancy, (float)cluster_new->wtrnum/par.numfr);
		fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>       % 8.4f %5.2f\n",
						ii, "CLT",
						cluster_new->center[0], cluster_new->center[1], cluster_new->center[2],
						"O", i+1, cluster_new->overall_occupancy, (float)cluster_new->wtrnum/par.numfr);
		ii++;

	}
	fclose(fo);

	//-----------------------OUTPUT cluster MOL2------------------------//

/*	FILE  *fo;
	fo = fopen( "Cluster_DBSCAN.mol2", "w");
	fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
	fprintf(fo, "% 5d % 5d     0   0     0\n", next_cluster-1, 0);
	fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
	for (j=0; j<num_clusters; j++)
	{
		fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  %1d <0>     % 8.4f\n", j, "WAT", save_cluster[j][0], save_cluster[j][1], save_cluster[j][2], "O", j+1,save_cluster[j][3]);
	}
	fclose(fo);


	mkdir("DBSCAN_mol2", S_IRWXU | S_IRWXG | S_IRWXO);
	chdir("DBSCAN_mol2");
	for (j=0; j<next_cluster-1; j++)
	{
		char fname[32]; 
		snprintf(fname, sizeof(char) * 32, "Cluster_DBSCAN_%i.mol2", j);
		fo = fopen( fname, "w");
		fprintf(fo, "@<TRIPOS>MOLECULE\n****\n");
		fprintf(fo, "% 5d % 5d     0   0     0\n", (int)save_cluster[j][4], 0);
		fprintf(fo, "SMALL\nUSER_CHARGES\n\n\n@<TRIPOS>ATOM\n");
	
		int counter = 0;
		for(i = 0; i < DATASET_SIZE; i++)
		{
			if (clusters[i] == j+1)
			{
				fprintf(fo, "% 7d %-4s    % 10.4f% 10.4f% 10.4f %-8s  1 <0>     % 8.4f\n", counter, "WAT", fgrd->GP_r[i][0], fgrd->GP_r[i][1], fgrd->GP_r[i][2], "O", fgrd->WaterOccupancy[i]);
				counter ++;
			}
		}
	//	printf("%f  %f  %f\n", save_cluster[j][0]/counter, save_cluster[j][1]/counter, save_cluster[j][2]/counter);
		fclose(fo);

	} 
	chdir("../");
*/

}



void expandCluster(float cluster_no, int num_npoints, int index, int DATASET_SIZE, int FEATURES, float **data, int *neigh_points, int *clusters, int *visited, int MIN_POINTS, int *clustered)
{
	clusters[index] = cluster_no;
	clustered[index] ++;
	data[index][0] = 0; 
	data[index][1] = 0; 
	data[index][2] = 0; 

	int i, count = 0;
	
	for(i = 0; i < num_npoints; i++)
	{
		if(!visited[neigh_points[i]])
		{
			visited[neigh_points[i]] = 1;
			
			count = regionQuery(num_npoints, neigh_points[i], DATASET_SIZE, FEATURES, data, neigh_points, MIN_POINTS);
			
			if(count >= MIN_POINTS)
			{
				num_npoints += count;
			}
		}
		
		if(!clusters[neigh_points[i]])
		{
			clusters[neigh_points[i]] = cluster_no;
			data[neigh_points[i]][0] = 0; 
			data[neigh_points[i]][1] = 0; 
			data[neigh_points[i]][2] = 0; 
			clustered[neigh_points[i]] ++;
		}
	}
}

int regionQuery(int start, int index, int DATASET_SIZE, int FEATURES, float **data, int *neigh_points, int MIN_POINTS)
{
	int i, j, count = 0;
	float distance, temp;

	//printf("function regionQuery\n");	
	for(i = 0; i < DATASET_SIZE; i++)
	{
		if(i != index)
		{
			distance = 0;
		
			for(j = 0; j < FEATURES; j++)
			{
				temp = data[i][j] - data[index][j];
			
				distance += temp * temp;
			}
		
			if(distance <= EPSILON)
			{
				//printf("distance %f\n", distance );
				neigh_points[start+count] = i;
				count++;
			}
		}
	}
	return count;
}



//=================================================================================================
// Bulk of DBSCAN Clutering
//=================================================================================================

void BulkCluster(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head)
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

	float 	cluster_found;
	float	cluster_bulk[5] = {0.25, 0.333, 0.5, 0.667, 0.75};
	int		gx, gy, gz, gnum;
	int		k2, l2, m2, num;
	float	x, y, z, xx, yy, zz, dd2;
	int		mi, i, ii;
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

	for ( ii = 0; ii < 5; ii++ )
	{
		cluster_found = cluster_bulk[ii];
		printf("Cluster_found = %f\n", cluster_found);

		largest_cluster_occupancy = -999.0;
		largest_cluster_avg_occupancy = -999.0;

		gx = round(grd->NumGP_r[0] * cluster_found);
		gy = round(grd->NumGP_r[1] * cluster_found);
		gz = round(grd->NumGP_r[2] * cluster_found);
		gnum = gx*grd->numGP2 + gy*grd->numGP1 + gz;
		printf("gx gy gz gnum %d %d %d %d\n", gx, gy, gz, gnum);

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

		printf("largest_cluster_occupancy: %f %f %d\n", largest_cluster_occupancy, largest_cluster_avg_occupancy, num_largest_cluster_members);

		//add into cluster list
		if(largest_cluster_occupancy > 0)
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
		printf("num_clusters %d\n", num_clusters);
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
	ii = 0;
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

