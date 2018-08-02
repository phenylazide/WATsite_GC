#include<stdio.h>
#include<stdlib.h>
#include "molecule.h"

#ifndef _DBSCAN_H
#define _DBSCAN_H

#define EPSILON 2
//#define MIN_POINTS 200

void Cluster_DB(fgrid *fgrd, ligand *lig, grid *grd, water *wtr, hcluster *cluster_head);
void dbScan(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head);
void BulkCluster(fgrid *fgrd, grid *grd, water *wtr, hcluster *cluster_head);
int regionQuery(int , int,  int , int , float **, int *, int);
void expandCluster(float , int , int,  int , int , float **, int *, int *, int *, int, int *);
int check(float *, int);

void OutputWaterIndex(hcluster *, int, water *);
void ClusterFilter_QT(hcluster *, int , water *);
void OutputWaterOccupancyGrid(hcluster *, grid *);

#endif
