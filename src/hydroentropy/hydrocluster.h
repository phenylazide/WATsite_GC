#include <sys/stat.h>
#include "molecule.h"
#include "util.h"
#include <unistd.h>
#include <time.h>

#ifndef _HYDROCLUSTER_H
#define _HYDROCLUSTER_H

void ReadLigand(char *, ligand *);
void DefineBindingPocket(ligand *, grid *);
void WaterOccupancy(char *, ligand *, grid *, water *);
void WaterOccupancy2(char *, ligand *, grid *, water *);
void WaterOccupancy3(char *, ligand *, grid *, fgrid *, water *);
void ClusterHydrogrid(ligand *, grid *, water *, hcluster *);
void ClusterHydrogrid_QT(ligand *, grid *, water *, hcluster *);
void QTClustering(grid *, hcluster *, water *);
void ClusterFilter_QT(hcluster *, int , water *);
void ClusterFilter1(hcluster *, float, water *);
//void ClusterFilter2(hcluster *, float, water *);
void OutputWaterIndex(hcluster *, int, water *);
void OutputWaterOccupancyGrid(hcluster *, grid *);
void BulkSolvent(grid *);

#endif
