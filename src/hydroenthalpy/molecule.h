
#ifndef _MOLECULE_H
#define _MOLECULE_H

typedef struct water
{
//	int		**atomlist;
//	int		**resnlist;
	int		*clusterlist;
	int		*framelist;
	int		numpresent;		//number of times this water present in all the clusters
} water;

typedef struct energy
{
//	int		**atomlist;
//	int		**resnlist;
	float	*enthalpy;
	float	*entropy;
	float	*free_energy;
	float	**coor;
	float	*radius;
	int		*occupancy;
	float		*occupancy_probability;
} energy;

typedef struct pharmacophore_grid
{
    int			numpoints;
    int			*clusternum;
    float		*free_energy;
    float		griddelta;
    int			numGP[3];
    int			numGP1, numGP2, numGP3;
    float		ZeroPoint[3];

	float		*WaterOccupancy;

} grid;
/*
typedef struct hydro_cluster
{
	int			num;
	float		**hydropoint;	//coor of grid point inside the cluster
	float		*hydroscore;	//occupancy for each point
	int			numpoints;

	float		maxoccupancy;	//maximum occupancy among the cluster points
	float		*occupancy;
	float		*free_energy;
}hcluster;*/

#endif
