#include <stdlib.h>
#include <stdio.h>

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#define kString 254
#define lString 508

#define Null 0.00001
#define Pi 				3.141592654
#define sqrtPi			1.772453851
#define sqrt_2			1.414213562
#define E 				2.718281828
// gas constant in cal/(K*mol)
#define R 				1.9858775		
// room temperature in kelvin
#define T0				298.15	
// degrees to radians		
#define	DtoR			0.0174532925199	
#define RtoD			57.2957795131


struct parameters
{
    char        ParameterFile[lString];
    char        LigandFile[lString];
    char		ProtFile[lString];
    char		outfolder[lString];
	float		griddelta;
	float		griddensity;
	float		FarestDist;
	int			number_clusters;
	float		cluster_mean_dist;
	float		maxdist;
	int			covdimension;
	float		size_water;
	int			numfr;
	int			TIP4P;
	int			occupancy_cutoff;
	char		indexfile[lString];
	int			numbins;
	float		NeatS;
	float		QT_occupancy_cutoff;
	int			dbscan_START;
	int			dbscan_END;
	int			cluster_method;
} par;

#endif
