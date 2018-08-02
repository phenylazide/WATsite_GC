#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#define kString 254
#define lString 508
#define Null 0.00001
#define Pi 				3.141592654
#define E 				2.718281828
// gas constant in cal/(K*mol)
#define R 				1.9858775		
// room temperature in kelvin
#define T0				298.15	
// degrees to radians		
#define	DtoR			0.0174532925199	
#define RtoD			57.2957795131

//enthalpy for bulk solvent water in kcal/mol (-8.13 from experimental data, -18.5 from Freisner paper)
#define BulkSolE		-8.13
//entropy for bulk solvent in cal/mol-k (16.42 from experimental data, 5.03 from Freisner paper)
#define BulkSolS		16.42

struct parameters
{
    char        waterindex[kString];

    char        wtrocupygrd[kString];
    char		entropy[kString];
    char		enthalpy_dir[kString];

	char		outfolder[kString];
	
	char		SAMGP0[kString];
	char		SAMGP1[kString];
	char		SAMGP2[kString];

	int			numwaters;
	int			numclusters;
	int			numframes;
	int			outputcluster;			//number of clusters output for docking (cutoff)
	
	float		Ecutoff, Scutoff, Ereward, TSreward;
	float		E_bulk;
	
	float		cluster_x0, cluster_y0, cluster_z0;		//zero point of cluster grid
} par;

#endif
