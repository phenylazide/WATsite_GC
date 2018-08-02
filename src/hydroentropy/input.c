#include "input.h"

//=================================================================================================
// ReadUsersParameters
//=================================================================================================
void ReadUsersParameters(int argc, char *argv[])
{
    char    filecomm[kString];
    int     j;
    // default: 
    sprintf(filecomm, "HydroEntro.bcf");
    sprintf(par.outfolder, "OUTPUT");

    // Read command line arguments
    // Standard name of command file
    
    for(j = 1; j < argc; j+=2)
    {
        if(!strcmp(argv[j], "-c"))
        {
            strcpy(filecomm, argv[j+1]);
        }
        else if(!strcmp(argv[j], "-O") || !strcmp(argv[j], "-o"))
        {
        	strcpy(par.outfolder, argv[j+1]);
        }
        else if(!strcmp(argv[j], "-H") || !strcmp(argv[j], "-h"))
        {
            printf("Usuage: hydroentropy \n  -c       HydroEntro.bcf \n  -o       Output folder \n");
            exit(0);
        }
    }

	//read command file
	printf("command file %s\n", filecomm);
	ReadCommandFile(filecomm);

    // Overwrite parameters with external flags
    for(j = 1; j < argc; j+=2)
    {
        if(!strcmp(argv[j], "-f"))
        {
//        strcpy(par->filepdb, argv[j+1]);
        }
    }
}

//=================================================================================================
// ReadCommandFile
//=================================================================================================
void ReadCommandFile(char *filename)
{
    FILE    *fp, *fi;
    char    line[lString];
	
	fp = fopen(filename, "r");
    if (!fp)
    {
		printf("Command file %s is missing.\n", filename);
		exit(0);
    }
    else
    {
		// input of ligands
		while(!feof(fp))
		{
	        fgets(line, lString, fp);
			if(strstr(line, "$Input_Ligand:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%s", par.LigandFile);
		printf("ligand file: %s\n", par.LigandFile);

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, lString, fp);
			if(strstr(line, "$Input_MDsnapshots:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%s", par.ProtFile);
		printf("protein file: %s\n", par.ProtFile);
        fgets(line, lString, fp);
		sscanf(line, "%d", &par.numfr);
		printf("number of frames: %d\n", par.numfr);
        fgets(line, lString, fp);
		sscanf(line, "%d", &par.TIP4P);
		printf("water model: %d (0=SPC, 1=TIP4P)\n", par.TIP4P);
		fgets(line, lString, fp);
		sscanf(line, "%s", par.indexfile);
		printf("index file: %s\n", par.indexfile);

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$Grid_Parameters:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.size_water));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.griddelta));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.FarestDist));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.griddensity));
		printf("griddelta %f\n", par.griddelta);

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, lString, fp);
			if(strstr(line, "$Hydro_Cluster:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%d", &(par.number_clusters));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.cluster_mean_dist));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.maxdist));
        fgets(line, lString, fp);
		sscanf(line, "%f", &(par.QT_occupancy_cutoff));

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, lString, fp);
			if(strstr(line, "$Cov_parameters:"))
				break;
		}
        fgets(line, lString, fp);
		sscanf(line, "%d", &(par.covdimension));

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$Entropy:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.numbins));
		fgets(line, lString, fp);
		sscanf(line, "%f", &(par.NeatS));

		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$DBSCAN:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.dbscan_START));
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.dbscan_END));



		rewind(fp);
		while(!feof(fp))
		{
	        fgets(line, kString, fp);
			if(strstr(line, "$Cluster_method:"))
				break;
		}
		fgets(line, lString, fp);
		sscanf(line, "%d", &(par.cluster_method));

	}
	fclose(fp);
}
		
