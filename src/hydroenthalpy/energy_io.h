
#ifndef _ENERGY_IO_H
#define _ENERGY_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"
#include "molecule.h"
void ReadUsersParameters(int, char **);
void ReadCommandFile(char *);
water* ReadWaterIndexFile(void);
void ReadEnthalpyFiles(water *, energy *);
void ReadEntropyFile(energy *);
//void ReadWaterOccupancyGrid(energy *, grid *);
//void ReadSAMGP(grid *);
//int* GetGridPoint(float *, grid *, int);
void OutputClusterInfo(energy *);
#endif
