#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <sys/stat.h>
#include "molecule.h"
#include "util.h"
#include <unistd.h>

#ifndef _ENTROPY_H
#define _ENTROPY_H
void CalcEntropy(water *, hcluster *);
void EulerAngle(water *, hcluster *);
void CovarianceMatrix6x6(water *, hcluster *);
void CovarianceMatrix3x3(water *, hcluster *);
void ProbabilityDistributionFunction(water *, hcluster *);
void SimpsonsIntegral(water *, hcluster *);
#endif
