/*
	'ParSplit' - Identification of the most parsimonious split(s) of a phylogenetic tree with regards of character values at tips.
	'Fix' - Fixing data to be read by ParSplit

    Copyright (C) 2017 Gilles DIDIER

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#ifndef ContinuousParsimonyF
#define ContinuousParsimonyF
#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Piecewise.h"
#include "StateTree.h"



typedef double TypeFunctTime(double);



typedef struct RESULT {
	int sizePar;
	TypeFraction *boundPar;
	double *val;
} TypeResult;

typedef struct COEFF_RESULT_FUNC {
	double a, b;
} TypeCoeffResultFunc;

typedef struct RESULT_FUNC {
	int sizePar;
	TypeFraction *boundPar;
	TypeCoeffResultFunc *coeff;
} TypeResultFunc;

TypeResult *getResult(TypeTree *tree, TypeFunctTime phi);
void fillNormFH(TypeTree *tree, TypeFunctTime phi, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_);
void getAllCut(TypeTree *tree, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_, int **tnmin, int **ttype, double **tmin);
void getMinCut(TypeTree *tree, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_, int *nmin, int *type, double *);
TypeResult getResultFromF(TypePiecewise *f, int unknown);
TypeResult getResultFromF_(TypePiecewise *h, TypeResult *a, int unknown);
TypeResultFunc getMinimumCost(TypePiecewise *F);
void freeResult(TypeResult r);
/*print result as flat list*/	
void fprintResultList(FILE *f, TypeResult *res, TypeTree *tree);
/*print result as flat list*/	
void fprintResultList(FILE *f, TypeResult *res, TypeTree *tree);
void fprintResult(FILE *f, TypeResult *r);
void fprintResultDebug(FILE *f, TypeResult *r);
void fprintStateTreeResult(FILE *f, TypeResult *res, TypeTree *tree);
void sprintResult(char *s, TypeResult *r);
void sprintResultDebug(char *s,  TypeResult *r);
void fprintResultDebug(FILE *f, TypeResult *r);
void fprintResultFuncDebug(FILE *f, TypeResultFunc *r);
void compactResult(TypeResult *r);
void compactResultFunc(TypeResultFunc *r);

#endif
