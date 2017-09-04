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




#ifndef PiecewiseF
#define PiecewiseF
#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Fraction.h"


typedef struct COEFF {
	double a, b, c, d;
} TypeCoeff;

typedef struct COEFF_SPECIAL {
	double a1, a2, b, c1, c2, d;
} TypeCoeffSpecial;


typedef struct FRACTION_INTER {
	TypeFraction inf, sup;
} TypeFractionInter;

/*a 'column' of a piecewise bilinear function*/
/*by convention, a special case occurs when sizeVal = -1: boundVal[0] contains a value v and 
 * the function is such that f(g,v)=coeff[0].c*g+coeff[0].d and f(g,x) = infty for a x != v*/
typedef struct COLUMN {
	int sizeVal;
	TypeCoeff *coeff;
	double *boundVal;
} TypeColumn;

typedef struct COLUMN_SPECIAL {
	int sizeVal;
	TypeCoeffSpecial *coeff;
	double *boundVal;
} TypeColumnSpecial;

/*Type for storin a piecewise bilinear function*/
typedef struct PIECEWISE_FUNCTION {
	int sizePar;
	TypeFraction *boundPar;
	TypeColumn *col;
} TypePiecewise;

int isConvex(TypePiecewise *F);
double intPiecewise(TypePiecewise *f);
double sumColumnSpecial(TypeColumnSpecial *col, double xi, double xs, double yi, double ys);
void sumPiecewise(TypePiecewise **f, int size, TypePiecewise *sum);
void treatColumSpecial(FILE *f, TypeColumnSpecial *col, double xi, double xs, double yi, double ys);
void fprintColumSpecialFig(FILE *f, TypeColumnSpecial *col, double xi, double xs, double yi, double ys);
void getColumnSpecial(TypeColumn *colA, TypeColumn *colB, TypeColumnSpecial *col);
void freeCol(TypeColumn c);
void freePiecewise(TypePiecewise f);
/*print result as flat list*/	
void fprintFunc(FILE *f, TypePiecewise *func);
void fprintColSpeDebug(FILE *f, TypeColumnSpecial *c);
void fprintFuncGnuplot(FILE *f, TypePiecewise *func);
void fprintFuncDebug(FILE *f, TypePiecewise *func);
void fprintColDebug(FILE *f, TypeColumn *c);
void fprintfsDebugX(FILE *f, TypeCoeff c);
void fprintfsSpecialDebugX(FILE *f, TypeCoeffSpecial c);
void fprintColDebugX(FILE *f, TypeColumn *c);
void fprintColSpecialDebugX(FILE *f, TypeColumnSpecial *c);
#endif
