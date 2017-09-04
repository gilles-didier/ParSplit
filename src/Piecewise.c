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




#include <limits.h>
#include <string.h>
#include <math.h>
#include "Floating.h"
#include "Piecewise.h"

#define RADIUS 0.5
#define XOFFSET 0.0
#define YOFFSET 0.0
#define XOFF 0.2
#define YOFF -0.4
#define SIZEBUF 10
#define FZERO {0,1}
#define FINFTY {1,0}
#define FNINFTY {-1,0}
#define INFTY 1E20
#define SIZE_BUFFER_CHAR 2000
#define WIDTH "stringWidth"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define FALSE 0
#define TRUE 1
#define EPSI 0.00000000001
static int isConvexCol(TypeColumn *col, double x, double y);
static double sumGen(double bi, double bs, double di, double ds, double ui, double us, TypeCoeffSpecial co, double v);

double sumGen(double bs, double be, double us, double ue, double ds, double de, TypeCoeffSpecial co, double v) {
	double O, P, Q, R;
	if(bs == be)
		return 0.;
	O = (de-ds)/(be-bs);
	P = ds-O*bs;
	Q = (ue-us)/(be-bs);
	R = us-Q*bs;
	return 
		0.5*(((pow(be, 3.)-pow(bs, 3.))/3.)*(2.*(co.a1*v+co.c1)*(Q-O)+(co.a2*v+co.c2)*(pow(Q, 2.)-pow(O, 2.)))+
		(pow(be, 2.)-pow(bs, 2.))*((co.a1*v+co.c1)*(R-P)+(co.b*v+co.d)*(Q-O)+(co.a2*v+co.c2)*(Q*R-O*P))+
		(be-bs)*(2.*(co.b*v+co.d)*(R-P)+(co.a2*v+co.c2)*(pow(R, 2.)-pow(P, 2.))));
}

void freeCol(TypeColumn c) {
	if(c.boundVal != NULL)
		free((void*)c.boundVal);
	if(c.coeff != NULL)
		free((void*)c.coeff);
}

void freePiecewise(TypePiecewise f) {
	if(f.boundPar != NULL)
		free((void*)f.boundPar);
	if(f.col != NULL) {
		int i;
		for(i=0; i<=f.sizePar; i++)
			freeCol(f.col[i]);
		free((void*)f.col);
	}
}

int isConvexCol(TypeColumn *col, double x, double y) {
	int i;
	for(i=0; i<col->sizeVal && 
		col->coeff[i].a*x+col->coeff[i].b-(col->coeff[i+1].a*x+col->coeff[i+1].b)<=EPSI &&
		col->coeff[i].a*y+col->coeff[i].b-(col->coeff[i+1].a*y+col->coeff[i+1].b)<=EPSI; i++)
		;
if(i<col->sizeVal)
printf("i %d %.2le %.2le (a%.2le b%.2le -- a%.2le b%.2le)\n", i, col->coeff[i].a*x+col->coeff[i].b-(col->coeff[i+1].a*x+col->coeff[i+1].b), col->coeff[i].a*y+col->coeff[i].b-(col->coeff[i+1].a*y+col->coeff[i+1].b), col->coeff[i].a, col->coeff[i].b, col->coeff[i+1].a, col->coeff[i+1].b);
	return i>=col->sizeVal;
}
		
int isConvex(TypePiecewise *F) {
	int p;
	for(p=1; p<F->sizePar && isConvexCol(&(F->col[p]), F->boundPar[p-1].num/F->boundPar[p-1].den, F->boundPar[p].num/F->boundPar[p].den); p++)
		;
if(p<F->sizePar)
printf("p %d (%.2le-%.2le)\n", p, F->boundPar[p-1].num/F->boundPar[p-1].den, F->boundPar[p].num/F->boundPar[p].den);
	return p>=F->sizePar;
}

void fprintCS(FILE *f, TypeCoeffSpecial c) {
	fprintf(f, "%.2lfx+%.2lfy+%.2lf\n", c.a1, c.a2, c.b);
}

double intPiecewise(TypePiecewise *f) {
	int i;
	double sum;
	sum = 0.;
	for(i=0; i<=f->sizePar; i++) {
		int j, jmin, jmax;
		double boundP, boundF;
		if(i>0)
			boundP = f->boundPar[i-1].num/f->boundPar[i-1].den;
		else			
			boundP = 0.;
		if(i<f->sizePar)
			boundF = f->boundPar[i].num/f->boundPar[i].den;
		else
			boundF = 1.;
		for(jmin=0; jmin<=f->col[i].sizeVal && cmpDouble(f->col[i].coeff[jmin].b/-f->col[i].coeff[jmin].a, boundP) <= 0; jmin++);
		for(jmax=f->col[i].sizeVal; jmax>0 && cmpDouble(f->col[i].coeff[jmax].b/-f->col[i].coeff[jmax].a, boundF) >= 0; jmax--);
		for(j=jmin; j<=jmax; j++) {
			double boundC = f->col[i].coeff[j].b/-f->col[i].coeff[j].a;
			sum += (f->col[i].coeff[j].a*f->col[i].boundVal[j-1]+f->col[i].coeff[j].c)*(pow(boundC, 2.)-pow(boundP, 2.))/2.+(f->col[i].coeff[j].b*f->col[i].boundVal[j-1]+f->col[i].coeff[j].d)*(boundC-boundP);
			boundP = boundC;
		}
		sum += (f->col[i].coeff[j].a*f->col[i].boundVal[j-1]+f->col[i].coeff[j].c)*(pow(boundF, 2.)-pow(boundP, 2.))/2.+(f->col[i].coeff[j].b*f->col[i].boundVal[j-1]+f->col[i].coeff[j].d)*(boundF-boundP);
	}
	return sum;
}

double sumColumnSpecial(TypeColumnSpecial *col, double xi, double xs, double yi, double ys) {
	int i, ind, *indice, last;
	double *xstart, *ystart, *xend, *yend, sum = 0.;
	indice = (int*) malloc((col->sizeVal+1)*sizeof(int));
	xstart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	ystart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	xend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	yend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	ind = 0;
	for(i=0; i<=col->sizeVal && cmpDouble(col->coeff[i].a1*xi+col->coeff[i].b+col->coeff[i].a2*yi, 0.) <= 0 && cmpDouble(col->coeff[i].a1*xs+col->coeff[i].b+col->coeff[i].a2*yi, 0.) <= 0; i++)
		;
	for(; i<=col->sizeVal && (cmpDouble(col->coeff[i].a1*xi+col->coeff[i].b+col->coeff[i].a2*ys, 0.) < 0 || cmpDouble(col->coeff[i].a1*xs+col->coeff[i].b+col->coeff[i].a2*ys, 0.) < 0); i++) {
		if(col->coeff[i].a2 != 0) {
			indice[ind] = i;
			ystart[ind] = -(col->coeff[i].a1*xi+col->coeff[i].b)/col->coeff[i].a2;
			if(ystart[ind]>=yi && ystart[ind]<=ys) {
				xstart[ind] = xi;
			} else {
				if(col->coeff[i].a1 != 0) {
					if(ystart[ind]<yi) {
						xstart[ind] = -(col->coeff[i].a2*yi+col->coeff[i].b)/col->coeff[i].a1;
						if(xstart[ind]>=xi && xstart[ind]<=xs)
							ystart[ind] = yi;
						else {
							xstart[ind] = HUGE_VAL;
							ystart[ind] = HUGE_VAL;
						}
					} else { /* ystart[ind]>ys */
						xstart[ind] = -(col->coeff[i].a2*ys+col->coeff[i].b)/col->coeff[i].a1;
						if(xstart[ind]>=xi && xstart[ind]<=xs)
							ystart[ind] = ys;
						else {
							xstart[ind] = HUGE_VAL;
							ystart[ind] = HUGE_VAL;
						}
					}
				} else {
					xstart[ind] = HUGE_VAL;
					ystart[ind] = HUGE_VAL;
				}
			}
			yend[ind] = -(col->coeff[i].a1*xs+col->coeff[i].b)/col->coeff[i].a2;
			if(yend[ind]>=yi && yend[ind]<=ys) {
				xend[ind] = xs;
			} else {
				if(col->coeff[i].a1 != 0.) {
					if(yend[ind]<yi) {
						xend[ind] = -(col->coeff[i].a2*yi+col->coeff[i].b)/col->coeff[i].a1;
						if(xend[ind]>=xi && xend[ind]<=xs)
							yend[ind] = yi;
						else {
							xend[ind] = HUGE_VAL;
							yend[ind] = HUGE_VAL;
						}
					} else { /* yend[ind]>ys */
						xend[ind] = -(col->coeff[i].a2*ys+col->coeff[i].b)/col->coeff[i].a1;
						if(xend[ind]>=xi && xend[ind]<=xs)
							yend[ind] = ys;
						else {
							xend[ind] = HUGE_VAL;
							yend[ind] = HUGE_VAL;
						}
					}
				} else {
					xend[ind] = HUGE_VAL;
					yend[ind] = HUGE_VAL;
				}
			}
			if((xstart[ind] != HUGE_VAL && ystart[ind] != HUGE_VAL) && (xend[ind] != HUGE_VAL && yend[ind] != HUGE_VAL))
				ind++;
		} else {
			fprintf(stderr, "Execution error : null a2 coefficient not handled\n");
			printf("\nColumn (%.2lf, %.2lf) (%.2lf, %.2lf)\n", xi, xs, yi, ys);
			fprintColSpeDebug(stdout, col);
			exit(1);
		}
	}
	last = i;
	sum = 0.;
	double pxstart, pxend, pystart, pyend;
	pxstart = xi; pxend = xs; pystart = yi; pyend = yi;
	for(i=0; i<ind; i++) {
		if(xstart[i] < pxstart) {
			sum += sumGen(xi, xstart[i], ystart[i], ystart[i], pystart, pystart, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
			if(xend[i] < pxstart) {
				sum += sumGen(xstart[i], xend[i], ystart[i], yend[i], pystart, pystart, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(xend[i], pxstart, yend[i], yend[i], pystart, pystart, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(pxstart, pxend, yend[i], yend[i], pystart, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(pxend, xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
			} else {
				double tmp1 = ystart[i]+((pxstart-xstart[i])*(yend[i]-ystart[i]))/(xend[i]-xstart[i]);
				sum += sumGen(xstart[i], pxstart, ystart[i], tmp1, pystart, pystart, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				if(xend[i] < pxend) {
					double tmp2 = pystart+((xend[i]-pxstart)*(pyend-pystart))/(pxend-pxstart);
					sum += sumGen(pxstart, xend[i], tmp1, yend[i], pystart, tmp2, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(xend[i], pxend, yend[i], yend[i], tmp2, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(pxend, xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				} else {
					double tmp2 = ystart[i]+((pxend-xstart[i])*(yend[i]-ystart[i]))/(xend[i]-xstart[i]);
					sum += sumGen(pxstart, pxend, tmp1, tmp2, pystart, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(pxend, xend[i], tmp2, yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(xend[i], xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				}
			}
		} else {
			sum += sumGen(xi, pxstart, ystart[i], ystart[i], pystart, pystart, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
			if(pxend < xstart[i]) {
				sum += sumGen(pxstart, pxend, ystart[i], ystart[i], pystart, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(pxend, xstart[i], ystart[i], ystart[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(xstart[i], xend[i], ystart[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				sum += sumGen(xend[i], xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
			} else {
				double tmp1 = pystart+((xstart[i]-pxstart)*(pyend-pystart))/(pxend-pxstart);
				sum += sumGen(pxstart, xstart[i], ystart[i], ystart[i], pystart, tmp1, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				if(xend[i] < pxend) {
					double tmp2 = pystart+((xend[i]-pxstart)*(pyend-pystart))/(pxend-pxstart);
					sum += sumGen(xstart[i], xend[i], ystart[i], yend[i], tmp1, tmp2, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(xend[i], pxend, yend[i], yend[i], tmp2, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(pxend, xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				} else {
					double tmp2 = ystart[i]+((pxend-xstart[i])*(yend[i]-ystart[i]))/(xend[i]-xstart[i]);
					sum += sumGen(xstart[i], pxend, ystart[i], tmp2, tmp1, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(pxend, xend[i], tmp2, yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
					sum += sumGen(xend[i], xs, yend[i], yend[i], pyend, pyend, col->coeff[indice[i]], col->boundVal[indice[i]-1]);
				}
			}
		}
		pxstart = xstart[i]; pxend = xend[i]; pystart = ystart[i]; pyend = yend[i];
	}
	sum += sumGen(xi, pxstart, ys, ys, pystart, pystart, col->coeff[last], col->boundVal[last-1]);
	sum += sumGen(pxstart, pxend, ys, ys, pystart, pyend, col->coeff[last], col->boundVal[last-1]);
	sum += sumGen(pxend, xs, ys, ys, pyend, pyend, col->coeff[last], col->boundVal[last-1]);
	free((void*)xstart);
	free((void*)ystart);
	free((void*)xend);
	free((void*)yend);
	free((void*)indice);
	return sum;
}

void treatColumSpecial(FILE *f, TypeColumnSpecial *col, double xi, double xs, double yi, double ys) {
	int i;
	double *xstart, *ystart, *xend, *yend;
	size_t  *istart, *iend;
	istart = (size_t*) malloc((col->sizeVal+1)*sizeof(size_t));
	iend = (size_t*) malloc((col->sizeVal+1)*sizeof(size_t));
	xstart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	ystart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	xend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	yend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	for(i=0; i<=col->sizeVal; i++) {
		if(col->coeff[i].a2 != 0) {
			ystart[i] = -(col->coeff[i].a1*xi+col->coeff[i].b)/col->coeff[i].a2;
			yend[i] = -(col->coeff[i].a1*xs+col->coeff[i].b)/col->coeff[i].a2;
		} else {
			ystart[i] = HUGE_VAL;
			yend[i] = HUGE_VAL;
		}
	}
	istart = qsortTable(ystart, col->sizeVal+1, sizeof(double), compareDouble);
	iend = qsortTable(yend, col->sizeVal+1, sizeof(double), compareDouble);
	for(i=0; i<=col->sizeVal; i++)
		if(istart[i] != iend[i] && ystart[istart[i]] != ystart[iend[i]] && yend[istart[i]] != yend[iend[i]] && ystart[istart[i]]>=yi && ystart[istart[i]]<=ys) {
			fprintf(f, "XXXXX rank %d start %ld %.2le end %ld %.2le\n", i, (long)istart[i], ystart[istart[i]], (long)iend[i], yend[iend[i]]);
			fprintf(f, "start  %.2le %.2le end %.2le %.2le (%.2le, %.2le)\n", ystart[istart[i]], ystart[iend[i]], yend[istart[i]], yend[iend[i]], ystart[istart[i]]-ystart[iend[i]], yend[istart[i]]-yend[iend[i]]);
			fprintColumSpecialFig(f, col, xi, xs, yi, ys);
		}
	free((void*)xstart);
	free((void*)ystart);
	free((void*)xend);
	free((void*)yend);
	free((void*)istart);
	free((void*)iend);
}

#define SIZE_PICTURE 10
#define OUT_OF_PICTURE HUGE_VAL

void fprintColumSpecialFig(FILE *f, TypeColumnSpecial *col, double xi, double xs, double yi, double ys) {
	int i, ind, *indice;
	double scalex, scaley, *xstart, *ystart, *xend, *yend;
//	scale = SIZE_PICTURE/utils_MAX(xs-xi, ys-yi);
	scalex = SIZE_PICTURE/(xs-xi);
	scaley = SIZE_PICTURE/(ys-yi);
	indice = (int*) malloc((col->sizeVal+1)*sizeof(int));
	xstart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	ystart = (double*) malloc((col->sizeVal+1)*sizeof(double));
	xend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	yend = (double*) malloc((col->sizeVal+1)*sizeof(double));
	ind = 0;
	for(i=0; i<=col->sizeVal; i++) {
		if(col->coeff[i].a2 != 0 && col->coeff[i].a1 != 0) {
			indice[ind] = i;
			ystart[ind] = -(col->coeff[i].a1*xi+col->coeff[i].b)/col->coeff[i].a2;
			if(ystart[ind]>=yi && ystart[ind]<=ys) {
				xstart[ind] = xi;
			} else {
				xstart[ind] = -(col->coeff[i].a2*yi+col->coeff[i].b)/col->coeff[i].a1;
				if(xstart[ind]>=xi && xstart[ind]<=xs)
					ystart[ind] = yi;
				else {
					xstart[ind] = OUT_OF_PICTURE;
					ystart[ind] = OUT_OF_PICTURE;
				}
			}
			yend[ind] = -(col->coeff[i].a1*xs+col->coeff[i].b)/col->coeff[i].a2;
			if(yend[ind]>=yi && yend[ind]<=ys) {
				xend[ind] = xs;
			} else {
				xend[ind] = -(col->coeff[i].a2*ys+col->coeff[i].b)/col->coeff[i].a1;
				if(xend[ind]>=xi && xend[ind]<=xs)
					yend[ind] = ys;
				else {
					xend[ind] = OUT_OF_PICTURE;
					yend[ind] = OUT_OF_PICTURE;
				}
			}
			if(xstart[ind] != OUT_OF_PICTURE && ystart[ind] != OUT_OF_PICTURE && xend[ind] != OUT_OF_PICTURE && yend[ind] != OUT_OF_PICTURE)
				ind++;
		}
	}
	if(ind>0) {
		fprintf(f, "\\begin{tikzpicture}\n");
		fprintf(f, "\\node[text=black,anchor=east] at (%.2lf,%.2lf) {%.2lf};\n", 0., 0., yi);
		fprintf(f, "\\node[text=black,anchor=east] at (%.2lf,%.2lf) {%.2lf};\n", 0., (ys-yi)*scaley, ys);
		fprintf(f, "\\node[text=black,anchor=north] at (%.2lf,%.2lf) {%.2lf};\n", 0., 0., xi);
		fprintf(f, "\\node[text=black,anchor=north] at (%.2lf,%.2lf) {%.2lf};\n", (xs-xi)*scalex, 0., xs);
		fprintf(f, "\\draw (0, 0) -- (%.2lf, 0) -- (%.2lf, %.2lf) -- (0, %.2lf) -- cycle;\n", (xs-xi)*scalex, (xs-xi)*scalex, (ys-yi)*scaley, (ys-yi)*scaley);
		for(i=0; i<ind; i++)
			fprintf(f, "\\draw (%.2lf, %.2lf) -- (%.2lf, %.2lf);\n\\node[text=black,anchor=west] at (%.2lf,%.2lf) {%d};\n", (xstart[i]-xi)*scalex, (ystart[i]-yi)*scaley, (xend[i]-xi)*scalex, (yend[i]-yi)*scaley, (xend[i]-xi)*scalex, (yend[i]-yi)*scaley, indice[i]);
		fprintf(f, "\\end{tikzpicture}\n");
	}
	free((void*)xstart);
	free((void*)ystart);
	free((void*)xend);
	free((void*)yend);
	free((void*)indice);
}

void getColumnSpecial(TypeColumn *colA, TypeColumn *colB, TypeColumnSpecial *col) {
	int curVal_f, curVal_g;
	double minVal;

	col->coeff = (TypeCoeffSpecial*) malloc((colA->sizeVal+colB->sizeVal+1)*sizeof(TypeCoeffSpecial));
	col->boundVal = (double*) malloc((colA->sizeVal+colB->sizeVal)*sizeof(double));
	col->sizeVal  = 0;
	curVal_f = 0;
	curVal_g = 0;
	do {
		col->coeff[col->sizeVal].a1 = colA->coeff[curVal_f].a;
		col->coeff[col->sizeVal].a2 = colB->coeff[curVal_g].a;
		col->coeff[col->sizeVal].b = colA->coeff[curVal_f].b+colB->coeff[curVal_g].b;
		col->coeff[col->sizeVal].c1 = colA->coeff[curVal_f].c;
		col->coeff[col->sizeVal].c2 = colB->coeff[curVal_g].c;
		col->coeff[col->sizeVal].d = colA->coeff[curVal_f].d+colB->coeff[curVal_g].d;
		minVal = INFTY;
		if(curVal_f<colA->sizeVal && colA->boundVal[curVal_f]<minVal)
			minVal = colA->boundVal[curVal_f];
		if(curVal_g<colB->sizeVal && colB->boundVal[curVal_g]<minVal)
			minVal = colB->boundVal[curVal_g];
		if(curVal_f<colA->sizeVal && colA->boundVal[curVal_f]<=minVal)
			curVal_f++;
		if(curVal_g<colB->sizeVal && colB->boundVal[curVal_g]<=minVal)
			curVal_g++;
		if(minVal<INFTY) {
			col->boundVal[col->sizeVal] = minVal;
			col->sizeVal++;
		}
	} while(minVal<INFTY);
	col->coeff = (TypeCoeffSpecial*) realloc((void*)col->coeff, (col->sizeVal+1)*sizeof(TypeCoeffSpecial));
	if(col->sizeVal>0)
		col->boundVal = (double*) realloc((void*)col->boundVal, col->sizeVal*sizeof(double));
	else {
		free((void*)col->boundVal);
		col->boundVal = NULL;
	}
}

void sumPiecewise(TypePiecewise **f, int size, TypePiecewise *sum) {
	int k, sizePar;
	TypeFraction minPar;
	int *curPar, *curVal;

	sizePar = 0;
	for(k=0; k<size; k++)
		sizePar += f[k]->sizePar;
	curPar = (int*) malloc(size*sizeof(int));
	curVal = (int*) malloc(size*sizeof(int));
	sum->sizePar = 0;
	sum->boundPar = (TypeFraction*) malloc(sizePar*sizeof(TypeFraction));
	sum->col = (TypeColumn*) malloc((sizePar+1)*sizeof(TypeColumn));
	for(k=0; k<size; k++)
		curPar[k] = 0;	
	do {
		int j;
		for(j=0; j<size && f[j]->col[curPar[j]].sizeVal>=0; j++)
			;
		if(j<size) {
			sum->col[sum->sizePar].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
			sum->col[sum->sizePar].boundVal = (double*) malloc(sizeof(double));
			sum->col[sum->sizePar].sizeVal = -1;
			sum->col[sum->sizePar].boundVal = (double*) malloc(sizeof(double));
			sum->col[sum->sizePar].coeff[0].a = 0;
			sum->col[sum->sizePar].coeff[0].b = 0;
			sum->col[sum->sizePar].coeff[0].c = 0;
			sum->col[sum->sizePar].coeff[0].d = 0;
			for(k=0; k<size; k++) {
				if(f[k]->col[curPar[k]].sizeVal >= 0) {
					int l;
					for(l=0; l<f[k]->col[curPar[k]].sizeVal && f[j]->col[curPar[j]].boundVal[0]>f[k]->col[curPar[k]].boundVal[l]; l++)
						;
					sum->col[sum->sizePar].coeff[0].c += f[k]->col[curPar[k]].coeff[l].c+f[k]->col[curPar[k]].coeff[l].a*f[j]->col[curPar[j]].boundVal[0];
					sum->col[sum->sizePar].coeff[0].d += f[k]->col[curPar[k]].coeff[l].d+f[k]->col[curPar[k]].coeff[l].b*f[j]->col[curPar[j]].boundVal[0];
				} else {
					if(f[j]->col[curPar[j]].boundVal[0] != f[k]->col[curPar[k]].boundVal[0]) {
						fprintf(stderr, "Execution error: summing incompatible degenerate functions (%.2le\t%.2le)\n", f[j]->col[curPar[j]].boundVal[0], f[k]->col[curPar[k]].boundVal[0]);
						exit(1);
					}
					sum->col[sum->sizePar].coeff[0].c += f[k]->col[curPar[k]].coeff[0].c;
					sum->col[sum->sizePar].coeff[0].d += f[k]->col[curPar[k]].coeff[0].d;
				}
			}
		} else {
			double minVal;
			int sizeVal = 0;
			for(k=0; k<size; k++)
				sizeVal += f[k]->col[curPar[k]].sizeVal;
			sum->col[sum->sizePar].coeff = (TypeCoeff*) malloc((sizeVal+1)*sizeof(TypeCoeff));
			sum->col[sum->sizePar].boundVal = (double*) malloc(sizeVal*sizeof(double));
			sum->col[sum->sizePar].sizeVal = 0;
			for(k=0; k<size; k++)
				curVal[k] = 0;
			do {
				sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].a = 0;
				sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].b = 0;
				sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].c = 0;
				sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].d = 0;
				for(k=0; k<size; k++) {
						sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].a += f[k]->col[curPar[k]].coeff[curVal[k]].a;
						sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].b += f[k]->col[curPar[k]].coeff[curVal[k]].b;
						sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].c += f[k]->col[curPar[k]].coeff[curVal[k]].c;
						sum->col[sum->sizePar].coeff[sum->col[sum->sizePar].sizeVal].d += f[k]->col[curPar[k]].coeff[curVal[k]].d;
					}
					minVal = INFTY;
					for(k=0; k<size; k++)
						if(curVal[k]<f[k]->col[curPar[k]].sizeVal && f[k]->col[curPar[k]].boundVal[curVal[k]]<minVal)
							minVal = f[k]->col[curPar[k]].boundVal[curVal[k]];
					for(k=0; k<size; k++)
						if(curVal[k]<f[k]->col[curPar[k]].sizeVal && f[k]->col[curPar[k]].boundVal[curVal[k]]<=minVal)
							curVal[k]++;
					if(minVal<INFTY) {
						sum->col[sum->sizePar].boundVal[sum->col[sum->sizePar].sizeVal] = minVal;
						sum->col[sum->sizePar].sizeVal++;
					}
			} while(minVal<INFTY);
			sum->col[sum->sizePar].coeff = (TypeCoeff*) realloc((void*)sum->col[sum->sizePar].coeff, (sum->col[sum->sizePar].sizeVal+1)*sizeof(TypeCoeff));
			if(sum->col[sum->sizePar].sizeVal>0)
				sum->col[sum->sizePar].boundVal = (double*) realloc((void*)sum->col[sum->sizePar].boundVal, sum->col[sum->sizePar].sizeVal*sizeof(double));
			else {
				free((void*)sum->col[sum->sizePar].boundVal);
				sum->col[sum->sizePar].boundVal = NULL;
			}
		}
		minPar = fract(1,0);
		for(k=0; k<size ; k++)
			if(curPar[k]<f[k]->sizePar && cmpFract(f[k]->boundPar[curPar[k]],minPar)<0)
				minPar = f[k]->boundPar[curPar[k]];
		if(cmpFract(minPar, fract(1,0)) < 0) {
			sum->boundPar[sum->sizePar] = minPar;
			sum->sizePar++;
			for(k=0; k<size; k++)
				if(curPar[k]<f[k]->sizePar && cmpFract(f[k]->boundPar[curPar[k]],minPar)<=0)
					curPar[k]++;
		}
	} while(cmpFract(minPar, fract(1,0))<0);
	free(curPar);
	free(curVal);
}

TypeColumn *cloneColumn(TypeColumn *c) {
    int i;
    TypeColumn *clonedColumn = (TypeColumn *) malloc(sizeof(TypeColumn));
    clonedColumn->sizeVal  = c->sizeVal;
    clonedColumn->boundVal  = (double*) malloc(clonedColumn->sizeVal*sizeof(double));
    clonedColumn->coeff = (TypeCoeff*) malloc((clonedColumn->sizeVal+1)*sizeof(TypeCoeff));
    for(i=0; i<clonedColumn->sizeVal; i++) {
        clonedColumn->boundVal[i] = c->boundVal[i];
        clonedColumn->coeff[i].a  = c->coeff[i].a;
        clonedColumn->coeff[i].b  = c->coeff[i].b;
        clonedColumn->coeff[i].c  = c->coeff[i].c;
        clonedColumn->coeff[i].d  = c->coeff[i].d;
    }
    return clonedColumn;
}



void fprintFuncGnuplot(FILE *f, TypePiecewise *func) {
	int p, v;
	fprintf(f, "splot ");
	if(func->sizePar>0) {
		fprintf(f, "(x<(%lf/%lf))?(", func->boundPar[0].num, func->boundPar[0].den);
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[0], func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[v], func->col[0].coeff[v].a, func->col[0].coeff[v].b, func->col[0].coeff[v].c, func->col[0].coeff[v].d);
			fprintf(f, ":%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[func->col[0].sizeVal].a, func->col[0].coeff[func->col[0].sizeVal].b, func->col[0].coeff[func->col[0].sizeVal].c, func->col[0].coeff[func->col[0].sizeVal].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
		}
		fprintf(f, ")");
		for(p=1; p<func->sizePar; p++) {
			fprintf(f, ":((x<(%lf/%lf))?", func->boundPar[p].num, func->boundPar[p].den);
			if(func->col[0].sizeVal>0) {
				fprintf(f, "(y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[0], func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
				for(v=1; v<func->col[p].sizeVal; v++)
					fprintf(f, ":((y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[v], func->col[p].coeff[v].a, func->col[p].coeff[v].b, func->col[p].coeff[v].c, func->col[p].coeff[v].d);
				fprintf(f, ":%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[func->col[p].sizeVal].a, func->col[p].coeff[func->col[p].sizeVal].b, func->col[p].coeff[func->col[p].sizeVal].c, func->col[p].coeff[func->col[p].sizeVal].d);
				for(v=1; v<func->col[p].sizeVal; v++)
					fprintf(f, ")");
			} else {
				fprintf(f, "%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
			}
		}
		fprintf(f, ":(");
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[0], func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
			for(v=1; v<func->col[p].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].boundVal[v], func->col[p].coeff[v].a, func->col[p].coeff[v].b, func->col[p].coeff[v].c, func->col[p].coeff[v].d);
			fprintf(f, ":%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[func->col[p].sizeVal].a, func->col[p].coeff[func->col[p].sizeVal].b, func->col[p].coeff[func->col[p].sizeVal].c, func->col[p].coeff[func->col[p].sizeVal].d);
			for(v=1; v<func->col[p].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[p].coeff[0].a, func->col[p].coeff[0].b, func->col[p].coeff[0].c, func->col[p].coeff[0].d);
		}
		fprintf(f, ")");
		for(p=1; p<func->sizePar; p++)
			fprintf(f, ")");
	} else {
		if(func->col[0].sizeVal>0) {
			fprintf(f, "(y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[0], func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ":((y<(%lf))?%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].boundVal[v], func->col[0].coeff[v].a, func->col[0].coeff[v].b, func->col[0].coeff[v].c, func->col[0].coeff[v].d);
			fprintf(f, ":%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[func->col[0].sizeVal].a, func->col[0].coeff[func->col[0].sizeVal].b, func->col[0].coeff[func->col[0].sizeVal].c, func->col[0].coeff[func->col[0].sizeVal].d);
			for(v=1; v<func->col[0].sizeVal; v++)
				fprintf(f, ")");
		} else {
			fprintf(f, "%lf*x*y+%lf*y+(%lf)*x+(%lf)", func->col[0].coeff[0].a, func->col[0].coeff[0].b, func->col[0].coeff[0].c, func->col[0].coeff[0].d);
		}
	}
}

void fprintfl(FILE *f, double a, double b) {
	if(a>0.0)
		fprintf(f, "%.1lf\\gamma", a);
	if(a>0.0 && b>0.0)
		fprintf(f, "+");
	if(b>0.0)
		fprintf(f, "%.1lf", b);
}

void fprintfrac(FILE *f, int n, int d) {
	if(d==0) {
		fprintf(f, " \\infty ");
	} else {
		if(d==1 || n ==0) {
			fprintf(f, " %d ", n);
		} else {
			fprintf(f, " \\frac{%d}{%d} ", n, d);
		}
	}
}

void fprintfs(FILE *f, TypeCoeff c) {
	if(c.a>0 || c.b>0) {
		if(c.a>0 && c.b>0)
			fprintf(f, "(");
		if(c.a>0) {
			if(c.a!=1)
				fprintf(f, "%.2lf", c.a);
			fprintf(f, "\\gamma");
		}
		if(c.a>0 && c.b>0)
			fprintf(f, "+");
		if(c.b>0)
			fprintf(f, "%.2lf", c.b);
		if(c.a>0 && c.b>0)
			fprintf(f, ")");
		fprintf(f, " x");
	}
	if(c.c>0. || c.d>0.) {
			fprintf(f, "+");
		if(c.c>0.0) {
			if(c.c!=1.)
				fprintf(f, "%.2lf", c.c);
			fprintf(f, "\\gamma");
		}
		if(c.c>0.0 && c.d>0.0)
			fprintf(f, "+");
		if(c.d!=0.0)
			fprintf(f, "%.2lf", c.d);
	}
}

void fprintFunc(FILE *f, TypePiecewise *func) {
	int p, v;
	fprintf(f, "$$\\begin{array}{");
	for(p=0; p<func->sizePar; p++)
		fprintf(f, "cc");
	fprintf(f, "c");
	fprintf(f, "}\n");
	for(p=0; p<func->sizePar; p++) {
		fprintf(f, "& ");
		fprintfrac(f, func->boundPar[p].num, func->boundPar[p].den);
		fprintf(f, "& ");
	}
	fprintf(f, "\\\\ \n");
	fprintf(f, "\\begin{array}{rc}\n");
	fprintf(f, "& ");
	fprintfs(f, func->col[0].coeff[0]);
	fprintf(f, "\\\\\n ");
	for(v=0; v<func->col[0].sizeVal; v++) {
		fprintf(f, "%.2lf & \\raisedrule[0.15em]{0.5pt}\\\\\n", func->col[0].boundVal[v]);
		fprintf(f, "& ");
		fprintfs(f, func->col[0].coeff[v+1]);
		fprintf(f, "\\\\\n ");
	}
	fprintf(f, "\\end{array} \n");
	for(p=1; p<=func->sizePar; p++) {
		fprintf(f, "& | &");
		fprintf(f, "\\begin{array}{rc}\n");
		fprintf(f, "& ");
		fprintfs(f, func->col[p].coeff[0]);
		fprintf(f, "\\\\\n ");
		for(v=0; v<func->col[p].sizeVal; v++) {
			fprintf(f, "%.2lf & \\raisedrule[0.15em]{0.5pt}\\\\\n", func->col[p].boundVal[v]);
		fprintf(f, "& ");
		fprintfs(f, func->col[p].coeff[v+1]);
		fprintf(f, "\\\\\n ");
		}
		fprintf(f, "\\end{array} \n");
	}
	fprintf(f, "\\end{array}$$");
}

void fprintptDebug(FILE *f, TypeCoeff c, double x) {
	fprintf(f, "(%.2lf, %.2lf)", -c.a*x+c.c, c.b*x+c.d);
}
void fprintfCoeff(FILE *f, TypeCoeff coeff) {
    if (coeff.a > 0)
        fprintf(f, "%.2lf g.x", coeff.a);
    else if (coeff.a == 0)
        fprintf(f, " 0 g.x");
    fprintf(f, " + ");
    fprintf(f, "%.2lf x", coeff.b);
    fprintf(f, " + ");
    fprintf(f, "%.2lf g", coeff.c);
    fprintf(f, " + ");
    fprintf(f, "%.2lf\n", coeff.d);
    
}

void fprintfsDebug(FILE *f, TypeCoeff c) {
	if(c.a!=0 || c.b!=0) {
		if(c.a!=0 && c.b!=0)
			fprintf(f, "(");
		if(c.a!=0) {
			if(c.a!=1)
				fprintf(f, "%.2lf", c.a);
			fprintf(f, "g");
		}
		if(c.a!=0 && c.b>0)
			fprintf(f, "+");
		if(c.b!=0 && (c.b != 1 || c.a != 0))
			fprintf(f, "%.2lf", c.b);
		if(c.a!=0 && c.b!=0)
			fprintf(f, ")");
		fprintf(f, "x");
	}
	if(c.c!=0. || c.d!=0.) {
		if(c.c>0.0) {
			fprintf(f, "+");
			if(c.c!=1.)
				fprintf(f, "%.2lf", c.c);
			fprintf(f, "g");
		} else {
			if(c.c<0.0) {
				if(c.c!=-1.) {
					fprintf(f, "%.2lfg", c.c);
				} else {
					fprintf(f, "-g");
				}
			}
		}
		if(c.d>0.0) {
			fprintf(f, "+");
			fprintf(f, "%.2lf", c.d);
		} else {
			if(c.d<0.0) {
				fprintf(f, "%.2lf", c.d);
			}
		}
	}
}

void fprintfsDebugX(FILE *f, TypeCoeff c) {
	if(c.a!=0 || c.b!=0) {
		if(c.a!=0) {
			if(c.a!=1)
				fprintf(f, "%.2lf", c.a);
			fprintf(f, "x");
		}
		if(c.a!=0 && c.b>0)
			fprintf(f, "+");
		if(c.b!=0 && (c.b != 1 || c.a != 0))
			fprintf(f, "%.2lf", c.b);
	}
}

void fprintfsSpecialDebugX(FILE *f, TypeCoeffSpecial c) {
	if(c.a1!=0 || c.a2!=0 ||c.b!=0) {
		if(c.a1!=0) {
			if(c.a1!=1)
				fprintf(f, "%.2lf", c.a1);
			fprintf(f, "x");
		}
		if(c.a1!=0 && c.a2>0)
			fprintf(f, "+");
		if(c.a2!=0) {
			if(c.a2 != 1)
				fprintf(f, "%.2lf", c.a2);
			fprintf(f, "y");
		}
		if((c.a1!=0 || c.a2!=0) && c.b>0)
			fprintf(f, "+");
		if(c.b!=0)
			fprintf(f, "%.2lf", c.b);

	}
}

void fprintfsSpeDebug(FILE *f, TypeCoeffSpecial c) {
	if(c.a1!=0 || c.a2!=0 ||c.b!=0) {
		fprintf(f, "(");
		if(c.a1!=0) {
			if(c.a1!=1)
				fprintf(f, "%.2lf", c.a1);
			fprintf(f, "a");
		}
		if(c.a1!=0 && c.a2>0)
			fprintf(f, "+");
		if(c.a2!=0) {
			if(c.a2!=1)
				fprintf(f, "%.2lf", c.a2);
			fprintf(f, "a'");
		}
		if((c.a1!=0 || c.a2!=0) && c.b>0)
			fprintf(f, "+");
		if(c.b!=0 && (c.b != 1 || c.a1 != 0 || c.a2 != 0))
			fprintf(f, "%.2lf", c.b);
		fprintf(f, ")");
		fprintf(f, "x");
	}
	if(c.c1!=0. || c.c2!=0. || c.d!=0.) {
		if(c.c1>0.0) {
			fprintf(f, "+");
			if(c.c1!=1.)
				fprintf(f, "%.2lf", c.c1);
			fprintf(f, "a");
		} else {
			if(c.c1<0.0) {
				if(c.c1!=-1.) {
					fprintf(f, "%.2lfg", c.c1);
				} else {
					fprintf(f, "-a'");
				}
			}
		}
		if(c.c2>0.0) {
			fprintf(f, "+");
			if(c.c2!=1.)
				fprintf(f, "%.2lf", c.c2);
			fprintf(f, "a'");
		} else {
			if(c.c2<0.0) {
				if(c.c2!=-1.) {
					fprintf(f, "%.2lfg", c.c2);
				} else {
					fprintf(f, "-a'");
				}
			}
		}
		if(c.d>0.0) {
			fprintf(f, "+");
			fprintf(f, "%.2lf", c.d);
		} else {
			if(c.d<0.0) {
				fprintf(f, "%.2lf", c.d);
			}
		}
	}
}

void fprintColSpeDebug(FILE *f, TypeColumnSpecial *c) {
	int v;
	fprintfsSpecialDebugX(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n:::::: %.2lf ::::::", c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsSpeDebug(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}

void fprintColDebugX(FILE *f, TypeColumn *c) {
	int v;
	fprintfsDebugX(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n:::::: %.2lf ::::::", c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsDebugX(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}

void fprintColSpecialDebugX(FILE *f, TypeColumnSpecial *c) {
	int v;
	fprintfsSpecialDebugX(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n:::::: %.2lf ::::::", c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsSpecialDebugX(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}

void fprintColDebug(FILE *f, TypeColumn *c) {
	int v;
	fprintfsDebug(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n:::::: %.2lf ::::::", c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsDebug(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}

void fprintColSpecialDebug(FILE *f, TypeColumnSpecial *c) {
	int v;
	fprintfsSpecialDebugX(f, c->coeff[0]);
	for(v=1; v<=c->sizeVal; v++) {
		fprintf(f, "\n:::::: %.2lf ::::::", c->boundVal[v-1]);
		fprintf(f,"\n");
		fprintfsSpecialDebugX(f, c->coeff[v]);
	}
	fprintf(f,"\n");
}


void fprintFuncDebug(FILE *f, TypePiecewise *func) {
	int p;
	fprintColDebug(f, &(func->col[0]));
	for(p=1; p<=func->sizePar; p++) {
		fprintf(f, "\n*********%.2le*********\n\n", func->boundPar[p-1].num/func->boundPar[p-1].den);
		fprintColDebug(f,  &(func->col[p]));
	}
}
