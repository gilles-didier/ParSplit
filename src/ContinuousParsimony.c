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
#include "ContinuousParsimony.h"
#include "Floating.h"

#define RADIUS 0.5
#define XOFFSET 0.0
#define YOFFSET 0.0
#define XOFF 0.2
#define YOFF -0.4
#define SIZEBUF 10
#define FZERO {0,1}
#define FINFTY {1,0}
#define FNINFTY {-1,0}
#define EPSI 1E-10
#define SIZE_BUFFER_CHAR 2000
#define WIDTH "stringWidth"
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define FALSE 0
#define TRUE 1

static void getColumnF(TypeColumn *colF, TypeColumn *colF_, int jp, int jm, double timeFact);
static void getColumnH(TypeColumn *colH, TypeColumn *colH_, int jp, int jm, double timeFact);
static void fillF(int n, TypePiecewise *F, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi);
static void fillH(int n, TypePiecewise *H, TypePiecewise *H_, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi);
static void FtoF_(TypePiecewise *f, TypePiecewise *h, double timeFact, double val);
static void HtoH_(TypePiecewise *f, TypePiecewise *h, double timeFact, double val);
static TypeResult *fillResult(TypePiecewise *F, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi);
static void fillResultRec(int n, TypeResult *res, TypePiecewise *funcFH, TypeTree *tree, TypeFunctTime phi);
void transformPieceWise(TypePiecewise *F);
static int sgnLin(double a, double b, TypeFraction f);
int compareResultFunc(TypeResultFunc *r, TypeResultFunc *s);


int sgnLin(double a, double b, TypeFraction f) {
	if(a == 0.) {
		if(b>0)
			return 1;
		if(b<0)
			return -1;
		return 0;
	}
	if(f.den>0) {
		if(a*f.num+b*f.den>0)
			return 1;
		if(a*f.num+b*f.den<0)
			return -1;
		return 0;
	}
	if(f.den<0) {
		if(a*f.num+b*f.den>0)
			return -1;
		if(a*f.num+b*f.den<0)
			return 1;
		return 0;
	}
	if(a*f.num>0)
		return 1;
	if(a*f.num<0)
		return -1;
	return 1;
}


void transformPieceWise(TypePiecewise *F) {
	int i, j;
	for(i=0; i<F->sizePar; i++)
		F->boundPar[i].den = F->boundPar[i].den+F->boundPar[i].num;
	for(i=0; i<=F->sizePar; i++)
		for(j=0; j<=F->col[i].sizeVal; j++) {
			F->col[i].coeff[j].a -= F->col[i].coeff[j].b;
			F->col[i].coeff[j].c -= F->col[i].coeff[j].d;
		}
}

TypeResult *getResult(TypeTree *tree, TypeFunctTime phi) {
	TypePiecewise *F, *F_;
	TypeResult *res;
	int n;
	F = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
	F_ = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
	fillF(tree->root, F, F_, tree, phi);
	res = fillResult(F, F_,tree, phi);
	for(n=0; n<tree->size; n++) {
		freePiecewise(F[n]);
		if(n != tree->root)
			freePiecewise(F_[n]);
	}
	free((void*)F);
	free((void*)F_);
	return res;
}

void fillNormFH(TypeTree *tree, TypeFunctTime phi, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_) {
	int n;
	fillF(tree->root, F, F_, tree, phi);
	H_[tree->root].sizePar = 0;
	H_[tree->root].boundPar = NULL;
	H_[tree->root].col = (TypeColumn*) malloc(sizeof(TypeColumn));
	H_[tree->root].col[0].sizeVal = 0;
	H_[tree->root].col[0].boundVal = NULL;
	H_[tree->root].col[0].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
	H_[tree->root].col[0].coeff[0].a = 0;
	H_[tree->root].col[0].coeff[0].b = 0;
	H_[tree->root].col[0].coeff[0].c = 0.;
	H_[tree->root].col[0].coeff[0].d = 0.;
	fillH(tree->root, H, H_,F_, tree, phi);
	for(n=0; n<tree->size; n++) {
		transformPieceWise(&(F[n]));
		if(n != tree->root) {
			transformPieceWise(&(F_[n]));
			transformPieceWise(&(H[n]));
		}
		if(((double*)tree->info)[n] == UNKNOWN) {
			transformPieceWise(&(H_[n]));
		}
	}
}

void getAllCut(TypeTree *tree, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_, int **tnmin, int **ttype, double **tmin) {
	int n, size = 0, *xnmin, *xtype;
	size_t *index;
	double min, *xmin;
	FILE *fo;
	if(!(fo = fopen("report.txt", "w"))) {
		fprintf(stderr, "Error while writing 'report.txt'\n");
		exit(1);
	}
	min = intPiecewise(&(F[tree->root]));
	xnmin = (int*) malloc(tree->size*sizeof(int));
	xtype = (int*) malloc(tree->size*sizeof(int));
	xmin = (double*) malloc(tree->size*sizeof(double));;
	for(n=0; n<tree->size; n++) {
		fprintf(fo, "%d", n);
		if(tree->name != NULL && tree->name[n] != NULL)
			fprintf(fo, "\t%s", tree->name[n]);
		else
			fprintf(fo, "\t");
		if(n != tree->root && ((double*)tree->info)[n] == UNKNOWN) {
			int i, j;
			double sum = 0.;
			for(i=0; i<=H_[n].sizePar; i++)
				for(j=0; j<=F[n].sizePar; j++) {
					TypeColumnSpecial col;
					double xi, xs, yi, ys;
					if(i>0)
						xi = H_[n].boundPar[i-1].num/H_[n].boundPar[i-1].den;
					else
						xi = 0;
					if(i<H_[n].sizePar)
						xs = H_[n].boundPar[i].num/H_[n].boundPar[i].den;
					else
						xs = 1;
					if(j>0)
						yi = F[n].boundPar[j-1].num/F[n].boundPar[j-1].den;
					else
						yi = 0;
					if(j<F[n].sizePar)
						ys = F[n].boundPar[j].num/F[n].boundPar[j].den;
					else
						ys = 1;
					getColumnSpecial(&(H_[n].col[i]), &(F[n].col[j]), &col);
					sum += sumColumnSpecial(&col, xi, xs, yi, ys);
					free((void*)col.coeff);
					free((void*)col.boundVal);
				}
			fprintf(fo, "\t%.4lf", sum);
			if(cmpDouble(sum, min) == -1) {
				(xnmin)[size] = n;
				(xmin)[size] = sum;
				(xtype)[size] = 0;
				size++;
			}
		} else {
			if(n == tree->root && ((double*)tree->info)[n] == UNKNOWN)
				fprintf(fo, "\t%.4lf", intPiecewise(&(F[tree->root])));
			else
				fprintf(fo, "\t");
		}
		if(n != tree->root) {
			int i, j;
			double sum = 0.;
			for(i=0; i<=H[n].sizePar; i++)
				for(j=0; j<=F_[n].sizePar; j++) {
					TypeColumnSpecial col;
					double xi, xs, yi, ys;
					if(i>0)
						xi = H[n].boundPar[i-1].num/H[n].boundPar[i-1].den;
					else
						xi = 0;
					if(i<H[n].sizePar)
						xs = H[n].boundPar[i].num/H[n].boundPar[i].den;
					else
						xs = 1;
					if(j>0)
						yi = F_[n].boundPar[j-1].num/F_[n].boundPar[j-1].den;
					else
						yi = 0;
					if(j<F_[n].sizePar)
						ys = F_[n].boundPar[j].num/F_[n].boundPar[j].den;
					else
						ys = 1;
					getColumnSpecial(&(H[n].col[i]), &(F_[n].col[j]), &col);
					sum += sumColumnSpecial(&col, xi, xs, yi, ys);
					free((void*)col.coeff);
					free((void*)col.boundVal);
				}
			fprintf(fo, "\t%.4lf\n", sum);
			if(cmpDouble(sum, min) == -1) {
				(xnmin)[size] = n;
				(xmin)[size] = sum;
				(xtype)[size] = 1;
				size++;
			}
		} else
			fprintf(fo, "\t\n");

	}
	fclose(fo);
	index = qsortTable(xmin, size, sizeof(double), compareDouble);
	*tnmin = (int*) malloc((size+1)*sizeof(int));
	*ttype = (int*) malloc(size*sizeof(int));
	*tmin = (double*) malloc(size*sizeof(double));;
	for(n=0; n<size; n++) {
		(*tnmin)[n] = xnmin[index[n]];
		(*tmin)[n] = xmin[index[n]];
		(*ttype)[n] = xtype[index[n]];
	}
	(*tnmin)[size] = NOSUCH;
	free((void*)xnmin);
	free((void*)xmin);
	free((void*)xtype);
	free((void*)index);
}

void getMinCut(TypeTree *tree, TypePiecewise *F, TypePiecewise *F_, TypePiecewise *H, TypePiecewise *H_, int *nmin, int *type, double *min) {
	int n;
	*min = intPiecewise(&(F[tree->root]));
	*nmin = NOSUCH;
	*type = 0;
	FILE *fo;
	if(!(fo = fopen("report.txt", "w"))) {
		fprintf(stderr, "Error while writing 'report.txt'\n");
		exit(1);
	}
	for(n=0; n<tree->size; n++) {
		fprintf(fo, "%d", n);
		if(tree->name != NULL && tree->name[n] != NULL)
			fprintf(fo, "\t%s", tree->name[n]);
		else
			fprintf(fo, "\t");
		if(n != tree->root && ((double*)tree->info)[n] == UNKNOWN) {
			int i, j;
			double sum = 0.;
			for(i=0; i<=H_[n].sizePar; i++)
				for(j=0; j<=F[n].sizePar; j++) {
					TypeColumnSpecial col;
					double xi, xs, yi, ys;
					if(i>0)
						xi = H_[n].boundPar[i-1].num/H_[n].boundPar[i-1].den;
					else
						xi = 0;
					if(i<H_[n].sizePar)
						xs = H_[n].boundPar[i].num/H_[n].boundPar[i].den;
					else
						xs = 1;
					if(j>0)
						yi = F[n].boundPar[j-1].num/F[n].boundPar[j-1].den;
					else
						yi = 0;
					if(j<F[n].sizePar)
						ys = F[n].boundPar[j].num/F[n].boundPar[j].den;
					else
						ys = 1;
					getColumnSpecial(&(H_[n].col[i]), &(F[n].col[j]), &col);
					sum += sumColumnSpecial(&col, xi, xs, yi, ys);
					free((void*)col.coeff);
					free((void*)col.boundVal);
				}
			fprintf(fo, "\t%.4lf", sum);
			if(sum<*min) {
				*nmin = n;
				*min = sum;
				*type = 0;
			}
		} else {
			if(n == tree->root && ((double*)tree->info)[n] == UNKNOWN)
				fprintf(fo, "\t%.4lf", intPiecewise(&(F[tree->root])));
			else
				fprintf(fo, "\t");
		}
		if(n != tree->root) {
			int i, j;
			double sum = 0.;
			for(i=0; i<=H[n].sizePar; i++)
				for(j=0; j<=F_[n].sizePar; j++) {
					TypeColumnSpecial col;
					double xi, xs, yi, ys;
					if(i>0)
						xi = H[n].boundPar[i-1].num/H[n].boundPar[i-1].den;
					else
						xi = 0;
					if(i<H[n].sizePar)
						xs = H[n].boundPar[i].num/H[n].boundPar[i].den;
					else
						xs = 1;
					if(j>0)
						yi = F_[n].boundPar[j-1].num/F_[n].boundPar[j-1].den;
					else
						yi = 0;
					if(j<F_[n].sizePar)
						ys = F_[n].boundPar[j].num/F_[n].boundPar[j].den;
					else
						ys = 1;
					getColumnSpecial(&(H[n].col[i]), &(F_[n].col[j]), &col);
					sum += sumColumnSpecial(&col, xi, xs, yi, ys);
					free((void*)col.coeff);
					free((void*)col.boundVal);
				}
			fprintf(fo, "\t%.4lf\n", sum);
			if(sum<*min) {
				*nmin = n;
				*min = sum;
				*type = 1;
			}
		} else
			fprintf(fo, "\t\n");
	}
	fclose(fo);
}


void fillF(int n, TypePiecewise *F, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi) {
	if(tree->node[n].child != NOSUCH) {
		int c, nchild, k, sizePar;
		TypePiecewise **list;
		nchild = 0;
		for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
			fillF(c, F, F_, tree, phi);
			nchild++;
		}
		list = (TypePiecewise**) malloc(nchild*sizeof(TypePiecewise*));
		k = 0;
		sizePar = 0;
		for(c=tree->node[n].child; c!= NOSUCH; c=tree->node[c].sibling) {
			FtoF_(&(F[c]), &(F_[c]), phi(tree->time[c]), ((double*)tree->info)[n]);
			sizePar += F_[c].sizePar;
			list[k++] = &(F_[c]);
		}
		sumPiecewise(list, nchild, &(F[n]));
		free((void*)list);
	} else {
		if(isUnknown(n, tree)) {
			F[n].sizePar = 0;
			F[n].boundPar = NULL;
			F[n].col = (TypeColumn*) malloc(sizeof(TypeColumn));
			F[n].col[0].sizeVal = 0;
			F[n].col[0].boundVal = NULL;
			F[n].col[0].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
			F[n].col[0].coeff[0].a = 0.;
			F[n].col[0].coeff[0].b = 0.;
			F[n].col[0].coeff[0].c = 0.;
			F[n].col[0].coeff[0].d = 0.;			
		} else {
			F[n].sizePar = 0;
			F[n].boundPar = NULL;
			F[n].col = (TypeColumn*) malloc(sizeof(TypeColumn));
			F[n].col[0].sizeVal = -1;
			F[n].col[0].boundVal = (double*) malloc(sizeof(double));
			F[n].col[0].boundVal[0] = ((double*)tree->info)[n];
			F[n].col[0].coeff = (TypeCoeff*) malloc(sizeof(TypeCoeff));
			F[n].col[0].coeff[0].a = 0.;
			F[n].col[0].coeff[0].b = 0.;
			F[n].col[0].coeff[0].c = 0.;
			F[n].col[0].coeff[0].d = 0.;			
		}
	}
}

void FtoF_(TypePiecewise *F, TypePiecewise *F_, double timeFact, double val) {
	int sizeParBuf, i;
	TypeFraction boundp, boundm;
	if(val == UNKNOWN) {
		sizeParBuf = F->sizePar;
		for(i=0; i<=F->sizePar; i++)
			if(F->col[i].sizeVal>0)
				sizeParBuf += 2*F->col[i].sizeVal;
		F_->sizePar = 0;
		F_->boundPar = (TypeFraction*) malloc(sizeParBuf*sizeof(TypeFraction));
		F_->col = (TypeColumn*) malloc((sizeParBuf+1)*sizeof(TypeColumn));
		for(i=0; i<=F->sizePar; i++) {
			TypeFraction finf, fsup;
			int jp, jm, sjp, sjm, ejp, ejm;
			if(i>0)
				finf = F->boundPar[i-1];
			else 
				finf = fract(0,1);			
			if(i<F->sizePar)
				fsup = F->boundPar[i];
			else
				fsup = fract(1,0);
			if(F->col[i].sizeVal >= 0) {
				for(sjp=0; sjp<=F->col[i].sizeVal && (sgnLin(F->col[i].coeff[sjp].a+timeFact, F->col[i].coeff[sjp].b,finf)<=0 && sgnLin(F->col[i].coeff[sjp].a+timeFact, F->col[i].coeff[sjp].b,fsup)<=0); sjp++)
					;
				for(ejp=sjp; ejp<=F->col[i].sizeVal && sgnLin(F->col[i].coeff[ejp].a+timeFact, F->col[i].coeff[ejp].b,finf)*sgnLin(F->col[i].coeff[ejp].a+timeFact, F->col[i].coeff[ejp].b,fsup)<0; ejp++)
					;
				for(ejm=F->col[i].sizeVal+1; ejm>0 && (sgnLin(F->col[i].coeff[ejm-1].a, F->col[i].coeff[ejm-1].b-timeFact,finf)>=0 && sgnLin(F->col[i].coeff[ejm-1].a, F->col[i].coeff[ejm-1].b-timeFact,fsup)>=0); ejm--)
					;
				for(sjm=ejm; sjm>0 && sgnLin(F->col[i].coeff[sjm-1].a, F->col[i].coeff[sjm-1].b-timeFact,finf)*sgnLin(F->col[i].coeff[sjm-1].a, F->col[i].coeff[sjm-1].b-timeFact,fsup)<0; sjm--)
					;
				jp = sjp;
				jm = sjm;
				getColumnF(&(F->col[i]), &(F_->col[F_->sizePar]), jp, jm, timeFact);
				while(jp<ejp || jm<ejm) {
					if(jp<ejp)
						boundp = fract(F->col[i].coeff[jp].b, -F->col[i].coeff[jp].a-timeFact);
					else
						boundp = fsup;
					if(jm<ejm)
						boundm = fract(F->col[i].coeff[jm].b-timeFact, -F->col[i].coeff[jm].a);
					else
						boundm = fsup;
					if(F_->sizePar >= sizeParBuf) {
						sizeParBuf += F->sizePar+2;
						F_->boundPar = (TypeFraction*) realloc((void*) F_->boundPar, sizeParBuf*sizeof(TypeFraction));
						F_->col = (TypeColumn*) realloc((void*) F_->col, (sizeParBuf+1)*sizeof(TypeColumn));
					}
					if(cmpFract(boundp, boundm)<=0)
						F_->boundPar[F_->sizePar] = boundp;
					else
						F_->boundPar[F_->sizePar] = boundm;
					if(cmpFract(boundp, F_->boundPar[F_->sizePar])<=0)
						jp++;
					if(cmpFract(boundm, F_->boundPar[F_->sizePar])<=0)
						jm++;
					F_->sizePar++;
					getColumnF(&(F->col[i]), &(F_->col[F_->sizePar]), jp, jm, timeFact);
				}
				if(i<F->sizePar) {
					F_->boundPar[F_->sizePar] = F->boundPar[i];
					F_->sizePar++; 
				}
			} else {
				F_->sizePar = F->sizePar;
				for(i=0; i<=F->sizePar; i++) {
					F_->col[i].sizeVal = 1;
					F_->col[i].boundVal = (double*) malloc(sizeof(double));
					F_->col[i].boundVal[0] = F->col[i].boundVal[0];
					F_->col[i].coeff = (TypeCoeff*) malloc(2*sizeof(TypeCoeff));
					F_->col[i].coeff[0].a = -timeFact;
					F_->col[i].coeff[0].b = 0;
					F_->col[i].coeff[0].c = F->col[i].boundVal[0]*(timeFact+F->col[i].coeff[0].a)+F->col[i].coeff[0].c;
					F_->col[i].coeff[0].d = F->col[i].boundVal[0]*F->col[i].coeff[0].b+F->col[i].coeff[0].d;
					F_->col[i].coeff[1].a = 0;
					F_->col[i].coeff[1].b = timeFact;
					F_->col[i].coeff[1].c = F->col[i].boundVal[0]*F->col[i].coeff[0].a+F->col[i].coeff[0].c;
					F_->col[i].coeff[1].d = F->col[i].boundVal[0]*(-timeFact+F->col[i].coeff[0].b)+F->col[i].coeff[0].d;
					if(i<F->sizePar)
						F_->boundPar[i] = F->boundPar[i];
				}
			}
		}
		F_->boundPar = (TypeFraction*) realloc((void*) F_->boundPar, F_->sizePar*sizeof(TypeFraction));
		F_->col = (TypeColumn*) realloc((void*) F_->col, (F_->sizePar+1)*sizeof(TypeColumn));
	} else {
		fprintf(stderr, "Critical error: known internal nodes are not handled in this version\n");
		exit(1);
	}
}

void getColumnF(TypeColumn *colF, TypeColumn *colF_, int jp, int jm, double timeFact) {
	int j;
	colF_->sizeVal = jm-jp-1;
	if(jp > 0)
		colF_->sizeVal++;
	if(jm <= colF->sizeVal)
		colF_->sizeVal++;
	colF_->boundVal = (double*) malloc(colF_->sizeVal*sizeof(double));
	colF_->coeff = (TypeCoeff*) malloc((colF_->sizeVal+1)*sizeof(TypeCoeff));
	colF_->sizeVal = 0;
	if(jp > 0) {
		colF_->coeff[colF_->sizeVal].a = -timeFact;
		colF_->coeff[colF_->sizeVal].b = 0.;
		colF_->coeff[colF_->sizeVal].c = colF->boundVal[jp-1]*(timeFact+colF->coeff[jp].a)+colF->coeff[jp].c;
		colF_->coeff[colF_->sizeVal].d = colF->boundVal[jp-1]*colF->coeff[jp].b+colF->coeff[jp].d;
		colF_->boundVal[colF_->sizeVal++] = colF->boundVal[jp-1];
	}
	if(jm <= colF->sizeVal) {
		for(j=jp; j<jm; j++) {
			colF_->coeff[colF_->sizeVal] = colF->coeff[j];
			colF_->boundVal[colF_->sizeVal++] = colF->boundVal[j];
		}
		colF_->coeff[colF_->sizeVal].a = 0.;
		colF_->coeff[colF_->sizeVal].b = timeFact;
		colF_->coeff[colF_->sizeVal].c = colF->boundVal[jm-1]*colF->coeff[jm].a+colF->coeff[jm].c;
		colF_->coeff[colF_->sizeVal].d = colF->boundVal[jm-1]*(colF->coeff[jm].b-timeFact)+colF->coeff[jm].d;
	} else {
		for(j=jp; j<colF->sizeVal; j++) {
			colF_->coeff[colF_->sizeVal] = colF->coeff[j];
			colF_->boundVal[colF_->sizeVal++] = colF->boundVal[j];
		}
		colF_->coeff[colF_->sizeVal] = colF->coeff[colF->sizeVal];
	}
}

void fillH(int n, TypePiecewise *H, TypePiecewise *H_, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi) {
	if(tree->node[n].child != NOSUCH) {
		if(((double*)tree->info)[n] == UNKNOWN) {
			int c, nchild, k;
			TypePiecewise **list;
			nchild = 0;
			for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling)
				nchild++;
			list = (TypePiecewise**) malloc(nchild*sizeof(TypePiecewise*));
			for(c=tree->node[n].child,k=0; c != NOSUCH; c=tree->node[c].sibling,k++)
				list[k] = &(F_[c]);
			for(c=tree->node[n].child,k=0; c!= NOSUCH; c=tree->node[c].sibling,k++) {
				list[k] = &(H_[n]);
				sumPiecewise(list, nchild, &(H[c]));
				HtoH_(&(H[c]), &(H_[c]), phi(tree->time[c]), ((double*)tree->info)[c]);
				list[k] = &(F_[c]);
				fillH(c, H, H_, F_, tree, phi);
			}
			free((void*)list);
		} else {
			fprintf(stderr, "Critical error: known internal nodes are not handled in this version\n");
			exit(1);
		}
	}
}

void HtoH_(TypePiecewise *H, TypePiecewise *H_, double timeFact, double val) {
	int sizeParBuf, i;
	TypeFraction boundp, boundm;
	if(val == UNKNOWN) {
		sizeParBuf = H->sizePar+1;
		for(i=0; i<=H->sizePar; i++)
			sizeParBuf += 2*H->col[i].sizeVal;
		H_->sizePar = 0;
		H_->boundPar = (TypeFraction*) malloc(sizeParBuf*sizeof(TypeFraction));
		H_->col = (TypeColumn*) malloc((sizeParBuf+1)*sizeof(TypeColumn));
		for(i=0; i<=H->sizePar; i++) {
			TypeFraction finf, fsup;
			int jp, jm, sjp, sjm, ejp, ejm;
			if(i>0)
				finf = H->boundPar[i-1];
			else 
				finf = fract(0,1);	
			if(i<H->sizePar)
				fsup = H->boundPar[i];
			else
				fsup = fract(1,0);
			if(H->col[i].sizeVal >= 0) {
				for(sjp=0; sjp<=H->col[i].sizeVal && (sgnLin(H->col[i].coeff[sjp].a, H->col[i].coeff[sjp].b+timeFact,finf)<=0 && sgnLin(H->col[i].coeff[sjp].a, H->col[i].coeff[sjp].b+timeFact,fsup)<=0); sjp++)
					;
				for(ejp=sjp; ejp<=H->col[i].sizeVal && sgnLin(H->col[i].coeff[ejp].a, H->col[i].coeff[ejp].b+timeFact,finf)*sgnLin(H->col[i].coeff[ejp].a, H->col[i].coeff[ejp].b+timeFact,fsup)<0; ejp++)
					;
				for(ejm=H->col[i].sizeVal+1; ejm>0 && (sgnLin(H->col[i].coeff[ejm-1].a-timeFact, H->col[i].coeff[ejm-1].b,finf)>=0 && sgnLin(H->col[i].coeff[ejm-1].a-timeFact, H->col[i].coeff[ejm-1].b,fsup)>=0); ejm--)
					;
				for(sjm=ejm; sjm>0 && sgnLin(H->col[i].coeff[sjm-1].a-timeFact, H->col[i].coeff[sjm-1].b,finf)*sgnLin(H->col[i].coeff[sjm-1].a-timeFact, H->col[i].coeff[sjm-1].b,fsup)<0; sjm--)
					;
				jp = sjp;
				jm = sjm;
				getColumnH(&(H->col[i]), &(H_->col[H_->sizePar]), jp, jm, timeFact);
				while(jp<ejp || jm<ejm) {
					if(jp<ejp)
						boundp = fract(H->col[i].coeff[jp].b+timeFact, -H->col[i].coeff[jp].a);
					else
						boundp = fsup;
					if(jm<ejm)
						boundm = fract(H->col[i].coeff[jm].b, -H->col[i].coeff[jm].a+timeFact);
					else
						boundm = fsup;
					if(H_->sizePar >= sizeParBuf) {
						sizeParBuf += H->sizePar+2;
						H_->boundPar = (TypeFraction*) realloc((void*) H_->boundPar, sizeParBuf*sizeof(TypeFraction));
						H_->col = (TypeColumn*) realloc((void*) H_->col, (sizeParBuf+1)*sizeof(TypeColumn));
					}
					if(cmpFract(boundp, boundm)<=0)
						H_->boundPar[H_->sizePar] = boundp;
					else
						H_->boundPar[H_->sizePar] = boundm;
					if(cmpFract(boundp, H_->boundPar[H_->sizePar])<=0)
						jp++;
					if(cmpFract(boundm, H_->boundPar[H_->sizePar])<=0)
						jm++;
					H_->sizePar++;
					getColumnH(&(H->col[i]), &(H_->col[H_->sizePar]), jp, jm, timeFact);
				}
				if(i<H->sizePar) {
					H_->boundPar[H_->sizePar] = H->boundPar[i];
					H_->sizePar++; 
				}
			} else {
				fprintf(stderr, "Critical error: known internal nodes are not handled in this version\n");
				exit(1);
			}
		}
		H_->boundPar = (TypeFraction*) realloc((void*) H_->boundPar, H_->sizePar*sizeof(TypeFraction));
		H_->col = (TypeColumn*) realloc((void*) H_->col, (H_->sizePar+1)*sizeof(TypeColumn));
	} else {
		H_->sizePar = H->sizePar;
		H_->boundPar = (TypeFraction*) malloc(H_->sizePar*sizeof(TypeFraction));
		H_->col = (TypeColumn*) malloc((H_->sizePar+1)*sizeof(TypeColumn));
		for(i=0; i<=H_->sizePar; i++) {
			H_->col[i].sizeVal = 1;
			H_->col[i].boundVal = (double*) malloc(sizeof(double));
			H_->col[i].boundVal[0] = H->col[i].boundVal[0];
			H_->col[i].coeff = (TypeCoeff*) malloc(2*sizeof(TypeCoeff));
			H_->col[i].coeff[0].a = 0;
			H_->col[i].coeff[0].b = -timeFact;
			H_->col[i].coeff[0].c =  H->col[i].coeff[0].c - H->col[i].coeff[0].a * H->col[i].boundVal[0];
			H_->col[i].coeff[0].d =  H->col[i].coeff[0].d + H->col[i].boundVal[0] * (H->col[i].coeff[0].b + timeFact);
			
			H_->col[i].coeff[1].a = timeFact;
			H_->col[i].coeff[1].b = 0;
			H_->col[i].coeff[1].c = H->col[i].coeff[0].c - H->col[i].boundVal[0] * (H->col[i].coeff[0].a + timeFact);
			H_->col[i].coeff[1].d =  H->col[i].boundVal[0] * H->col[i].coeff[0].b + H->col[i].coeff[0].d;
			if(i<H_->sizePar)
				H_->boundPar[i] = H->boundPar[i];
		}
	}
}


void getColumnH(TypeColumn *colH, TypeColumn *colH_, int jp, int jm, double timeFact) {
	int j;
	colH_->sizeVal = jm-jp-1;
 	if(jp > 0)
		colH_->sizeVal++;
	if(jm <= colH->sizeVal)
		colH_->sizeVal++;
	colH_->boundVal = (double*) malloc(colH_->sizeVal*sizeof(double));
	colH_->coeff = (TypeCoeff*) malloc((colH_->sizeVal+1)*sizeof(TypeCoeff));
	colH_->sizeVal = 0;
	if(jp > 0) {
		colH_->coeff[colH_->sizeVal].a = 0.;
		colH_->coeff[colH_->sizeVal].b = -timeFact;
		colH_->coeff[colH_->sizeVal].c = colH->coeff[jp].c+colH->coeff[jp].a*colH->boundVal[jp-1];
		colH_->coeff[colH_->sizeVal].d = colH->coeff[jp].d+colH->boundVal[jp-1]*(colH->coeff[jp].b+timeFact);
		colH_->boundVal[colH_->sizeVal++] = colH->boundVal[jp-1];
	}
	if(jm <= colH->sizeVal) {
		for(j=jp; j<jm; j++) {
			colH_->coeff[colH_->sizeVal] = colH->coeff[j];
			colH_->boundVal[colH_->sizeVal++] = colH->boundVal[j];
		}
		colH_->coeff[colH_->sizeVal].a = timeFact;
		colH_->coeff[colH_->sizeVal].b = 0.;
		colH_->coeff[colH_->sizeVal].c = colH->coeff[jm].c+colH->boundVal[jm-1]*(colH->coeff[jm].a-timeFact);
		colH_->coeff[colH_->sizeVal].d = colH->coeff[jm].d+colH->boundVal[jm-1]*colH->coeff[jm].b;
	} else {
		for(j=jp; j<colH->sizeVal; j++) {
			colH_->coeff[colH_->sizeVal] = colH->coeff[j];
			colH_->boundVal[colH_->sizeVal++] = colH->boundVal[j];
		}
		colH_->coeff[colH_->sizeVal] = colH->coeff[colH->sizeVal];
	}
}

TypeResult getResultFromF_(TypePiecewise *h, TypeResult *a, int unknown) {
	TypeResult r;
	if(unknown) {
		TypeFraction minPar;
		int ih, ia, sizeBuf = SIZEBUF;
		r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
		r.val = (double*) malloc((sizeBuf+1)*sizeof(double));
		r.sizePar = 0;
		ih = 0; ia = 0;
		do {
			if(r.sizePar >= sizeBuf) {
				sizeBuf += SIZEBUF;
				r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, sizeBuf*sizeof(TypeFraction));
				r.val = (double*) realloc((void*) r.val, (sizeBuf+1)*sizeof(double));
			}
			if(h->col[ih].sizeVal == 0 || (a->val[ia] >= h->col[ih].boundVal[0] && a->val[ia] <= h->col[ih].boundVal[h->col[ih].sizeVal-1])) {
				r.val[r.sizePar] = a->val[ia];
			} else {
				if(a->val[ia] < h->col[ih].boundVal[0]) {
					r.val[r.sizePar] = h->col[ih].boundVal[0];
				} else {
					r.val[r.sizePar] = h->col[ih].boundVal[h->col[ih].sizeVal-1];
				}
			}
			minPar = fract(1,0);
			if(ih<h->sizePar && cmpFract(h->boundPar[ih],minPar)<0)
				minPar = h->boundPar[ih];
			if(ia<a->sizePar && cmpFract(a->boundPar[ia],minPar)<0)
				minPar = a->boundPar[ia];
			if(ih<h->sizePar && cmpFract(h->boundPar[ih],minPar)<=0)
				ih++;
			if(ia<a->sizePar && cmpFract(a->boundPar[ia],minPar)<=0)
				ia++;
			if(cmpFract(fract(1,0), minPar) != 0)
				r.boundPar[r.sizePar++] = minPar;
		} while(cmpFract(fract(1,0),minPar) != 0);		
	} else {
		r.sizePar = 0;
		r.boundPar = NULL;
		r.val = (double*) malloc(sizeof(double));
		r.val[0] = h->col[0].boundVal[0];
	}
	compactResult(&r);
	return r;	
}

void fillResultRec(int n, TypeResult *res, TypePiecewise *funcFH, TypeTree *tree, TypeFunctTime phi) {
	int c;
	for(c=tree->node[n].child; c >= 0; c=tree->node[c].sibling) {
		res[c] = getResultFromF_(&(funcFH[c]), &(res[n]), isUnknown(c, tree));
		fillResultRec(c, res, funcFH, tree, phi);
	}
}


void compactResult(TypeResult *r) {
	int i, ind = 0;
	for(i=0; i<r->sizePar; i++) {
		if(r->val[i+1] != r->val[i]) {
			r->val[ind] = r->val[i];
			r->boundPar[ind++] = r->boundPar[i];
			cannonizeFract(&(r->boundPar[ind-1]));
		}
	}
	r->val[ind] = r->val[r->sizePar];
	r->sizePar = ind;
	r->boundPar = (TypeFraction*) realloc((void*)r->boundPar, r->sizePar*sizeof(TypeFraction));
	r->val = (double*) realloc((void*)r->val, (r->sizePar+1)*sizeof(double));
}

void compactResultFunc(TypeResultFunc *r) {
	int i, ind = 0;
	for(i=0; i<r->sizePar; i++) {
		if(fabs(r->coeff[i+1].a-r->coeff[i].a)>EPSI || fabs(r->coeff[i+1].b-r->coeff[i].b)>EPSI) {
			r->coeff[ind] = r->coeff[i];
			r->boundPar[ind++] = r->boundPar[i];
		}
	}
	r->coeff[ind] = r->coeff[r->sizePar];
	r->sizePar = ind;
	r->boundPar = (TypeFraction*) realloc((void*)r->boundPar, r->sizePar*sizeof(TypeFraction));
	r->coeff = (TypeCoeffResultFunc*) realloc((void*)r->coeff, (r->sizePar+1)*sizeof(TypeCoeffResultFunc));
}

int compareResultFunc(TypeResultFunc *r, TypeResultFunc *s) {
	int i;
	if(r->sizePar != s->sizePar) {
		return 0;
	}
	for(i=0; i<r->sizePar; i++)
		if(fabs(r->coeff[i].a-s->coeff[i].a)>EPSI || fabs(r->coeff[i].b-s->coeff[i].b)>EPSI) {
			return 0;
		}
	return 1;
}


void freeResult(TypeResult r) {
	if(r.boundPar != NULL)
		free((void*)r.boundPar);
	if(r.val != NULL)
		free((void*)r.val);
}


TypeResult *fillResult(TypePiecewise *F, TypePiecewise *F_, TypeTree *tree, TypeFunctTime phi) {
	TypeResult *res;
	res = (TypeResult*) malloc(tree->size*sizeof(TypeResult));
	res[tree->root] = getResultFromF(&(F[tree->root]), isUnknown(tree->root, tree));
	fillResultRec(tree->root, res, F_, tree, phi);
	return res;
}


TypeResult getResultFromF(TypePiecewise *f, int unknown) {
	TypeResult r;
	if(unknown) {
		int i, sizeBuf;
		sizeBuf = f->sizePar+1;
		for(i=0; i<=f->sizePar; i++)
			sizeBuf += f->col[i].sizeVal;
		r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
		r.val = (double*) malloc((sizeBuf+1)*sizeof(double));
		r.sizePar = 0;
		for(i=0; i<=f->sizePar; i++) {
			TypeFraction finf, fsup;
			int j, jmin, jmax;
			if(i>0)
				finf = f->boundPar[i-1];
			else 
				finf = fract(0,1);			
			if(i<f->sizePar)
				fsup = f->boundPar[i];
			else
				fsup = fract(1,0);
			for(jmin=0; jmin<=f->col[i].sizeVal && (finf.num*f->col[i].coeff[jmin].a+finf.den*f->col[i].coeff[jmin].b)<=0; jmin++);
			for(jmax=f->col[i].sizeVal+1; jmax>=0 && (fsup.num*f->col[i].coeff[jmax-1].a+fsup.den*f->col[i].coeff[jmax-1].b)>=0; jmax--);
			for(j=jmin; j<=jmax; j++) {			
				r.val[r.sizePar] = f->col[i].boundVal[j-1];
				if(j<jmax)
					r.boundPar[r.sizePar++] = fract(f->col[i].coeff[j].b, -f->col[i].coeff[j].a);
			}
			if(i<f->sizePar) 
				r.boundPar[r.sizePar++] = f->boundPar[i];
		}
		r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, r.sizePar*sizeof(TypeFraction));
		r.val = (double*) realloc((void*) r.val, (r.sizePar+1)*sizeof(double));
		compactResult(&r);
	} else {
		r.boundPar = NULL;
		r.val = (double*) malloc(sizeof(double));
		r.sizePar = 0;
		r.val[0] = f->col[0].boundVal[0];
	}
	return r;
}

TypeResultFunc getMinimumCost(TypePiecewise *F) {
	TypeResultFunc r;
	int i, sizeBuf;
	sizeBuf = F->sizePar+1;
	for(i=0; i<=F->sizePar; i++)
		if(F->col[i].sizeVal>0)
			sizeBuf += F->col[i].sizeVal;
	r.boundPar = (TypeFraction*) malloc(sizeBuf*sizeof(TypeFraction));
	r.coeff = (TypeCoeffResultFunc*) malloc((sizeBuf+1)*sizeof(TypeCoeffResultFunc));
	r.sizePar = 0;
	for(i=0; i<=F->sizePar; i++) {
		TypeFraction finf, fsup;
		int j, jmin, jmax;
		if(i>0)
			finf = F->boundPar[i-1];
		else 
			finf = fract(0,1);			
		if(i<F->sizePar)
			fsup = F->boundPar[i];
		else
			fsup = fract(1,0);
		for(jmin=1; jmin<=F->col[i].sizeVal && ((finf.num==0 && F->col[i].coeff[jmin].a<0 && F->col[i].coeff[jmin].b==0.) || sgnLin(F->col[i].coeff[jmin].a, F->col[i].coeff[jmin].b,finf)<0) && sgnLin(F->col[i].coeff[jmin].a, F->col[i].coeff[jmin].b,fsup)<0; jmin++)
			;
			
		for(jmax=F->col[i].sizeVal; jmax>0 && (sgnLin(F->col[i].coeff[jmax-1].a, F->col[i].coeff[jmax-1].b,finf)>=0 && sgnLin(F->col[i].coeff[jmax-1].a, F->col[i].coeff[jmax-1].b,fsup)>=0); jmax--)
			;
		for(j=jmin; j<=jmax; j++) {
			r.coeff[r.sizePar].a = F->col[i].coeff[j].c+F->col[i].coeff[j].a*F->col[i].boundVal[j-1];
			r.coeff[r.sizePar].b = F->col[i].coeff[j].d+F->col[i].coeff[j].b*F->col[i].boundVal[j-1];
			if(j<jmax)
				r.boundPar[r.sizePar++] = fract(F->col[i].coeff[j].b, -F->col[i].coeff[j].a);
		}
		if(i<F->sizePar) 
			r.boundPar[r.sizePar++] = F->boundPar[i];
	}
	r.boundPar = (TypeFraction*) realloc((void*) r.boundPar, r.sizePar*sizeof(TypeFraction));
	r.coeff = (TypeCoeffResultFunc*) realloc((void*) r.coeff, (r.sizePar+1)*sizeof(TypeCoeffResultFunc));
	compactResultFunc(&r);
	return r;
}

void fprintResultFuncDebug(FILE *f, TypeResultFunc *r) {
	int i;
	fprintf(f, "%.2leg+%.2le\n", r->coeff[0].a, r->coeff[0].b);
	for(i=1; i<=r->sizePar; i++)
		fprintf(f, "******* %.2le *******\n%.2leg+%.2le\n", r->boundPar[i-1].num/r->boundPar[i-1].den, r->coeff[i].a, r->coeff[i].b);
}

void fprintStateTreeResult(FILE *f, TypeResult *res, TypeTree *tree) {
	int n;
	char **myName, **myComment, **tmpComment;
	
	myName = (char**) malloc(tree->size*sizeof(char*));
	myComment = (char**) malloc(tree->size*sizeof(char*));
	tmpComment = tree->comment;
	for(n=0; n<tree->size; n++) {
		char buffer[500000];
		if(isUnknown(n, tree)) {
			myName[n] = (char*) malloc((strlen(UNKNOWN_STRING)+1)*sizeof(char));
			strcpy(myName[n], UNKNOWN_STRING);
		} else {
			sprintf(buffer, "%.3lf", ((double*)tree->info)[n]);
			myName[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
			strcpy(myName[n], buffer);
		}
		sprintResult(buffer, &(res[n]));
		myComment[n] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(myComment[n], buffer);
	}
	tmpComment = tree->comment;
	tree->comment = myComment;
	fprintTree(f, tree, display_name);
	tree->comment = tmpComment;
	for(n=0; n<tree->size; n++) {
		free((void*)myName[n]);
		free((void*)myComment[n]);
	}
	free((void*)myName);
	free((void*)myComment);
}


void sprintResult(char *s, TypeResult *r) {
	int i;
	char *tmp;
	tmp = s;
	tmp += sprintf(tmp, "\\hbox{$\\scriptstyle\\{%.3lf", r->val[0]);
	for(i=1; i<=r->sizePar; i++)
		tmp += sprintf(tmp, "\\stackrel{\\frac{%lf}{%lf}}{,}%.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
	tmp += sprintf(tmp, "\\}$}");
} 

void sprintResultDebug(char *s,  TypeResult *r) {
	int i;
	char *tmp;
	tmp = s;
	tmp += sprintf(tmp, "%.3lf", r->val[0]);
	for(i=1; i<r->sizePar; i++)
			tmp += sprintf(tmp, " >%lf/%lf< %.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
}

void fprintResultDebug(FILE *f, TypeResult *r) {
	int i;
	fprintf(f, "%.3lf", r->val[0]);
	for(i=1; i<=r->sizePar; i++)
		fprintf(f, " >%lf/%lf< %.3lf", r->boundPar[i-1].num, r->boundPar[i-1].den, r->val[i]);
}

TypeFraction *mergeFractList(TypeFraction *l, int *size, TypeResult *res) {
	int i, j, ind;
	TypeFraction *new;
	new = (TypeFraction*) malloc((*size+res->sizePar)*sizeof(TypeFraction));
	i = 0; j = 0; ind = 0;
	while(i<*size || j<res->sizePar) {
		TypeFraction min = fract(1,0);
		if(i<*size && cmpFract(l[i], min)<=0)
			min = l[i];
		if(j<res->sizePar && cmpFract(res->boundPar[j], min)<=0)
			min = res->boundPar[j];
		if(i<*size && cmpFract(min, l[i])<=0)
			i++;
		if(j<res->sizePar && cmpFract(min, res->boundPar[j])<=0)
			j++;
		new[ind++] = min;
	}
	*size = ind;
	new = (TypeFraction*) realloc((void*)new, ind*sizeof(TypeFraction));
	return new;
}



	
/*print result as flat list*/	
void fprintResultList(FILE *f, TypeResult *res, TypeTree *tree) {
	int n;
	if(tree->size<=0)
		return;
	for(n=0; n<tree->size; n++) {
		fprintIdentTimeComment(f, n, tree, display_name);
		fprintf(f, "\t");
		fprintResultDebug(f, &(res[n]));
		fprintf(f, "\n");
	}
}

TypeListDouble *getErrors(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree) {
	int i, k, size;
	TypeFraction *tot;
	TypeListDouble *err;
	double min;
	err = (TypeListDouble*) malloc(sizeof(TypeListDouble));
	tot = (TypeFraction*) malloc(res[sta->index[0]].sizePar*sizeof(TypeFraction));
	for(size=0; size<res[sta->index[0]].sizePar; size++)
		tot[size] = res[sta->index[0]].boundPar[size];
	for(k=1; k<sta->size; k++) {
		TypeFraction *tmp;
		tmp = mergeFractList(tot, &size, &(res[sta->index[k]]));
		free((void*) tot);
		tot = tmp;
	}
	err->size = size+1;
	err->val = (double*) malloc(err->size*sizeof(double));

	for(i=0; i<size; i++)
		err->val[i] = 0.;
	for(k=0; k<sta->size; k++) {
		int i, j;
		j = 0;
		for(i=0; i<=size; i++) {
			err->val[i] += fabs(sta->state[k]-res[sta->index[k]].val[j]);
			if(cmpFract(res[sta->index[k]].boundPar[j], tot[i]) <= 0)
				j++;
		}
	}
	min = err->val[0];
	for(i=1; i<=size; i++)
		if(err->val[i]<min) {
			min = err->val[i];
		}
	return err;
}
		
TypeFractionInter getBestInter(TypeResult *res, TypeIndexStateList *sta, TypeTree *tree) {
	int i, k, imin, size;
	TypeFraction *tot;
	TypeListDouble *err;
	TypeFractionInter inter;
	double min;
	err = (TypeListDouble*) malloc(sizeof(TypeListDouble));
	tot = (TypeFraction*) malloc(res[sta->index[0]].sizePar*sizeof(TypeFraction));
	for(size=0; size<res[sta->index[0]].sizePar; size++)
		tot[size] = res[sta->index[0]].boundPar[size];
	for(k=1; k<sta->size; k++) {
		TypeFraction *tmp;
		tmp = mergeFractList(tot, &size, &(res[sta->index[k]]));
		free((void*) tot);
		tot = tmp;
	}
	err->size = size+1;
	err->val = (double*) malloc(err->size*sizeof(double));

	for(i=0; i<size; i++)
		err->val[i] = 0.;
	for(k=0; k<sta->size; k++) {
		int i, j;
		j = 0;
		for(i=0; i<=size; i++) {
			err->val[i] += fabs(sta->state[k]-res[sta->index[k]].val[j]);
			if(cmpFract(res[sta->index[k]].boundPar[j], tot[i]) <= 0)
				j++;
		}
	}
	min = err->val[0]; imin = 0;
	for(i=1; i<=size; i++)
		if(err->val[i]<min) {
			min = err->val[i];
			imin = i;
		}
	if(imin > 0)
		inter.inf = tot[imin-1];
	else
		inter.inf = fract(0,1);
	if(imin < size)
		inter.sup = tot[imin];
	else
		inter.sup = fract(1,0);
	free((void*) tot);
	return inter;
}

