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




#include <math.h>

#include "Fraction.h"

#define EPSILON 0.0000000000001
static int sign(double x);


TypeFraction fract(double num, double den) {
	TypeFraction f;
	f.num = num;
	f.den = den;
	return f;
}

/*double fractToPercent(TypeFraction a) {
	if(a.num == 0 && a.den == 0)
		return 50.;
	if(a.num<a.den)
		return (50.*a.num)/a.den;
	else
		return 1-(50.*a.den)/a.num;
}
*/

double fractToPercent(TypeFraction a) {
	if(a.num == 0 && a.den == 0)
		return 50.;
	return 100.-(100.*(((double)a.den))/((double)a.num + a.den));
}

int sign(double x) {
	return (x>0) - (x<0);
}

int cmpFract(TypeFraction a, TypeFraction b) {
	if(a.den == 0.) {
		if(b.den != 0) {
			if(a.num>0)
				return 1;
			if(a.num<0)
				return -1;
			if(b.num/b.den<1)
				return 1;
			if(b.num/b.den>1)
				return -1;
			return 0;
		} else {
			if(a.num>0) {
				if(b.num>0)
					return 0;
				return 1;
			}
			if(a.num<0) {
				if(b.num<0)
					return 0;
				return -1;
			}
			if(b.num>0)
				return -1;
			if(b.num<0)
				return 1;
			return 0;
		}
	}
	if(b.den == 0.) {
		if(b.num>0)
			return -1;
		if(b.num<0)
			return 1;
		if(a.num/a.den<1)
			return -1;
		if(a.num/a.den>1)
			return 1;
		return 0;
	}
	if(fabs(a.num*b.den-a.den*b.num)<EPSILON)
		return 0;
	if((sign(a.den)*sign(b.den)) < 0) {
		if(a.num*b.den>a.den*b.num)
			return -1;
		if(a.num*b.den<a.den*b.num)
			return 1;
		return 0;
	}
	if(a.num*b.den>a.den*b.num)
		return 1;
	if(a.num*b.den<a.den*b.num)
		return -1;
	return 0;
}

void cannonizeFract(TypeFraction *fr) {
	return;
/*	int a, b;
	if(fr->den == 0 || fr->num == 0)
		return;
	if(fr->den<fr->num) {
		a = fr->num;
		b = fr->den;
	} else {
		a = fr->den;
		b = fr->num;
	}
	do {
		int r;
		r = a%b;
		a = b;
		b = r;
	} while(b != 0);
	fr->den /= a;
	fr->num /= a;*/
}

void fprintfracDebug(FILE *f, double n, double d) {
	if(d==0) {
		fprintf(f, "inf");
	} else {
		if(d==1 || n ==0) {
			fprintf(f, "%.2lf", n);
		} else {
			fprintf(f, "%.2lf/%.2lf", n, d);
		}
	}
}

void fprintFrac(FILE *f, TypeFraction frac) {
	if(frac.den == 0) {
		if(frac.num > 0)
			fprintf(f, "+inf");
		if(frac.num < 0)
			fprintf(f, "-inf");
		if(frac.num == 0)
			fprintf(f, "?");
		return;
	}
	if(frac.num == 0 || frac.den == 1)
		fprintf(f, "%.2lf", frac.num);
	else
		fprintf(f, "%.2lf/%.2lf", frac.num, frac.den);
}

void fprintFracLatex(FILE *f, TypeFraction frac) {
	if(frac.den == 0) {
		if(frac.num > 0)
			fprintf(f, "+\\infty");
		if(frac.num < 0)
			fprintf(f, "-\\infty");
		if(frac.num == 0)
			fprintf(f, "?");
		return;
	}
	if(frac.num == 0 || frac.den == 1)
		fprintf(f, "%.2lf", frac.num);
	else
		fprintf(f, "\\frac{%.2lf}{%.2lf}", frac.num, frac.den);
}	

void sprintFrac(char *s, TypeFraction frac) {
	char *tmp = s;
	if(frac.den == 0) {
		if(frac.num > 0)
			tmp += sprintf(tmp, "+inf");
		if(frac.num < 0)
			tmp += sprintf(tmp, "-inf");
		if(frac.num == 0)
			tmp += sprintf(tmp, "?");
		return;
	}

	if(frac.num == 0 || frac.den == 1)
		tmp += sprintf(tmp, "%.2lf", frac.num);
	else
		tmp += sprintf(tmp, "%.2lf/%.2lf", frac.num, frac.den);
}

void sprintFracLatex(char *s, TypeFraction frac) {
	char *tmp = s;
	if(frac.den == 0) {
		if(frac.num > 0)
			tmp += sprintf(tmp, "+\\infty");
		if(frac.num < 0)
			tmp += printf(tmp, "-\\infty");
		if(frac.num == 0)
			tmp += printf(tmp, "?");
		return;
	}
	if(frac.den == 1)
		tmp += printf(tmp, "%.2lf", frac.num);
	else
		tmp += printf(tmp, "\\frac{%.2lf}{%.2lf}", frac.num, frac.den);
}	
