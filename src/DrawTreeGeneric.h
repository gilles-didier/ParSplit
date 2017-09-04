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




#ifndef DrawTreeGenericF
#define DrawTreeGenericF
#include <stdio.h>
#include "Tree.h"

typedef struct RGB_COLOR {
	double red, green, blue;
} TypeRGB;

typedef struct PARAM_DRAW_TREE_GENERIC {
	double scale, leafSep, leafCur, labelSep, ycenter, ydec, radius, roffset, xoffset, yoffset, height, width, labelWidth, tickLength, vmin, vmax, vscale;
	double tmin, tmax;
	TypeRGB start, end, curgb;
	void *info;
} TypeParamDrawTreeGeneric;

typedef struct FUNCT_DRAW_TREE_GENERIC {
	void (*startStd)(char *, double, double, TypeParamDrawTreeGeneric *), (*start)(char *, TypeTree *, TypeParamDrawTreeGeneric *), (*end)(TypeParamDrawTreeGeneric *), (*drawText)(double, double, char*, char*, TypeParamDrawTreeGeneric*), 
	(*drawLine)(double, double, double, double, TypeParamDrawTreeGeneric*), (*drawDottedLine)(double, double, double, double, TypeParamDrawTreeGeneric*), (*drawTextAngle)(double, double, double, char*, char*, TypeParamDrawTreeGeneric*), 
	(*fillWedge)(TypeRGB, double, double, double, double, TypeParamDrawTreeGeneric*), 
	(*drawWedge)(double, double, double, double, TypeParamDrawTreeGeneric*), (*fillGradient)(TypeRGB, TypeRGB, double, double, double, double, TypeParamDrawTreeGeneric*);
} TypeFunctDrawTreeGeneric;

typedef struct INFO_DRAW_TREE_GENERIC {
	TypeParamDrawTreeGeneric param;
	TypeFunctDrawTreeGeneric funct;
} TypeInfoDrawTreeGeneric;

void fillTime(int n, double tanc, TypeTree *tree, double *min, double *max, int *dmax);
void fillBounds(int n, double tmin, double tmax, TypeTree *tree, double *min, double *max, int *dmax);
void fillUnknownTimes(double tmin, double tmax, TypeTree *tree);

void drawScaleGeneric(double x, double y, TypeInfoDrawTreeGeneric *info);
void drawTreeFileGeneric(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info);
void drawTreeFileGenericDebug(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info);
#endif
