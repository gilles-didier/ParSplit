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




#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include "Utils.h"
#include "DrawTreeGeneric.h"

#define NB_MIN 6
#define OFFSET 10
#define TICK_LENGTH 10
#define FONT_NAME "Helvetica"
#define FONT_SIZE 9.
#define LABEL_SEP 10.
#define MAX_STRING_SIZE 500
#define STANDARD_WIDTH 1000
#define CHAR_WIDTH 10

static double drawNodeGeneric(int n, int parent, TypeTree *tree,  TypeInfoDrawTreeGeneric *param);
static void drawTreeGeneric(TypeTree *tree, TypeInfoDrawTreeGeneric *param);

void fillTime(int n, double tanc, TypeTree *tree, double *min, double *max, int *dmax) {
	int c;
	if(tree->time[n] == NO_TIME) {
		double tmp = utils_MAX(tanc, min[n]);
		if(max[n]<tanc)
			printf("\n\nProblem %.2lf %.2lf\n%d\n", max[n], tanc, n);
		if(tree->node[n].child>=0)
			tree->time[n] = tmp+(max[n]-tmp)/((double)(2+dmax[n]));
		else
			tree->time[n] = max[n];
	}
	for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling)
		fillTime(c, tree->time[n], tree, min, max, dmax);
}

void fillBounds(int n, double tmin, double tmax, TypeTree *tree, double *min, double *max, int *dmax) {
	int c;
	if(tree->time[n] != NO_TIME) {
		min[n] = tree->time[n];
	} else {
		min[n] = tmin;
	}
	for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling)
		fillBounds(c, min[n], tmax, tree, min, max, dmax);
	if(tree->time[n] != NO_TIME) {
		max[n] = tree->time[n];
		dmax[n] = 0;
	} else {
		if(tree->node[n].child==NOSUCH) {
			max[n] = tmax;
			dmax[n] = 0;
		} else {	
			max[n] = tmax+1;
			dmax[n] = 0;
			for(c=tree->node[n].child; c!=NOSUCH; c=tree->node[c].sibling) {
				if((max[c])<(max[n])) {
					max[n] = max[c];
					dmax[n] = dmax[c]+1;
				}
			}
		}
	}
}

void fillUnknownTimes(double tmin, double tmax, TypeTree *tree) {
	int *dmax;
	double *min, *max;
	
	min = (double*) malloc(tree->size*sizeof(double));
	max = (double*) malloc(tree->size*sizeof(double));
	dmax = (int*) malloc(tree->size*sizeof(int));
	fillBounds(tree->root, tmin, tmax, tree, min, max, dmax);
	fillTime(tree->root, tmin, tree, min, max, dmax);
	free((void*)min);
	free((void*)max);
	free((void*)dmax);
}

void drawScaleGeneric(double x, double y, TypeInfoDrawTreeGeneric *info) {
	int flag = 0;
	double start, end, step, cur, width;
	width = info->param.width-info->param.labelWidth-info->param.xoffset-info->param.labelSep;
	info->funct.drawLine(x, y, x+width, y, &(info->param));
	if((info->param.tmax-info->param.tmin) <= 0.)
		return;
	step = pow(10., floor(log10(info->param.tmax-info->param.tmin)));
	if((info->param.tmax-info->param.tmin)/step<NB_MIN)
		step /= 2.;
	flag = step<1.;
	start = step*ceil(info->param.tmin/step);
	end = step*floor(info->param.tmax/step);
	for(cur=start; cur<=end; cur += step) {
		char tmp[500];
		info->funct.drawLine(x+(cur-info->param.tmin)*info->param.scale, y, x+(cur-info->param.tmin)*info->param.scale, y+info->param.tickLength, &(info->param));
		if(flag)
			sprintf(tmp, "%.1lf", cur);
		else
			sprintf(tmp, "%.0lf", cur);
		info->funct.drawText(x+(cur-info->param.tmin)*info->param.scale, y+info->param.tickLength, tmp, "t", &(info->param));
	}
}

void drawTreeFileGeneric(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	int i;
	double *timeSave;
	
	timeSave = tree->time;
	tree->time = (double*) malloc(tree->size*sizeof(double));
	if(tree->time != NULL) {
		int hasUnknown = 0;
		for(i=0; i<tree->size; i++) {
			tree->time[i] = timeSave[i];
			if(timeSave[i] == NO_TIME)
				hasUnknown = 1;
		}
		if(hasUnknown)
			fillUnknownTimes(info->param.tmin, info->param.tmax,  tree);
	} else {
		for(i=0; i<tree->size; i++)
			tree->time[i] = 1.;
		bltoabsTime(tree);
	}
	info->funct.start(filename, tree, &(info->param));
	drawTreeGeneric(tree, info);
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	tree->time = timeSave;
}

void drawTreeFileGenericDebug(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	int i;
	double *timeSave;
	char **nameSave;
	
	timeSave = tree->time;
	tree->time = (double*) malloc(tree->size*sizeof(double));
	if(tree->time != NULL) {
		int hasUnknown = 0;
		for(i=0; i<tree->size; i++) {
			tree->time[i] = timeSave[i];
			if(timeSave[i] == NO_TIME)
				hasUnknown = 1;
		}
		if(hasUnknown)
			fillUnknownTimes(info->param.tmin, info->param.tmax,  tree);
	} else {
		for(i=0; i<tree->size; i++)
			tree->time[i] = 1.;
		bltoabsTime(tree);
	}
	nameSave = tree->name;
	tree->name = (char**) malloc(tree->size*sizeof(char*));
	for(i=0; i<tree->size; i++) {
		char buffer[200];
		if(nameSave != NULL && nameSave[i] != NULL)
			sprintf(buffer, "%s/%d", nameSave[i], i);
		else
			sprintf(buffer, "%d", i);
		tree->name[i] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(tree->name[i], buffer);
	}
	info->funct.start(filename, tree, &(info->param));
	drawTreeGeneric(tree, info);
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	tree->time = timeSave;
}


void drawTreeGeneric(TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	int tmp;
	double min, max, y;
	if(tree->size<=0)
		return;
	if((tmp = tree->node[tree->root].child) >= 0) {
		min = drawNodeGeneric(tmp, tree->root, tree, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling) {
			max = drawNodeGeneric(tmp, tree->root, tree, info);
		}
	} else {
		max = info->param.leafCur+info->param.leafSep/2.;
		min = max;
	}
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, min, (tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, max, &(info->param));
	y = (min+max)/2;
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, info->param.xoffset, info->param.yoffset+y, &(info->param));
	if(tree->name != NULL && tree->name[tree->root] != NULL)
		info->funct.drawText((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[tree->root], "l", &(info->param));
}

double drawNodeGeneric(int n, int parent, TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	double min, max, y;
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		min = drawNodeGeneric(tmp, n, tree, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling)
			max = drawNodeGeneric(tmp, n, tree, info);
		y = (min+max)/2;
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+min, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+max, &(info->param));
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[n], "l", &(info->param));
	} else {
		info->param.leafCur += info->param.leafSep;
		y = info->param.leafCur;
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, info->param.yoffset+y, tree->name[n], "l", &(info->param));
	}
	info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[parent]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	return y;
}
















