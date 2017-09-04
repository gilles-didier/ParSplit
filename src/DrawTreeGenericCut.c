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
#include "DrawTreeGenericCut.h"

#define NB_MIN 6
#define OFFSET 10
#define TICK_LENGTH 10
#define FONT_NAME "Helvetica"
#define FONT_SIZE 9.
#define LABEL_SEP 10.
#define MAX_STRING_SIZE 500
#define STANDARD_WIDTH 1000
#define CHAR_WIDTH 10
#define NONE -1

static void fillName(TypeTree *tree, char **nameSave);
static double drawNodeGeneric(int n, int parent, int *type, TypeTree *tree, TypeInfoDrawTreeGeneric *info);
static double drawNodeGenericAll(int n, int parent, int *type, int *changeA, int *changeB, double *level, TypeTree *tree, TypeInfoDrawTreeGeneric *info);
static void fillTypeMode(int n, int *type, TypeTree *tree);
static double convex(double s, double e, double v);
static TypeRGB getRGB(TypeRGB start, TypeRGB end, double v);

void fillName(TypeTree *tree, char **nameSave) {
	int i;
	char buffer[300], *tmp;
	tree->name = (char**) malloc(tree->size*sizeof(char*));
	for(i=0; i<tree->size; i++) {
		buffer[0] = '\0';
		tmp = buffer;
		if(nameSave != NULL && nameSave[i] != NULL)
			tmp += sprintf(buffer, "%s", nameSave[i]);
		if(tree->info != NULL &&  ((double*)tree->info)[i] != UNKNOWN)
			tmp += sprintf(tmp, " (%.1lf m)", ((double*)tree->info)[i]);
		tree->name[i] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(tree->name[i], buffer);
	}
}

double convex(double s, double e, double v) {
	return (1.-v)*s+v*e;
}

TypeRGB getRGB(TypeRGB start, TypeRGB end, double v) {
	TypeRGB res = {.red = convex(start.red, end.red, v), .green = convex(start.green, end.green, v), .blue = convex(start.blue, end.blue, v)};
	return res;
}

void fillTypeMode(int n, int *type, TypeTree *tree) {
	int c;
	if(n == NOSUCH)
		return;
	for(c=tree->node[n].child; c != NOSUCH; c=tree->node[c].sibling) {
		if(type[n] != NONE && type[c] == NONE)
			type[c] = type[n];
		fillTypeMode(c, type, tree);
	}
}

void drawTreeFileGenericCutAll(char *filename, TypeTree *tree, int *node, int *mode, double *min, TypeInfoDrawTreeGeneric *info) {
	int i, *type, *changeA, *changeB;
	double y;
	double *timeSave;
	char **nameSave;
	TypeRGB rgbSave;
	rgbSave = info->param.curgb;
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
	fillName(tree, nameSave);
	changeA = (int*) malloc(tree->size*sizeof(int));
	changeB = (int*) malloc(tree->size*sizeof(int));
	for(i=0; i<tree->size; i++) {
		changeA[i] = 0;
		changeB[i] = 0;
	}
	for(i=0; node[i] != NOSUCH; i++) {
		if(mode[i] == 0)
			changeA[node[i]] = i+1;
		else
			changeB[node[i]] = i+1;
	}
	type = (int*) malloc(tree->size*sizeof(int));
	for(i=0; i<tree->size; i++)
		type[i] = NONE;
	for(i=0; node[i] != NOSUCH; i++)
		if(mode[i] == 0) {
			int c;
			for(c=tree->node[node[i]].child; c != NOSUCH; c=tree->node[c].sibling)
				type[c] = i;
		} else
			type[node[i]] = i;
	fillTypeMode(tree->root, type, tree);
	info->funct.start(filename, tree, &(info->param));
	y = drawNodeGenericAll(tree->root, tree->root, type, changeA, changeB, min, tree, info);
	if(type[tree->root] == 1)
		info->param.curgb = info->param.end;
	else
		info->param.curgb = info->param.start;
	info->funct.drawLine(info->param.xoffset, info->param.yoffset+y, (tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	info->param.curgb = rgbSave;
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	free((void*) type);
	free((void*) changeA);
	free((void*) changeB);
	tree->time = timeSave;
	if(tree->name != NULL) {
		for(i=0; i<tree->size; i++)
			if(tree->name[i] != NULL)
				free((void*) tree->name[i]);
		free((void*) tree->name);
	}
	tree->name = nameSave;
}

void drawTreeFileGenericCut(char *filename, TypeTree *tree, int node, int mode, TypeInfoDrawTreeGeneric *info) {
	int i, *type;
	double y;
	double *timeSave;
	char **nameSave;
	TypeRGB rgbSave;
	rgbSave = info->param.curgb;
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
	fillName(tree, nameSave);
	type = (int*) malloc(tree->size*sizeof(int));
	for(i=0; i<tree->size; i++)
		type[i] = NONE;
	if(node != NOSUCH) {
		if(mode == 0) {
			int c;
			for(c=tree->node[node].child; c != NOSUCH; c=tree->node[c].sibling)
				type[c] = 0;
		} else
			type[node] = 0;
	}
	fillTypeMode(tree->root, type, tree);
	info->funct.start(filename, tree, &(info->param));
	y = drawNodeGeneric(tree->root, tree->root, type, tree, info);
	if(type[tree->root] == 1)
		info->param.curgb = info->param.end;
	else
		info->param.curgb = info->param.start;
	info->funct.drawLine(info->param.xoffset, info->param.yoffset+y, (tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	info->param.curgb = rgbSave;
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	free((void*) type);
	tree->time = timeSave;
	if(tree->name != NULL) {
		for(i=0; i<tree->size; i++)
			if(tree->name[i] != NULL)
				free((void*) tree->name[i]);
		free((void*) tree->name);
	}
	tree->name = nameSave;
}

void drawTreeFileGenericCutDebug(char *filename, TypeTree *tree, int node, int mode, TypeInfoDrawTreeGeneric *info) {
	int i, *type;
	double *timeSave;
	char **nameSave;
	TypeRGB rgbSave;
	rgbSave = info->param.curgb;
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
		if(tree->info != NULL) {
			if(nameSave != NULL && nameSave[i] != NULL)
				sprintf(buffer, "%s/%d (%.1lf m)", nameSave[i], i, ((double*)tree->info)[i]);
			else
				sprintf(buffer, "%d (%.1lf m)", i, ((double*)tree->info)[i]);
		} else {
			if(nameSave != NULL && nameSave[i] != NULL)
				sprintf(buffer, "%s/%d", nameSave[i], i);
			else
				sprintf(buffer, "%d", i);
		}
		tree->name[i] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(tree->name[i], buffer);
	}
	type = (int*) malloc(tree->size*sizeof(int));
	for(i=0; i<tree->size; i++)
		type[i] = 0;
	if(mode == 0) {
		int c;
		for(c=tree->node[node].child; c != NOSUCH; c=tree->node[c].sibling)
			type[c] = 1;
	} else
		type[node] = 1;
	fillTypeMode(tree->root, type, tree);
	info->funct.start(filename, tree, &(info->param));
	drawNodeGeneric(tree->root, tree->root, type, tree, info);
	info->param.curgb = rgbSave;
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	free((void*) tree->time);
	free((void*) type);
	tree->time = timeSave;
	if(tree->name != NULL) {
		for(i=0; i<tree->size; i++)
			if(tree->name[i] != NULL)
				free((void*) tree->name[i]);
		free((void*) tree->name);
	}
	tree->name = nameSave;
}

double drawNodeGeneric(int n, int parent, int *type, TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	double min, max, y;
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		min = drawNodeGeneric(tmp, n, type, tree, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling)
			max = drawNodeGeneric(tmp, n, type, tree, info);
		y = (min+max)/2;
		if(type[tree->node[n].child] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = info->param.end;
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+min, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
		if(type[tree->node[tree->node[n].child].sibling] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = info->param.end;
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+max, &(info->param));
		if(type[n] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = info->param.end;
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[n], "l", &(info->param));
	} else {
		info->param.leafCur += info->param.leafSep;
		y = info->param.leafCur;
		if(type[n] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = info->param.end;
		if(tree->name && tree->name[n]) {
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, info->param.yoffset+y, tree->name[n], "l", &(info->param));
		}
	}
	if(type[n] == NONE)
		info->param.curgb = info->param.start;
	else
		info->param.curgb = info->param.end;
	info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[parent]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	return y;

}

double drawNodeGenericAll(int n, int parent, int *type, int *changeA, int *changeB, double *level, TypeTree *tree, TypeInfoDrawTreeGeneric *info) {
	double min, max, y;
	char buffer[10];
	TypeRGB black = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	if(tree->node[n].child >= 0) {
		int tmp = tree->node[n].child;
		min = drawNodeGenericAll(tmp, n, type, changeA, changeB, level, tree, info);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling)
			max = drawNodeGenericAll(tmp, n, type, changeA, changeB, level, tree, info);
		y = (min+max)/2;
		if(type[tree->node[n].child] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = getRGB(info->param.start, info->param.end, level[type[tree->node[n].child]]);
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+min, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
		if(type[tree->node[tree->node[n].child].sibling] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = getRGB(info->param.start, info->param.end, level[type[tree->node[tree->node[n].child].sibling]]);
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+max, &(info->param));
		if(type[n] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = getRGB(info->param.start, info->param.end, level[type[n]]);
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[n], "l", &(info->param));
	} else {
		info->param.leafCur += info->param.leafSep;
		y = info->param.leafCur;
		if(type[n] == NONE)
			info->param.curgb = info->param.start;
		else
			info->param.curgb = getRGB(info->param.start, info->param.end, level[type[n]]);
		if(tree->name && tree->name[n]) {
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, info->param.yoffset+y, tree->name[n], "l", &(info->param));
		}
	}
	if(type[n] == NONE)
		info->param.curgb = info->param.start;
	else
		info->param.curgb = getRGB(info->param.start, info->param.end, level[type[n]]);
	info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[parent]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	info->param.curgb = black;
	if(changeA[n] >= 1) {
		info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, "X", "c", &(info->param));
		sprintf(buffer, "%d", changeA[n]);
		info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.tickLength, info->param.yoffset+y, buffer, "l", &(info->param));
	}
	if(changeB[n] >= 1) {
		double tmp = (3.*tree->time[parent]+tree->time[n])/4.;
		info->funct.drawText((tmp-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, "X", "c", &(info->param));
		sprintf(buffer, "%d", changeB[n]);
		info->funct.drawText((tmp-info->param.tmin)*info->param.scale+info->param.xoffset+2*info->param.tickLength, info->param.yoffset-2*info->param.tickLength+y, buffer, "lb", &(info->param));
	}
	return y;
}

