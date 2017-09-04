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




#ifndef DrawTreeTikzF
#define DrawTreeTikzF
#include "Tree.h"
#include "DrawTreeGeneric.h"

#ifdef __cplusplus
extern "C" {
#endif



void drawTextTikz(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawTextAngleTikz(double x0, double y0, double a, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawDottedLineTikz(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void setFunctTikz(TypeFunctDrawTreeGeneric *funct);
void startTikz(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endTikz(TypeParamDrawTreeGeneric *param);
char *sprintRGBTikz(char *buffer, TypeRGB rgb);
double getMaxLeafLabelWidthTikz(TypeTree *tree);

void startTikzStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startTikz(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endTikz(TypeParamDrawTreeGeneric *param);
void drawTextTikz(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void fillWedgeTikz(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void drawWedgeTikz(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void fillGradientTikz(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);

#ifdef __cplusplus
}
#endif

#endif
