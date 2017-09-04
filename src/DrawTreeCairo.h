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




#ifndef DrawTreeCairoF
#define DrawTreeCairoF
#include "Tree.h"
#include "DrawTreeGeneric.h"

void drawTextCairo(double x, double y, char *text, char *al, TypeParamDrawTreeGeneric *param);
void drawTextAngleCairo(double x, double y, double a, char *text, char *al, TypeParamDrawTreeGeneric *param);
void drawDottedLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void fillWedgeCairo(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void drawWedgeCairo(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void fillGradientCairo(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void setFunctSVG(TypeFunctDrawTreeGeneric *funct);
void startSVGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startSVG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPDF(TypeFunctDrawTreeGeneric *funct);
void startPDFStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPDF(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPS(TypeFunctDrawTreeGeneric *funct);
void startPSStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPS(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPNG(TypeFunctDrawTreeGeneric *funct);
void startPNGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPNG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endPNG(TypeParamDrawTreeGeneric *param);
void endCairo(TypeParamDrawTreeGeneric *param);

#endif
