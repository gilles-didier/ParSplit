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
#include <string.h>
#include <math.h>
#include <time.h>
#include "Utils.h"
#include "Tree.h"
#include "StateTree.h"
#include "ContinuousParsimony.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawTreeGenericCut.h"

#define OUTPUT "output"
#define OMEGA "\u03a9"
#define gamma "\u0263"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "\nNAME\n\tParSplit - Identification of the most parsimonious split(s) of a phylogenetic tree with regards of character values at tips.\n\t\nSYNOPSIS\n\tParSplit [OPTIONS] <input Tree File> <input State File> <output File>\n\nDESCRIPTION\n\tIdentify the most parsimonious split(s) of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards of character values at tips contained in <input State File> (which must be in csv format) and output the results in text format in the file <output File> and, when one sets option '-f', in graphic format in the file <output File> with the corresponding extension (i.e .pdf, .png etc.)\n\n\tOptions are\n\t-t <type>\n\t\tset the function applied on the branch length of the tree(which is 'inverse' by default): \n\t\t\t-t i -> inverse\n\t\t\t-t u -> identity\n\t-z <input Tree File>\n\t\toutput the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit\n\t-n <number>\n\t\tset the max number of splits which will be displayed\n\t-f <number>\n\t\tset the graphic format of the output (option is required if one wants a graphic output)\n\t\t\t-f 1 -> pdf\n\t\t\t-f 2 -> postscript\n\t\t\t-f 3 -> png\n\t\t\t-f 4 -> svg\n\t\t\t-f 5 -> LaTeX (psTricks)\n\t\t\t-f 6 -> LaTeX (TikZ)\n\t-c <r1> <g1> <b1> <r2> <g2> <b2>\n\t\tset the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])\n\t-h\n\t\tdisplay help\n"

double fUn(double t) {
	return 1;
}

double fInv(double t) {
	return 1/t;
}

int main(int argc, char **argv) {
	char *inputFileNameTree, *inputFileNameState, *outputFileName, option[256], type = 'i', format = '0';
	int i, number = -1;
	FILE *fi, *fo, *fs;
	TypeInfoDrawTreeGeneric info;
	
	info.param.start = (TypeRGB) {.red = 0., .green = 0., .blue = 1.};
	info.param.end = (TypeRGB) {.red = 1., .green = 0., .blue = 0.};
	info.param.curgb = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &type) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -t");
		}
		if(option['n']) {
			option['n'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &number) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -t");
		}
		if(option['z']) {
			option['z'] = 0;
			if((i+1)<argc)
				inputFileNameTree = argv[++i];
			else
				exitProg(ErrorArgument, "a file name is required after option -x");
			if((fi = fopen(inputFileNameTree, "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				printTreeDebug(stdout, tree->root, tree, tree->name);
				bltoabsTime(tree);
				reorderTreeSize(tree);
				if(tree->minTime == NO_TIME || tree->minTime == 0.)
					tree->minTime = tree->time[tree->root]*0.9;
				if(tree->maxTime == NO_TIME) {
					int n;
					tree->maxTime = 0.;
					for(n=0; n<tree->size; n++)
						if(tree->time[n]>tree->maxTime)
							tree->maxTime = tree->time[n];
				}
				setFunctPDF(&(info.funct));
				drawTreeFileGenericDebug("tree.pdf", tree, &info);
			} else {
				fprintf(stderr, "Error while reading %s.\n", inputFileNameTree);
				exit(1);
			}
			exit(0);
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.red)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.green)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.start.blue)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.red)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.green)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(info.param.end.blue)) == 1)
				i++;
			else
				exitProg(ErrorArgument, "6 numbers are expected after -c");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -f");
		}
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exitProg(ExitOk, NULL);
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameState = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing character values of tips in csv format\n");
		exit(1);
	}
	if(i<argc) {
		outputFileName = argv[i++];
	} else {
		outputFileName = "out.txt";
	}
	if((fi = fopen(inputFileNameTree, "r"))) {
		int n, *tchange, *tmode;
		double *tmin;
		TypeTree *tree;
		TypePiecewise *F, *F_, *H, *H_;
		tree = readTree(fi);
		fclose(fi);
		fixTreeTime(tree);
		fillState(tree);
		if(strlen(inputFileNameState)>0 && (fs = fopen(inputFileNameState, "r"))) {
			setState(fs, tree);
			fclose(fs);
		} else {
			fprintf(stderr, "Error while reading %s.\n", inputFileNameState);
			exit(1);
		}
		F = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
		F_ = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
		H = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
		H_ = (TypePiecewise*) malloc(tree->size*sizeof(TypePiecewise));
		switch(type) {
			case 'i':
				fillNormFH(tree, fInv, F, F_, H, H_);
				break;
			case 'u':
			default:
				fillNormFH(tree, fUn, F, F_, H, H_);
		}
		getAllCut(tree, F, F_, H, H_, &tchange, &tmode, &tmin);
		if((fo = fopen(outputFileName, "w"))) {
			for(i=0; tchange[i] != NOSUCH; i++) {
				fprintf(fo, "%d", i+1);
				fprintf(fo, "\t%.4lf", tmin[i]);
				switch(tmode[i]) {
					case 1:
						fprintf(fo, "\tB");
						break;
					case 0:
					default:
						fprintf(fo, "\tA");
						break;
				}
				fprintf(fo, "\t%d", tchange[i]);
				fprintf(fo, "\n");
			}
			fprintf(fo, "no split\t%.4lf\n", intPiecewise(&(F[tree->root])));
			fclose(fo);
		} else {
			fprintf(stderr, "Error while writing %s.\n", outputFileName);
			exit(1);
		}
		if(number>=0) {
			for(i=0; i<number && tchange[i] != NOSUCH; i++)
				;
			tchange[i] = NOSUCH;
		}	
		if(format != '0') {
			char *tmp, outputFileNameG[SIZE_BUFFER_CHAR];
			double *level;
			strcpy(outputFileNameG, outputFileName);
			if((tmp = strrchr(outputFileNameG, '.')) != NULL)
				tmp[0] = '\0';
			bltoabsTime(tree);
			reorderTreeSize(tree);
			if(tree->minTime == NO_TIME || tree->minTime == 0.)
				tree->minTime = tree->time[tree->root]*0.9;
			if(tree->maxTime == NO_TIME) {
				int n;
				tree->maxTime = 0.;
				for(n=0; n<tree->size; n++)
					if(tree->time[n]>tree->maxTime)
						tree->maxTime = tree->time[n];
			}
			switch(format) {
				case '1':
					strcat(outputFileNameG, ".pdf");
					setFunctPDF(&(info.funct));
					break;
				case '2':
					strcat(outputFileNameG, ".ps");
					setFunctPS(&(info.funct));
					break;
				case '3':
					strcat(outputFileNameG, ".png");
					setFunctPNG(&(info.funct));
					break;
				case '4':
					strcat(outputFileNameG, ".svg");
					setFunctSVG(&(info.funct));
					break;
				case '5':
					strcat(outputFileNameG, ".tex");
					setFunctPSTricks(&(info.funct));
					break;
				case '6':
					strcat(outputFileNameG, ".tex");
					setFunctTikz(&(info.funct));
					break;
				default:
					strcat(outputFileNameG, ".pdf");
					setFunctPDF(&(info.funct));
			}
			if(tchange[0] != NOSUCH) {
				double min, max;
				int imin;
				min = tmin[0];
				imin = 0;
				max = intPiecewise(&(F[tree->root]));
				for(i=1; tchange[i] != NOSUCH; i++)
					if(tmin[i]<min) {
						imin = i;
						min = tmin[i];
					}
				if(imin != 0) {
					int itmp;
					double dtmp;
					itmp = tchange[0]; tchange[0] = tchange[imin]; tchange[imin] = itmp;
					itmp = tmode[0]; tmode[0] = tmode[imin]; tmode[imin] = itmp;
					dtmp = tmin[0]; tmin[0] = tmin[imin]; tmin[imin] = dtmp;
				}
				level = (double*) malloc(i*sizeof(double));
				if(max == min)
					min--;
				for(i=0; tchange[i] != NOSUCH; i++)
					level[i] = (max-tmin[i])/(max-min);
			} else {
				level = (double*) malloc(sizeof(double));
				level[0] = 1;
			}
			drawTreeFileGenericCutAll(outputFileNameG, tree, tchange, tmode, level, &info);
			free((void*)level);
			free((void*)tchange);
			free((void*)tmode);
			free((void*)tmin);
		}
		for(n=0; n<tree->size; n++) {
			if(n != tree->root)
				freePiecewise(F_[n]);
			freePiecewise(F[n]);
			if(n != tree->root)
				freePiecewise(H[n]);
			freePiecewise(H_[n]);
		}
		free((void*)F);
		free((void*)F_);
		free((void*)H);
		free((void*)H_);
		if(tree->info != NULL) {
			free((void*)tree->info);
			tree->info = NULL;
		}
		freeTree(tree);
	} else {
		fprintf(stderr, "Error while reading %s.\n", inputFileNameTree);
		exit(1);
	}
	return 0;
}

