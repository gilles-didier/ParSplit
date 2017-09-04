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

#define OUTPUT "output"

#define MIN_TIME 0.001

#define SIZE_BUFFER_CHAR 300
#define HELPMESSAGE "\nNAME\n\tFix - Fixing data to be used by ParSplit.\n\t\nSYNOPSIS\n\tFix [OPTIONS] <input Tree File> <input State File>\n\nDESCRIPTION\n\tMake sure that all the tip idents of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) are also idents of character values in the file <input State File> (which must be in csv format) and reciprocally. It outputs two files with '_fixed' containing the fixed tree and character values.\n\n\tOptions are\n\t-h\n\t\tdisplay help\n"

void fixName(char **name) {
    int i, ind;
    if((*name) == NULL)
		return;
    for(i=0; (*name)[i] != '\0'; i++)
        if((*name)[i] == '_')
            (*name)[i] = ' ';	
    for(i=0; (*name)[i] != '\0' && (*name)[i]==' '; i++)
		;
	ind = 0;
    for(; (*name)[i] != '\0'; i++)
		(*name)[ind++] = (*name)[i];
	ind--;
    for(; ind>=0 && (*name)[ind] == ' '; ind--)
		;
	(*name)[++ind] = '\0';
	(*name) = (char*) realloc((void*)(*name), ind*sizeof(char));
}

int main(int argc, char **argv) {		
	char *inputFileNameTree, *inputFileNameState, outputFileNameTree[SIZE_BUFFER_CHAR], outputFileNameState[SIZE_BUFFER_CHAR],
	option[256];
	int i;
	FILE *fi, *fo, *fs;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['h']) {
			printf("%s\n", HELPMESSAGE);
			exit(0);
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
	if((fi = fopen(inputFileNameTree, "r"))  && (fs = fopen(inputFileNameState, "r"))) {
		TypeTree *tree, *tree1, *tree2;
		TypeNameStateList *list;
		TypeLexiTree *dict;
		char *ctmp;
		int n, size;

		tree = readTree(fi);
		fclose(fi);
		list = readStateList(fs);
		fclose(fs);
		dict = newLexiTree();
		for(i=0; i<list->size; i++)
			if(addWordLexi(list->name[i], i, dict)>=0)
                fprintf(stderr, "Warning! duplicate identifier '%s' in states\n", list->name[i]);
		tree1 = pruneLeavesFromDict(tree, dict);
		freeTree(tree);
		tree2 = fixBinary(tree1);
		freeTree(tree1);
		freeLexiTree(dict);
		dict = newLexiTree();
		for(n=0; n<tree2->size; n++)
			if(addWordLexi(tree2->name[n], n, dict)>=0)
                fprintf(stderr, "Warning! duplicate identifier '%s' in tree\n", tree2->name[n]);
		size = 0;
		for(i=0; i<list->size; i++) {
			if(findWordLexi(list->name[i], dict) != -1) {
				list->name[size] = list->name[i];
				list->state[size] = list->state[i];
				size++;
			} else {
				if(list->name[i] != NULL)
					free((void*)list->name[i]);
			}
		}
		list->size = size;
		freeLexiTree(dict);
		for(i=0; i<list->size; i++)
			fixName(&(list->name[i]));
		for(n=0; n<tree->size; n++)
			fixName(&(tree->name[n]));
		ctmp = strrchr(inputFileNameTree, '.');
		if(ctmp != NULL)
			ctmp[0] = '\0';
		sprintf(outputFileNameTree, "%s_fixed.newick", inputFileNameTree);
		ctmp = strrchr(inputFileNameState, '.');
		if(ctmp != NULL)
			ctmp[0] = '\0';
		sprintf(outputFileNameState, "%s_fixed.csv", inputFileNameState);
		if((fo = fopen(outputFileNameTree, "w"))) {
			fprintTree(fo, tree2, display_time_name);
			fclose(fo);
		} else {
			fprintf(stderr, "Error while writing %s.\n", outputFileNameTree);
			exit(1);
		}
		if((fo = fopen(outputFileNameState, "w"))) {
			fprintStateList(fo, list);
			fclose(fo);
		} else{
			fprintf(stderr, "Error while writing %s.\n", outputFileNameState);
			exit(1);
		}		freeTree(tree2);
		freeNameStateList(list);
	} else {
		fprintf(stderr, "Error while reading %s or %s.\n", inputFileNameTree, inputFileNameState);
		exit(1);
	}
	return 0;
}

