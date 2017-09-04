ParSplit

Two software are provided
 - 'ParSplit'
	Identification of the most parsimonious split(s) of a phylogenetic tree with regards of character values at tips.
 - 'Fix'
	Fixing data to be read by ParSplit

type
	> make ParSplit
	in a console opened on the src directory for compiling the software.

Directory "src" contains the C sources of the software and Directory "data" contains two toy examples.
'tree.newick' is a file containing a phylogenetic tree in newick format
'statesA.csv' and 'statesB.csv' are two csv files containing two different set of character values for the tips of the tree in 'tree.newick'(lines of the form '<node ident> <state value>')

If the software 'ParSplit' and the example are in the same file hierarchy as in the repository typing in a console opened in the 'src' directory :

./ParSplit -t i -f 1 -s ../example/tree.newick ../example/statesA.csv resultsA.txt

will output two files:

	'resultsA.txt': a text file which contains the list of the parsimonious split of 'tree.newick' with regard to 'statesA.csv' with the corresponding cost and the no-split cost.

	'resultsA.pdf': is a pdf file resulting for option '-f 1' and displaying the figure of the corresponding splits

A complete description of the options allowed is given below.

------------
| ParSplit |
------------

--------------------------
REQUIREMENT

	The software needs the Cairo library.

--------------------------
COMPILING

	Just type
	> make ParSplit
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'ParSplit' identifies the most parsimonious split of a phylogenetic with regards of character values at tips.


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	ParSplit - Identification of the most parsimonious split(s) of a phylogenetic tree with regards of character values at tips.
	
SYNOPSIS
	ParSplit [OPTIONS] <input Tree File> <input State File> <output File>

DESCRIPTION
	Identify the most parsimonious split(s) of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) with regards of character values at tips contained in <input State File> (which must be in csv format) and output the results in text format in the file <output File> and, when one sets option '-f', in graphic format in the file <output File> with the corresponding extension (i.e .pdf, .png etc.)

	Options are
	-t <type>
		set the function applied on the branch length of the tree (which is 'inverse' by default): 
			-t i -> inverse
			-t u -> identity
	-z <input Tree File>
		output the tree in 'text' format in the console and save it in pdf format as 'tree.pdf' with the internal idents of the nodes (for debug purposes) next exit
	-n <number>
		set the max number of splits which will be displayed
	-f <number>
		set the graphic format of the output (option is required if one wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> LaTeX (psTricks)
			-f 6 -> LaTeX (TikZ)
	-c <r1> <g1> <b1> <r2> <g2> <b2>
		set the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])
	-h
		display help

--------------------------


-------
| Fix |
-------

--------------------------
REQUIREMENT

	None

--------------------------
COMPILING

	Just type
	> make Fix
	in a console opened on the directory containing the source files to build the binary.

--------------------------
DESCRIPTION

	'Fix' fix data to be used by ParSplit.


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	Fix - Fixing data to be used by ParSplit.
	
SYNOPSIS
	Fix [OPTIONS] <input Tree File> <input State File>

DESCRIPTION
	Make sure that all the tip idents of the phylogenetic tree contained in <input Tree File> (which must be in Newick format) are also idents of character values in the file <input State File> (which must be in csv format) and reciprocally. It outputs two files with '_fixed' containing the fixed tree and character values.

	Options are
	-h
		display help

--------------------------
