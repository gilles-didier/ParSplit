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
#include <float.h>

#include "Floating.h"


double epsilon = 0.0000001;


int cmpDouble(double a, double b) {
	if(fabs(a - b) < epsilon) {
		return 0;
	} else {
		if(a>b)
			return 1;
		else
			return -1;
	}
}