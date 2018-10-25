/*

mm.c		Joe Mietus		Dec 19 2001
_______________________________________________________________________________
mm: Print minimum and maximum of data
Copyright (C) 2002 Joe Mietus

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place - Suite 330, Boston, MA 02111-1307, USA.

You may contact the author by e-mail (joe@physionet.org) or postal mail
(MIT Room E25-505A, Cambridge, MA 02139 USA).  For updates to this software,
please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________

Print minimum and maximum of data

Usage : mm

*/

#include <stdio.h>

main(argc, argv)
int argc;
char *argv[];
{
	double n, min, max;

	if (argc != 1) {
	    fprintf(stderr, "Usage : %s\n", argv[0]);
	    fprintf(stderr, " Print minimum and maximum of data\n");
	    exit(1);
	}

	scanf("%lf", &min);
        max = min;
	while (scanf("%lf", &n) == 1) {
		if(n < min)
			min = n;
		if(n > max)
			max = n;
	}
	printf("%g %g\n", min, max);
}
