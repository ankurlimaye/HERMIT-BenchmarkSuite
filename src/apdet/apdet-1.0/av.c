/*
av.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
av: Reads stdin and outputs n, av, and sd of data
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

Usage : av

Reads stdin and outputs n, av, and sd of data

*/


#include <stdio.h>
#include <math.h>

main(argc, argv)
int argc;
char *argv[];
{
	int cnt;
	double n, sum, sum2;
	double av, var, sd;

	if (argc != 1) {
            usage(argv[0]);
            exit(1);
	}

	cnt = 0;
	sum = sum2 = 0.;

	while (scanf("%lf", &n) == 1) {
	    cnt++;
	    sum += n;
	    sum2 += n*n;
	}

	av = sum/cnt;
	var = (sum2 - sum*sum/cnt) / (cnt-1);
	sd = sqrt(var);

	printf("%d %g %g\n", cnt, av, sd);
}


usage(prog)
char *prog;
{
    fprintf(stderr, "Usage : %s\n\n", prog);
    fprintf(stderr, " Reads stdin and outputs n, av, and sd of data\n");
}
