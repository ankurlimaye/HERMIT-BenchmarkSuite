/*
filt.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
filt: moving average signal filter
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

filt filt hwin [options]

Reads 2 columns from stdin and filters outliers by
deleting those intervals outside of `filt' range of
the average within a window hwin' distance on either
side of the current interval.

options are :
  [-x min max] : exclude date outside min - max
*/

#include <stdio.h>
#include <stdlib.h>

#define MAXWIN 1024


main(argc, argv)
int argc;
char * argv[];
{
    int hwin, win, i, j, k;
    int xflag;
    double filt, filtmin, filtmax, min, max, mtmp, sum, av;
    static double x[MAXWIN+1], y[MAXWIN+1], lastx, lasty;
    FILE *fp;


    xflag = 0;
    sum = 0.0;

    if (argc < 3) {
        usage(argv[0]);
	exit(1);
    }

    i = 1;

    if ((filt = atof(argv[i])) <= 0) {
            fprintf(stderr, "%s: filt must be greater than 0\n",
                    argv[0]);
            exit(2);
    }
    if ((hwin = atoi(argv[++i])) <= 1) {
            fprintf(stderr, "%s: hwin must be integer greater than 1\n",
                    argv[0]);
            exit(2);
    }
    if ((win = 2*hwin) > MAXWIN) {
            fprintf(stderr, "%s: max hwin = %d\n", argv[0], MAXWIN/2);
            exit(2);
    }

    for ( ; ++i < argc && *argv[i] == '-'; ) {
        switch(argv[i][1]) {
            case 'x': if (i+3 > argc) {
	                  usage();
			  exit(1);
		      }
	              min = atof(argv[++i]);
	              max = atof(argv[++i]);
                      if (min > max) {
		          mtmp = min;
			  min = max;
			  max = mtmp;
                      }
	              xflag = 1;
                      break;
            default:  usage(argv[0]);
                      exit(1);
        }
    }

    for (i=0; i<=win; i++) {
        if (scanf("%lf %lf", &x[i], &y[i]) != 2) {
	    fprintf(stderr, "%s : insufficient data\n", argv[0]);
            exit(2);
	}
	if (xflag) {
            while (xflag && (y[i]<min || y[i]>max)) {
                if (scanf("%lf %lf", &x[i], &y[i]) != 2) {
	            fprintf(stderr, "%s : insufficient data\n", argv[0]);
	            exit(2);
	        }
            }
        }
    }

    i = 0;
    j = hwin;

    for (i=0; i<=win; i++) {
        sum += y[i];
    }
    sum -= y[j];
    av = sum/win;
    sum += y[j] - y[0];

    filtmax = filtmin = filt * av;

    if (y[j] <= av+filtmax && y[j] >= av-filtmin) {
        printf("%.9g %g\n", x[j], y[j]);
	lastx = x[j];
	lasty = y[j];
    }

    i = 0;

    while (scanf("%lf %lf", &x[i], &y[i]) == 2) {
	if (xflag) {
            while (y[i]<min || y[i]>max) {
                if (scanf("%lf %lf", &x[i], &y[i]) != 2) {
		    exit(0);
	        }
            }
	}

        if (++j > win)
            j = 0;

	sum += y[i] - y[j];
	av = sum/win;
	if (++i > win)
            i = 0;
	sum += y[j] - y[i];

        filtmax = filtmin = filt * av;

        if (y[j] <= av+filtmax && y[j] >= av-filtmin) {
            printf("%.9g %g\n", x[j], y[j]);
	    lastx = x[j];
	    lasty = y[j];
        }
    }

    exit(0);
}


usage(prog)
char *prog;
{
    fprintf(stderr, "Usage : %s filt hwin [options]\n\n", prog);
    fprintf(stderr, " Reads 2 columns from stdin and filters outliers by\n");
    fprintf(stderr, " deleting those intervals outside of `filt' range of\n");
    fprintf(stderr, " the average within a window hwin' distance on either\n");
    fprintf(stderr, " side of the current interval.\n\n");
    fprintf(stderr, " options are :\n");
    fprintf(stderr, " [-x min max] : exclude date outside min - max\n");
}
