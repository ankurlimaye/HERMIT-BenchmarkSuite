/*
htmedfilt.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
htmedfilt: Separately median filter Hilbert amplitudes and frequencies
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

Separately median filter Hilbert amplitudes and frequencies

Usage : htmedfilt win

*/

#include <stdio.h>

#define MAXWIN 2000

static double time[MAXWIN], x[MAXWIN], y[MAXWIN], sx[MAXWIN], sy[MAXWIN];


main(argc, argv)
int argc;
char *argv[];
{
    int i, j, k, win, hwin;
    int dblcmp();

    if (argc != 2) {
        usage(argv[0]);
	exit(1);
    }

    if ((win = atoi(argv[1])) < 3) {
        usage(argv[0]);
        exit(1);
    }
    if (win > MAXWIN) {
	fprintf(stderr, "%s : maximum window size is MAXWIN\n", argv[0]);
	exit(2);
    }

    for (i=0; i<win && scanf("%lf %lf %lf", &time[i], &x[i], &y[i]) == 3; i++)
        ;

    if (i < win) {
        fprintf(stderr, "%s : not enough points in input\n", argv[0]);
	exit(2);
    }
    if (++i >= win)
	i = 0;
    j = hwin = win/2 -1;

    for (k=0; k<win; k++) {
        sx[k] = x[k];
        sy[k] = y[k];
    }

    qsort(sx, win, sizeof(double), dblcmp);
    qsort(sy, win, sizeof(double), dblcmp);
    printf("%g %g %g\n", time[j], sx[hwin], sy[hwin]);

    while (scanf("%lf %lf %lf", &time[i], &x[i], &y[i]) == 3) {

	if (++i >= win)
	    i = 0;
	if (++j >= win)
	    j = 0;

        for (k=0; k<win; k++) {
            sx[k] = x[k];
            sy[k] = y[k];
        }

        qsort(sx, win, sizeof(double), dblcmp);
        qsort(sy, win, sizeof(double), dblcmp);
	printf("%g %g %g\n", time[j], sx[hwin], sy[hwin]);

    }

}


int dblcmp(y1, y2)
double *y1, *y2;
{

    if (*y1 < *y2)
        return(-1);
    else if (*y1 > *y2)
        return(1);
    else
        return(0);

}


usage(prog)
char * prog;
{
	fprintf(stderr, "Usage : %s  win\n\n", prog);
	fprintf(stderr, " Separately median filter Hilbert amplitudes and\n");
	fprintf(stderr, " frequencies with a sliding window hwin points wide.\n");

}
