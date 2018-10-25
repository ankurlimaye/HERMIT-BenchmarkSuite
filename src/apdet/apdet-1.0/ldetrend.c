/*

ldetrend.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
ldetrend: Locally detrend signals 
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

Reads 2 columns from stdin and locally detrends
by subtracting a least squares fitted line
over a sliding window 2*hwin+1 points wide.

Usage: ldetrend hwin

*/


#include <stdio.h>

#define MAXDAT 262144


main(argc, argv)
int argc;
char *argv[];
{
    int i, n, hwin, win;
    double a, b, sumx, sumy, sumxy, sumx2;
    static double x[MAXDAT+1], y[MAXDAT+1];

    if (argc != 2) {
        usage(argv[0]);
	exit(1);
    }

    if ((hwin = atoi(argv[1])) < 1) {
        fprintf(stderr, "%s : hwin must be greater than 0\n", argv[0]);
	exit(1);
    }
    win = 2*hwin+1;

    sumx = sumy = sumxy = sumx2 = 0.;

    for (n=0; n<MAXDAT && scanf("%lf %lf", &x[n], &y[n]) == 2; n++) {
        sumx += x[n];
        sumy += y[n];
        sumxy += x[n]*y[n];
	sumx2 += x[n]*x[n];
    }
    if (n == MAXDAT && scanf("%lf %lf", &x[i], &y[i]) == 2) {
        fprintf(stderr, "%s : input exceeds buffer length, ", argv[0]);
        fprintf(stderr, "truncating to %d points\n", MAXDAT);
    }

    if (win >= n) {

        b = (sumxy - sumx*sumy/n) / (sumx2 - sumx*sumx/n);
        a = sumy/n - b*sumx/n;

        for (i=0; i<n; i++)
            printf("%g %g\n", x[i], y[i] - (a + b*x[i]));
    }

    else {

        sumx = sumy = sumxy = sumx2 = 0.;

        for (i=0; i<win; i++) {
            sumx += x[i];
            sumy += y[i];
	    sumxy += x[i]*y[i];
	    sumx2 += x[i]*x[i];
        }

        b = (sumxy - sumx*sumy/win) / (sumx2 - sumx*sumx/win);
        a = sumy/win - b*sumx/win;

        for (i=0; i<=hwin; i++)
            printf("%g %g\n", x[i], y[i] - (a + b*x[i]));

        for (i=win ; i<n; i++) {
            sumx += x[i]-x[i-win];
	    sumy += y[i]-y[i-win];
	    sumxy += x[i]*y[i]-x[i-win]*y[i-win];
	    sumx2 += x[i]*x[i]-x[i-win]*x[i-win];

	    b = (sumxy - sumx*sumy/win) / (sumx2 - sumx*sumx/win);
	    a = sumy/win - b*sumx/win;

            printf("%g %g\n", x[i-hwin], y[i-hwin] - (a + b*x[i-hwin]));
	}

        for (i=i-hwin; i<n; i++)
            printf("%g %g\n", x[i], y[i] - (a + b*x[i]));
    }
}


usage(prog)
char *prog;
{
    fprintf(stderr, "Usage : %s hwin\n\n", prog);
    fprintf(stderr, " Reads 2 columns from stdin and locally detrends\n");
    fprintf(stderr, " by subtracting a least squares fitted line\n");
    fprintf(stderr, " over a sliding window 2*hwin+1 points wide.\n");
}
