/*

smooth.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
smooth: Smooth two column stdin with a sliding window average.
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

Smooth two column stdin with a sliding window average.

Usage : smooth window
*/

#include <stdlib.h>
#include <stdio.h>

#define MAXWIN    7200

main(argc, argv
)
int argc;
char *argv[];
{
int i, win;
static double x[MAXWIN], y[MAXWIN];
double sumx, sumy;

if (argc < 2) {
usage(argv[0]);
exit(1);
}

if ((
win = atoi(argv[1])
) <= 1) {
fprintf(stderr,
"%s : maximum window size must be greater than 1\n", argv[0]);
exit(2);
}
else if (win > MAXWIN) {
fprintf(stderr,
"%s : maximum window size is %d\n", argv[0], MAXWIN);
exit(2);
}

sumx = sumy = 0.0;

for (
i = 0;
i<win && scanf("%lf %lf", &x[i], &y[i]) == 2; i++) {
sumx += x[i];
sumy += y[i];
}

if (i < win) {
fprintf(stderr,
"%s : not enough points in input\n", argv[0]);
exit(2);
}

printf("%g %g\n", sumx/win, sumy/win);

i = 0;
sumx -= x[i];
sumy -= y[i];

while (scanf("%lf %lf", &x[i], &y[i]) == 2) {

sumx += x[i];
sumy += y[i];

printf("%g %g\n", sumx/win, sumy/win);

if (++i >= win)
i = 0;
sumx -= x[i];
sumy -= y[i];

}

}

usage(prog)
char *prog;
{
fprintf(stderr,
"Usage : %s window\n\n", prog);
fprintf(stderr,
" Smooth two column stdin with  sliding window average.\n");
}
