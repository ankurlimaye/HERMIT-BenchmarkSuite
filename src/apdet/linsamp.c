/*

linsamp.c		Joe Mietus		Dec 19 2001
-------------------------------------------------------------------------------
linsamp: Resample stdin by linear interpolation, output evenly sampled data
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

Resample stdin by linear interpolation, output evenly sampled data

Usage : linsamp dt

dt : sampling interval

*/

#include <stdio.h>
#include <stdlib.h>

main(argc, argv
)
int argc;
char *argv[];
{

int i;
double dx, x0, x[2], y[2];
double intrpolt();

if (argc  != 2 ) {
usage(argv[0]);
exit(1);
}

if ((
dx = atof(argv[1])
) <= 0) {
fprintf(stderr,
"%s : dt must be greater than 0\n", argv[0]);
exit(1);
}

for (
i = 0;
i<2; i++)
if (scanf("%lf %lf", &x[i], &y[i]) < 2)
exit(1);

printf("%g %g\n", x[0], y[0]);
x0 = x[0] + dx;

for (;; ) {

while (x0 > x[1]) {
x[0] = x[1];
y[0] = y[1];

if (scanf("%lf %lf", &x[1], &y[1]) < 2 )
exit(0);       /* Normal exit from program */

}

printf("%g %g\n", x0,
intrpolt(x0, x, y
));

x0 +=
dx;
}
}


double intrpolt(x0, x, y)
    double x0, *x, *y;
{
  double a, b;

  b = (y[1] - y[0]) / (x[1] - x[0]);
  a = y[0] - b * x[0];

  return (b * x0 + a);
}

usage(prog)
char *prog;
{

fprintf(stderr,
"Usage : %s dt\n\n", prog);
fprintf(stderr,
" Resample stdin by linear interpolation\n");
fprintf(stderr,
" Output evenly sampled data\n\n");
fprintf(stderr,
" dt : sampling interval\n");

}
