/*
htavsd.c		Joe Mietus		Dec 19 2001

-------------------------------------------------------------------------------
htavsd: compute Hilbert transform averages, standard deviations and threshold 
crossings for amplitude and frequency over moving window
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


get Hilbert transform averages, standard deviations and threshold crossings
for amplitude and frequency over moving window of length `win'

Usage : htavsd incr win ampthres

Reads stdin of time in secs, amp and freq
outputs time, av, sd and % time within threshold limits
for both amp and freq 

incr : output for every 'incr' time steps
win : window length
ampthres : min amp

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define MAXFREQ 0.06

int main(int argc, char *argv[]) {
  int i, j, k, win, start, incr, ydet, zdet;
  double ampthres;
  double sumy, sumz, sumyy, sumzz, avy, avz, sdy, sdz;
  double *x, *y, *z;
  double atof();

  if (argc < 4) {
    usage(argv[0]);
    exit(1);
  }

  if ((incr = strtim(argv[1])) <= 0) {
    fprintf(stderr, "%0 : incr must be greater than 0\n", argv[0]);
    exit(1);
  }

  if ((win = strtim(argv[2])) <= 0) {
    fprintf(stderr, "%0 : win must be greater than 0\n", argv[0]);
    exit(1);
  }

  if ((ampthres = atof(argv[3])) <= 0) {
    fprintf(stderr, "%0 : ampthres must be greater than 0\n", argv[0]);
    exit(1);
  }

  if (incr > win) {
    fprintf(stderr, "%s : incr must be <= win\n", argv[0]);
    exit(1);
  }

  if ((x = (double *) calloc(win, sizeof(double))) == NULL || (y = (double *) calloc(win, sizeof(double))) == NULL || (z = (double *) calloc(win, sizeof(double))) == NULL) {
    fprintf(stderr, "%s : insufficient memory\n", argv[0]);
    exit(1);
  }

  ydet = zdet = 0;
  sumy = sumz = sumyy = sumzz = 0.;

  if (scanf("%lf %lf %lf", &x[0], &y[0], &z[0]) != 3) {
    fprintf(stderr, "%s : no data read\n", argv[0]);
    exit(2);
  }

  start = x[0];

  sumy = y[0];
  sumz = z[0];
  sumyy = y[0] * y[0];
  sumzz = z[0] * z[0];

  if (y[0] >= ampthres) {
    ydet++;
  }

  if (z[0] <= MAXFREQ) {
    zdet++;
  }

  for (i = 1 ; i < win && scanf("%lf %lf %lf", &x[i], &y[i], &z[i]) == 3 ; i++) {
    sumy += y[i];
    sumz += z[i];
    sumyy += y[i] * y[i];
    sumzz += z[i] * z[i];

    if (y[i] >= ampthres) {
      ydet++;
    }

    if (z[i] <= MAXFREQ) {
      zdet++;
    }
  }

  if (i < win) {
    fprintf(stderr, "not enough points in input\n");
    exit(2);
  }

  avy = sumy / win;
  avz = sumz / win;
  sdy = sqrt((sumyy - sumy * sumy / win) / (win - 1));
  sdz = sqrt((sumzz - sumz * sumz / win) / (win - 1));

  printf("%02d:%02d:%02d ", start / 3600, (start % 3600) / 60, start % 60);
  printf("%f %f %f %f %f %f\n", avy, sdy, ((double) ydet) / win, avz, sdz, ((double) zdet) / win);

  if (incr == win) {
    sumy = sumz = sumyy = sumzz = 0.0;
    ydet = zdet = 0;
  } else {
    for (j = 0 ; j < incr ; j++) {
      if (j >= win) {
        j = 0;
      }

      sumy -= y[j];
      sumz -= z[j];
      sumyy -= y[j] * y[j];
      sumzz -= z[j] * z[j];
      if (y[j] >= ampthres) {
        ydet--;
      }
      if (z[j] <= MAXFREQ) {
        zdet--;
      }
    }
  }

  i = 0;
  start += incr;

  while (scanf("%lf %lf %lf", &x[i], &y[i], &z[i]) == 3) {

    while (x[i] < start) {
      if (scanf("%lf %lf %lf", &x[i], &y[i], &z[i]) != 3) {
        exit(0);
      }
    }

    sumy += y[i];
    sumz += z[i];
    sumyy += y[i] * y[i];
    sumzz += z[i] * z[i];
    if (y[i] >= ampthres)
      ydet++;
    if (z[i] <= MAXFREQ)
      zdet++;

    if (++i >= win)
      i = 0;

    for (
        j = 1 ;
        j < incr &&
            j < win ;
        i++, j++) {

      if (i >= win)
        i = 0;

      if (scanf("%lf %lf %lf", &x[i], &y[i], &z[i]) != 3)
        exit(0);

      sumy += y[i];
      sumz += z[i];
      sumyy += y[i] * y[i];
      sumzz += z[i] * z[i];

      if (y[i] >= ampthres)
        ydet++;
      if (z[i] <= MAXFREQ)
        zdet++;
    }

    if (i >= win)
      i = 0;

    avy = sumy / win;
    avz = sumz / win;
    sdy = sqrt((sumyy - sumy * sumy / win) / (win - 1));
    sdz = sqrt((sumzz - sumz * sumz / win) / (win - 1));

    printf("%02d:%02d:%02d ", start / 3600, (start % 3600) / 60, start % 60);
    printf("%f %f %f %f %f %f\n",
           avy, sdy, ((double) ydet) / win,
           avz, sdz, ((double) zdet) / win);

    if (incr == win) {
      sumy = sumz = sumyy = sumzz = 0.0;
      ydet = zdet = 0;
    } else {
      for (
          j = 0, k = j + i ;
          j < incr ;
          j++, k++) {
        if (k >= win)
          k = 0;
        sumy -= y[k];
        sumz -= z[k];
        sumyy -= y[k] * y[k];
        sumzz -= z[k] * z[k];
        if (y[k] >= ampthres)
          ydet--;
        if (z[k] <= MAXFREQ)
          zdet--;
      }
    }

    start +=
        incr;

  }
}

int strtim(buf)        /* convert string in [[HH:]MM:]SS format to seconds */
    char *buf;
{
  int x, y, z;

  switch (sscanf(buf, "%d:%d:%d", &x, &y, &z)) {
    case 1: return (x);
    case 2: return (60 * x + y);
    case 3: return (3600 * x + 60 * y + z);
    default: return (-1);
  }
}

usage(prog)
char *prog;
{
fprintf(stderr,
"Usage : %s incr win ampthres\n\n", prog);
fprintf(stderr,
" Reads stdin of time in secs, amp and freq\n");
fprintf(stderr,
" outputs time, av, sd, and");
fprintf(stderr,
" % time within threshold limits for both amp and freq\n\n");
fprintf(stderr,
"incr] : output for every 'incr' time steps\n");
fprintf(stderr,
"win : window length\n");
fprintf(stderr,
"ampthres : minumum amplitude threshold\n");
}
