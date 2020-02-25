// TODO: Notice stuff, etc.
/* file: activity.c	G. Moody	2 April 1992

-------------------------------------------------------------------------------
activity: Estimate activity level from heart rate signal
Copyright (C) 2002 George B. Moody

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

You may contact the author by e-mail (george@mit.edu) or postal mail
(MIT Room E25-505A, Cambridge, MA 02139 USA).  For updates to this software,
please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________

This program derives an "activity index" from a time series of instantaneous
heart rate measurements, such as can be produced by 'tach'.  'tach' is included
in the WFDB Software Package;  for details, see
    http://www.physionet.org/physiotools/wag/tach-1.htm

For example:
    tach -r RECORD -a ANNOTATOR -Vs > fileName.ts
    activity -i fileName.ts

Each value of the activity index is derived from 'len' values in the input
heart rate time series;  by default, 5 minutes of input data are used to
produce each output value.  The input windows overlap by 50%, so that the
interval between output values is half of that specified by 'len', or 2.5
minutes by default.  Other values of 'len' can be specified on the command
line, as in:
    tach -r RECORD -a ANNOTATOR | activity 240
which would yield outputs at 1-minute intervals, based on 2-minute windows.

Use activity's '-m' option to find and output only the interval for which
the activity index is minimum.

The activity index is based on mean heart rate, total power of the observed
heart rate signal, and a heart rate stationarity index.  For details, see
"ECG-based indices of physical activity", pp. 403-406, Computers in Cardiology
1992.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define DEFLEN  600  /* 5 minutes at 2 samples/sec */

int main(int argc, char *argv[]) {

  char buf[80];
  double *t, *hr, activity, meanhr, meanhr0, meanhr1, p, tpower, stationarity;
  double acmin = -1.0, hrmin, stmin, tpmin, tmin0, tmin1;
  int i = 0, len = DEFLEN;
  long tt = 0L;
  FILE *in_file, *out_file;
  bool mflag = false, inputFile = false;

  for (int j = 1 ; j < argc ; ++j) {
    if (strcmp(argv[j], "-r") == 0) {
      inputFile = true;
      if ((in_file = fopen(argv[(j + 1)], "r")) == NULL) {
        inputFile = false;
      }
    } else if (strcmp(argv[j], "-len") == 0) {
      len = atoi(argv[(j + 1)]);
    } else if (strcmp(argv[j], "-m") == 0) {
      mflag = true;
    }
  }

  if(! inputFile) {
    printf("Incorrect input file! Using test input: \"test-100.ts\" \n");
    in_file = fopen("test-100.ts", "r");
  }

  if ((t = (double *)malloc(len * sizeof(double))) == NULL || (hr = (double *)malloc(len * sizeof(double))) == NULL) {
    printf("%s: Insufficient memory\n", argv[0]);
    return -1;
  }

  out_file = fopen("activity-out.txt", "w");
  while (fgets(buf, 80, in_file)) {
    int n = sscanf(buf, "%lf %lf", &t[i], &hr[i]);

    if (n == 0) {
      continue;	/* skip empty lines in input */
    }

    if (n == 1) {
      hr[i] = t[i];
      t[i] = tt / 2.0;
    }

    ++tt;

    if (++i >= len) {
      /* hr buffer full -- emit output and reset */
      meanhr0 = meanhr1 = tpower = 0.0;

      for (i = 0 ; i < len / 2 ; i++) {
        meanhr0 += hr[i];
      }

      for ( ; i < len ; i++) {
        meanhr1 += hr[i];
      }

      meanhr0 /= len/2;
      meanhr1 /= len/2;
      meanhr = (meanhr0 + meanhr1) / 2;
      stationarity = meanhr0 - meanhr1;

      if (stationarity < 0) {
        stationarity = -stationarity;
      }

      for (i = 0 ; i < len ; i++) {
        p = hr[i] - meanhr;
        tpower += p*p;
      }

      tpower /= len;

      if (tpower > 100.) {
        tpower = 100.;
      }

      activity = sqrt(((meanhr - 40.) * (meanhr - 40.)) + (10. * stationarity * stationarity) + (100. * tpower));

      if (meanhr < 25.) {
        /* penalty for unbelievably low heart rates */
        activity += 25. - meanhr;
      }

      if (! mflag) {
        fprintf(out_file, "%g \t %g \t %g \t %g \t %g \n", t[(len / 4) - 1], meanhr, tpower, stationarity, activity);
        fprintf(out_file, "%g \t %g \t %g \t %g \t %g \n", t[(3 * (len / 4)) - 1], meanhr, tpower, stationarity, activity);
      } else if ((activity < acmin) || (acmin < 0.)) {
        acmin = activity;
        hrmin = meanhr;
        stmin = stationarity;
        tpmin = tpower;
        tmin0 = t[0];
        tmin1 = t[len-1];
      }

      for (i = 0 ; i < (len / 2) ; i++) {
        hr[i] = hr[i + (len / 2)];
        t[i] = t[i + (len / 2)];
      }
    }
  }

  if (mflag) {
    fprintf(out_file, "%g \t %g \t %g \t %g \t %g \t %g \n", tmin0, tmin1, hrmin, tpmin, stmin, acmin);
  }

  fclose(in_file);
  fclose(out_file);

  return 0;
}