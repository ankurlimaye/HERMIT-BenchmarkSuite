/*
detruns.c		Joe Mietus		Dec 19 2001 

-------------------------------------------------------------------------------
detruns: List detection runs
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

List detection runs

Usage: detruns incr win minlen

Reads one column of data (hour)
and outputs start and end times of run detection
windows and sum total detection time

incr : output increment in hh:mm:ss
win : window length in hh:mm:ss
minlen : min length of run in hh:mm:ss

*/

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  char hour[9], *timstr();
  int i, runflag, runflag0;
  long incr, win, minlen;
  long time, lasttime, runstart, runstart0, runend0, outtime, sum, strtim();

  sum = 0;

  if (argc < 4) {
    usage(argv[0]);
    exit(1);
  }

  incr = strtim(argv[1]);

  if (incr < 1) {
    fprintf(stderr, "%s : incr must be greater than 0\n", argv[0]);
    exit(2);
  }

  win = strtim(argv[2]);

  if (win < 1) {
    fprintf(stderr, "%s : win must be greater than 0\n", argv[0]);
    exit(2);
  }

  minlen = strtim(argv[3]);

  if (minlen < 1) {
    fprintf(stderr, "%s : minlen must be greater than 0\n", argv[0]);
    exit(2);
  }

  if (minlen < win) {
    minlen = win;
  }

  runflag = runflag0 = 0;

  if (scanf("%s", hour) != 1) {
    fprintf(stderr, "%s: incorrectly formatted data : ", argv[0]);
    exit(2);
  }

  runstart = lasttime = strtim(hour);

  while (scanf("%s", hour) == 1) {
    time = strtim(hour);

    if (time - lasttime != incr) {
      if (lasttime - runstart + win >= minlen) {
        if (!runflag0) {
          runflag0 = 1;
          runstart0 = runstart;
          runend0 = lasttime;
        } else if (minlen > win && runstart <= runend0 + win) {
          runend0 = lasttime;
        } else {
          printf("%s ", timstr(runstart0));
          printf("%s\n", timstr(runend0 + win));
          sum += runend0 - runstart0 + win;

          runstart0 = runstart;
          runend0 = lasttime;
        }
      }
      runstart = time;
    }

    lasttime = time;
  }

  if (lasttime - runstart + win >= minlen) {
    runflag = 1;
  }

  if (runflag0) {
    if (runflag && minlen > win && runstart <= runend0 + win) {
      runend0 = lasttime;
      runflag = 0;
    }
    printf("%s ", timstr(runstart0));
    printf("%s\n", timstr(runend0 + win));
    sum += runend0 - runstart0 + win;
  }

  if (runflag) {
    printf("%s ", timstr(runstart));
    printf("%s\n", timstr(lasttime + win));
    sum += lasttime - runstart + win;
  }

  printf("%s\n", timstr(sum));

  return 0;
}

long strtim(char *buf) {        /* convert string in [[HH:]MM:]SS format to seconds */
  long x, y, z;

  switch (sscanf(buf, "%ld:%ld:%ld", &x, &y, &z)) {
    case 1: return (x);
    case 2: return (60 * x + y);
    case 3: return (3600 * x + 60 * y + z);
    default: return (-1L);
  }
}

char *timstr(long time) {
  /* convert seconds to [HH:]MM:SS format */
  int hours, minutes, seconds;
  static char buf[9];

  hours = time / 3600L;
  time -= (long) hours * 3600;
  minutes = time / 60;
  seconds = time - minutes * 60;
  sprintf(buf, "%02d:%02d:%02d", hours, minutes, seconds);
  return (buf);
}

void usage(char *prog) {
  fprintf(stderr, "Usage : %s incr win minlen\n\n", prog);
  fprintf(stderr, " Reads one column of data (hour)\n");
  fprintf(stderr, " and outputs start and end times of run detection\n");
  fprintf(stderr, " windows and sum total detection time.\n\n");
  fprintf(stderr, " incr : output increment in hh:mm:ss\n");
  fprintf(stderr, " win : window length in hh:mm:ss\n");
  fprintf(stderr, " minlen : min length of run in hh:mm:ss\n");
}
