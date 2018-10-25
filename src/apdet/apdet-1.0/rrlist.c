/* rrlist.c		Joe Mietus		Apr 1 2011 */

/*

rrlist.c		Joe Mietus		Dec 19 2001
_______________________________________________________________________________
rrlist: List R-R intervals
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

You may contact the author by e-mail (joe at physionet dot org) or postal mail
(MIT Room E25-505A, Cambridge, MA 02139 USA).  For updates to this software,
please visit PhysioNet (http://www.physionet.org/).
_______________________________________________________________________________

List R-R intervals

Usage: rrlist annotator tape [options]

options are :

    [-f start] : begin at time 'start'
    [-t end] : end at time 'end'
    [-l length] : output for duration 'length'
    [-h] : output time in hours in first column
    [-m] : output time in minutes in first column
    [-s] : output time in sec in first column
    [-a annotation] : list only intervals between consecutive annotations

*/

#include <stdio.h>
#include <wfdb/wfdb.h>
#include <wfdb/ecgmap.h>


main(argc, argv)	
int argc;
char *argv[];
{
    int i, annotation;
    int hflag, mflag, sflag, aflag;
    long start, end, length;
    double sps;
    struct WFDB_anninfo ai[1];
    struct WFDB_ann annot[2];

    if (argc < 3) {
        usage(argv[0]);
	exit(1);
    }

    hflag = mflag = sflag = aflag = 0;
    start = end = length = 0L;

    sps = sampfreq(argv[2]);
    ai[0].name = argv[1];
    ai[0].stat = WFDB_READ;
    if (annopen(argv[2], ai, 1) < 0)
       exit(2);

    for (i=2; ++i < argc && *argv[i] == '-'; ) {
        switch(argv[i][1]) {
	    case 'f': start = strtim(argv[++i]);
	              break;
	    case 't': end = strtim(argv[++i]);
	              break;
	    case 'l': length = strtim(argv[++i]);
	              break;
	    case 'h': hflag = 1;
	              break;
	    case 'm': mflag = 1;
	              break;
	    case 's': sflag = 1;
	              break;
	    case 'a': annotation = strann(argv[++i]);
	              aflag = 1;
	              break;
	    default:  usage(argv[0]);
	              exit(1);
        }
    }
    if (end == 0L && length != 0L)
        end = start + length;

    if (iannsettime(start) < 0)
        exit(2);

    start = 0L;

    while (getann(0, &annot[0]) >= 0) {
        if (isqrs(annot[0].anntyp))
        if ((!aflag && isqrs(annot[0].anntyp)) ||
            (aflag && annot[0].anntyp == annotation))
            break;
    }

    while (getann(0, &annot[1]) >= 0 && 
           (annot[1].time <= end || end == 0L)) {
	if ((!aflag && isqrs(annot[1].anntyp)) ||
            (aflag && annot[0].anntyp == annotation && 
             annot[1].anntyp == annotation)) {

	    if (hflag)
		printf("%.8f ", (annot[1].time-start)/(sps*3600));
	    else if (mflag)
		printf("%.6f ", (annot[1].time-start)/(sps*60));
	    else if (sflag)
		printf("%.3f ", (annot[1].time-start)/sps);
	    printf("%.3f %s\n", (annot[1].time - annot[0].time)/sps,
		                annstr(annot[1].anntyp));
	}
        if (isqrs(annot[1].anntyp)) {
	    annot[0] = annot[1];
	}
    }

    exit(0);
}


usage(prog)
char *prog;
{
    fprintf(stderr, "Usage : %s annotator tape [options]\n\n", prog);
    fprintf(stderr, " options are :\n");
    fprintf(stderr, " [-f start] : begin at time 'start'\n");
    fprintf(stderr, " [-t end] : end at time 'end'\n");
    fprintf(stderr, " [-l length] : output for duration 'length'\n");
    fprintf(stderr, " [-h] : output time in hours in first column\n");
    fprintf(stderr, " [-m] : output time in minutes in first column\n");
    fprintf(stderr, " [-s] : output time in seconds in first column\n");
    fprintf(stderr, " [-a annotation] : ");
    fprintf(stderr, " list only intervals between consecutive annotations\n");
}
