
------------------------------
Install wfdb
------------------------------

1. Install pre-requisites

	gcc --version
	curl-config --version
	ls /usr/include/expat.h

If these commands work, go to step 2, else:

	apt-get install gcc libcurl4-openssl-dev libexpat1-dev

2. Configure, install and test package:

	cd ../sqrs_wabp
	tar xfvz wfdb.tar.gz
	cd wfdb-10.m.n
	./configure
	make install
	make check

------------------------------
Compile activity.c File
------------------------------

	cd activity
	g++ -o activity activity.c

------------------------------
Run perf
------------------------------

Run perf on 'tach -r record -a annotator | activity [-m] [len]'
	where record, annotator can be found in wfdb database; default len value = 600

------------------------------
Source
------------------------------

	https://www.physionet.org/physiotools/activity/
	https://www.physionet.org/physiotools/activity/activity.pdf
