
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

	tar xfvz wfdb.tar.gz
	cd wfdb-10.m.n
	./configure
	make install
	make check

------------------------------
Run perf
------------------------------

Run perf on 'sqrs' and 'wabp' in app folder

------------------------------
Source
------------------------------

PhysioNet:
	https://physionet.org/physiotools/wfdb-linux-quick-start.shtml
