
------------------------------
Install wfdb (as in sqrs_Wabp)
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
Install apdet
------------------------------

1. Unpack, and install apdet

	cd ../apdet
	tar xfvz apdet.src.tar.gz
	cd apdet-1.0
	make all
	make install

2. Download the Sleep-apnea Database

	wget -r -np http://www.physionet.org/physiobank/database/apnea-ecg/

3. Extract the database in ../apdet/apdet-1.0/apnea-ecg folder

------------------------------
Run perf
------------------------------

Run perf on 'get_apdet apnea-ecg/a03'

------------------------------
Source
------------------------------

PhysioNet:
	https://physionet.org/physiotools/wfdb-linux-quick-start.shtml

Apdet:
	https://www.physionet.org/physiotools/apdet/

Sleep apnea database:
	T Penzel, GB Moody, RG Mark, AL Goldberger, JH Peter. The Apnea-ECG Database. Computers in Cardiology 2000;27:255-258.
	https://www.physionet.org/physiobank/database/apnea-ecg/
	
	
