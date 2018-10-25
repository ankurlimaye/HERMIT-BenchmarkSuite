
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
Install HRV
------------------------------

1. Unpack, and install apdet

	cd ../hrv
	tar xfvz HRV.src.tar.gz
	cd HRV
	make install

2. Download the BIDMC Congestive Heart Failure Database:
ea
	wget -r -np http://www.physionet.org/physiobank/database/chfdb/

3. Extract the database in ../apdet/apdet-1.0/chfdb folder

------------------------------
Run perf
------------------------------

Run perf on 'get_hrv -f "0.2 20 -x 0.4 2.0" -p "20 50" chfdb/chf03 ecg'

------------------------------
Source
------------------------------

PhysioNet:
	https://physionet.org/physiotools/wfdb-linux-quick-start.shtml

HRV:
	https://www.physionet.org/tutorials/hrv-toolkit/	
	
BIDMC Congestive Heart Failure Database:
	Baim DS, Colucci WS, Monrad ES, Smith HS, Wright RF, Lanoue A, Gauthier DF, Ransil BJ, Grossman W, Braunwald E. Survival of patients with severe congestive heart failure treated with oral milrinone. J American College of Cardiology 1986 Mar; 7(3):661-670.
	https://www.physionet.org/physiobank/database/chfdb/
