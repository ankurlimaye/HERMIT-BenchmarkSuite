
------------------------------
Install OpenCV
------------------------------

1.	Pre-requisities

	sudo apt-get install build-essential
	sudo apt-get install cmake git libgtk2.0-dev pkg-config libavcodec-dev libavformat-dev libswscale-dev
	sudo apt-get install python-dev python-numpy libtbb2 libtbb-dev libjpeg-dev libpng-dev libtiff-dev libjasper-dev libdc1394-22-dev

2.	Getting the OpenCV from Git Repository

	cd ~/<_my_working_directory>
	git clone https://github.com/opencv/opencv.git
	git clone https://github.com/opencv/opencv_contrib.git

3.	Building OpenCV from source

	cd ~/opencv
	mkdir build
	cd build
	cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX=/usr/local ..
	make -j7

------------------------------
Install imghist
------------------------------

	cd imghist
	g++ imghist.cpp -o imghist.o

------------------------------
Run perf
------------------------------

------------------------------
Source
------------------------------

Bibek Subedi
	http://www.programming-techniques.com/2013/01/histogram-equalization-using-c-image.html
