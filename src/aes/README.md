
------------------------------
Install cryptopp
------------------------------

	cd cryptopp
	make static cryptest.exe
	./cryptest.exe v
	./cryptest.exe tv
	sudo make install PREFIX=/usr/local

------------------------------
Compile Example File
------------------------------

	g++ -o aesExampleBin aesExample.cpp -lcryptopp

------------------------------
Run perf
------------------------------

Run perf on 'aesExampleBin'

------------------------------
Source
------------------------------

	https://www.cryptopp.com/
