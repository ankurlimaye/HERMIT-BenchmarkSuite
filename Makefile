
curDir = $(shell pwd)
srcDir = $(curDir)/src
binDir = $(curDir)/bin

export srcDir
export binDir

kernels = activity aes apdet hrv imghist iradon kmeans lzw sqrs wabp

.PHONY : all
all : $(kernels)

% :
	mkdir -p $(binDir)
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)

.PHONY : clean
clean :
	rm -rf $(binDir)