
curDir = $(shell pwd)
srcDir = $(curDir)/src
binDir = $(curDir)/bin
inputsDir = $(curDir)/inputs

export srcDir
export binDir

#kernels = activity aes apdet hrv imghist iradon kmeans lzw sqrs wabp
kernels = activity sqrs wabp imghist

.PHONY : all
all : $(kernels)

% :
	mkdir -p $(binDir)/$@
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)/$@
	cp $(inputsDir)/$@/* $(binDir)/$@/

.PHONY : clean
clean :
	rm -rf $(binDir)