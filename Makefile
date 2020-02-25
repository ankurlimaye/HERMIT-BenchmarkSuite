
curDir = $(shell pwd)
srcDir = $(curDir)/src
binDir = $(curDir)/bin
scriptDir = $(curDir)/scripts
inputsDir = $(curDir)/inputs

export srcDir
export binDir

#kernels = activity aes apdet hrv imghist iradon kmeans lzw sqrs wabp
kernels = aes

.PHONY : all
all : $(kernels)

% :
	mkdir -p $(binDir)/$@
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)/$@
	cp $(inputsDir)/$@/* $(binDir)/$@/
	cp $(scriptDir)/run-$@.sh $(binDir)/$@/

.PHONY : clean
clean :
	rm -rf $(binDir)