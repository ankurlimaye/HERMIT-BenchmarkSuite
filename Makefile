
curDir = $(shell pwd)
srcDir = $(curDir)/src
binDir = $(curDir)/bin
scriptDir = $(curDir)/scripts
inputsDir = $(curDir)/inputs

export srcDir
export binDir

kernels = activity apdet hrv imghist iradon kmeans sqrs wabp
suppKernels = aes lzw

.PHONY : all
all : $(kernels)

aes :
	mkdir -p $(binDir)/$@
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)/$@
	cp $(inputsDir)/$@/* $(binDir)/$@/
	cp $(scriptDir)/run-$@.sh $(binDir)/$@/

lzw :
	mkdir -p $(binDir)/$@
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)/$@
	cp $(inputsDir)/$@/* $(binDir)/$@/
	cp $(scriptDir)/run-$@.sh $(binDir)/$@/

$(kernels) : aes lzw
	mkdir -p $(binDir)/$@
	cd $(srcDir)/$@; make
	mv $(srcDir)/$@/$@ $(binDir)/$@
	cp $(binDir)/aes/aes $(binDir)/$@/
	cp $(binDir)/lzw/lzw $(binDir)/$@/
	cp $(inputsDir)/$@/* $(binDir)/$@/
	cp $(scriptDir)/run-$@.sh $(binDir)/$@/


.PHONY : clean
clean :
	rm -rf $(binDir)