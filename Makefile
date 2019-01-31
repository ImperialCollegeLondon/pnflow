
# makefile for top directory

all:
	mkdir -p lib
	mkdir -p bin
	mkdir -p include
	(cd thirdparty && make -j$(nproc))
	(cd src && ./AllMake)


mgw:
	mkdir -p lib
	mkdir -p bin
	mkdir -p include
	(cd thirdparty && make -j$(nproc) OPT=.exe)
	(cd src && ./AllMakeMinGW)


clean:
	(cd thirdparty && make clean)
	(cd src && ./AllClean)

distclean:
	(cd thirdparty && make clean)
	(cd src && ./AllDistClean)
	@echo "\n**********************************************************\nWarning:\n"
	@echo " deleting log.make, bin/, include/, lib/, share/ and build/"
	@echo "\n Press Ctrl+c to stop "
	@echo "\n**********************************************************\n"
	sleep 6
	rm -rf   src/../bin  src/../include  src/../lib  src/../share  src/../build  log.make
