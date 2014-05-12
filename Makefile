UNIXPLATFORM=$(shell uname -s)
CPBA_ROOT=/work/mthomas/cpba
PAB_ROOT=/opt/cpba_nels

all: library executables

ifeq ($(strip $(UNIXPLATFORM)),Darwin)
	MAKE=make
endif

ifeq ($(strip $(UNIXPLATFORM)),SunOS)
	MAKE=gmake
endif

library:
	$(MAKE) -C src/core
	@echo "Library build complete"

executables:
	$(MAKE) -C src/exe
	@echo "Executable build complete"
	
clean:
	$(MAKE) -C src/core clean
	$(MAKE) -C src/exe clean

sync-data:
	scp -r InputData pinnbeta@pinnbeta:cpba/trunk

sync:
	scp Makefile pinnbeta@pinnbeta:cpba/trunk
	scp src/core/*.h src/core/*.cpp src/core/Makefile pinnbeta@pinnbeta:cpba/trunk/src/core
	scp src/exe/*.cpp src/exe/Makefile pinnbeta@pinnbeta:cpba/trunk/src/exe
	scp src/python/*.py pinnbeta@pinnbeta:cpba/trunk/src/python
	scp src/pinnacle/*.Script pinnbeta@pinnbeta:cpba/trunk/src/pinnacle

sync-pabulum:
	ssh root@pabulum mkdir -p ${PAB_ROOT}/src/core
	ssh root@pabulum mkdir -p ${PAB_ROOT}/src/exe
	ssh root@pabulum mkdir -p ${PAB_ROOT}/src/python
	ssh root@pabulum mkdir -p ${PAB_ROOT}/src/pinnacle
	ssh root@pabulum mkdir -p ${PAB_ROOT}/src/examples
	
	scp Makefile root@pabulum:${PAB_ROOT}
	scp -r InputData root@pabulum:${PAB_ROOT}
	scp src/core/*.h src/core/*.cpp src/core/Makefile root@pabulum:${PAB_ROOT}/src/core
	scp src/exe/*.cpp src/exe/Makefile root@pabulum:${PAB_ROOT}/src/exe
	scp src/python/*.py root@pabulum:${PAB_ROOT}/src/python
	scp src/pinnacle/*.Script root@pabulum:${PAB_ROOT}/src/pinnacle
	scp src/examples/*.sh root@pabulum:${PAB_ROOT}/src/examples


sync-tezpur:
	ssh mthomas@tezpur.hpc.lsu.edu mkdir -p ${CPBA_ROOT}/src/core
	ssh mthomas@tezpur.hpc.lsu.edu mkdir -p ${CPBA_ROOT}/src/exe
	ssh mthomas@tezpur.hpc.lsu.edu mkdir -p ${CPBA_ROOT}/src/python
	ssh mthomas@tezpur.hpc.lsu.edu mkdir -p ${CPBA_ROOT}/src/pinnacle
	ssh mthomas@tezpur.hpc.lsu.edu mkdir -p ${CPBA_ROOT}/src/examples
	
	scp Makefile mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}
	scp -r InputData mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}
	scp src/core/*.h src/core/*.cpp src/core/Makefile mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}/src/core
	scp src/exe/*.cpp src/exe/Makefile mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}/src/exe
	scp src/python/*.py mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}/src/python
	scp src/pinnacle/*.Script mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}/src/pinnacle
	scp src/examples/*.sh mthomas@tezpur.hpc.lsu.edu:${CPBA_ROOT}/src/examples
	
