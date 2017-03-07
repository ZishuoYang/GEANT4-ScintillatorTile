#!/bin/csh
set SCRAM_ARCH = slc6_amd64_gcc491
set MYLOCATION = `pwd`
echo $MYLOCATION

set path = ($path /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/cmake/2.8.10-cms/bin)
set path = ($path /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/expat/2.0.1-cms/)
set path = ($path /data/users/eno/geant4.10.01.p02-install/share/Geant4-10.1.2/data)

cd /home/eno/CMSSW_7_4_0
cmsenv

set MYGEANT = /data/users/eno/geant4.10.01.p02-install
set MYGEANT2 = $MYGEANT/lib64/Geant4-10.1.2
set MYGEANT3 = /data/users/eno/geant4.10.01.p02

cd $MYGEANT/bin
pwd
source geant4.csh
cd $MYLOCATION

