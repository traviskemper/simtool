#!/bin/bash

rm *.mod *.x *~

#ifort  main.f -o mdt.x
gfortran  main.f -o mdt.x
#f95  main.f -o mdt.x
#/programs/intel/Compiler/11.1/073/bin/intel64/ifort  main.f -o polmol.x
#/programs/intel/Compiler/11.1/073/bin/intel64/
#ifort -check all -traceback  main.f -o ch.x

