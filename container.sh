#!/bin/bash

#source runimage.sh
source setup.sh
#make
#./bin/RunppTestAna -i /gpfs01/star/pwg/elayavalli/ppRun12Datapicos/ppJP2Run12/sum8.root -intype pico -c JetTree -trig ppJP2 -o Results/sum8_20pt30.root -N -1 -pj 20 30 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999 -geantnum 0
./bin/RunppTestAna -i /gpfs01/star/pwg/svomich/Jets_pp_2012_Raghav/ppRun12Datapicos/ppJP2Run12/$1 -intype pico -c JetTree -trig ppJP2 -o /gpfs01/star/pwg/svomich/Jets_pp_2012_Raghav/workDir/out/out_$1 -N -1 -pj 10 40 -pc 0.2 30 -lja antikt -ec 1 -R 0.4 -hadcorr 0.9999999 -geantnum 0
