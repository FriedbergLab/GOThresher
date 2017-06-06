#!/bin/bash

#pip install debias

debias_prep -i go.obo
debias -cprot 100 -i goa_yeast.gaf goa_dicty.gaf -a C -WCTHRESHp 2 -recal 1 -odir output0
debias -i goa_yeast.gaf goa_dicty.gaf -a C P -PLTHRESHp 30 -e EXPEC IBA -single 1 -odir output1
debias -cattn 1000 -i goa_yeast.gaf goa_dicty.gaf -a C P -einv COMPEC -pref testing -selrefinv Reactome -odir output2

rm -rf ./output* ./data
#pip uninstall debias

