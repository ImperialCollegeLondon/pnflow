
import os, sys;  "msRoot" in os.environ or sys.exit("Retry after `source .../src/bashrc`");  from msrc import *; DbgMsg('*** main, ignore ***')



# pnflow does not support single-pore network
#Paint sphere  I J K  R  InV OutV(ignored in Paint)
nErrs=0
cmds=''
with open("spack4.mhd", 'w') as f1:	
	f1.write('DimSize = 40 40 40\n')
	f1.write('replaceRange 0 255 0\n')
	f1.write('reset dx  1   1   1\n\n')

	for kk in range(0,3,1):
		for jj in range(0,3,1):
			for ii in range(0,3,1):
				f1.write("Paint sphere  "+str(ii*20)+" "+str(jj*20)+" "+str(kk*20)+"     11   1\n")
	f1.write('\nreset dx 2e-6  2e-6  2e-6\n')



with open("CAs.tsv", 'w') as f1:
	f1.write('''
0
0
30
30
30
30
30
40
30
40
'''); f1.close();

with open("input_pnflow.dat", 'w') as f1:
	f1.write('''// -*- C -*-
NETWORK: F spack4;
TITLE:       spack4Net;
cycle1:    0.      5.0E+05     0.05          T            T;
cycle2:    1.     -5.0E+05     0.05          T            T;
cycle3:    0.      5.0E+05     0.05          T            T;
cycle1_BC:   T     F       T     T       DP   1.  1.;
cycle2_BC:   T     F       T     T       DP   1.  1.;
cycle3_BC:   T     F       T     T       DP   1.  1.;
CALC_BOX:  0. 1.;
INIT_CONT_ANG:   1   0   10  -0.2    -3.   rand   0.;
//# EQUIL_CON_ANG:   4   30   50   -0.2   -3.   rand   0.;
READ_ALTR_CA:  CAs.tsv pore 4 90;
WaterOilInterface 0.03;
//#         xR  nTheta   init   Drainage  Imbibition  corners   all steps
visualize:  1.     10      T        F          F         T         F ;
visuaLight                  T        T          T         F         T ;
writeStatistics T;
''');

runSh('.', "pnextract spack4.mhd");
runSh('.', "pnflow    input_pnflow.dat");
runSh('.', "voxelImageProcess   spack4.mhd dump_spack4.raw");
if fileFloatDiffersFrom("voxelImageProcess.log","totalPorosity:",0.32,0.05): nErrs+=1
if fileFloatDiffersFrom("pnflow.log","Absolute permeability:",1.5e-12,0.05): nErrs+=1


exit(nErrs)
