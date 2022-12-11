

#define HASH_ENDS_LINE //backward compatibility, deprecated


#include "FlowDomain.h"


#ifndef __DATE__
 #define __DATE__  "2021"
#endif


using namespace std;


void usageCNF(string exename, int detailed)  {

	std::cout<<"Pore Network Flow Simulation Code, pnflow  version  " << RELEASE_DATE << std::endl;
	if(detailed)  {
		std::cout<<"For more information, please visit Imperial College pore-scale modelling website:"<<std::endl
			 <<"https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling"<<std::endl
			 <<"or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk"<<endl;
			std::cout<<"\nUsage:"<< std::endl;
			std::cout<<" Run a flow simulation:"<< std::endl;
			std::cout<<"  "<<exename<< std::endl;
			std::cout<<"  "<<exename<<"  input_pnflow.dat\n"<< std::endl;
			std::cout<<" Generate a sample input file:"<< std::endl;
			std::cout<<"  "<<exename<<"  -g  input_pnflow.dat\n"<< std::endl;
			std::cout<<" Generate a sample input file with extra/experimental keywords:"<< std::endl;
			std::cout<<"  "<<exename<<"  -ge input_pnflow.dat\n"<< std::endl;
	}
}


void pnflowQD(InputFile input)  {

	if(input.getOr("invert_X", false))  { /// invert the flow direction
		int icycl=0;
		istringstream cycleBCData;
		while(input.giv("cycle"+_s(++icycl)+"_BC", cycleBCData))  {
			string str; double Dps[2]; char cc[] = {' ',' ',' ',' '};
			cycleBCData >>cc[0]>>cc[1]>>cc[2]>>cc[3] >>str>>Dps[0]>>Dps[1] ;
			ostringstream cycleDataNew;
			cycleDataNew <<cc[1]<<' '<<cc[0]<<' '<<cc[3]<<' '<<cc[2]<<"   "<<str<<"  "<<-Dps[0]<<"  "<<-Dps[1];
			input.set("cycle"+_s(icycl)+"_BC", cycleDataNew.str());
		}
		istringstream calcBox;
		if(input.giv("CALC_BOX", calcBox))  {
			double frac[2]={0.,0.};
			calcBox >>frac[0]>>frac[1];
			ostringstream boxDataNew;
			boxDataNew <<1.-frac[1]<<"  "<<1.-frac[0];
			input.set("CALC_BOX", boxDataNew.str());
		}
	}


	FlowDomain flowsim(input);


	int ic = 1;
	for(string cStr=input.kwrd("cycle"+_s(ic),0); !cStr.empty(); cStr=input.kwrd("cycle"+_s(++ic),0))  {
		const int dir=1-2*(ic&1);
		bool wantRelPerm = true;
		bool wantResIdx = true;
		double endSw(0.5*(1+dir)), endPc(-1e6*dir), delSw(0.05);
		bool entreL(true), entreR(false), exitL(true), exitR(true);


		double                             NcO_; // not used
		double                             NcW_; // not used
		char                               BC_type_[2]; // not used
		BC_type_[0]='D'; BC_type_[1]='P'; NcW_=1.; NcO_=1.;
		{
			istringstream instrim;
			instrim.str(cStr);
			instrim >> endSw >> endPc >> delSw ;

			wantRelPerm=readBoolOr("T",instrim);
			wantResIdx=readBoolOr("T",instrim);
			input.checkEndOfData(instrim,"cycle"+_s(ic));

			string BCStr=input.kwrd("cycle"+_s(ic)+"_BC",0);
			if (BCStr.empty()) BCStr = input.kwrd("BC",0);
			if (!BCStr.empty())  {
				instrim.clear();   instrim.str(BCStr);
				if (instrim.good())  entreL=readBoolOr("T",instrim);
				if (instrim.good())  entreR=readBoolOr("F",instrim);
				if (instrim.good())  exitL =readBoolOr("T",instrim);
				if (instrim.good())  exitR =readBoolOr("T",instrim);
				instrim >> BC_type_[0] >> BC_type_[1] >> NcW_ >> NcO_;
				input.checkEndOfData(instrim,"cycle"+_s(ic)+"_BC");
				input.Assert((BC_type_[0]=='D'||BC_type_[0]=='N')  && (BC_type_[1]=='V'||BC_type_[1]=='P'),"cycle"+_s(ic)+"_BC", " First argument after BC should be DP, DV, NP or NV");
			}
		}

		if(ic%2)
			flowsim.Drainage(endSw, endPc, delSw, wantRelPerm, wantResIdx, entreL, entreR, exitL, exitR);
		else
			flowsim.Imbibition(endSw, endPc, delSw, wantRelPerm, wantResIdx, entreL, entreR, exitL, exitR);
	}

	(cout<<"\n **** upscaled "+input.outputName()+", network "+input.kwrd("networkFile")+" ****\n").flush();
}


#ifndef MAIN
thread_local std::ofstream    outD;  //!< alias to mstream::dbgFile
int main(int argc, char *argv[])  {

 try
 {
	//int myid(0), nprocs(1);  MPI_Init(&argc, &argv);  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	cout<<"\n Network Model Code version 2 alpha, built: "<<__DATE__ << endl;

	string arg1;
	if (argc>1)  arg1 = argv[1];
	else  {  usageCNF(argv[0],1); cout<< "\nPlease input data file: ";   cin >> arg1;  }
	if(arg1=="-h" || arg1=="?")  {  usageCNF(argv[0],1); return 0; }

	InputFile inFile(arg1,true);
	inFile.renameKeys("ALTR_CONT_ANG","EQUIL_CON_ANG");

	pnflowQD(inFile);

	//MPI_Finalize();
 }
 catch (std::exception &exc) {  std::cerr << "\n\n Exception on processing: \n" << exc.what() << "Aborting! \n"  << endl;	return 1; }
 catch (...)                 {  std::cerr << "\n\n Unknown exception! \n Aborting! \n" << endl;	return 1; }
 return 0;
}

#endif
