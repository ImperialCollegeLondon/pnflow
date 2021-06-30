#ifdef WIN32
#pragma warning(disable:4786)
#endif

#define HASH_ENDS_LINE //backward compatibility, deprecated


#include "FlowDomain.h"
//#include "mpi.h"

#ifndef __DATE__
 #define __DATE__  "2016"
#endif 


using namespace std;


void usageCNF(string exename, int detailed)
{

	std::cout<<"Pore Network Flow Simulation Code, cnflow  version  " << RELEASE_DATE << std::endl;
	if(detailed)
	{
		std::cout<<"For more information, please visit Imperial College pore-scale modelling website:"<<std::endl
			 <<"https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling"<<std::endl
			 <<"or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk"<<endl;
			std::cout<<"\nUsage:"<< std::endl;
			std::cout<<" Run a flow simulation:"<< std::endl;
			std::cout<<"  "<<exename<< std::endl;
			std::cout<<"  "<<exename<<"  input_cnflow.dat\n"<< std::endl;
			std::cout<<" Generate a sample input file:"<< std::endl;
			std::cout<<"  "<<exename<<"  -g  input_cnflow.dat\n"<< std::endl;
			std::cout<<" Generate a sample input file with extra/experimental keywords:"<< std::endl;
			std::cout<<"  "<<exename<<"  -ge input_cnflow.dat\n"<< std::endl;
	}
}


int cnflowQD(InputFile input, unsigned randShift, bool invert)
{ // sync_cn_gnflowQD
	//char dir='X';
	//if(input.lookup("direction",dir) && dir!='X' && input.outputName().back()!='Y'&& input.outputName().back()!='Z'))
		//input.setTitle(input.outputName().append(1,toupper(dir)));

	if(invert^input.getOr(false,"invert_X"))
	{ /// invert the flow direction
		int icycl=0;
		istringstream cycleBCData;
		while(input.getData(cycleBCData, "cycle"+_s(++icycl)+"_BC"))
		{
			{
				string str; double Dps[2]; char cc[] = {' ',' ',' ',' '};
				cycleBCData >>cc[0]>>cc[1]>>cc[2]>>cc[3] >>str>>Dps[0]>>Dps[1] ;
				ostringstream cycleDataNew;
				cycleDataNew <<cc[1]<<' '<<cc[0]<<' '<<cc[3]<<' '<<cc[2]<<"   "<<str<<"  "<<-Dps[0]<<"  "<<-Dps[1];
				input.setKeyword("cycle"+_s(icycl)+"_BC", cycleDataNew.str());
			}
		}
		istringstream calcBox;
		if(input.getData(calcBox, "CALC_BOX"))
		{
			double frac[2]={0.0,0.0};
			calcBox >>frac[0]>>frac[1];
			ostringstream boxDataNew;
			boxDataNew <<1.0-frac[1]<<"  "<<1.0-frac[0];
			input.setKeyword("CALC_BOX", boxDataNew.str());
		}
	}


	 unsigned int randSeed = input.getOr((unsigned)time( NULL ),"RAND_SEED")+randShift; ///. to generate the same randoms for all scanning runs don't move to inside of scanning loop
	//input.setKeyword("RAND_SEED", _s(randSeed));


	FlowDomain flowsim(input,randSeed);


	int ic = 1;
	for(string cStr=input.keyvals("cycle"+_s(ic),0); !cStr.empty(); cStr=input.keyvals("cycle"+_s(++ic),0))
	{
		const double dir=1.0-2.0*(ic%2);
		bool wantRelPerm = true;
		bool wantResIdx = true;
		double endSw(0.5*(1.0+dir)), endPc(-1.0e06*dir), delSw(0.05);
		bool entreL(true), entreR(false), exitL(true), exitR(true);


		double                             NcO_; // not used
		double                             NcW_; // not used
		char                               BC_type_[2]; // not used
		BC_type_[0]='D'; BC_type_[1]='P'; NcW_=1.0; NcO_=1.0;
		{
			istringstream instrim;
			instrim.str(cStr);
			instrim >> endSw >> endPc >> delSw ;

			wantRelPerm=readBoolOr("T",instrim);
			wantResIdx=readBoolOr("T",instrim);
			input.checkEndOfData(instrim,"cycle"+_s(ic));

			string BCStr=input.keyvals("cycle"+_s(ic)+"_BC",0);
			if (BCStr.empty()) BCStr = input.keyvals("BC",0);
			if (!BCStr.empty())
			{
				instrim.clear();   instrim.str(BCStr);
				if (instrim.good())  entreL=readBoolOr("T",instrim);
				if (instrim.good())  entreR=readBoolOr("F",instrim);
				if (instrim.good())  exitL =readBoolOr("T",instrim);
				if (instrim.good())  exitR =readBoolOr("T",instrim);
				instrim >> BC_type_[0] >> BC_type_[1] >> NcW_ >> NcO_;
				input.checkEndOfData(instrim,"cycle"+_s(ic)+"_BC");
				input.Assert((BC_type_[0]=='D'||BC_type_[0]=='N')  && (BC_type_[1]=='V'||BC_type_[1]=='P'),"cycle"+_s(ic)+"_BC", " First argument after BC should be DP, DV, NP or NV");
			}

			//out_<<"Stop@ Sw: "<<endSw<<"  Pc: "<<endPc<<endl
				 //<<"BC Entry:" <<entreL <<" "<< entreR<<" Exit: " <<exitL <<" "<< exitR<<endl
				 //<<BC_type_[0]<<BC_type_[1]<<"  o: "<< NcW_ <<" w: "<< NcO_
				 //<<endl;
		}
		//endPc /= comn_.watoil().sigma();

		if(ic%2)
			flowsim.Drainage(endSw, endPc, delSw, wantRelPerm, wantResIdx, entreL, entreR, exitL, exitR);
		else
			flowsim.Imbibition(endSw, endPc, delSw, wantRelPerm, wantResIdx, entreL, entreR, exitL, exitR);
	}

	(cout<<"\n **** upscaled "+input.outputName()+", network "+input.keyvals("networkFile")+" ****\n").flush();
	return 0;
}


#ifndef MAIN
#include "globals.cpp"
thread_local OnDemandStream    outD;  ///. alias to OnDemandStream::dbgFile  //thread_local
int main(int argc, char *argv[])
{

 try
 {
	//int myid(0), nprocs(1);  MPI_Init(&argc, &argv);  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	cout<<"\n Network Model Code version 2 alpha, built: "<<__DATE__ << endl;

	string arg1;
	if (argc > 1)                                 arg1 = argv[1];
	else  {  cout<< "Please input data file: ";   cin >> arg1;  }

	InputFile input(arg1,true);
	input.renameKeys("ALTR_CONT_ANG","EQUIL_CON_ANG");

	cnflowQD(input, 0,false);

   //MPI_Finalize();
 }
 catch (std::exception &exc) {  std::cerr << "\n\n Exception on processing: \n" << exc.what() << "Aborting! \n"  << endl;	return 1; }
 catch (...)                 {  std::cerr << "\n\n Unknown exception! \n Aborting! \n" << endl;	return 1; }
 return 0;
}

#endif

