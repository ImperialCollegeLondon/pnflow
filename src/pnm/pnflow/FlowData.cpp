#include "FlowData.h"
#if defined SVG
#include "PcKrRI_plotter.h"
#endif


#include "FlowDomain.h"


using namespace std;
void getWetIndices(double& amottWat, double& amottOil, double& USBMIndx, const stvec<stvec<stvec<double> > >& KcSwKrsQsss)
{
	int idrain=3,   iimbib=3-1;

	double amottImbib[] = {(*KcSwKrsQsss[iimbib].begin())[1],1.0, (*KcSwKrsQsss[iimbib].rbegin())[1]};
	double amottDrain[] = {(*KcSwKrsQsss[idrain].begin())[1],0.0, (*KcSwKrsQsss[idrain].rbegin())[1]};

	d_assert((*KcSwKrsQsss[idrain].begin())[0]<0.0);
	d_assert((*KcSwKrsQsss[idrain].rbegin())[0]>0.0);

	d_assert((*KcSwKrsQsss[iimbib].begin())[0]>0.0);
	d_assert((*KcSwKrsQsss[iimbib].rbegin())[0]<0.0);


	double areaDrain=0.0, areaImbib=0.0;

	for(size_t i = 1; i < KcSwKrsQsss[idrain].size(); ++i)
	{
		const stvec<double>& pt1(KcSwKrsQsss[idrain][i]), pt2(KcSwKrsQsss[idrain][i-1]); ///. Sw - Pc
		areaDrain += (pt2[1]-pt1[1]) * (max(pt2[0],0.0)+max(pt1[0],0.0))/2.0;
		if(pt1[0]==0.0 || pt1[0]*pt2[0] < 0.0)
			amottDrain[1] = pt1[1] + (pt2[1]-pt1[1])/(pt2[0]-pt1[0]) * (0.0-pt1[0]);
	}
	for(size_t i = 1; i < KcSwKrsQsss[iimbib].size(); ++i)
	{
		const vector<double>& pt1(KcSwKrsQsss[iimbib][i-1]), pt2(KcSwKrsQsss[iimbib][i]);  ///. Sw - Pc
		areaImbib -= (pt2[1]-pt1[1]) * (min(pt2[0],0.0)+min(pt1[0],0.0))/2.0;
		if(pt1[0]==0.0 || pt1[0]*pt2[0] < 0.0)
			amottImbib[1] = pt1[1] + (pt2[1]-pt1[1])/(pt2[0]-pt1[0]) * (0.0-pt1[0]);
	}

	amottWat = (amottImbib[1]-amottImbib[0])/(amottImbib[2]-amottImbib[0]);
	amottOil = (amottDrain[1]-amottDrain[0])/(amottDrain[2]-amottDrain[0]);
	USBMIndx = (areaDrain==0.0  ?  -std::numeric_limits<double>::infinity()  :  areaImbib==0.0  ?  std::numeric_limits<double>::infinity()  :  (log10(areaDrain/areaImbib)));
}




void FlowDomain::writeResultData(bool wantRelPerm, bool wantResIdx)
{

	ostringstream resFileName;
	resFileName << title_;

	bool matlabFormat(false), excelFormat(false), MCPMode(true);
	input_.resFormat(matlabFormat, excelFormat, MCPMode);


	//formatResults(matlabFormat, excelFormat, MCPMode);
	stvec< ststr >                          results;
	{
		bool wantRelPerm(resultWaterFlowRate_.size() > 0);
		bool wantResIdx(resultResistivityIdx_.size() > 0);

		string whiteSpace;
		if(matlabFormat) whiteSpace = ", ";
		else if(MCPMode) whiteSpace = "\t";
		else if(excelFormat) whiteSpace = ",";
		else whiteSpace = "	";


		for (size_t ii = 0; ii < resultWaterSat_.size(); ++ii)
		{
			ostringstream out;
			out.flags(ios::showpoint);
			out.flags(ios::fixed);

			out << resultWaterSat_[ii] << whiteSpace;
			out.flags(ios::scientific);
			out << setw(15) << resultCappPress_[ii] << whiteSpace;


					if(wantRelPerm)
					{
						double krw(0.0), kro(0.0);
						if(resultWaterFlowRate_[ii].first > 0.0)
						{
							if(useAvrPrsAsBdr_)
								krw = resultWaterFlowRate_[ii].first*singlePhaseDprs_ / (singlePhaseWaterQ_*resultWaterFlowRate_[ii].second);
							else
								krw = resultWaterFlowRate_[ii].first / singlePhaseWaterQ_;
						}
						if(resultOilFlowRate_[ii].first > 0.0)
						{
							if(useAvrPrsAsBdr_)
								kro = resultOilFlowRate_[ii].first*singlePhaseDprs_ / (singlePhaseOilQ_*resultOilFlowRate_[ii].second);
							else
								kro = resultOilFlowRate_[ii].first / singlePhaseOilQ_;
						}
						out << setw(15) << krw << whiteSpace << setw(15) << kro << whiteSpace;
					}
					else if(MCPMode)
						out_ << "Error microporosity format requires calc_perm to be true (See SAT_CONTROL keyword)";

					if(wantResIdx)
					{
						out << setw(15) << resultResistivityIdx_[ii] << whiteSpace;
					}
					else if(MCPMode)
						out << 1 << whiteSpace;

					if(reportMaterialBal_ && !MCPMode)
						out << setw(15) << resultWaterMass_[ii] << whiteSpace
							<< setw(15) << resultOilMass_[ii] << whiteSpace;


			if(matlabFormat) out << "; ...";

			results.push_back(out.str());
		}
	}


	ostringstream out;
	string legend;
	int icy = comn_.dispCycle();

	if(icy==0 && MCPMode)
	{
		string mcpFilename=title_+"_upscaled.tsv";
		ofstream of(mcpFilename);

		of<<endl<<"\nhomogeneous: \t"<<title_<<";"<<endl<<endl;

		double absPermeability = (singlePhaseWaterQ_ * water_.viscosity() * box_.x * (solverBoxEnd_-solverBoxStart_))
			/ (box_.y * box_.z * (inletSolverPrs_ - outletSolverPrs_));
		double formationFactor = ((box_.y * box_.z) * (inletSolverPrs_ - outletSolverPrs_))
			/ (water_.resistivity() * (singlePhaseCurrent_+1.0e-200) * box_.x * (solverBoxEnd_-solverBoxStart_));


		of<<endl<<title_<<"_permeability: \t" <<absPermeability<<";"<<endl;

		of<<endl<<title_<<"_porosity: \t"   <<flowVolume_ / satBoxVolume_<<";"<<endl;

		of<<endl<<title_<<"_formationfactor: \t" <<formationFactor<<";"<<endl<<endl;
	}

	if(comn_.injectant() == &oil_)
	{
		legend = matlabFormat ? "draincycle_" : " drainage, cycle  ";
		resFileName << "_cycle"<<dispCycle()<<"_drain";
	}
	else
	{
		legend = matlabFormat ? "imb_": " imbcycle ";
		resFileName << "_cycle"<<dispCycle()<<"_imb";
	}

	if(!MCPMode)
	{
		char seperator = excelFormat ? ',': ' ';
		if(!matlabFormat)
		{
			out << title_<<"  "<<legend<<"\n";
			out << (results.size()) << seperator << "data points" << endl;
			out << "Sw,			Pc(Pa)			";

			if(wantRelPerm)		out << "   ,Krw,		  	   ,Kro			   ";
			if(wantResIdx)		out << "   ,RI		   	  CC ";

			if(reportMaterialBal_)
				out << seperator << "Mass_w (kg) " << legend << icy << ")"
					<< seperator << "Mass_o (kg) " << legend << icy << ")";
			out << endl;
		}
		else
		{
			out << "DataLegend = {'Sw " << legend << icy << "' 'Pc (Pa) " << legend << icy << "'";
			if(wantRelPerm) out << "'Krw " << legend << icy << "' 'Kro " << legend << icy << "'";
			if(wantResIdx)  out << "'I " << legend << icy << "'";
			if(reportMaterialBal_) out << "'Mass_w (kg) " << legend << icy << "' 'Mass_o (kg) " << legend << icy << "'";
			out << "};" << endl << "Res_" << legend << icy << " = [";
		}
	}

	if(matlabFormat || excelFormat || !MCPMode)
	{
		string extension(".out");
		if(matlabFormat) extension = ".m";
		else if(excelFormat) extension = ".csv";

		ofstream of(resFileName.str()+extension);
		if(!of)	{ cerr << "\n\n *** Error: Could not open " << (resFileName.str()+extension) << " for writing ***\n\n"  << endl;		exit(-1);	}

		of << out.str();
		for(size_t j = 0; j < results.size(); ++j)		of << results[j] << endl;
		if(matlabFormat) of << "];" << endl;
	}




	double amottWat(0), amottOil(0), USBMIndx(0);
	if(icy == 3)
	{ getWetIndices(amottWat,amottOil,USBMIndx,comn_.KcSwKrsQsss_);

		out_<< "================== Wettability State ================="
			<< "\n Amott water index, Iw:            " << amottWat
			<< "\n Amott oil index, Io:              " << amottOil
			<< "\n Amott wettability index, I=Iw-Io: " << amottWat-amottOil
			<< "\n USBM  wettability index:           " << USBMIndx << endl;
	}


	if(MCPMode)
	{
		const string mcpFilename=title_+"_upscaled.tsv";
		ofstream of(mcpFilename.c_str(), std::ios_base::app);
		ensure(of,"can not open "+mcpFilename,2);

		of<<"\n\n\n"<<title_<<"_SwPcKrwKroRI_cycle"<<icy << ":\t  //"<< legend  << endl;
		of<< "//Sw     \t   Pc(Pa)   \t   Krw      \t   Kro      \t   RI"<<endl;

		for(size_t j = 0; j < results.size(); ++j)	of << results[j] << endl;
		of << "\n\n" << endl;

		if(icy == 3)
		{
			of  << "wettability indices: [" <<endl;
			of  << " AmottI= "<< setw(12) <<amottWat-amottOil<<"\t;\t AmottIw= "<< setw(12) <<amottWat<<"\t;\t AmottIo= "<< setw(12) <<amottOil<<"\t;\t USBM="<< setw(12) << USBMIndx << "\t;];" <<endl;
		}

		of.close();
	}




	#ifdef SVG
	if(MCPMode)
	{
		out_<<"\n\n Writting SVG ";
		#ifdef BOOST_SVPLOT_HPP
		out_<<input_.name()<<" ... >> "<<title_+"_upscaled.svg";
		plotSwPcKrRI(comn_.KcSwKrsQsss_,  title_+"_upscaled.svg", title_, amottWat, amottOil, USBMIndx);
		#else
		out_<<" : no ";
		#endif
	}
	#endif
	out_<<"\n\n";


	resultWaterFlowRate_.clear();
	resultOilFlowRate_.clear();
	resultWaterSat_.clear();
	resultCappPress_.clear();
	//resultBoundingPc_.clear();
	resultResistivityIdx_.clear();
	resultWaterMass_.clear();
	resultOilMass_.clear();
}




