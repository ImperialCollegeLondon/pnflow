#include "FlowData.h"
#ifdef SVG
#include "PcKrRI_plotter.h"
#endif
#include <iomanip>




using namespace std;
void getWetIndices(double& amottWat, double& amottOil, double& USBMIndx, const vector<vector<array<double,8> > >& SwKcKrsQsss)  {
	int idrain=3,   iimbib=3-1;

	double amottImbib[] = {(*SwKcKrsQsss[iimbib].begin())[ISw],1., (*SwKcKrsQsss[iimbib].rbegin())[ISw]};
	double amottDrain[] = {(*SwKcKrsQsss[idrain].begin())[ISw],0., (*SwKcKrsQsss[idrain].rbegin())[ISw]};

	dbgAsrt((*SwKcKrsQsss[idrain].begin())[IPc]<0.);
	dbgAsrt((*SwKcKrsQsss[idrain].rbegin())[IPc]>0.);

	dbgAsrt((*SwKcKrsQsss[iimbib].begin())[IPc]>0.);
	dbgAsrt((*SwKcKrsQsss[iimbib].rbegin())[IPc]<0.);


	double areaDrain=0., areaImbib=0.;

	for(size_t i = 1; i<SwKcKrsQsss[idrain].size(); ++i)  {
		const auto& pt1(SwKcKrsQsss[idrain][i]), pt2(SwKcKrsQsss[idrain][i-1]); ///. Sw - Pc
		areaDrain += (pt2[ISw]-pt1[ISw]) * (max(pt2[IPc],0.)+max(pt1[IPc],0.))/2.;
		if(pt1[IPc]==0. || pt1[IPc]*pt2[IPc] < 0.)
			amottDrain[ISw] = pt1[ISw] + (pt2[ISw]-pt1[ISw])/(pt2[IPc]-pt1[IPc]) * (0.-pt1[IPc]);
	}
	for(size_t i = 1; i<SwKcKrsQsss[iimbib].size(); ++i)  {
		const auto& pt1(SwKcKrsQsss[iimbib][i-1]), pt2(SwKcKrsQsss[iimbib][i]);  ///. Sw - Pc
		areaImbib -= (pt2[ISw]-pt1[ISw]) * (min(pt2[IPc],0.)+min(pt1[IPc],0.))/2.;
		if(pt1[IPc]==0. || pt1[IPc]*pt2[IPc] < 0.)
			amottImbib[ISw] = pt1[ISw] + (pt2[ISw]-pt1[ISw])/(pt2[IPc]-pt1[IPc]) * (0.-pt1[IPc]);
	}

	amottWat = (amottImbib[1]-amottImbib[0])/(amottImbib[2]-amottImbib[0]);
	amottOil = (amottDrain[1]-amottDrain[0])/(amottDrain[2]-amottDrain[0]);
	USBMIndx = (areaDrain==0.  ?  -std::numeric_limits<double>::infinity()  :  areaImbib==0.  ?  std::numeric_limits<double>::infinity()  :  (log10(areaDrain/areaImbib)));
}


void GNMData::writeResultData(bool wantRelPerm, bool wantResIdx, const std::string& nam) const
{



	unsigned int nMinNKrLines(46);
	string resForm("svg");
	std::istringstream resformat;
	if (input_.giv("RES_FORMAT", resformat))  resformat>>resForm>>nMinNKrLines;

	char format = tolower(resForm.front());

	int icy = (*this).dispCycle();
	double amottWat(0), amottOil(0), USBMIndx(0);
	if(icy == 3)  getWetIndices(amottWat,amottOil,USBMIndx,(*this).KrQsss_);



	if(format=='s')  {
		#ifdef SVG
		out_<<"\n\nWriting "<<nam+"_upscaled.svg";
		plotSwPcKrRI(KrQsss_,  nam+"_upscaled.svg", nam, amottWat, amottOil, USBMIndx);
		out_<<"\n\n";
		return;
		#endif
		out_<<" : svplot not supported, switching to tsv"; format='t';
	}



	if (icy==0)  {  if(format=='t')  {
			ofstream of(nam+"_upscaled.tsv");   ensure(of,"can not open "+nam+"_upscaled.tsv",2);

			of<<"//-*- C++ -*-\t upscaled results\n\n";
			of<<"\n\nporousRock:  \t"<<nam<<"\n"<<endl;

			ensure((*this).KrQsss_.size(), "SwKcKrsQsss should have been initialized before writeResultData");
			of<<nam<<"_porosity:        \t"<<(*this).KrQsss_[0][0][ISw]<<"\t;"<<endl;
			of<<nam<<"_permeability:    \t"<<(*this).KrQsss_[0][0][2]<<"\t;"<<endl;
			of<<nam<<"_formationfactor: \t"<<(*this).KrQsss_[0][0][4]<<"\t;"<<endl;
	}	}
	else {

		vector<string> results;
		string wspace;
		if(format=='m') wspace = ", ";       // Matlab
		else if(format=='e') wspace = ",";   // Excel -> csv
		else { wspace = "\t "; format='t'; }    // tsv
		const vector<array<double,8> >& KcSwKrsQss = (*this).KrQsss_[(*this).dispCycle()];
		for(size_t istp=0; istp < KcSwKrsQss.size(); ++istp)	{ //SYNC39756397846
			ostringstream out;
			out.flags(ios::showpoint);
			out.flags(ios::fixed);

			out << min(max(-1., KcSwKrsQss[istp][ISw]),2.)<< wspace;
			out.flags(ios::scientific);
			out << setw(14) << KcSwKrsQss[istp][IPc]*(*this).watoil().sigma() << wspace;


			if(wantRelPerm)
				 out << setw(14) << min(max(-1., KcSwKrsQss[istp][2]),2.) << wspace << setw(14) << min(max(-1., KcSwKrsQss[istp][3]),2.) << wspace;

			if(wantResIdx)       out << setw(14) << min(max(0., KcSwKrsQss[istp][4]),1e32) << wspace;
			else if (format=='t')  out  << 1 << wspace;


			if(format=='m') out << "; ...";

			results.push_back(out.str());
		}



		if (format=='m' || format=='e')  {
			ostringstream out;
			string legend;
			std::string resFileName = nam;

			if((*this).isDrainage())  {
			  legend = format=='m' ? "draincycle_" : " drainage  ";
			  resFileName += "_cycle"+_s(icy)+"_drain";
			}
			else
			{
			  legend = format=='m' ? "imb_": " imbibition ";
			  resFileName += "_cycle"+_s(icy)+"_imb";
			}
			resFileName+=(format=='m')? ".m" : ".csv";

			if(format=='e')  {
				out << input_.getOr("networkFile", nam) <<"  "<<legend<<"\n";
				out << (results.size()) << ", data points" << endl;
				out << "Sw,            Pc(Pa)            ";
				if(wantRelPerm)	out << ",Krw,               ,Kro               ";
				if(wantResIdx)				out << ",RI           ";
				out << endl;
			}
			else if(format=='m')  {
				out << "DataLegend = {'Sw " << legend << icy << "' 'Pc (Pa) " << legend << icy << "'";
				if(wantRelPerm) out << "'Krw " << legend << icy << "' 'Kro " << legend << icy << "'";
				if(wantResIdx)  out << "'I " << legend << icy << "'";
				out << "};" << endl << "Res_" << legend << icy << " = [";
			}

			ofstream of(resFileName);
			ensure(of,  " Could not open " + resFileName + " for writing ",2);

			of << out.str();
			for(size_t j=0; j < results.size(); ++j)	of << results[j] << endl;
			if(format=='m') of << "];" << endl;
		}




		if(icy == 3 && informative())
			out_<< "================== Wettability State ================="
				<< "\n Amott water index, Iw:            " << amottWat
				<< "\n Amott oil index, Io:              " << amottOil
				<< "\n Amott wettability index, I=Iw-Io: " << amottWat-amottOil
				<< "\n USBM  wettability index:           " << USBMIndx << endl;


		if(format=='t')  {
			ofstream of(nam+"_upscaled.tsv", std::ios_base::app);   ensure(of,"can not open "+nam+"_upscaled.tsv",2);

			of<<"\n"<<nam<<"_SwPcKrwKroRI_cycle"<<icy << ": \t// "<< ((*this).isDrainage() ? "drainage" : "imbibition")<<" \t cycle_"<<icy<<endl;
			of << "//Sw          \t Pc(Pa)        \t Krw           \t Kro            \t  RI"<<endl;

			size_t krLine=0;
			for(; krLine < results.size(); ++krLine) 	of << results[krLine] << "\n";
			of << "\n" << endl; ++krLine;
			for(; krLine < nMinNKrLines; ++krLine)   	of << "\n";

			if(icy == 3)  {
				of  << "wettability indices: [" <<endl;
				of  << " AmottI= "<< setw(12) <<amottWat-amottOil<<"\t;\t AmottIw= "<< setw(12) <<amottWat<<"\t;\t AmottIo= "<< setw(12) <<amottOil<<"\t;\t USBM="<< setw(12) << USBMIndx << "\t;];" <<endl;
			}
		}
	}


	out_<<"\n\n";

}


