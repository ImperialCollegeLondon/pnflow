
/*---------------------------------------------------------------------------*\
Developed by (2015): Ali Q Raeini  email: a.q.raeini@gmail.com
\*---------------------------------------------------------------------------*/


#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <valarray>


#include "Element.h"

#include "elem_Model.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "../compareFuncs.h"






template<class Elemnt,  class Elemp>
double accumulate(std::function<double(double, double)> operatorFunc, const std::vector< Elemnt const * >& elems,  double (Elemp::*radFunc)() const,  double result=0.)  {	for(const Elemnt* v : elems)	result = operatorFunc(result,(v->*radFunc)());	return result; }

template<class IterT,  class Elemp>
double accumulate(std::function<double(double, double)> operatorFunc, IterT iter, const IterT last,  double (Elemp::*radFunc)() const,  double result=0.)  {	for(;iter<last;++iter)	result = operatorFunc(result,((*iter)->*radFunc)());	return result; }

template<typename T>
std::valarray<std::valarray<T> > transpose(const std::valarray<std::valarray<T> >& vecvec)  {
	if(!vecvec.size()) return std::valarray<std::valarray<T> >();

	std::valarray<std::valarray<T> > trans(std::valarray<T>(0.,vecvec.size()), vecvec[0].size());
	for (size_t i=0; i<vecvec[0].size();++i)  {
		for (size_t j=0; j<vecvec.size();++j) trans[i][j] = vecvec[j][i] ;
	}
	return trans;
}
template<typename T>
std::ostream & operator << (std::ostream & outstream, const std::valarray<T>& vec)  {
	if(vec.size() && vec.size()<10)  for (auto v : vec) outstream << v << ' ';
	else                             for (auto v : vec) outstream << v << '\n';
	return outstream;
}

#include "netStatistics.h"

using namespace std;
//const static double PI = 3.14159265358979;


void printRadiusStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {
		cout<<"\n\n//           \t arithmetic \t volume-weighted\n";

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadPore = elms[elms.size()/2]->model()->RRR();
		double RadAvgPore=0.; ///VolWeighted pore shape factor
		double VWRadAvgPore=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgPore+= elms[i]->flowVolume()*elms[i]->model()->RRR();
			RadAvgPore+= elms[i]->model()->RRR();
		}
		RadAvgPore=RadAvgPore/(nBpPors_-nBSs_);
		VWRadAvgPore=VWRadAvgPore/sumFlowVolume;
		//cout<<"nPores_:\t"<<nPores_<<endl;
		cout<<"RadAvgPore:  \t "<<RadAvgPore<<" \t "<<VWRadAvgPore<<endl;
	}


	{
		double nTrots=elemans.size()-nBpPors_;
		vector< Elem const * > elms(elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadThroat = elms[elms.size()/2]->model()->RRR();
		double RadAvgThroat=0.; ///VolWeighted pore shape factor
		double VWRadAvgThroat=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgThroat+= elms[i]->flowVolume()*elms[i]->model()->RRR();
			RadAvgThroat+= elms[i]->model()->RRR();
		}
		RadAvgThroat=RadAvgThroat/nTrots;
		VWRadAvgThroat=VWRadAvgThroat/sumFlowVolume;
		//cout<<"numThroats:\t"<<nTrots<<endl;
		//cout<<"RadAvgThroat:\t"<<RadAvgThroat<<endl;
		//cout<<"VWRadAvgThroat:\t"<<VWRadAvgThroat<<endl;
		cout<<"RadAvgThroat:\t "<<RadAvgThroat<<"\t "<<VWRadAvgThroat<<endl;
	}



	{
		double numElems_=elemans.size()-nBSs_;
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		elms.insert(elms.end(), elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadElem = elms[elms.size()/2]->model()->RRR();
		double RadAvgElem=0.; ///VolWeighted pore shape factor
		double VWRadAvgElem=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgElem+= elms[i]->flowVolume()*elms[i]->model()->RRR();
			RadAvgElem+= elms[i]->model()->RRR();
		}
		RadAvgElem=RadAvgElem/numElems_;
		VWRadAvgElem=VWRadAvgElem/sumFlowVolume;
		//cout<<"numElems_:\t"<<numElems_<<endl;
		//cout<<"RadAvgElem:\t"<<RadAvgElem<<endl;
		//cout<<"VWRadAvgElem:\t"<<VWRadAvgElem<<endl;
		cout<<"RadAvgElem:  \t "<<RadAvgElem<<"\t "<<VWRadAvgElem<<endl;
	}
}



void printDistanceMapStatistics( const std::vector<Elem const *>&  elemans, int nBSs_, int nBpPors_)  {

		int nsteps = 32;

//int ErOrr; //trote->length() is assumed to be total length CHECK needed
		///. max radius:
		double maxRad = accumulate((double const & (*) (double const &, double const &))(std::max<double>), elemans, &Elem::RRR)+1e-15;
		double dr = maxRad/nsteps*(1.+1e-14);
		cout<<"dr: "<<dr<< "   maxRad: "<<maxRad<<endl;

		std::valarray<std::valarray<double> > distribG(std::valarray<double>(0., nsteps),2);
		for (int i=0; i<nsteps; ++i)	distribG[0][i] = 0.+dr/2+i*dr;
		for(size_t i = nBpPors_; i < elemans.size(); ++i) if(elemans[i]->nCncts()==2)  {
			const Throat* trote = dynamic_cast<const Throat*>(elemans[i]);

			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elemans[i]->model());
			double rt=trote->RRR();
			for (double rr=0.5*dr; rr<rt; rr+=dr)  {
				if(shyp0)
				 for(int ic=0;ic<shyp0->numCorners();++ic)
					  distribG[1][rr/dr] += 2.*((rt-rr)/tan(shyp0->cornerHalfAngles(ic)))*max(0.,(trote->throatLength())); //int Warn = -trote->poreLength(0)-trote->poreLength(1)  /// dA
				else
					 distribG[1][rr/dr] += 2.*PI*((rt-rr))*(trote->throatLength()); //int Warn = -trote->poreLength(0)-trote->poreLength(1)  /// dA
			}

			rt=elemans[i]->neib(0)->RRR();
			const Polygon*  shyp1 = dynamic_cast<const Polygon*>(elemans[i]->neib(0)->model());
			for (double rr=0.5*dr; rr<rt; rr+=dr)  {
				if(shyp1)
				for(int ic=0;ic<shyp1->numCorners();++ic)  {
					//double dA =  dA;
					distribG[1][rr/dr] += 2.*((rt-rr)/tan(shyp1->cornerHalfAngles(ic)))*(trote->poreLength(0));
				}
				else
					distribG[1][rr/dr] += 2.*PI*((rt-rr))*(trote->poreLength(0));  /// dA
			}
			rt=elemans[i]->neib(1)->RRR();
			const Polygon*  shyp2 = dynamic_cast<const Polygon*>(elemans[i]->neib(1)->model());
			for (double rr=0.5*dr; rr<rt; rr+=dr)  {
				if(shyp2)
				for(int ic=0;ic<shyp2->numCorners();++ic)  {
					//double dA = dA;
					distribG[1][rr/dr] += 2.*((rt-rr)/tan(shyp2->cornerHalfAngles(ic)))*(trote->poreLength(1));
				}
				else
					distribG[1][rr/dr] += 2.*PI*((rt-rr))*(trote->poreLength(1));  /// dA
			}

		}

		cout<<"\n\n"<< "// cnm distance map distribution from G, "<<"maxRad "<<maxRad<<"dr "<<dr<<endl;
		cout<<"distanceMapDistribution: //distance\tfrequency"<<endl;
		cout<<transpose(distribG);

}




void printCornerAngStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.; ///VolWeighted pore shape factor
		double VWhAng=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=1e-32;
		double sumhAng=1e-32;
		for(size_t i = 0; i < elms.size(); ++i)  {
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)  {
				double rt=shyp0->RRR();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngPore:    \t "<<hAng<<"\t "<<VWhAng<<endl;
	}


	{
		vector< Elem const * > elms(elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.; ///VolWeighted pore shape factor
		double VWhAng=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=1e-32;
		double sumhAng=1e-32;
		for(size_t i = 0; i < elms.size(); ++i)  {
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)  {
				double rt=shyp0->RRR();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngThroat:   \t "<<hAng<<"\t "<<VWhAng<<endl;
	}



	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		elms.insert(elms.end(), elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.; ///VolWeighted pore shape factor
		double VWhAng=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=1e-32;
		double sumhAng=1e-32;
		for(size_t i = 0; i < elms.size(); ++i)  {
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)  {
				double rt=shyp0->RRR();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngElem:  \t "<<hAng<<"\t "<<VWhAng<<endl;
	}
}




void printCornerNumStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double Ncor=0.; ///VolWeighted pore shape factor
		double VWNcor=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcor+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			Ncor+= elms[i]->model()->numCorners();
		}
		Ncor=Ncor/(nBpPors_-nBSs_);
		VWNcor=VWNcor/sumFlowVolume;
		cout<<"NcorPore:    \t "<<Ncor<<"\t "<<VWNcor<<endl;
	}


	{
		double nTrots=elemans.size()-nBpPors_;
		vector< Elem const * > elms(elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double NcorThroat=0.; ///VolWeighted pore shape factor
		double VWNcorThroat=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcorThroat+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			NcorThroat+= elms[i]->model()->numCorners();
		}
		NcorThroat=NcorThroat/nTrots;
		VWNcorThroat=VWNcorThroat/sumFlowVolume;
		cout<<"NcorThroat:  \t "<<NcorThroat<<"\t "<<VWNcorThroat<<endl;
	}



	{
		double numElems_=elemans.size()-nBSs_;
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		elms.insert(elms.end(), elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double NcorElem=0.; ///VolWeighted pore shape factor
		double VWNcorElem=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcorElem+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			NcorElem+= elms[i]->model()->numCorners();
		}
		NcorElem=NcorElem/numElems_;
		VWNcorElem=VWNcorElem/sumFlowVolume;
		cout<<"NcorElem:     \t "<<NcorElem<<"\t "<<VWNcorElem<<endl;
	}
}


void printShapeFactorStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGPore = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgPore=0.; ///VolWeighted pore shape factor
		double VWGAvgPore=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgPore+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgPore+= elms[i]->model()->shapeFactor();
		}
		GAvgPore=GAvgPore/(nBpPors_-nBSs_);
		VWGAvgPore=VWGAvgPore/sumFlowVolume;
		//cout<<"nPores_:\t"<<nPores_<<endl;
		cout<<"GAvgPore:    \t "<<GAvgPore<<"\t "<<VWGAvgPore<<endl;
	}


	{
		double nTrots=elemans.size()-nBpPors_;
		vector< Elem const * > elms(elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGThroat = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgThroat=0.; ///VolWeighted pore shape factor
		double VWGAvgThroat=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgThroat+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgThroat+= elms[i]->model()->shapeFactor();
		}
		GAvgThroat=GAvgThroat/nTrots;
		VWGAvgThroat=VWGAvgThroat/sumFlowVolume;
		//cout<<"nTrots:\t"<<nTrots<<endl;
		//cout<<"GAvgThroat:\t"<<GAvgThroat<<endl;
		//cout<<"VWGAvgThroat:\t"<<VWGAvgThroat<<endl;
		cout<<"GAvgThroat:  \t "<<GAvgThroat<<"\t "<<VWGAvgThroat<<endl;
	}



	{
		double numElems_=elemans.size()-nBSs_;
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		elms.insert(elms.end(), elemans.begin()+nBpPors_, elemans.end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGElem = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgElem=0.; ///VolWeighted pore shape factor
		double VWGAvgElem=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgElem+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgElem+= elms[i]->model()->shapeFactor();
		}
		GAvgElem=GAvgElem/numElems_;
		VWGAvgElem=VWGAvgElem/sumFlowVolume;
		//cout<<"numElems_:\t"<<numElems_<<endl;
		//cout<<"GAvgElem:\t"<<GAvgElem<<endl;
		//cout<<"VWGAvgElem:\t"<<VWGAvgElem<<endl;
		cout<<"GAvgElem:     \t "<<GAvgElem<<"\t "<<VWGAvgElem<<endl;
	}
}


void printAspectRatioStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		//sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianAspRatioPore = elms[elms.size()/2]->model()->shapeFactor();
		double AspRatioAvgPore=0.; ///VolWeighted pore shape factor
		double VWAspRatioAvgPore=0.; ///VolWeighted pore shape factor
		double VWAWAspRatioAvgPore=0.; ///VolWeighted pore shape factor
		double VWQWAspRatioAvgPore=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			const Pore* pore = dynamic_cast< Pore const * >(elms[i]);
			//if (pore->nCncts()>0)
			double AvgAspectRatio=0.;
			double AWAvgAspectRatio=0.,sumAW=1e-31;
			double QWAvgAspectRatio=0.,sumQW=1e-31;
			//double maxAspectRatio=0.;
			double poreR=pore->model()->RRR();
			for(int conn = 0; conn < pore->nCncts(); ++conn)  {
				const ElemModel* shape = pore->neib(conn)->model();
				double AspectRatio= shape->RRR()/poreR;
				AvgAspectRatio+= AspectRatio;
				sumAW+=shape->area();
				 AWAvgAspectRatio+= shape->area()*AspectRatio;
				sumQW+=shape->SPConductance(shape->area(),1.);
				 QWAvgAspectRatio+=shape->SPConductance(shape->area(),1.)*AspectRatio;
				//maxAspectRatio=max(maxAspectRatio,AspectRatio);
			}
				AvgAspectRatio/=pore->nCncts()+1e-16;
				AWAvgAspectRatio/=sumAW;
				QWAvgAspectRatio/=sumQW;
			sumFlowVolume+=elms[i]->flowVolume();
			VWAWAspRatioAvgPore+= elms[i]->flowVolume()*AWAvgAspectRatio;
			VWQWAspRatioAvgPore+= elms[i]->flowVolume()*QWAvgAspectRatio;
			VWAspRatioAvgPore+= elms[i]->flowVolume()*AvgAspectRatio;
			AspRatioAvgPore+= AvgAspectRatio;
		}
		AspRatioAvgPore=AspRatioAvgPore/(nBpPors_-nBSs_);
		VWAspRatioAvgPore=VWAspRatioAvgPore/sumFlowVolume;
		VWAWAspRatioAvgPore=VWAWAspRatioAvgPore/sumFlowVolume;
		VWQWAspRatioAvgPore=VWQWAspRatioAvgPore/sumFlowVolume;
		//cout<<"nPores_:\t"<<nPores_<<endl;
		cout<<"//               \tArithmatic\t  VW    \t  VWAW   \t  VWQW   \n";
		cout<<"AspRatioAvg: \t "
			<<AspRatioAvgPore<<"\t "
			<<VWAspRatioAvgPore<<"\t "
			<<VWAWAspRatioAvgPore<<"\t "
			<<VWQWAspRatioAvgPore<<endl;
	}


}




void printCoordinaNumStatistics( const std::vector<Elem const*>&  elemans, int nBSs_, int nBpPors_)  {

	{
		vector< Elem const * > elms(elemans.begin()+nBSs_, elemans.begin()+nBpPors_);
		//sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianCoordNPore = elms[elms.size()/2]->model()->shapeFactor();
		double CoordNAvgPore=0.; ///VolWeighted pore shape factor
		double VWCoordNAvgPore=0.; ///VolWeighted pore shape factor
		double VWAWCoordNAvgPore=0.; ///VolWeighted pore shape factor
		double VWQWCoordNAvgPore=0.; ///VolWeighted pore shape factor
		double sumFlowVolume=0.;
		for(size_t i = 0; i < elms.size(); ++i)  {
			const Pore* pore = dynamic_cast< Pore const * >(elms[i]);
			//if (pore->nCncts()>0)
			double AvgCoordN=0.;
			double AWAvgCoordN=0.,sumAW=0.;
			double QWAvgCoordN=0.,sumQW=0.;
			//double maxCoordN=0.;
			//double poreR=pore->model()->RRR();
			for(int conn = 0; conn < pore->nCncts(); ++conn)  {
				const ElemModel* shape = pore->neib(conn)->model();
				AvgCoordN+= 1;
				sumAW+=shape->area();
				 AWAvgCoordN+= shape->area()*shape->area();
				sumQW+=shape->SPConductance(shape->area(),1.);
				 QWAvgCoordN+=shape->SPConductance(shape->area(),1.)*shape->SPConductance(shape->area(),1.);
				//maxCoordN=max(maxCoordN,CoordN);
			}
				AvgCoordN/=1;
				AWAvgCoordN=sumAW*sumAW/(AWAvgCoordN+1e-63);
				QWAvgCoordN=sumQW*sumQW/(QWAvgCoordN+1e-63);
			sumFlowVolume+=elms[i]->flowVolume();
			VWAWCoordNAvgPore+= elms[i]->flowVolume()*AWAvgCoordN;
			VWQWCoordNAvgPore+= elms[i]->flowVolume()*QWAvgCoordN;
			VWCoordNAvgPore+= elms[i]->flowVolume()*AvgCoordN;
			CoordNAvgPore+= AvgCoordN;
		}
		CoordNAvgPore=CoordNAvgPore/(nBpPors_-nBSs_);
		VWCoordNAvgPore=VWCoordNAvgPore/sumFlowVolume;
		VWAWCoordNAvgPore=VWAWCoordNAvgPore/sumFlowVolume;
		VWQWCoordNAvgPore=VWQWCoordNAvgPore/sumFlowVolume;
		//cout<<"nPores_:\t"<<nPores_<<endl;
		cout<<"CoordNAvgPore:  \t "
			<<CoordNAvgPore<<"\t "
			<<VWCoordNAvgPore<<"\t "
			<<VWAWCoordNAvgPore<<"\t "
			<<VWQWCoordNAvgPore<<endl;
	}


}
