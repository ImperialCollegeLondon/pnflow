
/*---------------------------------------------------------------------------*\
Developed by (2015): Ali Q Raeini  email: a.qaseminejad-raeini09@imperial.ac.uk
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
//#include "Element.h"
//#include "typses.h"
//#include "fluid.h"

#include "Element.h"

#include "elem_Model.h"
#include "polygon.h"
//#include "elem_hetroPorous.h"
//#include "elem_porous.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "../compareFuncs.h"

template<class Elemnt,  class Elemp>
double accumulate(std::function<double(double, double)> operatorFunc, const std::vector< Elemnt const * >& elems,  double (Elemp::*radFunc)() const,  double result=0.0)
{	for(const Elemnt* v : elems)	result = operatorFunc(result,(v->*radFunc)());	return result; }

template<class IterT,  class Elemp>
double accumulate(std::function<double(double, double)> operatorFunc, IterT iter, const IterT last,  double (Elemp::*radFunc)() const,  double result=0.0)
{	for(;iter<last;++iter)	result = operatorFunc(result,((*iter)->*radFunc)());	return result; }

template<typename T>
std::valarray<std::valarray<T> > transpose(const std::valarray<std::valarray<T> >& vecvec)
{
	if(!vecvec.size()) return std::valarray<std::valarray<T> >();
	
	std::valarray<std::valarray<T> > trans(std::valarray<T>(0.0,vecvec.size()), vecvec[0].size());
	for (size_t i=0; i<vecvec[0].size();++i)
	{
		for (size_t j=0; j<vecvec.size();++j) trans[i][j] = vecvec[j][i] ;
	}
	return trans;
}
template<typename T>
std::ostream & operator << (std::ostream & outstream, const std::valarray<T>& vec)
{
	if(vec.size() && vec.size()<10)  for (auto v : vec) outstream << v << ' ';
	else                             for (auto v : vec) outstream << v << '\n';
	return outstream;
}

#include "netStatistics.h"

using namespace std;
const static double PI = 3.14159265358979;

void printDistanceMapStatistics( const std::vector<Element const *>&  m_rockLattices, int m_numPores)
{

		int nsteps = 32;

		///. max radius:
		double maxRad = accumulate((double const & (*) (double const &, double const &))(std::max<double>), m_rockLattices, &Element::RRR)+1.0e-15;
		double dr = maxRad/nsteps*(1.0+1.0e-14);
		cout<<"dr "<<dr<<endl;
		cout<<"maxRad "<<maxRad<<endl;

		std::valarray<std::valarray<double> > distribG(std::valarray<double>(0.0, nsteps),2);
		for (int i=0; i<nsteps; ++i)	distribG[0][i] = 0.0+dr/2+i*dr;
		for(size_t i = m_numPores+2; i < m_rockLattices.size(); ++i)
		{
			const Throat* trote = dynamic_cast<const Throat*>(m_rockLattices[i]);
			
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(m_rockLattices[i]->model());
			double rt=trote->RRR();
			for (double rr=0.5*dr; rr<rt; rr+=dr)
			{	
				if(shyp0)
				 for(int ic=0;ic<shyp0->numCorners();++ic)
					distribG[1][rr/dr] += 2.0*((rt-rr)/tan(shyp0->cornerHalfAngles(ic)))*max(0.0,(trote->length()-trote->poreLength(trote->connection(0))-trote->poreLength(trote->connection(1))));  /// dA
				else 
					distribG[1][rr/dr] += 2.0*PI*((rt-rr))*(trote->length()-trote->poreLength(trote->connection(0))-trote->poreLength(trote->connection(1)));  /// dA
			}

			rt=m_rockLattices[i]->connection(0)->RRR();
			const Polygon*  shyp1 = dynamic_cast<const Polygon*>(m_rockLattices[i]->connection(0)->model());
			for (double rr=0.5*dr; rr<rt; rr+=dr)
			{	
				if(shyp1)
				for(int ic=0;ic<shyp1->numCorners();++ic)
				{
					//double dA =  dA;
					distribG[1][rr/dr] += 2.0*((rt-rr)/tan(shyp1->cornerHalfAngles(ic)))*(trote->poreLength(trote->connection(0)));
				}
				else 
					distribG[1][rr/dr] += 2.0*PI*((rt-rr))*(trote->poreLength(trote->connection(0)));  /// dA
			}
			rt=m_rockLattices[i]->connection(1)->RRR();
			const Polygon*  shyp2 = dynamic_cast<const Polygon*>(m_rockLattices[i]->connection(1)->model());
			for (double rr=0.5*dr; rr<rt; rr+=dr)
			{	
				if(shyp2)
				for(int ic=0;ic<shyp2->numCorners();++ic)
				{
					//double dA = dA;
					distribG[1][rr/dr] += 2.0*((rt-rr)/tan(shyp2->cornerHalfAngles(ic)))*(trote->poreLength(trote->connection(1)));
				}
				else 
					distribG[1][rr/dr] += 2.0*PI*((rt-rr))*(trote->poreLength(trote->connection(1)));  /// dA
			}

		}
		
		cout<<"\n\n"<< " distance map distribution from G"<<endl;
		cout<<"distance frequency "<<endl;
		cout<<transpose(distribG);

}




void printCornerAngStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.0; ///VolWeighted pore shape factor
		double VWhAng=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=1.0e-32;
		double sumhAng=1.0e-32;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)
			 {
				double rt=shyp0->radius();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.0;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngPore   VWhAngPore:  \t "<<hAng<<"\t "<<VWhAng<<endl;
	}


	{
		vector< Element const * > elms((*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.0; ///VolWeighted pore shape factor
		double VWhAng=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=1.0e-32;
		double sumhAng=1.0e-32;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)
			 {
				double rt=shyp0->radius();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.0;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngThroat VWhAngThroat:\t "<<hAng<<"\t "<<VWhAng<<endl;
	}



	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		elms.insert(elms.end(), (*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double hAng=0.0; ///VolWeighted pore shape factor
		double VWhAng=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=1.0e-32;
		double sumhAng=1.0e-32;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			const Polygon* shyp0 = dynamic_cast<const Polygon*>(elms[i]->model());
			if(shyp0)
			 for(int ic=0;ic<shyp0->numCorners();++ic)
			 {
				double rt=shyp0->radius();
				double vol = (rt*rt/tan(shyp0->cornerHalfAngles(ic)))/shyp0->area()*elms[i]->flowVolume();
				VWhAng += vol*shyp0->cornerHalfAngles(ic);  /// dA
				sumFlowVolume += vol;  /// dA
				hAng += shyp0->cornerHalfAngles(ic);  /// dA
				sumhAng += 1.0;  /// dA
			 }
		}
		hAng=hAng/sumhAng;
		VWhAng=VWhAng/sumFlowVolume;
		cout<<"hAngElem   VWhAngElem:  \t "<<hAng<<"\t "<<VWhAng<<endl;
	}
}




void printCornerNumStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double Ncor=0.0; ///VolWeighted pore shape factor
		double VWNcor=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcor+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			Ncor+= elms[i]->model()->numCorners();
		}
		Ncor=Ncor/m_numPores;
		VWNcor=VWNcor/sumFlowVolume;
		cout<<"NcorPore   VWNcorPore:  \t "<<Ncor<<"\t "<<VWNcor<<endl;
	}


	{
		double m_numThroats=(*rockLattices).size()-2-m_numPores;
		vector< Element const * > elms((*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double NcorThroat=0.0; ///VolWeighted pore shape factor
		double VWNcorThroat=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcorThroat+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			NcorThroat+= elms[i]->model()->numCorners();
		}
		NcorThroat=NcorThroat/m_numThroats;
		VWNcorThroat=VWNcorThroat/sumFlowVolume;
		cout<<"NcorThroat VWNcorThroat:\t "<<NcorThroat<<"\t "<<VWNcorThroat<<endl;
	}



	{
		double m_numElems=(*rockLattices).size()-2;
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		elms.insert(elms.end(), (*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		double NcorElem=0.0; ///VolWeighted pore shape factor
		double VWNcorElem=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWNcorElem+= elms[i]->flowVolume()*elms[i]->model()->numCorners();
			NcorElem+= elms[i]->model()->numCorners();
		}
		NcorElem=NcorElem/m_numElems;
		VWNcorElem=VWNcorElem/sumFlowVolume;
		cout<<"NcorElem   VWNcorElem:  \t "<<NcorElem<<"\t "<<VWNcorElem<<endl;
	}
}


void printShapeFactorStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGPore = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgPore=0.0; ///VolWeighted pore shape factor
		double VWGAvgPore=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgPore+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgPore+= elms[i]->model()->shapeFactor();
		}
		GAvgPore=GAvgPore/m_numPores;
		VWGAvgPore=VWGAvgPore/sumFlowVolume;
		//cout<<"m_numPores:\t"<<m_numPores<<endl;
		cout<<"GAvgPore   VWGAvgPore:  \t "<<GAvgPore<<"\t "<<VWGAvgPore<<endl;
	}


	{
		double m_numThroats=(*rockLattices).size()-2-m_numPores;
		vector< Element const * > elms((*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGThroat = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgThroat=0.0; ///VolWeighted pore shape factor
		double VWGAvgThroat=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgThroat+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgThroat+= elms[i]->model()->shapeFactor();
		}
		GAvgThroat=GAvgThroat/m_numThroats;
		VWGAvgThroat=VWGAvgThroat/sumFlowVolume;
		//cout<<"m_numThroats:\t"<<m_numThroats<<endl;
		//cout<<"GAvgThroat:\t"<<GAvgThroat<<endl;
		//cout<<"VWGAvgThroat:\t"<<VWGAvgThroat<<endl;
		cout<<"GAvgThroat VWGAvgThroat:\t "<<GAvgThroat<<"\t "<<VWGAvgThroat<<endl;
	}



	{
		double m_numElems=(*rockLattices).size()-2;
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		elms.insert(elms.end(), (*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianGElem = elms[elms.size()/2]->model()->shapeFactor();
		double GAvgElem=0.0; ///VolWeighted pore shape factor
		double VWGAvgElem=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWGAvgElem+= elms[i]->flowVolume()*elms[i]->model()->shapeFactor();
			GAvgElem+= elms[i]->model()->shapeFactor();
		}
		GAvgElem=GAvgElem/m_numElems;
		VWGAvgElem=VWGAvgElem/sumFlowVolume;
		//cout<<"m_numElems:\t"<<m_numElems<<endl;
		//cout<<"GAvgElem:\t"<<GAvgElem<<endl;
		//cout<<"VWGAvgElem:\t"<<VWGAvgElem<<endl;
		cout<<"GAvgElem   VWGAvgElem:  \t "<<GAvgElem<<"\t "<<VWGAvgElem<<endl;
	}      
}

void printRadiusStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadPore = elms[elms.size()/2]->model()->radius();
		double RadAvgPore=0.0; ///VolWeighted pore shape factor
		double VWRadAvgPore=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgPore+= elms[i]->flowVolume()*elms[i]->model()->radius();
			RadAvgPore+= elms[i]->model()->radius();
		}
		RadAvgPore=RadAvgPore/m_numPores;
		VWRadAvgPore=VWRadAvgPore/sumFlowVolume;
		//cout<<"m_numPores:\t"<<m_numPores<<endl;
		cout<<"RadAvgPore   VWRadAvgPore:  \t "<<RadAvgPore<<"\t "<<VWRadAvgPore<<endl;
	}


	{
		double m_numThroats=(*rockLattices).size()-2-m_numPores;
		vector< Element const * > elms((*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadThroat = elms[elms.size()/2]->model()->radius();
		double RadAvgThroat=0.0; ///VolWeighted pore shape factor
		double VWRadAvgThroat=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgThroat+= elms[i]->flowVolume()*elms[i]->model()->radius();
			RadAvgThroat+= elms[i]->model()->radius();
		}
		RadAvgThroat=RadAvgThroat/m_numThroats;
		VWRadAvgThroat=VWRadAvgThroat/sumFlowVolume;
		//cout<<"m_numThroats:\t"<<m_numThroats<<endl;
		//cout<<"RadAvgThroat:\t"<<RadAvgThroat<<endl;
		//cout<<"VWRadAvgThroat:\t"<<VWRadAvgThroat<<endl;
		cout<<"RadAvgThroat VWRadAvgThroat:\t "<<RadAvgThroat<<"\t "<<VWRadAvgThroat<<endl;
	}



	{
		double m_numElems=(*rockLattices).size()-2;
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		elms.insert(elms.end(), (*rockLattices).begin()+2+m_numPores, (*rockLattices).end());
		sort(elms.begin(), elms.end(), ElemRadCmpRed());
		//double medianRadElem = elms[elms.size()/2]->model()->radius();
		double RadAvgElem=0.0; ///VolWeighted pore shape factor
		double VWRadAvgElem=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			sumFlowVolume+=elms[i]->flowVolume();
			VWRadAvgElem+= elms[i]->flowVolume()*elms[i]->model()->radius();
			RadAvgElem+= elms[i]->model()->radius();
		}
		RadAvgElem=RadAvgElem/m_numElems;
		VWRadAvgElem=VWRadAvgElem/sumFlowVolume;
		//cout<<"m_numElems:\t"<<m_numElems<<endl;
		//cout<<"RadAvgElem:\t"<<RadAvgElem<<endl;
		//cout<<"VWRadAvgElem:\t"<<VWRadAvgElem<<endl;
		cout<<"RadAvgElem   VWRadAvgElem:  \t "<<RadAvgElem<<"\t "<<VWRadAvgElem<<endl;
	}      
}




void printAspectRatioStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		//sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianAspRatioPore = elms[elms.size()/2]->model()->shapeFactor();
		double AspRatioAvgPore=0.0; ///VolWeighted pore shape factor
		double VWAspRatioAvgPore=0.0; ///VolWeighted pore shape factor
		double VWAWAspRatioAvgPore=0.0; ///VolWeighted pore shape factor
		double VWQWAspRatioAvgPore=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			const Pore* pore = dynamic_cast< Pore const * >(elms[i]); 
			//if (pore->connectionNum()>0)
			double AvgAspectRatio=0.0;
			double AWAvgAspectRatio=0.0,sumAW=1.0e-31;
			double QWAvgAspectRatio=0.0,sumQW=1.0e-31;
			//double maxAspectRatio=0.0;
			double poreR=pore->model()->radius();
			for(int conn = 0; conn < pore->connectionNum(); ++conn)
			{
				const ElemModel* shape = pore->connection(conn)->model();
				double AspectRatio= shape->radius()/poreR;
				//Element *nextPore = currPore->get Connection Prop(conn, conductance, deltaGrav, fluid, resistivitySolve);
				AvgAspectRatio+= AspectRatio;
			    sumAW+=shape->area();  
			     AWAvgAspectRatio+= shape->area()*AspectRatio; 
			    sumQW+=shape->SPConductance(shape->area(),1.0);  
			     QWAvgAspectRatio+=shape->SPConductance(shape->area(),1.0)*AspectRatio;
			    //maxAspectRatio=max(maxAspectRatio,AspectRatio);
			}
				AvgAspectRatio/=pore->connectionNum()+1.0e-16;
			    AWAvgAspectRatio/=sumAW; 
			    QWAvgAspectRatio/=sumQW;
			sumFlowVolume+=elms[i]->flowVolume();
			VWAWAspRatioAvgPore+= elms[i]->flowVolume()*AWAvgAspectRatio;
			VWQWAspRatioAvgPore+= elms[i]->flowVolume()*QWAvgAspectRatio;
			VWAspRatioAvgPore+= elms[i]->flowVolume()*AvgAspectRatio;
			AspRatioAvgPore+= AvgAspectRatio;
		}
		AspRatioAvgPore=AspRatioAvgPore/m_numPores;
		VWAspRatioAvgPore=VWAspRatioAvgPore/sumFlowVolume;
		VWAWAspRatioAvgPore=VWAWAspRatioAvgPore/sumFlowVolume;
		VWQWAspRatioAvgPore=VWQWAspRatioAvgPore/sumFlowVolume;
		//cout<<"m_numPores:\t"<<m_numPores<<endl;
		cout<<"AspRatioAvgPore:  Arithmatic   VW   VWAW    VWQW :  \t "
			<<AspRatioAvgPore<<"\t "
			<<VWAspRatioAvgPore<<"\t "
			<<VWAWAspRatioAvgPore<<"\t "
			<<VWQWAspRatioAvgPore<<endl;
	}

   
}




void printCoordinaNumStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores)
{

	{
		vector< Element const * > elms((*rockLattices).begin()+1, (*rockLattices).begin()+1+m_numPores);
		//sort(elms.begin(), elms.end(), ElemGCmpRed());
		//double medianCoordNumPore = elms[elms.size()/2]->model()->shapeFactor();
		double CoordNumAvgPore=0.0; ///VolWeighted pore shape factor
		double VWCoordNumAvgPore=0.0; ///VolWeighted pore shape factor
		double VWAWCoordNumAvgPore=0.0; ///VolWeighted pore shape factor
		double VWQWCoordNumAvgPore=0.0; ///VolWeighted pore shape factor
		double sumFlowVolume=0.0;
		for(size_t i = 0; i < elms.size(); ++i)
		{
			const Pore* pore = dynamic_cast< Pore const * >(elms[i]); 
			//if (pore->connectionNum()>0)
			double AvgCoordNum=0.0;
			double AWAvgCoordNum=0.0,sumAW=0.0;
			double QWAvgCoordNum=0.0,sumQW=0.0;
			//double maxCoordNum=0.0;
			//double poreR=pore->model()->radius();
			for(int conn = 0; conn < pore->connectionNum(); ++conn)
			{
				const ElemModel* shape = pore->connection(conn)->model();
				//Element *nextPore = currPore->get Connection Prop(conn, conductance, deltaGrav, fluid, resistivitySolve);
				AvgCoordNum+= 1;
			    sumAW+=shape->area();  
			     AWAvgCoordNum+= shape->area()*shape->area(); 
			    sumQW+=shape->SPConductance(shape->area(),1.0);  
			     QWAvgCoordNum+=shape->SPConductance(shape->area(),1.0)*shape->SPConductance(shape->area(),1.0);
			    //maxCoordNum=max(maxCoordNum,CoordNum);
			}
				AvgCoordNum/=1;
			    AWAvgCoordNum=sumAW*sumAW/(AWAvgCoordNum+1.0e-63); 
			    QWAvgCoordNum=sumQW*sumQW/(QWAvgCoordNum+1.0e-63); 
			sumFlowVolume+=elms[i]->flowVolume();
			VWAWCoordNumAvgPore+= elms[i]->flowVolume()*AWAvgCoordNum;
			VWQWCoordNumAvgPore+= elms[i]->flowVolume()*QWAvgCoordNum;
			VWCoordNumAvgPore+= elms[i]->flowVolume()*AvgCoordNum;
			CoordNumAvgPore+= AvgCoordNum;
		}
		CoordNumAvgPore=CoordNumAvgPore/m_numPores;
		VWCoordNumAvgPore=VWCoordNumAvgPore/sumFlowVolume;
		VWAWCoordNumAvgPore=VWAWCoordNumAvgPore/sumFlowVolume;
		VWQWCoordNumAvgPore=VWQWCoordNumAvgPore/sumFlowVolume;
		//cout<<"m_numPores:\t"<<m_numPores<<endl;
		cout<<"CoordNumAvgPore:  Arithmatic   VW   VWAW    VWQW :  \t "
			<<CoordNumAvgPore<<"\t "
			<<VWCoordNumAvgPore<<"\t "
			<<VWAWCoordNumAvgPore<<"\t "
			<<VWQWCoordNumAvgPore<<endl;
	}

   
}



