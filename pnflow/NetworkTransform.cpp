#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
//#include <cassert>
//#include <ctime>
//#include <set>
#include <vector>
//#include <algorithm>
//#include <functional>
//#include <map>
using namespace std;

//#include "f2c.h"
//#include "sortedEvents.h"
//#include "threeSome.h"
#include "inputData.h"
//#include "node.h"


#include "Element.h"
#include "elem_Model.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"
#include "NetworkTransform.h"




void NetworkTransform::getModif(	vector<Pore*>& pores2BAltered, vector<Throat*>& throats2BAltered, const InputFile & input, istringstream& ins, const vector<Elem*> & elemans, size_t nBpPors)  {

	bool volBased(false), totalFraction(false), ffClust(true);
	double fraction(1.);
	string  spatialDistrib("rand");
	WeibulParam  clustDs; clustDs.correlation = "rand";


	if(ins.good())  {

		if(ins.good()) ins >> fraction;
		if(ins.good())  {	char vB('x');    ins >>vB;  	volBased = (vB=='T'||vB=='t'||vB=='V'||vB=='v');
			if (!(vB=='T'||vB=='t'||vB=='V'||vB=='v'||vB=='F'||vB=='f'||vB=='N'||vB=='n'))	{ out_<<" Error Wrong choice for \"volume-/number- based\" fraction, expected T/t/V/v or F/f/N/n"<<endl; exit(-1); }
		}
		if(ins.good()) 
		{ 	char tF('x');  ins >> tF;  totalFraction = (tF=='T'||tF=='t'); 
			if (!(tF=='T'||tF=='t'||tF=='F'||tF=='f'||tF=='O'||tF=='o'))	{ out_<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}
		}

		if(ins.good()) ins >> spatialDistrib;
		if (spatialDistrib[1] == 'o' || spatialDistrib[1] == 'O')  {
			clustDs.read(ins);
			clustDs.maxV+=0.999999;
			if (ins.fail()) out_ <<"Error: cant read cluster size "<<endl;
			input.Assert(clustDs.minV>1 && clustDs.maxV>=clustDs.minV, "Error: wrong cluster length");
			if (ins.good()) 
			{	char oInW('T'); ins >> oInW;		ffClust = (oInW == 'Y' || oInW == 'y' || oInW == 'T' || oInW == 't');
				if (!(oInW=='T'||oInW=='t'||oInW=='F'||oInW=='f'||oInW=='O'||oInW=='o'))	{ out_<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}
			}
		}
		else clustDs.correlation = spatialDistrib;
	}


	vector<Pore*> poresToSeed;poresToSeed.reserve(nBpPors);
	vector<int> clustDiams;
	{
		for(size_t i = 0; i < nBpPors; ++i)
			if (/*elemans[i]->exists(Ff) &&*/ dynamic_cast<Pore*>(elemans[i]))  poresToSeed.push_back(dynamic_cast<Pore*>(elemans[i]));


		Correlate<Pore, double(Elem::*)()const, &Elem::RRR>()(poresToSeed, clustDs.correlation, *this);
		clustDiams = randfield<int>(poresToSeed.size(),clustDs, *this);

	}

	pores2BAltered.reserve(poresToSeed.size());
	throats2BAltered.reserve(poresToSeed.size()*3);


	double ffInvadedVol(0.), TotalVol(0.), slctedVol(0.);
	int totElem(0);

	for(size_t i = 1; i < elemans.size(); ++i)  {
		//if (dynamic_cast<VoidElem*>(elemans[i]))
		if (i!=size_t(nBpPors-1))  {
			//if(elemans[i]->exists(Ff))
			{
				if(i < size_t(nBpPors)) ++totElem;                        // Only count pores
				ffInvadedVol += elemans[i]->flowVolume();
			}
			TotalVol += elemans[i]->flowVolume();
		}
	}

	double targetFrac = fraction;
	if (totalFraction)  targetFrac = min(fraction*TotalVol/ffInvadedVol, 1.);
	if (ffInvadedVol/TotalVol>1e-12)  {

		vector<int> clusterIndx(elemans.size()+2, 0);

		int numFracWetted(0), clusterIdx(1);
		if(spatialDistrib[1] == 'o' || spatialDistrib[1] == 'O')    // Spatial correlation approach
		{
			bool clct_more(true);


			targetFrac = ffClust ? targetFrac: 1.-targetFrac;

			while(clct_more)  {
				set< Elem * > frontier;
				Pore* elem = poresToSeed.back(); poresToSeed.pop_back();
				int   clustDiam = clustDiams.back(); clustDiams.pop_back();
					//int clustDiam = (*this).weibull(clustDiam1,clustDiam2, deltaClustD, etaClustD) +0.5;

				if (clusterIndx[elem->index()]==0)  {

					frontier.insert(elem);
					typedef set<Elem*>::iterator ItrSet;
					int frontExpansions(0);
					while(!frontier.empty() && clct_more && frontExpansions < clustDiam)  // Front expansions will go over thraots as well as pores, wheras
					{                                                                           // cluster diam is wrt pores  => Don't divide by 2
						set<Elem*> oldFrontier = frontier;
						frontier.clear();
						++frontExpansions;
						for(ItrSet itrElem = oldFrontier.begin(); itrElem != oldFrontier.end(); ++itrElem)  {

							if( clusterIndx[(*itrElem)->index()] != clusterIdx && clct_more)  {
								if(clusterIndx[(*itrElem)->index()] == 0)  {// && (*itrElem)->exists(Ff))
									clusterIndx[(*itrElem)->index()] = clusterIdx;                    /// Set flag
									slctedVol += (*itrElem)->flowVolume();
									if(dynamic_cast<Pore*>(*itrElem)) ++numFracWetted;   /// Only count pores

									clct_more = volBased ? slctedVol/ffInvadedVol < targetFrac :  double(numFracWetted)/totElem < targetFrac;

									if(!clct_more)			break;
								}

								for(int i = 0; i < (*itrElem)->nCncts(); ++i)
									if ( clusterIndx[(*itrElem)->neib(i)->index()] != clusterIdx
										  && (dynamic_cast<Throat*>((*itrElem)->neib(i)) || dynamic_cast<Pore*>((*itrElem)->neib(i))) )
										frontier.insert((*itrElem)->neib(i));
							}
						}
					}
					++clusterIdx;
				}
			}
			numFracWetted = 0;
			slctedVol = 0.;
			if(ffClust)  {
				for(size_t i = 0; i < elemans.size(); ++i)  {
					if( clusterIndx[elemans[i]->index()] > 0)  {// && elemans[i]->exists(Ff))
						if( i < nBpPors && dynamic_cast<Pore*>(elemans[i]))  {
							++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
							pores2BAltered.push_back(dynamic_cast<Pore*>(elemans[i]));
						}
						else if( dynamic_cast<Throat*>(elemans[i]))  {
							++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
							throats2BAltered.push_back(dynamic_cast<Throat*>(elemans[i]));
						}
					}
				}
			}
			else
			{
				for(size_t i = 0; i < elemans.size(); ++i)  {
					if( clusterIndx[elemans[i]->index()] == 0)  {// && elemans[i]->exists(Ff))
						if( i < nBpPors && dynamic_cast<Pore*>(elemans[i]))  {
							++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
							pores2BAltered.push_back(dynamic_cast<Pore*>(elemans[i]));
						}
						else if( dynamic_cast<Throat*>(elemans[i]))  {
							++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
							throats2BAltered.push_back(dynamic_cast<Throat*>(elemans[i]));
						}
					}
				}
			}
			out_<< "Number of correlated regions: " << clusterIdx-1                         << endl;

		}
		else ///. spatialDistrib == rand
		{
			//vector< Pore*> > toBeAltered;

			//for(int i = 2; i < nBpPors; ++i)
			//{
				//if(elemans[i]->exists(Ff) && dynamic_cast<Pore*>(elemans[i]))
				//{
					//pair<double, Pore*> rockEntry;
					//rockEntry.second = dynamic_cast<Pore*>(elemans[i]);
					//if(spatialDistrib[1] == 'a' || spatialDistrib[1] == 'A')        ///. rand // Random
					//{
						//rockEntry.first = (*this).rand01();
					//}
					//else if(spatialDistrib[0] == 'r' || spatialDistrib[0] == 'R')   ///. rMax/rMin // Radius
						//rockEntry.first = rockEntry.second->RRR();
					//else
					//{
						//out_ << endl << endl << "Error: Did not recognize fractional wetting model: " << spatialDistrib  << endl  << endl << endl;
						//exit(-1);
					//}
					//toBeAltered.push_back(rockEntry);
				//}
			//}
			//if(spatialDistrib[2] == 'a' || spatialDistrib[2] == 'A') ///. rM'a'x
				//sort(toBeAltered.begin(), toBeAltered.end(), FracWettInc());
			//else
				//sort(toBeAltered.begin(), toBeAltered.end(), FracWettDec()); ///. rMin/rand

			int targetNum(fraction*poresToSeed.size());

			while(!poresToSeed.empty() &&
				((volBased && slctedVol/ffInvadedVol < targetFrac) || (!volBased && numFracWetted < targetNum)))  {
				++numFracWetted;
				Pore* elem = poresToSeed.back(); poresToSeed.pop_back();
				slctedVol += elem->flowVolume();
				pores2BAltered.push_back(elem);
				int pind = elem->index();
				clusterIndx[pind]=1;
				for(int t = 0; t < elem->nCncts(); ++t)  {
					Throat* throat = dynamic_cast<Throat*>(elem->neib(t));
					if( throat && clusterIndx[throat->index()] == 0) // throat->exists(Ff) &&
					{
						double randNum = (*this).rand01();
						if(randNum < fraction)  {
							throats2BAltered.push_back(throat);
							clusterIndx[throat->index()]=1;
							slctedVol += throat->flowVolume();
						}
					}
				}
			}
		}
		out_
			<< "selected volume (fraction): " << slctedVol/ffInvadedVol   << endl
			<< "selected volume (of total net pore volume): " << slctedVol/TotalVol   << endl
			<< "Number of pores selected: " << numFracWetted                            << endl
			<< "================================================================="  << endl;

	}
	else out_
		<< "\nWarning:  nothing selected, too low volume fraction: " << ffInvadedVol/TotalVol             << endl
		<< "=================================================================\n"      << endl;




	//setContactAngles(pores2BAltered, throats2BAltered, minVal, maxVal, deltaExp, etaExp, cntctAngModel, modelTwoSepAng,  globalCorrelation, nBpPors, (*this));


}



void NetworkTransform::modify(const InputFile & input, vector<Elem*> & elemans, size_t nBpPors, double boxVolume)//, dbl3 boxSize
{
	///. modify network
	//IncreaseClay     0.3   0.3  -0.2  -3.   rand  ;
	//ConvertToClay    0.3   0.3  -0.2  -3.   rand  T  ;
	//ScaleRadius     0.9   1.1  -0.2  -3.   rand;
	//ScalePoreRadius     0.9   1.1  -0.2  -3.   rand;
	//ScaleThroatRadius     0.9   1.1  -0.2  -3.   rand;
	//MatchMICP     MICPCurve.txt  0.485  45;

	istringstream ins;
	if (input.giv("AddClay", ins))  {  out_<<"AddClay: Adding clay porosity: "<<endl;

		WeibulParam  wbdist(ins);
		input.Assert(!ins.fail(), "AddClay", "wrong data",true);

		vector<Throat*> throats;throats.reserve(elemans.size()-nBpPors);
		for(size_t i = nBpPors; i < elemans.size(); ++i)
			if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));


		vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
		getModif( pores2BAltered,  throats2BAltered,  input,  ins, elemans, nBpPors);

		{	vector<double> rands = randfield<double>(pores2BAltered.size(),wbdist, *this);
			Correlate<Pore, double(Elem::*)()const, &Elem::RRR>()(pores2BAltered, wbdist.correlation, *this);

			for(size_t i = 0; i < pores2BAltered.size(); ++i)  {
				Pore* elm = pores2BAltered[i];
				double clayFrac = rands[i];
				elm->adjustVolume(-1., clayFrac*elm->flowVolume()/(1.-clayFrac));
			}
		}

		{	vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
			Correlate<Throat, double(Elem::*)()const, &Elem::RRR>()(throats2BAltered, wbdist.correlation, *this);

			for(size_t i = 0; i < throats2BAltered.size(); ++i)  {
				Throat* elm = throats2BAltered[i];
				double clayFrac = rands[i];
				elm->adjustVolume(-1., clayFrac*elm->flowVolume()/(1.-clayFrac));
			}
		}
	}

	if (input.giv("FillWithClay", ins))  {
		out_<<"Filling void space by clay : "<<endl;
		WeibulParam  wbdist(ins);
		input.Assert(!ins.fail(), "FillWithClay", "wrong data",true);

		vector<Throat*> throats;throats.reserve(elemans.size()-nBpPors);
		for(size_t i = nBpPors; i < elemans.size(); ++i)
			if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));

		vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
		getModif( pores2BAltered,  throats2BAltered,  input,  ins, elemans, nBpPors);

		vector<double> clayporesRad = randfield<double>(pores2BAltered.size(),wbdist, *this);
		Correlate<Pore, double(Elem::*)()const, &Elem::RRR>()(pores2BAltered, wbdist.correlation, *this);

		vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
		Correlate<Throat, double(Elem::*)()const, &Elem::RRR>()(throats2BAltered, wbdist.correlation, *this);

		double totalFlowVolume(0.), totalClayVolSum(0.);
		for(size_t i = 0; i < elemans.size(); ++i)  {
			totalFlowVolume -= elemans[i]->flowVolume();
			totalClayVolSum -= elemans[i]->clayVolume();
		}
		//double selectedFlowVolume(0.), selectedClayVolSum(0.);
		//for(size_t i = 0; i < throats2BAltered.size(); ++i)
		//{
			//double totalFlowVolume = throats2BAltered[i]->flowVolume();
			//double totalClayVolSum = throats2BAltered[i]->clayVolume();
		//}

		for(size_t i = 0; i < throats2BAltered.size(); ++i)  {
			Throat* trot = throats2BAltered[i];
			double clayFrac = rands[i];
			double newNetVol = trot->flowVolume()*(1.-clayFrac);
			double newClayVol = trot->flowVolume()*(clayFrac) +trot->clayVolume();
			trot->adjustVolume(newNetVol, newClayVol);
			transform(trot, sqrt(1.-clayFrac));
			//for (int i=0;i<2;++i)
			//{
				//Pore* por = dynamic_cast<Pore*>(trot->neib(i));
				//double newNetVol = por->flowVolume()*(1.-clayFrac);
				//double newClayVol = por->flowVolume()*(clayFrac) +por->clayVolume();
				//por->adjustVolume(newNetVol, newClayVol, netVolSum, clayVolSum);
				//transform(por, 0.5+0.5*sqrt(1.-clayFrac));
			//}
		}

		for(size_t i = 0; i < pores2BAltered.size(); ++i)  {
			Pore* por = pores2BAltered[i];
			double clayFrac = clayporesRad[i];
			double newNetVol = por->flowVolume()*(1.-clayFrac);
			double newClayVol = por->flowVolume()*(clayFrac) +por->clayVolume();
			por->adjustVolume(newNetVol, newClayVol);
			transform(por, sqrt(1.-clayFrac));
		}

		for(size_t i = nBpPors; i < elemans.size(); ++i)
			if (dynamic_cast<Throat*>(elemans[i]))	  fixRadius(dynamic_cast<Throat*>(elemans[i]));

			//exit(-1);
	}

	double radiusScale;
	if (input.giv("scaleRadius", radiusScale))  {
		out_<<" scalling radius "<<endl;
		WeibulParam  wbdist(ins);
		input.Assert(!ins.fail(), "scaleRadius", "wrong data",true);

		vector<Throat*> throats;throats.reserve(elemans.size()-nBpPors);
		for(size_t i = nBpPors; i < elemans.size(); ++i)
			if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));

		vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
		getModif( pores2BAltered,  throats2BAltered,  input,  ins, elemans, nBpPors);

		//vector<double> claypores = randfield<double>(pores2BAltered.size(),minVal, maxVal, deltaExp, etaExp, *this);
		//Correlate<Pore, double(Elem::*)()const, &Elem::RRR>()(pores2BAltered, wbdist.correlation, *this);

		vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
		Correlate<Throat, double(Elem::*)()const, &Elem::RRR>()(throats2BAltered, wbdist.correlation, *this);

		for(size_t i = 0; i < throats2BAltered.size(); ++i)  {
			Throat* trot = throats2BAltered[i];
			double scaleFact = rands[i];
			transform(trot, scaleFact);
		}
	}
	if (input.giv("scaleRadiusMICP", radiusScale))  {
		out_<<"scaleMICP: ToBeImplemented"<<endl;
	}
}


inline double tableLookUpX(double yVal, const vector<dbl2>& lookupTable)  {
	if(yVal <= lookupTable.front().b)
		return lookupTable.front().a;
	else if(yVal >= lookupTable.back().b)
		return lookupTable.back().a;

	vector<dbl2>::const_iterator itr = lookupTable.begin()+1;
	while(itr->b < yVal && itr != lookupTable.end()) ++itr;
	assert(itr->b >= yVal && (itr-1)->b < yVal);
	double xOne((itr-1)->a), xTwo(itr->a), yOne((itr-1)->b), yTwo(itr->b);
	double logResX = log10(xOne)+log10(xTwo/xOne)*(yVal-yOne)/(yTwo-yOne);
	return pow(10., logResX);
}


inline double tableLookUpY(double xVal, const vector<dbl2>& lookupTable)  {
	vector<dbl2>::const_iterator itr = lookupTable.begin()+1;
	while(itr->a > xVal && itr != lookupTable.end()) ++itr;
	assert(itr->a <= xVal && (itr-1)->a > xVal);
	double xOne((itr-1)->a), xTwo(itr->a), yOne((itr-1)->b), yTwo(itr->b);
	return yOne+log10(xVal/xOne)*(yTwo-yOne)/log10(xTwo/xOne);
}



void setShapeFactFromFile(vector<Elem*>& elems, const string& gFileName,
							  double lowCutOffDiam, double highCutOffDiam)  {
	ifstream fin(gFileName);
	ensure (fin, "Unable to open G distribution file " + gFileName, -1);

	vector<dbl2> targetDist;        // g (decr) vs. numeric fraction (incr)
	dbl2 funcPt, oldPt(1., -1.);
	while(fin >> funcPt.a)  {
		fin >> funcPt.b;
		ensure(funcPt.b >= oldPt.b && funcPt.a <= oldPt.a, 
				"Fractional index should increase and value should decrease, monotonically.", -1 );
		if(funcPt.b != oldPt.b)	targetDist.push_back(funcPt);
		oldPt = funcPt;
	}

	if(lowCutOffDiam > targetDist.back().a && lowCutOffDiam < targetDist.front().a)  {
		double volCutoff = tableLookUpY(lowCutOffDiam, targetDist);
		while(!targetDist.empty() && targetDist.back().a < lowCutOffDiam)
			targetDist.pop_back();

		for(size_t i = 0; i < targetDist.size(); ++i)
			targetDist[i].b /= volCutoff;
		dbl2 endPoint(lowCutOffDiam, 1.);
		targetDist.push_back(endPoint);
	}

	if(highCutOffDiam > targetDist.back().a && highCutOffDiam < targetDist.front().a)  {
		double volCutoff = tableLookUpY(highCutOffDiam, targetDist);
		vector<dbl2> newTargetDist;
		dbl2 pt(highCutOffDiam, 0.);
		newTargetDist.push_back(pt);

		for(size_t i = 0; i < targetDist.size(); ++i)  {
			if(targetDist[i].a < highCutOffDiam)  {
				dbl2 pt(targetDist[i].a,
					(targetDist[i].b-volCutoff)/(1.-volCutoff));
				newTargetDist.push_back(pt);
			}
		}
		targetDist = newTargetDist;
	}

	if(targetDist.empty() || targetDist.front().b != 0. || targetDist.back().b != 1.)  {
		cerr << endl
			<< "================================================="   << endl
			<< "Error: Target distribution contains errors."         << endl
			<< " The cumulative volume should go from 0. to 1.."   << endl
			<< "================================================="   << endl;        exit(-1);
	}

	int totNum(0);
	for(size_t i = 0; i < elems.size(); ++i)  {
		if(elems[i]->model()->shapeFactor() <=  sqrt(3.)/36.)
			++totNum;
	}

	int runIdx(0);
	for(size_t j = 0; j < elems.size(); ++j)  {
		if(elems[j]->model()->shapeFactor() <=  sqrt(3.)/36.)  {
			++runIdx;
			double frac = static_cast< double >(runIdx)/static_cast< double >(totNum);
			double newG = tableLookUpX(frac, targetDist);
			elems[j]->ChModel()->setShapeFactor(newG);
		}
	}
}

void modifyShapeFactor(int model, const string& options, vector<Elem*>& elems,
							   const string& fileName, bool writeToFile, int numPts)  {
	istringstream in(options);

	if(model == 0)  {
		// Do nothing option
	}
	else if(model == 1)  {
		string distFile;
		double lowCutOffG, highCutOffG;
		in >> distFile >> lowCutOffG >> highCutOffG;
		setShapeFactFromFile(elems, distFile, lowCutOffG, highCutOffG);
	}
	else if(model == 3)  {
		double multFactor;
		in >> multFactor;
		for(size_t i = 0; i < elems.size(); ++i)  {
			double gFact = min(sqrt(3.)/36.-1.0E-6, elems[i]->model()->shapeFactor()*multFactor);
			elems[i]->ChModel()->setShapeFactor(gFact);
		}
	}
	else if(model == 4)     // Not a very good model this.....
	{                       // NEED TO CHANGE THIS
		double aConst;
		in >> aConst;
		sort(elems.begin(), elems.end(), ElemGCmpRed());
		double medianG = elems[elems.size()/2]->model()->shapeFactor();
		for(size_t i = 0; i < elems.size(); ++i)  {
			double gFact = aConst*(elems[i]->model()->shapeFactor()-medianG)+elems[i]->model()->shapeFactor();
			if(gFact > sqrt(3.)/36.) gFact = sqrt(3.)/36.-1.0E-6;
			if(gFact <=  0.) gFact = 0.0001;
			elems[i]->ChModel()->setShapeFactor(gFact);
		}
	}
	else if(model == 5)  {
		double minVal, maxVal, deltaExp, etaExp;
		in >> minVal >> maxVal >> deltaExp >> etaExp;
		if(maxVal > sqrt(3.)/36.)  {
			cout<< endl
				<< "==========================================================="    << endl
				<< "Warning: Maximum shape factor should be less than 0.048113 "     << endl
				<< "which represents a equilateral triangle."                       << endl
				<< "==========================================================="    << endl
				<< endl;
			maxVal = sqrt(3.)/36. - 1.0E-6;
		}

		for(size_t i = 0; i < elems.size(); ++i)  {
			elems[i]->ChModel()->setShapeFactor(weibull(minVal, maxVal, deltaExp, etaExp));
		}
	}
	else
	{
		cerr << "========================================== " << endl
			 << "Error: Uknown shape factor modifier option" << endl
			 << "========================================== " << endl;        exit(-1);
	}

	if(writeToFile)  {
		ofstream fout(fileName);
		fout.flags(ios::showpoint);
		fout.flags(ios::scientific);
		fout.precision(6);

		int totNum(0);
		for(size_t i = 0; i < elems.size(); ++i)  {
			if(elems[i]->model()->shapeFactor() <=  sqrt(3.)/36.)
				++totNum;
		}

		sort(elems.begin(), elems.end(), ElemGCmpRed());

		double step(0.), stepDelta(1./numPts);
		int runIdx(0);
		for(size_t j = 0; j < elems.size(); ++j)  {
			if(elems[j]->model()->shapeFactor() <=  sqrt(3.)/36.)  {
				++runIdx;
				double frac = static_cast< double >(runIdx)/static_cast< double >(totNum);
				if(frac > step || frac == 1.)  {
					fout << elems[j]->model()->shapeFactor() << ",   "  << frac << endl;
					step += stepDelta;
				}
			}
		}
	}
}

void modifyInscribedRadii(int model, const string& options, vector<Elem*>& elems,
								  const string& fileName, bool writeToFile, int numPts, bool forPores)  {   ///. MODIFY_RAD_DIST
	istringstream in(options);
	sort(elems.begin(), elems.end(), ElemRadCmpRed());

	cout<< "\nmodify R, model == "<< model << endl << endl;

	if(model == 0)  {
		// Do nothing option
	}
	else if(model == 1)  {
		string distFile;
		double lowCutOffDiam, highCutOffDiam;
		in >> distFile >> lowCutOffDiam >> highCutOffDiam;
		setRadiiFromFile(elems, distFile, lowCutOffDiam*1.0E-6, highCutOffDiam*1.0E-6, forPores);
	}
	else if(model == 2)  {
		string whichAspectRatioType;
		double aspectRatioMin, aspectRatioMax, delta, eta;
		in  >> aspectRatioMin >> aspectRatioMax >> delta >> eta;
		for(size_t i = 0; i < elems.size(); ++i)  {
			cout<<"not supported AQ"<<endl;
			//double aspectRatio(-1.);
			//if(aspectRatioMin > 0. && aspectRatioMax > 0.)
				//aspectRatio = weibull(aspectRatioMin, aspectRatioMax, delta, eta);
			//elems[i]->setRadiusFromAspectRatio(aspectRatio);
		}
	}
	else if(model == 3)  {
		double multFactor;
		in >> multFactor;
		cout<< "multFactor "<< multFactor << endl;

		for(size_t i = 0; i < elems.size(); ++i)  {
			double maxRadius(elems[i]->model()->RRR()*multFactor);
			if(forPores)  {
				for(int j = 0; j < elems[i]->nCncts(); ++j)  {
					maxRadius = max(maxRadius, elems[i]->neib(j)->model()->RRR());
				}
			}
			elems[i]->ChModel()->setRadius(maxRadius);
		}
	}
	else if(model == 4)  {
		double aConst;
		bool applyAboveMedian, applyBelowMedian;
		char above, below;
		in >> aConst >> above >> below;
		applyAboveMedian = (above == 't' || above == 'T');
		applyBelowMedian = (below == 't' || below == 'T');

		vector<dbl2> poreSizeDist;
		double totVolume(0.);
		for(size_t i = 0; i < elems.size(); ++i)  {
			totVolume += elems[i]->flowVolume();
			dbl2 entry(elems[i]->model()->RRR(), totVolume);
			poreSizeDist.push_back(entry);
		}
		double volBasedAvr = tableLookUpX(totVolume/2., poreSizeDist);
		for(size_t j = 0; j < elems.size(); ++j)  {
			double rad(elems[j]->model()->RRR());
			if((rad > volBasedAvr && applyAboveMedian) || (rad < volBasedAvr && applyBelowMedian))
				rad = pow(elems[j]->model()->RRR(), aConst)/pow(volBasedAvr, aConst-1.);
			elems[j]->ChModel()->setRadius(rad);
		}
	}
	else if(model == 5)  {
		double minVal, maxVal, deltaExp, etaExp;
		in >> minVal >> maxVal >> deltaExp >> etaExp;
		minVal *= 1.0E-6;
		maxVal *= 1.0E-6;
		for(size_t i = 0; i < elems.size(); ++i)  {
			elems[i]->ChModel()->setRadius(weibull(minVal, maxVal, deltaExp, etaExp));
		}
	}
	else if(model == 6)  {
		double a, b;
		in >> a >> b ;
		a *= 1.0E-6;
		b *= 1.0E-6;
		for(auto& elm:elems)
			elm->ChModel()->setRadius( elm->model()->RRR()-a*tanh( elm->model()->RRR()/b )+0.0001*a );
	}
	else
	{
		cerr << "==================================== " << endl
			 << "Error: Uknown radius modifier option" << endl
			 << "==================================== " << endl;        exit(-1);
	}

	if(writeToFile)  {
		///. write old distribution
		ofstream fout(fileName);
		fout.flags(ios::showpoint);
		fout.flags(ios::scientific);
		fout.precision(6);

		double cumVol(0.), totCumVol(0.);
		sort(elems.begin(), elems.end(), ElemRadCmpRed());
		for(auto& elm:elems)  totCumVol += elm->flowVolume();
		fout << "diameter       , cumulative volume"<< endl;
		double step(0.), stepDelta(1./numPts);
		for(auto& elm:elems)  {
			cumVol += elm->flowVolume();
			double poreVolFrac = cumVol/totCumVol;
			if(poreVolFrac > step || poreVolFrac == 1.)  {
				fout << setw(15) << elm->model()->RRR()*2.0E6 << "," << setw(15) << poreVolFrac << endl;
				step += stepDelta;
			}
		}

		///. write old distribution
		fout << endl;
		fout << "diameter       , pdf"<<endl;
		unsigned int elemFrequency(0);
		double  bin_radius(elems[0]->model()->RRR()*(1.-1./numPts));
		stepDelta = bin_radius/numPts;
		for(size_t j = 0; j < elems.size(); ++j)  {
			//double poreVolFrac = cumVol/totCumVol;
			if(elems[j]->model()->RRR() < bin_radius || j == elems.size()-1 )  {
				fout << setw(15) << (bin_radius+stepDelta/2.)*1.0E6 << "," << setw(15) << (double)elemFrequency/elems.size()/stepDelta/1.0E6 << endl;
				bin_radius -= stepDelta;
				elemFrequency = 0;
			}
			elemFrequency ++;
		}

	}
}

void setRadiiFromFile(vector<Elem*>& elems, const string& radiiFileName, double lowCutOffDiam, double highCutOffDiam, bool forPores)  {
	ifstream fin(radiiFileName);
	if (!fin)  {
		cerr << "==========================================================="       << endl
			 << "Error: Unable to open radius distribution file " << radiiFileName  << endl
			 << "==========================================================="       << endl;
		exit( -1 );
	}

	vector<dbl2> targetDist;        // diameter (decr) vs. pore vol fraction (incr)
	dbl2 funcPt, oldPt(1., -1.);
	while(fin >> funcPt.a)  {
		fin >> funcPt.b;
		funcPt.a *= 1.0E-6;                    // Convert from micro.m to m
		if(funcPt.b < oldPt.b || funcPt.a >= oldPt.a)  {
			cerr << "==========================================================="   << endl
				<< "Error: Fractional pore volume should increase monotonically."   << endl
				<< "==========================================================="    << endl;
			exit( -1 );
		}
		else if(funcPt.b != oldPt.b)  {
			targetDist.push_back(funcPt);
		}
		oldPt = funcPt;
	}

	if(lowCutOffDiam > targetDist.back().a && lowCutOffDiam < targetDist.front().a)  {
		double volCutoff = tableLookUpY(lowCutOffDiam, targetDist);
		while(!targetDist.empty() && targetDist.back().a < lowCutOffDiam)
			targetDist.pop_back();

		for(size_t i = 0; i < targetDist.size(); ++i)
			targetDist[i].b /= volCutoff;
		dbl2 endPoint(lowCutOffDiam, 1.);
		targetDist.push_back(endPoint);
	}

	if(highCutOffDiam > targetDist.back().a && highCutOffDiam < targetDist.front().a)  {
		double volCutoff = tableLookUpY(highCutOffDiam, targetDist);
		vector<dbl2> newTargetDist;
		dbl2 pt(highCutOffDiam, 0.);
		newTargetDist.push_back(pt);

		for(size_t i = 0; i < targetDist.size(); ++i)  {
			if(targetDist[i].a < highCutOffDiam)  {
				dbl2 pt(targetDist[i].a,
					(targetDist[i].b-volCutoff)/(1.-volCutoff));
				newTargetDist.push_back(pt);
			}
		}
		targetDist = newTargetDist;
	}

	if(targetDist.empty() || targetDist.front().b != 0. || targetDist.back().b != 1.)  {
		cerr << endl
			<< "============================================="        << endl
			<< "Error: Target distribution contains errors."          << endl
			<< "The cumulative volume should go from 0. to 1.."     << endl
			<< targetDist.size() << "  " << targetDist.back().b  << endl
			<< "============================================="    << endl;    exit(-1);
	}

	double totElemVol(0.), runningVol(0.);
	for(size_t i = 0; i < elems.size(); ++i)  {
		totElemVol += elems[i]->flowVolume();
	}

	for(size_t j = 0; j < elems.size(); ++j)  {
		runningVol += elems[j]->flowVolume();
		double poreVolFrac = runningVol/totElemVol;
		double maxRadius(tableLookUpX(poreVolFrac, targetDist)/2.);
		if(forPores)  {
			for(int k = 0; k < elems[j]->nCncts(); ++k)  {
				maxRadius = max(maxRadius, elems[j]->neib(k)->model()->RRR());
			}
		}
		elems[j]->ChModel()->setRadius(maxRadius);
	}
}










