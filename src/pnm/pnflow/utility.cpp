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
//#include "fluid.h"
//#include "apex.h"

#include "Element.h"
//#include "solver.h"

#include "elem_Model.h"
#include "polygon.h"
//#include "elem_hetroPorous.h"
//#include "elem_porous.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"

#include "utility.h"


void printInfo(Element & elem)
{
	cout<<" index: "<<elem.orenIndex() <<"lattice index: "<<elem.latticeIndex()<<endl;
}



void setShapeFactFromFile(vector<Element*>& elems, const string& gFileName,
                              double lowCutOffDiam, double highCutOffDiam)
{
    ifstream fin(gFileName.c_str());
    if (!fin)
    {
        cerr << "==========================================================="   << endl
             << "Error: Unable to open G distribution file " << gFileName       << endl
             << "==========================================================="   << endl;
        exit( -1 );
    }

    vector< pair< double, double > > targetDist;        // g (dec) vs. numeric fraction (inc)
    pair< double, double > funcPt, oldPt(1.0, -1.0);
    while(fin >> funcPt.first)
    {
        fin >> funcPt.second;
        if(funcPt.second < oldPt.second || funcPt.first > oldPt.first)
        {
            cerr << "==========================================================="   << endl
                << "Error: Fractional index should increase monotonically."         << endl
                << "==========================================================="    << endl;
            exit( -1 );
        }
        else if(funcPt.second != oldPt.second)
        {
            targetDist.push_back(funcPt);
        }
        oldPt = funcPt;
    }

    if(lowCutOffDiam > targetDist.back().first && lowCutOffDiam < targetDist.front().first)
    {
        double volCutoff = tableLookUpY(lowCutOffDiam, targetDist);
        while(!targetDist.empty() && targetDist.back().first < lowCutOffDiam)
            targetDist.pop_back();

        for(size_t i = 0; i < targetDist.size(); ++i)
            targetDist[i].second /= volCutoff;
        pair< double, double > endPoint(lowCutOffDiam, 1.0);
        targetDist.push_back(endPoint);
    }

    if(highCutOffDiam > targetDist.back().first && highCutOffDiam < targetDist.front().first)
    {
        double volCutoff = tableLookUpY(highCutOffDiam, targetDist);
        vector< pair< double, double > > newTargetDist;
        pair< double, double > pt(highCutOffDiam, 0.0);
        newTargetDist.push_back(pt);

        for(size_t i = 0; i < targetDist.size(); ++i)
        {
            if(targetDist[i].first < highCutOffDiam)
            {
                pair< double, double > pt(targetDist[i].first,
                    (targetDist[i].second-volCutoff)/(1.0-volCutoff));
                newTargetDist.push_back(pt);
            }
        }
        targetDist = newTargetDist;
    }

    if(targetDist.empty() || targetDist.front().second != 0.0 || targetDist.back().second != 1.0)
    {
        cerr << endl
            << "================================================="   << endl
            << "Error: Target distribution contains errors."         << endl
            << " The cumulative volume should go from 0.0 to 1.0."   << endl
            << "================================================="   << endl;        exit(-1);
    }

    int totNum(0);
    for(size_t i = 0; i < elems.size(); ++i)
    {
        if(elems[i]->model()->shapeFactor() <=  sqrt(3.0)/36.0)
            ++totNum;
    }

    int runIdx(0);
    for(size_t j = 0; j < elems.size(); ++j)
    {
        if(elems[j]->model()->shapeFactor() <=  sqrt(3.0)/36.0)
        {
            ++runIdx;
            double frac = static_cast< double >(runIdx)/static_cast< double >(totNum);
            double newG = tableLookUpX(frac, targetDist);
            elems[j]->ChModel()->setShapeFactor(newG);
        }
    }
}

void modifyShapeFactor(int model, const string& options, vector<Element*>& elems,
                               const string& fileName, bool writeToFile, int numPts)
{
    istringstream in(options);

    if(model == 0)
    {
        // Do nothing option
    }
    else if(model == 1)
    {
        string distFile;
        double lowCutOffG, highCutOffG;
        in >> distFile >> lowCutOffG >> highCutOffG;
        setShapeFactFromFile(elems, distFile, lowCutOffG, highCutOffG);
    }
    else if(model == 3)
    {
        double multFactor;
        in >> multFactor;
        for(size_t i = 0; i < elems.size(); ++i)
        {
            double gFact = min(sqrt(3.0)/36.0-1.0E-6, elems[i]->model()->shapeFactor()*multFactor);
            elems[i]->ChModel()->setShapeFactor(gFact);
        }
    }
    else if(model == 4)     // Not a very good model this.....
    {                       // NEED TO CHANGE THIS
        double aConst;
        in >> aConst;
        sort(elems.begin(), elems.end(), ElemGCmpRed());
        double medianG = elems[elems.size()/2]->model()->shapeFactor();
        for(size_t i = 0; i < elems.size(); ++i)
        {
            double gFact = aConst*(elems[i]->model()->shapeFactor()-medianG)+elems[i]->model()->shapeFactor();
            if(gFact > sqrt(3.0)/36.0) gFact = sqrt(3.0)/36.0-1.0E-6;
            if(gFact <=  0.0) gFact = 0.0001;
            elems[i]->ChModel()->setShapeFactor(gFact);
        }
    }
    else if(model == 5)
    {
        double minVal, maxVal, deltaExp, etaExp;
        in >> minVal >> maxVal >> deltaExp >> etaExp;
        if(maxVal > sqrt(3.0)/36.0)
        {
            cout<< endl
                << "==========================================================="    << endl
                << "Warning: Maximum shape factor should be less than 0.048113 "     << endl
                << "which represents a equilateral triangle."                       << endl
                << "==========================================================="    << endl
                << endl;
            maxVal = sqrt(3.0)/36.0 - 1.0E-6;
        }

        for(size_t i = 0; i < elems.size(); ++i)
        {
            elems[i]->ChModel()->setShapeFactor(weibull(minVal, maxVal, deltaExp, etaExp));
        }
    }
    else
    {
        cerr << "========================================== " << endl
             << "Error: Uknown shape factor modifier option" << endl
             << "========================================== " << endl;        exit(-1);
    }

    if(writeToFile)
    {
        ofstream fout(fileName.c_str());
        fout.flags(ios::showpoint);
        fout.flags(ios::scientific);
        fout.precision(6);

        int totNum(0);
        for(size_t i = 0; i < elems.size(); ++i)
        {
            if(elems[i]->model()->shapeFactor() <=  sqrt(3.0)/36.0)
                ++totNum;
        }

        sort(elems.begin(), elems.end(), ElemGCmpRed());

        double step(0.0), stepDelta(1.0/numPts);
        int runIdx(0);
        for(size_t j = 0; j < elems.size(); ++j)
        {
            if(elems[j]->model()->shapeFactor() <=  sqrt(3.0)/36.0)
            {
                ++runIdx;
                double frac = static_cast< double >(runIdx)/static_cast< double >(totNum);
                if(frac > step || frac == 1.0)
                {
                    fout << elems[j]->model()->shapeFactor() << ",   "  << frac << endl;
                    step += stepDelta;
                }
            }
        }
    }
}

void modifyInscribedRadii(int model, const string& options, vector<Element*>& elems,
                                  const string& fileName, bool writeToFile, int numPts, bool forPores)
{   ///. MODIFY_RAD_DIST
    istringstream in(options);
    sort(elems.begin(), elems.end(), ElemRadCmpRed());

	cout<< "\nmodify R, model == "<< model << endl << endl;

    if(model == 0)
    {
        // Do nothing option
    }
    else if(model == 1)
    {
        string distFile;
        double lowCutOffDiam, highCutOffDiam;
        in >> distFile >> lowCutOffDiam >> highCutOffDiam;
        setRadiiFromFile(elems, distFile, lowCutOffDiam*1.0E-6, highCutOffDiam*1.0E-6, forPores);
    }
    else if(model == 2)
    {
        string whichAspectRatioType;
        double aspectRatioMin, aspectRatioMax, delta, eta;
        in  >> aspectRatioMin >> aspectRatioMax >> delta >> eta;
        for(size_t i = 0; i < elems.size(); ++i)
        {
			cout<<"not supported AQ"<<endl;
            //double aspectRatio(-1.0);
            //if(aspectRatioMin > 0.0 && aspectRatioMax > 0.0)
                //aspectRatio = weibull(aspectRatioMin, aspectRatioMax, delta, eta);
            //elems[i]->setRadiusFromAspectRatio(aspectRatio);
        }
    }
    else if(model == 3)
    {
        double multFactor;
        in >> multFactor;
		cout<< "multFactor "<< multFactor << endl;

        for(size_t i = 0; i < elems.size(); ++i)
        {
            double maxRadius(elems[i]->model()->radius()*multFactor);
            if(forPores)
            {
                for(int j = 0; j < elems[i]->connectionNum(); ++j)
                {
                    maxRadius = max(maxRadius, elems[i]->connection(j)->model()->radius());
                }
            }
            elems[i]->ChModel()->setRadius(maxRadius);
        }
    }
    else if(model == 4)
    {
        double aConst;
        bool applyAboveMedian, applyBelowMedian;
        char above, below;
        in >> aConst >> above >> below;
        applyAboveMedian = (above == 't' || above == 'T');
        applyBelowMedian = (below == 't' || below == 'T');

        vector< pair< double, double > > poreSizeDist;
        double totVolume(0.0);
        for(size_t i = 0; i < elems.size(); ++i)
        {
            totVolume += elems[i]->flowVolume();
            pair< double, double > entry(elems[i]->model()->radius(), totVolume);
            poreSizeDist.push_back(entry);
        }
        double volBasedAvr = tableLookUpX(totVolume/2.0, poreSizeDist);
        for(size_t j = 0; j < elems.size(); ++j)
        {
            double rad(elems[j]->model()->radius());
            if((rad > volBasedAvr && applyAboveMedian) || (rad < volBasedAvr && applyBelowMedian))
                rad = pow(elems[j]->model()->radius(), aConst)/pow(volBasedAvr, aConst-1.0);
            elems[j]->ChModel()->setRadius(rad);
        }
    }
    else if(model == 5)
    {
        double minVal, maxVal, deltaExp, etaExp;
        in >> minVal >> maxVal >> deltaExp >> etaExp;
        minVal *= 1.0E-6;
        maxVal *= 1.0E-6;
        for(size_t i = 0; i < elems.size(); ++i)
        {
            elems[i]->ChModel()->setRadius(weibull(minVal, maxVal, deltaExp, etaExp));
        }
    }
    else if(model == 6)
    {
        double a, b;
        in >> a >> b ;
        a *= 1.0E-6;
        b *= 1.0E-6;
        for(size_t i = 0; i < elems.size(); ++i)
        {
            elems[i]->ChModel()->setRadius
            (
				elems[i]->model()->radius()-a*tanh( elems[i]->model()->radius()/b )+0.0001*a
			);
        }
    }
    else
    {
        cerr << "==================================== " << endl
             << "Error: Uknown radius modifier option" << endl
             << "==================================== " << endl;        exit(-1);
    }

    if(writeToFile)
    {
		///. write old distribution
        ofstream fout(fileName.c_str());
        fout.flags(ios::showpoint);
        fout.flags(ios::scientific);
        fout.precision(6);

        double cumVol(0.0), totCumVol(0.0);
        sort(elems.begin(), elems.end(), ElemRadCmpRed());
        for(size_t i = 0; i < elems.size(); ++i)
            totCumVol += elems[i]->flowVolume();
		fout << "diameter       , cumulative volume"<< endl;
        double step(0.0), stepDelta(1.0/numPts);
        for(size_t j = 0; j < elems.size(); ++j)
        {
            cumVol += elems[j]->flowVolume();
            double poreVolFrac = cumVol/totCumVol;
            if(poreVolFrac > step || poreVolFrac == 1.0)
            {
                fout << setw(15) << elems[j]->model()->radius()*2.0E6 << "," << setw(15) << poreVolFrac << endl;
                step += stepDelta;
            }
        }

		///. write old distribution
		fout << endl;
		fout << "diameter       , pdf"<<endl;
		unsigned int elemFrequency(0);
        double  bin_radius(elems[0]->model()->radius()*(1.0-1.0/numPts));
        stepDelta = bin_radius/numPts;
        for(size_t j = 0; j < elems.size(); ++j)
        {
            //double poreVolFrac = cumVol/totCumVol;
            if(elems[j]->model()->radius() < bin_radius || j == elems.size()-1 )
            {
                fout << setw(15) << (bin_radius+stepDelta/2.0)*1.0E6 << "," << setw(15) << (double)elemFrequency/elems.size()/stepDelta/1.0E6 << endl;
                bin_radius -= stepDelta;
                elemFrequency = 0;
            }
            elemFrequency ++;
        }

    }
}

void setRadiiFromFile(vector<Element*>& elems, const string& radiiFileName,
                              double lowCutOffDiam, double highCutOffDiam, bool forPores)
{
    ifstream fin(radiiFileName.c_str());
    if (!fin)
    {
        cerr << "==========================================================="       << endl
             << "Error: Unable to open radius distribution file " << radiiFileName  << endl
             << "==========================================================="       << endl;
        exit( -1 );
    }

    vector< pair< double, double > > targetDist;        // diameter (dec) vs. pore vol fraction (inc)
    pair< double, double > funcPt, oldPt(1.0, -1.0);
    while(fin >> funcPt.first)
    {
        fin >> funcPt.second;
        funcPt.first *= 1.0E-6;                    // Convert from micro.m to m
        if(funcPt.second < oldPt.second || funcPt.first >= oldPt.first)
        {
            cerr << "==========================================================="   << endl
                << "Error: Fractional pore volume should increase monotonically."   << endl
                << "==========================================================="    << endl;
            exit( -1 );
        }
        else if(funcPt.second != oldPt.second)
        {
            targetDist.push_back(funcPt);
        }
        oldPt = funcPt;
    }

    if(lowCutOffDiam > targetDist.back().first && lowCutOffDiam < targetDist.front().first)
    {
        double volCutoff = tableLookUpY(lowCutOffDiam, targetDist);
        while(!targetDist.empty() && targetDist.back().first < lowCutOffDiam)
            targetDist.pop_back();

        for(size_t i = 0; i < targetDist.size(); ++i)
            targetDist[i].second /= volCutoff;
        pair< double, double > endPoint(lowCutOffDiam, 1.0);
        targetDist.push_back(endPoint);
    }

    if(highCutOffDiam > targetDist.back().first && highCutOffDiam < targetDist.front().first)
    {
        double volCutoff = tableLookUpY(highCutOffDiam, targetDist);
        vector< pair< double, double > > newTargetDist;
        pair< double, double > pt(highCutOffDiam, 0.0);
        newTargetDist.push_back(pt);

        for(size_t i = 0; i < targetDist.size(); ++i)
        {
            if(targetDist[i].first < highCutOffDiam)
            {
                pair< double, double > pt(targetDist[i].first,
                    (targetDist[i].second-volCutoff)/(1.0-volCutoff));
                newTargetDist.push_back(pt);
            }
        }
        targetDist = newTargetDist;
    }

    if(targetDist.empty() || targetDist.front().second != 0.0 || targetDist.back().second != 1.0)
    {
        cerr << endl
            << "============================================="        << endl
            << "Error: Target distribution contains errors."          << endl
            << "The cumulative volume should go from 0.0 to 1.0."     << endl
            << targetDist.size() << "  " << targetDist.back().second  << endl
            << "============================================="    << endl;    exit(-1);
    }

    double totElemVol(0.0), runningVol(0.0);
    for(size_t i = 0; i < elems.size(); ++i)
    {
        totElemVol += elems[i]->flowVolume();
    }

    for(size_t j = 0; j < elems.size(); ++j)
    {
        runningVol += elems[j]->flowVolume();
        double poreVolFrac = runningVol/totElemVol;
        double maxRadius(tableLookUpX(poreVolFrac, targetDist)/2.0);
        if(forPores)
        {
            for(int k = 0; k < elems[j]->connectionNum(); ++k)
            {
                maxRadius = max(maxRadius, elems[j]->connection(k)->model()->radius());
            }
        }
        elems[j]->ChModel()->setRadius(maxRadius);
    }
}




void setContactAngles(vector<Element*>& pores, vector<Element*>& throats, double min, double max,
                              double delta, double eta,
								int wettingClass, double modelTwoSepAng, string  angDistScheme, bool attractToMarked )
{

    vector< double > conAngles(pores.size());
    for(size_t i = 0; i < conAngles.size(); ++i)
    {
        conAngles[i] = weibull(min, max, delta, eta);
    }

    sort(conAngles.begin(), conAngles.end(), greater< double >());
    if(angDistScheme[2] == 'a' || angDistScheme[2] == 'A')      // rMax
        sort(pores.begin(), pores.end(), ElemRadCmpRed());
    else if(angDistScheme[2] == 'i' || angDistScheme[2] == 'I') // rMin
        sort(pores.begin(), pores.end(), ElemRadCmpInc());
    else
        sort(pores.begin(), pores.end(), ElemRandCmp());            // Random

    for(size_t j = 0; j < pores.size(); ++j)
    {
        //double rad = pores[j]->model()->radius();
        //double conAng = conAngles[j];
        pores[j]->ChModel()->setContactAngle(conAngles[j], wettingClass, modelTwoSepAng);
    }

    for(size_t k = 0; k < throats.size(); ++k)
    {
        const Element* pOne = throats[k]->connection(0);
        const Element* pTwo = throats[k]->connection(1);
        if((pOne->model()->clusterIndex() == 0 && pOne->model()->clusterIndex() == 0) ||    // None or both are marked ->
            (!pOne->model()->clusterIndex() == 0 && !pOne->model()->clusterIndex() == 0))   // assign on 50/50 basis
        {
            if(double(rand())/double(RAND_MAX) > 0.5)
                throats[k]->ChModel()->setContactAngle(pOne->model()->conAngEquil(), wettingClass, modelTwoSepAng);
            else
                throats[k]->ChModel()->setContactAngle(pTwo->model()->conAngEquil(), wettingClass, modelTwoSepAng);
        }
        else if((pOne->model()->clusterIndex() == 0 && attractToMarked) || (pTwo->model()->clusterIndex() == 0 && !attractToMarked))
        {
            throats[k]->ChModel()->setContactAngle(pTwo->model()->conAngEquil(), wettingClass, modelTwoSepAng);
        }
        else
        {
            throats[k]->ChModel()->setContactAngle(pOne->model()->conAngEquil(), wettingClass, modelTwoSepAng);
        }
    }
}





void loadRockTypeData(const string& fileName, vector<std::string>& poreTypes)
{
    ifstream input(fileName.c_str());
	if (!input)    cerr << "Error: Unable to open file " << fileName  << endl;
	else 		cout<< "Reading file " << fileName; cout.flush();

	//char tmp[10000];
	for(size_t i = 0; i < poreTypes.size(); ++i)
	{
		//int rockType;
		getline(input, poreTypes[i]);//
	}
	if (input.fail()) cout<<":\n      error in reading file"<<fileName;
    cout<<".\n"; cout.flush();

}



