#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include <cassert>

using namespace std;


#include "fluid.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"
#include "inputData.h"

#include "Element.h"



Throat::Throat(const CommonData& common, int indx, dbl3 nod, double radius, double vol, double volumeClay,
			   double shapeFactor, double length, double lengthPore1, double lengthPore2, int rockType)
	: Elem(common, indx, nod, radius, vol, volumeClay, shapeFactor, 2, false, rockType), index_(indx),
	poreLength_{lengthPore1,lengthPore2}, length_(length)  {}










void Throat::prepare2()  {
	if(index()<10) outD<<"g:"<<neib(0)->index()<<"-"<<neib(1)->index()<<":"<<index()<<" R:"<<RRR()<<" nCor_"<<index()<<":"<<model()->numCorners()<<"\n";
	checkConnections();
}



double Throat::snapOfLongitCurvature() const
{
	double delRSqr = 0.5*(cnctions_[1]->model()->RRR()+cnctions_[0]->model()->RRR())-model()->RRR();
	if (delRSqr<0.) return 0.; ///. Errror
	delRSqr *= delRSqr;

	return 0.;
	//return 4.*std::sqrt(delRSqr/(delRSqr+lengthEff*lengthEff*0.025330296))/lengthEff; //lengthEff/lengthEff;
	///. return 2*sigma*sin(45)*sqrt(2.*delRSqr/(delRSqr+(L_poreToPore/2PI)^2)/(lengthEff/2.) //lengthEff/lengthEff;

}



const Pore* Throat::neighbouringPore(const Elem* callingPore) const
{
	ensure(cnctions_.size() == 2, " Throat connection deleted, set drainSinglets to true.", -1);
	if(cnctions_[0] == callingPore)  return static_cast<Pore*>(cnctions_[1]);
	else                                return static_cast<Pore*>(cnctions_[0]);
}





/// Does a throat cross any given interior plane
bool Throat::crossesPlaneAt(double location) const
{
	double ptOne(cnctions_[0]->node().x), ptTwo(cnctions_[1]->node().x);

	if(ptTwo < ptOne)  {   // Ensure that that the two points we pick up acyually are within the inteior of the model
		double tmp = ptTwo;   ptTwo = ptOne;   ptOne = tmp;
	}

	return ptOne < location && ptTwo >= location;
}


/**
// When sorting the connecting pores in the throats we have to make sure that we also sort
// the associated pore lengths ;-).
*/
void Throat::sortConnectingElems_DistToExit()  {
	Elem *oldFirst = cnctions_[0];
	sort(cnctions_.begin(), cnctions_.end(), DistToExitComparePores());

	assert(cnctions_.size() == 2);
	if(oldFirst != cnctions_[0])    // Pores got switched
	{
		double tmp(poreLength_[0]);   poreLength_[0] = poreLength_[1];  poreLength_[1] = tmp;
	}
}




/**

// The connecting pores are added to the throat. From the pores it is also
// possible to determine if the throat is inside the calculation box. If part of
// the connection (pore-throat-pore) is outside the calculation box that lenght
// is set to zero. This has the effect when solving the pressure field, no
// pressure loss occurs outside the box. However if the pore inside the box itself
// is on the boundary, pressure losses do occur within that pore. This constraint
// is needed in the case we're using the whole network for rel perm calculations.
// The first pores are usually located at x position 0.. These are still within
// the box. There has to be some length to the outlet to make the problem
// solveable, hence we allow pressure drop to occur within that pore.
*/
void Throat::addConnections(Elem* p1, Elem* p2, double inletBdr, double outletBdr, bool moveBoundary, bool setNod)  {
	if(setNod)  {
		node_=0.5*(p2->node() + p1->node());
		if(p1->isEntryOrExitRes())  {
			node_.y= p2->node().y;  node_.z= p2->node().z;
		}

		if(p2->isEntryOrExitRes())  {
			node_.y= p1->node().y;  node_.z= p1->node().z;
		}
	}

	double mid=(inletBdr+outletBdr)*0.5;

	if(p1->isEntryRes() || p1->isExitRes())       setGravityCorrection(p2->node());
	else if(p2->isEntryRes() || p2->isExitRes())  setGravityCorrection(p1->node());
	else                                          setGravityCorrection(p1->node(), p2->node());

	if(p1->isEntryRes() || p2->isEntryRes())  {
		isConnectedToEntry_ = true;
		p1->isConnectedToEntry(p1->isEntryRes());
		p2->isConnectedToEntry(p2->isEntryRes());
	}
	else if(p1->isExitRes() || p2->isExitRes())  {
		isConnectedToExit_ = true;
		p1->isConnectedToExit(p1->isExitRes());
		p2->isConnectedToExit(p2->isExitRes());
	}

	isInsideSolverBox_ = (p1->isInsideSolverBox() || p2->isInsideSolverBox());
	isInCalcBox_ = (p1->isInCalcBox() || p2->isInCalcBox());
	connectedToEntryOrExit_ = (p1->isEntryOrExitRes() || p2->isEntryOrExitRes());

	double oldPOneLen(poreLength_[0]), oldPTwoLen(poreLength_[1]), oldThrLen(length_);

	if(isInsideSolverBox_ && !p2->isInsideSolverBox())  {
		p2->setOnInletSlvrBdr(p2->node().x < mid);
		p2->setOnOutletSlvrBdr(p2->node().x > mid);

		//if(moveBoundary) // this is ttoo risky, not compatible with small throat lenghts and twsted throats
		//{
			//double scaleFact = (p2->node().x - p1->node().x) / (length_+poreLength_[0]+poreLength_[1]);
			//if(p2->isEntryOrExitRes()) scaleFact = (p2->node().x-p1->node().x) / fabs(p2->node().x-p1->node().x);     // We don't know position of exit and entry res.

			//double bdr = (p2->node().x < inletBdr) ? inletBdr: outletBdr;
			//double throatStart = p1->node().x + poreLength_[0]*scaleFact;
			//double throatEnd = throatStart + length_*scaleFact;


			if(p2->isEntryOrExitRes())                                          // Keep throat lengths if whole model is being used
				poreLength_[1] = 1e-300;                                              // for calculations
			//else if(throatEnd > inletBdr && throatEnd < outletBdr)                  // Both pore1 and throat are within the box
				//poreLength_[1] *= (bdr - throatEnd)/(poreLength_[1]*scaleFact);
			//else if(throatStart > inletBdr && throatStart < outletBdr)              // Onle pore 1 is fully within box
			//{
				//poreLength_[1] = 1e-300;
				//length_ *= (bdr - throatStart)/(length_*scaleFact);
			//}
			//else                                                                    // Pore 1 is only partially within box
			//{
				//poreLength_[1] = 1e-300;
				//length_ = 1e-300;
			//}
		//}
	}
	else if(isInsideSolverBox_ && !p1->isInsideSolverBox())             // Pore 1 is outside box
	{
		ensure(poreLength_[0]>0.);	ensure(poreLength_[1]>0.);
		p1->setOnInletSlvrBdr(p1->node().x < mid);
		p1->setOnOutletSlvrBdr(p1->node().x > mid);

		//if(moveBoundary)
		//{
			//double scaleFact = (p1->node().x - p2->node().x) / (length_+poreLength_[0]+poreLength_[1]);
			//if(p1->isEntryOrExitRes())       // We don't know position of exit and entry res.
				//scaleFact = (p1->node().x-p2->node().x) / fabs(p1->node().x-p2->node().x);

			//double bdr = (p1->node().x < inletBdr) ? inletBdr: outletBdr;
			//double throatStart = p2->node().x + poreLength_[1]*scaleFact;
			//double throatEnd = throatStart + length_*scaleFact;

			if(p1->isEntryOrExitRes())
				poreLength_[0] = 1e-300;
			//else if(throatEnd > inletBdr && throatEnd < outletBdr)               // Both pore 2 and throat are within the box
				//poreLength_[0] *= (bdr - throatEnd)/(poreLength_[0]*scaleFact);
			//else if(throatStart > inletBdr && throatStart < outletBdr)      // Only pore 2 is fully within box
			//{
				//poreLength_[0] = 1e-300;
				//length_ *= (bdr - throatStart)/(length_*scaleFact);
			//}
			//else                                                            // Pore 2 is only partially within box
			//{
				//poreLength_[0] = 1e-300;
				//length_ = 1e-300;
			//}
			//ensure(poreLength_[0]>0.);	ensure(poreLength_[1]>0.);
		//}
	}

	if(poreLength_[0] > 1.1*oldPOneLen || poreLength_[1] > 1.1*oldPTwoLen || length_ > 1.1*oldThrLen)  {
		cout<< endl
			<< "==============================================="    << endl
			<< "Warning: The new lengths for elements connected"    << endl
			<< "to the pressure boundary are larger than the"       << endl
			<< "original ones. The lengths should be smaller"       << endl
			<< "since we do not want pressure drops occuring"       << endl
			<< "outside the box across which we're calculating"     << endl
			<< "relative permeability."                             << endl
			<< "============================================== "     << endl
			<< endl;
	}
	ensure(poreLength_[0]>0.);	ensure(poreLength_[1]>0.);

	cnctions_.push_back(p1);         // Add pore connections
	cnctions_.push_back(p2);

	neiPores_[0]=static_cast<Pore*>(p1);
	neiPores_[1]=static_cast<Pore*>(p2);
}








void Throat::modifyLength(double scaleFactor)  {
	poreLength_[0] *= scaleFactor;
	poreLength_[1] *= scaleFactor;
	length_ *= scaleFactor;
}
