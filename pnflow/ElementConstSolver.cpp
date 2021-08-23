#include <vector>
#include <stack>
#include <cassert>
#include "fluid.h"
#include "Element.h"
#include "inputData.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"




using namespace std;

/// averages pore and throat conductivities
void Throat::calcR2(const Fluid& fluid)
{
	auto throat = this;

	//Throat *throat = dynamic_cast< Throat* >(cnctions_[conn]);
	//ensure(throat);

	//double lengthPore2, throatLength, throatLength(throat->length());
	const Elem * p1=neib(0);
	const Elem * p2=neib(1);
	double 	poreLength0=poreLength(0);
	double 	poreLength2=poreLength(1);
	double length_t = length_;
	//Elem *p2;
	//if(throat->neib(0) == this)
	//{
		//p2 = throat->neib(1);
		//lengthPore2 = throat->poreLength(1);
		//throatLength = throat->poreLength(0);
	//}
	//else
	//{
		//p2 = throat->neib(0);
		//lengthPore2 = throat->poreLength(0);
		//throatLength = throat->poreLength(1);
	//}


	bool bulkConnection(model()->conductsAny(fluid) && 
		(p1->model()->conductsAny(fluid) || p1->isEntryOrExitRes()) && 
		(p2->model()->conductsAny(fluid) || p2->isEntryOrExitRes()));
	 

	bool paralell( !fluid.isOil() && 
		throat->model()->disConectedCentreWCornerW() && bulkConnection &&
		p1->model()->disConectedCentreWCornerW() && p2->model()->disConectedCentreWCornerW());
 
	if (bulkConnection && connectedToNetwork_)
	{
		double flowResistance(0.);
		//bool neighbourToExit();

		if (fluid.ff()>OIL)
		{
			double cond1, throatCond, cond2;
			//double resistivity = fluid.resistivity();
			cond1 = p1->model()->electricalConductance();
			throatCond = throat->model()->electricalConductance();
			cond2 = p2->model()->electricalConductance();

			//cout<<throatCond <<" "<< throatCond <<" "<< cond2 <<" "<<p2->isEntryOrExitRes()<< endl;
			ensure(cond1 > 0. &&  throatCond > 0. && cond2 > 0.);

			flowResistance = (length_t/throatCond + poreLength0/cond1 + poreLength2/cond2);
			double A_clay = (throat->clayVolume() + p1->clayVolume() + p2->clayVolume()) / (poreLength0+poreLength2+length_t);
			double condd = 1. / flowResistance + A_clay / (poreLength0+poreLength2+length_t) / (comn_.clay().resistivity());
			flowResistance = 1. / condd;
 
		}
		else if(paralell)
		{
			double cond1, throatCond, cond2;
			ensure(fluid.ff()==WTR);
			cond1 = p1->model()->getWaterConductance(filmBlob,p1->isEntryOrExitRes());
			throatCond = throat->model()->getWaterConductance(filmBlob);
			cond2 = p2->model()->getWaterConductance(filmBlob, p2->isEntryOrExitRes());
			ensure(cond1*throatCond*cond2 != 0.);

			if(cond1*throatCond*cond2 > 0.)
				flowResistance = (length_t/throatCond + poreLength0/cond1 + poreLength2/cond2);
			else
				flowResistance = 1e32;

			double cond1Bulk = p1->model()->getWaterConductance(bulkBlob, p1->isEntryOrExitRes());
			double throatCondBulk = throat->model()->getWaterConductance(bulkBlob);
			double cond2Bulk = p2->model()->getWaterConductance(bulkBlob, p2->isEntryOrExitRes());
			ensure(cond1Bulk*throatCondBulk*cond2Bulk != 0.);

			double flowResistBulk = 1e32;

			if(cond1Bulk*throatCondBulk*cond2Bulk > 0.)
				flowResistBulk = (length_t/throatCondBulk+poreLength0/cond1Bulk+poreLength2/cond2Bulk);

			flowResistance = 1./(1./flowResistance + 1./flowResistBulk);
		}
		else
		{
			double cond1, throatCond, cond2;
			cond1 = p1->model()->getConductance(fluid,p1->isEntryOrExitRes());
			throatCond = throat->model()->getConductance(fluid);
			cond2 = p2->model()->getConductance(fluid, p2->isEntryOrExitRes());



			if(cond1*throatCond*cond2 <=  1e-150 )
			{	if(debugLevel>10)
				{
					if  ( fluid.ff()==WTR)
					cout
					//<<" bulkSolverFlag  " << p1->canWBulkBePassedToSolver()<<" "<< (throat->canWBulkBePassedToSolver())<<" "<< (p2->canWBulkBePassedToSolver()) << "\n"
					<<" rockType        "<<p1->rockIndex()<<" "<<throat->rockIndex()<<" "<<p2->rockIndex()<<"\n"; 
					else
					cout
					//<<" filmSolverFlag  " << canOBePassedToSolver(filmBlob) <<" "<<(throat->canOBePassedToSolver(filmBlob)) <<" "<< (p2->canOBePassedToSolver(filmBlob)) << "\n"
					//<<" bulkSolverFlag  " << p1->canOBePassedToSolver(/*bulkBlob*/)<<" "<< (throat->canOBePassedToSolver(/*bulkBlob*/))<<" "<< (p2->canOBePassedToSolver(/*bulkBlob*/)) << "\n"
					<<" rockType        "<<p1->rockIndex()<<" "<<throat->rockIndex()<<" "<<p2->rockIndex()<<"\n"; 

					cout<< "cond: pore1 = " << cond1 <<" throat = "<< throatCond <<" next = "<< cond2 <<",  area = "<<throat->model()->area()<< "\n\n";
					//cout<<p2->isTrappedWat(filmBlob)<<p2->isTrappedWat(bulkBlob);
				}
				flowResistance = 5.e49;
			}
			else
				flowResistance = (length_t/throatCond + poreLength0/cond1 + poreLength2/cond2);
		}

		ensure(flowResistance > 0.);



		throat->setPoreToPoreCond(fluid.ff(), 1. / flowResistance);

	}
	else
	{
		throat->setPoreToPoreCond(fluid.ff(), 0.);
	}

}














/// Passes back the prevous pressure determined when solving the laplace eq.
bool Pore::prevSolvrRes(const Fluid& fluid, double loc, double& res, double& flowRate) const
{
	res = WOISolver[fluid.ff()-1]._condOrP2;

	if ( res<-1e+12 )    res = 0.; 

	return true;
}


void Pore::setSolverPrs(const Fluid& fluid, double res) const
{
	ensure(iAmAPore());
	if ( res<-1e-12 && (slvrCnct(fluid.ff()).isPassed() || isOnInletSlvrBdr_ || isOnOutletSlvrBdr_) && ( res>-1.4e24 || res<-1.6e24 ))  ++nErrs; 
	WOISolver[fluid.ff()-1]._condOrP2 = res;
}



/// Passes back the prevous pressure determined when solving the laplace eq.
/// The routine for throats interpolates between the connecting pore pressures and
/// also calculates the flowrate between the pores.  Needed only when useAvgPressure is activated
bool Throat::prevSolvrRes(const Fluid& fluid, double loc, double& res, double& flowRate) const
{
	double res0(0.), res1(0.), tmp(0.), gravCorr0(0.), gravCorr1(0.);

	cnctions_[0]->prevSolvrRes(fluid, loc, res0, tmp);
	cnctions_[1]->prevSolvrRes(fluid, loc, res1, tmp);

	if(USE_GRAV_IN_KR && fluid.ff()!=ELEC)
	{
		gravCorr0 = rhogh(fluid.density(), cnctions_[0]->node());
		gravCorr1 = rhogh(fluid.density(), cnctions_[1]->node());
	}

	if(cnctions_[0]->isEntryOrExitRes())         // Datum correct in/outlet reservoirs. These are always assumed
	{                                                // to be at the same lavel as connecting pore.
		res0 += gravCorr1-gravCorr0;
		gravCorr0 = gravCorr1;
	}

	if(cnctions_[1]->isEntryOrExitRes())
	{
		res1 += gravCorr0-gravCorr1;
		gravCorr1 = gravCorr0;
		if (fluid.ff()==ELEC) res1 = res0;
	}
	res0 -= gravCorr0;                          // Pressures passed back are datum levelled
	res1 -= gravCorr1;


	double loc0(cnctions_[0]->node().x);
	double loc1(cnctions_[1]->node().x);

	ensure(poreToPoreCond(fluid.ff())>=0.);

	if(poreToPoreCond(fluid.ff()) > 0.)
		flowRate = (2*(loc1<loc0)-1) * (res1-res0) * poreToPoreCond(fluid.ff());
	else flowRate = 0.;

	if (abs(loc1-loc0) <= 1e-32)    res = res0;
	else               res = res0 + (res1-res0) * (loc-loc0) / (loc1-loc0);

	ensure((loc >= loc0 && loc <=  loc1) || (loc >= loc1 && loc <=  loc0));
	ensure(res >= min(res0, res1) && res <=  max(res0, res1));

	return poreToPoreCond(fluid.ff()) > 0.;
}






	inline const Elem* Elem::nextSuccSolvrOil(bool& outletFound) const
	{
		for(auto nei:cnctions_)
		{
			if(!nei->slvrCnct(OIL).isSearched() && nei->isInsideSolverBox() &&
				((nei->model()->conductsAnyOil() && nei->model()->oilConductance() > COND_CUT_OFF) || nei->connectedToEntryOrExit()) 
				)
			{
				return nei;
			}
			else if(nei->isOnOutletSlvrBdr() &&  ( nei->model()->conductsAnyOil() || nei->isEntryOrExitRes() ) )
			{
				outletFound = true;
			}
		}

		return NULL;
	}



/// Before the solver can solve the pressure field we need to ensure that there is pressure
/// communication. Previous trapping routines are not sufficient for this as region connected
/// through oil layers might have become trapped without having been identified. This is done
/// through a depth first traversal from the elements at the inlet boundary. All connected
/// elements (through the defined fluid) are marked. These marked elements can then be passed to
/// the solver. The routine is somewhat more involved than the trapping routine since we do not
/// want connections through in/outlet to be marked. What we want is only the region of nodes
/// that are connected both to the in and outlet.
bool Elem::connectedToOutlet(const Fluid& fluid) const
{
	if(!model_->conductsAny(fluid) && !isEntryOrExitRes()) return false;


	if(fluid.ff()&WTR)
		for(int throat = 0; throat < nCncts_; ++throat)
		{            
			if(markWaterElemForSolver(cnctions_[throat], fluid.ff())) return true;
		}
	else
		for(int throat = 0; throat < nCncts_; ++throat)
		{  
			//if(markOilElemForSolver(cnctions_[throat])) return true;
			const Elem* elem = cnctions_[throat];
			
			if( (!elem->model()->conductsAnyOil() &&	!elem->connectedToEntryOrExit()) 
				||	elem->slvrCnct(OIL).isSearched() || !elem->isInsideSolverBox() )
			{
				continue;
			}

			vector<const Elem*> trappingStorage;
			stack<const Elem*> elemStack;
			elemStack.push(elem);
			bool outletFound(false);

			while(elem)
			{
				elem->slvrCnct(OIL).setAllFlags();
				trappingStorage.push_back(elem);
				elem = elem->nextSuccSolvrOil(outletFound);

				while(!elem)
				{
					elemStack.pop();
					if(elemStack.empty()) break;
					elem = elemStack.top();
					elem = elem->nextSuccSolvrOil(outletFound);
				}

				if(elem) elemStack.push(elem);
			}

			if(!outletFound)  // If region was isolated it should not be passed to solver
				for(size_t k = 0; k < trappingStorage.size(); ++k)   
					trappingStorage[k]->slvrCnct(OIL).clearFlag();
			else
				return true;
			
			
		}

	return false;
}



//&& (nei->iAmAPore()||nei->poreToPoreCond(WTR)>COND_CUT_OFF) 
		inline const Elem* Elem::nextSuccSolvrWat(bool& outletFound, fluidf ff) const // this is also used for ELECtricity phase
		{
			for(auto nei:cnctions_)
			{
				if(!nei->slvrCnct(ff).isSearched() && nei->isInsideSolverBox() &&
					( (nei->model()->conductAnyWater()) || nei->connectedToEntryOrExit()))
				{
					return nei;
				}
				else if(nei->isOnOutletSlvrBdr() &&
					(nei->model()->conductAnyWater() || nei->isEntryOrExitRes()))
				{
					outletFound = true;
				}
			}

			return NULL;
		}

	bool Elem::markWaterElemForSolver(const Elem* elem, fluidf ff) const
	{
		if( (!elem->model()->conductAnyWater() && !elem->connectedToEntryOrExit())
			 || elem->slvrCnct(ff).isSearched() ||	!elem->isInsideSolverBox())
			return false;

		vector< const Elem*> trappingStorage;
		stack< const Elem*> elemStack;
		elemStack.push(elem);
		bool outletFound(false);

		while(elem)
		{
			elem->slvrCnct(ff).setAllFlags();
			trappingStorage.push_back(elem);
			elem = elem->nextSuccSolvrWat(outletFound,ff);

			while(!elem)
			{
				elemStack.pop();
				if(elemStack.empty()) break;
				elem = elemStack.top();
				elem = elem->nextSuccSolvrWat(outletFound,ff);
			}

			if(elem) elemStack.push(elem);
		}

		if(!outletFound)
			for(size_t k = 0; k < trappingStorage.size(); ++k)         // If region was isolated it should not be passed to solver
				trappingStorage[k]->slvrCnct(ff).clearFlag();

		return outletFound;
	}





