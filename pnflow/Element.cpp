#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include <vector>
#include <stack>
#include <cassert>

#include "fluid.h"
#include "inputData.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"
#include "CommonData.h"


using namespace std;

bool Elem::USE_GRAV_IN_KR = false;
double Elem::COND_CUT_OFF = 1e-29;
const double Elem::PI = acos(-1.);

int Elem::nErrs = 0;

double Elem::RRR() const {return model_->RRR();}

Elem::Elem(const CommonData& comn, int indx, dbl3 nod, double radius, double volume,
				 double volClay, double shapeFact, int connNum, bool amPore, int typ)
	 :  Apex(), iAmAPore_(amPore), rockIndex_(typ),comn_(comn), flowVolume_(volume), clayVolume_(volClay), nCncts_(connNum),
	trapIndexOil_(-1, 0.),
	trapIndexWatBulk_(-1, 0.),
	trapIndexWatFilm_(-1, 0.),
	numOilCentreFeederNeis_(0),
	numWatCentreFeederNeis_(connNum), index_(indx), node_(nod)  {
	waterSaturation_ = 1.;
	isInsideSolverBox_ = false;
	isInCalcBox_ = false;
	isInWatFloodVec_ = false;
	isInOilFloodVec_ = false;
	isOnInletSlvrBdr_ = false;
	isOnOutletSlvrBdr_ = false;
	isExitRes_ = false;
	isEntryRes_ = false;
	isInWatFloodVec_ = false;
	isInOilFloodVec_ = false;
	connectedToNetwork_ = false;
	isConnectedToExit_ = false;
	isConnectedToEntry_ = false;
	connectedToEntryOrExit_ = false;



	eventI_ = -2;



	if(shapeFact <=  sqrt(3.)/36.+0.00001)                   // Triangular:  0 >= G >= sqrt(3.)/36.
	{
		model_ = new Triangle(*this, comn_, radius, min(shapeFact, sqrt(3.)/36.-0.00005), nCncts_, typ);
		 comn.countTriangle();

	}
	else if(shapeFact < 0.07)                         // Square:      G == 1./16/0
	{
		model_ = new Square(*this, comn_, radius, nCncts_, typ);
		comn.countSquare();
	}
	else                                                // Circular:    G == 1./4.*PI
	{
		model_ = new Circle(*this, comn_, radius, nCncts_, typ);
		comn.countCircle();
	}





	///.       set parent
	setConnections(model_,-1);


}


bool Elem::convertToMicroPorosityForSven(bool entry1Exit0)  {
		return false;
}

int Elem::ffaz() const {return int(model_->containCOil())+1;} ///. Viz only
//double Elem::saturation() const
//{ return
	//(flowVolume_*waterSaturation()+
	//cnctions_[0]->waterSaturation()*cnctions_[0]->flowVolume()/cnctions_[0]->nCncts()+
	//cnctions_[1]->waterSaturation()*cnctions_[1]->flowVolume()/cnctions_[1]->nCncts() ) /
	//(flowVolume_+
	 //cnctions_[0]->flowVolume()/cnctions_[0]->nCncts()+
	 //cnctions_[1]->flowVolume()/cnctions_[1]->nCncts()
	//);
//}

Elem::~Elem(){delete model_;}


void Elem::checkConnections() const
{
	if(nCncts_ != int(cnctions_.size()))  {
		cout<< "========================================" << endl
			<< "Error: The connection number is incorrect" << endl
			<< "=========================================" << endl;				exit(-1);//ERROR_STATE = true;
	}

	for(int ij = 0; ij < nCncts_; ++ij)  {
		if(cnctions_[ij] == NULL)  {
			cout<< "==================================="         << endl
				<< "Error: missing network connecetions"          << endl
				<< "Interanl pointer remains NULL      "          << endl
				<< "Elem index: " << indexOren()               << endl
				<< "Elem type: " << typeid(*this).name()       << endl
				<< "Connection: " << ij                         << endl
				<< "Total connection number: " << nCncts_ << endl
				<< "==================================="          << endl;	//exit(-1);//ERROR_STATE = true;
		}

		Elem* nextElem = cnctions_[ij];

		bool connExist(false);
		for(int ii = 0; ii < nextElem->nCncts(); ++ii)  {
			connExist = (nextElem->neib(ii) == this);
			if (connExist) break;
		}

		if(!connExist)  {
			cout<< "==================================="                               << endl
				<< "Error: mismatch in network connecetions"                                << endl
				<< "Interanal links do not match up    "                                << endl
				<< typeid(*this).name() <<" index: " << indexOren()                     << endl
				<< "Connection: " << ij                                               << endl
				<< "Neighbouring " <<typeid(*nextElem).name() << " index: " << nextElem->indexOren()    << endl
				<< "Neighbouring " <<typeid(*nextElem).name() << "  number of connections: " << nextElem->nCncts()    << endl
				<< "===================================" << endl;			//exit(-1);
		}
	}
}



/**
// This is the wrapper function for the trapping routine. In the trapping routine the
// visited elements are kept track of by inserting them in a vector. If an exit is found
// all the elements in that storage is unmarked and the storage is cleared. If the element
// indeed is trapped the region is given an identifier and the storage is copied such that
// during secondary drainage when this region again might be connected they can all be
// quickly unmarked.
*/
void Elem::findMarkTrappedOilGanglia(double prs, vector<Elem*>& trappingStorageOil,
								   double& elapsed, TrappingCriteria criteria)  {
	if(!isEntryOrExitRes() && model_->conductsAnyOil())  {
		clock_t startTrapRoutine(clock());
		if(criteria == escapeToBoth)  {
		   if ( foundEscapePathOil_trapOtherwise(prs, trappingStorageOil, escapeToOutlet)) ;
		   else foundEscapePathOil_trapOtherwise(prs, trappingStorageOil, escapeToInlet)   ;
		}
		else
			 foundEscapePathOil_trapOtherwise(prs, trappingStorageOil, criteria);

		elapsed += (double)(clock() - startTrapRoutine) / CLOCKS_PER_SEC;
	}
}





/**
// Marking during the trapping routine is done by increasing the trapping index.
// If the outlet is subsequently found, it is reset down to 0.
*/
inline void Elem::trapOil(double prs)  {
	trapIndexOil_.first = comn_.newOilTrappingIndex();
	trapIndexOil_.second = prs;
}




	/**
	* During depth first graph traversal we need to retrive the next sucsessor
	* after an elemnt. This will be an element that is not marked and contains
	* the fluid in question (oil in this case)
	*/
	inline Elem* Elem::nextUntrappedOil(TrappingCriteria criteria)  {
		if(criteria == escapeToOutlet)  ///.  escapeToInlet does not make any difference
		{
			for(short i = 0; i < nCncts_; ++i)  {
				if(!cnctions_[i]->isEntryRes() &&
					!cnctions_[i]->isTrappedOil() &&
					cnctions_[i]->model()->conductsAnyOil())
					return cnctions_[i];
			}
		}
		else if(criteria == escapeToEither)  {
			for(short i = 0; i < nCncts_; ++i)  {
				if(// !cnctions_[i]->isEntryOrExitRes() &&
					!cnctions_[i]->isTrappedOil() &&
					cnctions_[i]->model()->conductsAnyOil())
					return cnctions_[i];
			}
		}
		else if(criteria == escapeToInlet)  {
			for(short i = nCncts_-1; i >= 0; --i)  {
				if(!cnctions_[i]->isExitRes() &&
					!cnctions_[i]->isTrappedOil() &&
					cnctions_[i]->model()->conductsAnyOil())///. assumes elem oil is connected to all connections
					return cnctions_[i];
			}
		}
		else
		{
			for(short i = 0; i < nCncts_; ++i)  {
				if(//!cnctions_[i]->isEntryOrExitRes() &&
					!cnctions_[i]->isTrappedOil() &&
					cnctions_[i]->model()->conductsAnyOil())
					return cnctions_[i];
			}
		}
		return NULL;
}

/**
 * finds the list of elements forming a disconnected ganglia and calls their trapOil(localPc), or returns true otherwise.
 * Per:
// This is the trapping routine. It is a depth first traversal of the network. To prevent stack
// overflow this routine is implemented without recursion. It could have been made general with
// respect to fluid, however that would have reduced efficiency due to *MANY* checks for fluid.
// This concern was strong enough to warrant separate implementations for each fluid
*/
bool Elem::foundEscapePathOil_trapOtherwise(double pc, vector<Elem*>& trappingStorage, TrappingCriteria criteria)  {

	Elem* elemPtr = this;
	stack<Elem*> elemStack;                      // Simulate recursive behavior using a stack
	//elemStack.push(this);
	double datumPc(pc+gravityCorrection());

	while(elemPtr)  {
		if(  (elemPtr->isConnectedToEntry() && (criteria == escapeToInlet || criteria == escapeToEither))
		  || (elemPtr->isConnectedToExit() && (criteria == escapeToOutlet || criteria == escapeToEither))  )  {
			for(size_t i = 0; i < trappingStorage.size(); ++i)  {
				trappingStorage[i]->unTrapOil();
			}
			trappingStorage.clear();
			return false;
		}

		double localPc(datumPc - elemPtr->gravityCorrection());
		ensure(!elemPtr->isEntryOrExitRes());
		elemPtr->trapOil(localPc);
		trappingStorage.push_back(elemPtr);
		elemStack.push(elemPtr);            // Simulate recursive descent


		do{
			elemPtr = elemStack.top()->nextUntrappedOil(criteria); ///.  criteria is just for speed
			if(!elemPtr)
				elemStack.pop();  // Simulate recursive return. Will unwind the stack
		}while(!elemPtr && !elemStack.empty()); // until a new possible branch is found

	}
	return true;                                       // Stack is empty. The region is trapped.
}




void Elem::findMarkTrappedWaterGanglia(double prs, FluidBlob startPt, vector< pair<Elem*,FluidBlob> >& trappingStorageWat,
								   double& elapsed, TrappingCriteria criteria)  {
	if(!isEntryOrExitRes() && model_->conductAnyWaterBlob(startPt))  {
		clock_t startTrapRoutine(clock());
		if(criteria == escapeToBoth)  {
			if  (foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, escapeToOutlet)) ;
			else foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, escapeToInlet );
		}
		else     foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, criteria);

		elapsed += (double)(clock() - startTrapRoutine) / CLOCKS_PER_SEC;
	}
}




	inline void Elem::trapWat(double prs, FluidBlob blob)  {
		if(blob == bulkBlob)  {
			trapIndexWatBulk_.first = comn_.newWatTrappingIndex();
			trapIndexWatBulk_.second = prs;
		}
		else
		{
			trapIndexWatFilm_.first = comn_.newWatTrappingIndex();
			trapIndexWatFilm_.second = prs;
		}
	}





	inline Elem* Elem::nextSuccessorWat(TrappingCriteria criteria, FluidBlob& blob)  {
		if(criteria == escapeToInlet)  {
			for(short i = nCncts_; i > 0; --i)  {
				if(!cnctions_[i-1]->isEntryOrExitRes() &&
					!cnctions_[i-1]->isTrappedWat(blob) &&
					cnctions_[i-1]->model()->conductAnyWaterBlob(blob))
					//(cnctions_[i-1]->model()->conductsAnyWaterBlob(blob) || cnctions_[i-1]->connectedToEntryOrExit()))
					return cnctions_[i-1];
			}
		}
		else
		{
			for(short i = 0; i < nCncts_; ++i)  {
				if(!cnctions_[i]->isEntryOrExitRes() &&
					!cnctions_[i]->isTrappedWat(blob) &&
					cnctions_[i]->model()->conductAnyWaterBlob(blob))
					//(cnctions_[i]->model()->conductsAnyWaterBlob(blob) || cnctions_[i]->connectedToEntryOrExit()))
					return cnctions_[i];
			}
		}

		FluidBlob blobbo(blob == filmBlob ? bulkBlob: filmBlob);
		if(!model_->disConectedCentreWCornerW() && model_->conductAnyWaterBlob(blobbo) && !isTrappedWat(blobbo))  {
			blob = blobbo;
			return this;
		}

		return NULL;
	}

inline bool Elem::foundEscapePathWat_trapOtherwise(double pc, FluidBlob startPt, vector< pair<Elem*,FluidBlob> >& trappingStorage,
										TrappingCriteria criteria)  {
	pair<Elem*, FluidBlob> elem(this, startPt);
	stack< pair<Elem*, FluidBlob> > elemStack;
	elemStack.push(elem);
	double datumPc(pc+gravityCorrection());

	while(elem.first)  {
		if((elem.first->isConnectedToEntry() && (criteria == escapeToInlet || criteria == escapeToEither)) ||
			(elem.first->isConnectedToExit() && (criteria == escapeToOutlet || criteria == escapeToEither)))  {
			for(size_t i = 0; i < trappingStorage.size(); ++i)
				trappingStorage[i].first->unTrapWat(trappingStorage[i].second);
			trappingStorage.clear();
			return false;
		}

		ensure(!elem.first->isEntryOrExitRes());
		double localPc(datumPc-elem.first->gravityCorrection());
		elem.first->trapWat(localPc, elem.second);
		trappingStorage.push_back(elem);
		ensure(elem.first);
		elem.first = elem.first->nextSuccessorWat(criteria, elem.second);

		while(!elem.first)  {
			elemStack.pop();
			if(elemStack.empty()) break;
			elem = elemStack.top();
			ensure(elem.first);
			elem.first = elem.first->nextSuccessorWat(criteria, elem.second);
		}

		if(elem.first) elemStack.push(elem);
	}
	return true;
}


/// It is possible that some parts of the network in fact is physically isolated from the
/// network due to diagenesis processes. These parts will remain water filled but are not
/// available for displacement. Identify these elements during initialization by doing a
/// graph traversal from the inlet.
void Elem::identifyConnectedPoreElems()  {
	setConnectedToNetwork(true);      // Do allow escapes through entry reservoir

	for(int throat = 0; throat < nCncts_; ++throat)  {
		if(!cnctions_[throat]->connectedToNetwork())                // Might have been marked previously
		{
			set< Elem * > frontier, oldFrontier;
			frontier.insert(cnctions_[throat]);

			while(!frontier.empty())  {
				oldFrontier = frontier;
				frontier.clear();

				for(auto elm:oldFrontier)  {
					if(!elm->connectedToNetwork())                   // Flag is not set
					{
						elm->setConnectedToNetwork(true);               // Set flag
						for(int i = 0; i < elm->nCncts(); ++i)
							frontier.insert(elm->neib(i));
					}
				}
			}
		}
	}
}



/// We have the option to not drain pores that are connected to the network through a
/// single throat. These pores (and throats) are removed completley from the network
/// initially. This is found to greatly improve solver convergence
int Elem::removeFromNetwork()  {
	int numRemoved(1);
	for(int i = 0; i < nCncts_; ++i)  {
		if(cnctions_[i]->connectedToNetwork())  {
			if(cnctions_[i]->nCncts() > 2)  cnctions_[i]->severConnection(this);
			else                                      numRemoved += cnctions_[i]->removeFromNetwork();
		}
	}
	cnctions_.clear(); nCncts_=0; connectedToNetwork_=false;
	return numRemoved;
}

/// Physically remove the single connetion to the network. This cluster is now isolated
void Elem::severConnection(Elem* elm)  {
	ensure(nCncts_ > 0);
	vector<Elem*>::iterator delm = find(cnctions_.begin(), cnctions_.end(), elm);
	if(*delm == elm)  cnctions_.erase(delm);
	else              cerr << endl << "oooops..." << endl;
	--nCncts_;
}













int Throat::indexOren() const { return index_-comn_.numPores()-1; }


double Elem::rhogh(double rho, dbl3 dh) const {return comn_.rhogh(rho,dh);}

void Elem::setGravityCorrection(const dbl3& nod)  {
	gravityCorrection_ = (comn_.water().density()-comn_.oil().density())*(comn_.gravConst()&nod);

		gravityCorrection_ = 0.;///. ERROR

}

void Elem::setGravityCorrection(const dbl3& nod1, const dbl3& nod2)  {
	gravityCorrection_ = (comn_.water().density()-comn_.oil().density())*(comn_.gravConst()&((nod1+nod2)*0.5));

		gravityCorrection_ = 0.;///. ERROR
}



/**
// The corner conductancees are computed using the suggested procedure by Patzek SPE59312
// Elements that are trapped have their capillary pressure frozen at the time they're
// identified as being trapped.
*/
double Elem::updateSat_calcR(double cappPrs)  {
	if(isEntryOrExitRes()) return 0.;

	waterSaturation_ = model_->calcR(cappPrs);

	return isInCalcBox_ ? waterSaturation_ * flowVolume_ + clayVolume_ : 0.;
}


/**
// As an elemnt is drained, its neighbours have an increase in number of elems
// containing oil.
*/
void Elem::fillElemCentreWithOilRemoveLayers()  {
	//ensure((comn_.water()));
	//double warninnng;
	if(model_->bulkFluid()->isOil()) return;


	if (!isInOilFloodVec_ && debugLevel>0) cout<<"Qo";
	isInOilFloodVec_ = false;
	model_->fillCentreWithOilRemoveLayers();

	for(int i = 0; i < nCncts_; ++i)  {
		cnctions_[i]->IncreaseNumOilCentreFeederNeis();
		if (!rockIndex())	cnctions_[i]->ReduceNumWatCentreFeederNeis();
	}

	///. rubbish
	if(model_->displacementType() == 'S')  {
		//model_->addFillingEventToHistory(2);
		eventI_ = -1;
	}
	else
	{
		//model_->addFillingEventToHistory(1);
		eventI_ = nCncts_ - numOilCentreFeederNeis();
	}
}

/**
// As an elemt is imbided we need to update the imbibition entry pressures of
// neighbouring elements.
*/
void Elem::fillElemCentreWithWaterCreateLayers(bool snapOffOverRide)  {
	if (debugLevel>0)  {
		if (!snapOffOverRide && !isInWatFloodVec_) cout<<"Qw";
		else if (snapOffOverRide && isInWatFloodVec_) cout<<"Qi";
	}

	isInWatFloodVec_ = false;

	ensure(model_->bulkFluid()->isOil(), "RefillingWithWater ");

	model_->fillCentreWithWaterCreateLayers(snapOffOverRide);

	for(int i = 0; i < nCncts_; ++i)  {
		cnctions_[i]->ReduceNumOilCentreFeederNeis(); ///.
		if (!rockIndex())	cnctions_[i]->IncreaseNumWatCentreFeederNeis();
	}

	if(model_->displacementType() == 'S')  {
		//model_->addFillingEventToHistory(4);
		eventI_ = -1;
	}
	else
	{
		//model_->addFillingEventToHistory(3);
		eventI_ = numOilCentreFeederNeis();
	}
}



bool Elem::canBeAddedToEventVec(const Fluid& injectant) const
{
	if(model_->canNOTReconfigure(injectant) || !connectedToNetwork_)
		return false;

   if(injectant.isOil())  {
		if (model_->containCOil() && rockIndex() == 0)
		{
			if (debugLevel>0) cout<<" sgsobqbd ";
			return false;//         || nonTrappedOilCentreNeighbour(injectant));
		}

		if( isInOilFloodVec_ || (trapIndexWatBulk_.first > -1) ) return false;

		if (model_->conductsAnyOil()) return true;//         || nonTrappedOilCentreNeighbour(injectant));


		for(int i = 0; i < nCncts_; ++i)  {
			if(cnctions_[i]->model()->conductCOil() && !cnctions_[i]->isTrappedOil()) /// second is extra to delete
				return true;
		}

		if(numOilCentreFeederNeis_ == 0) return false;

		return false;
	}
	else
	{
		if (!model_->conductsAnyOil() && rockIndex() == 0)
		{
			if (debugLevel>0) cout<<" sjswbqld ";
			return false;//         || nonTrappedOilCentreNeighbour(injectant));
		}

		if(isInWatFloodVec_ || (trapIndexOil_.first > -1) ) return false;

		if (model_->conductAnyWater() ) return true;//         || notTrapdWNeighbour(injectant));


		for(int i = 0; i < nCncts_; ++i)
			if(cnctions_[i]->model()->conductCWater())
				return true;


		//if(numWatCentreFeederNeis_ == 0) return false;

		return false;
	}

}





///.  rare funcs
bool Elem::addToLayerVec(const Fluid& injectant, Polygon* shyp, std::vector< int >& cornToAdd) const
{
	if( !shyp->hasOilLayer() || shyp->containCOil() || trapIndexOil_.first != -1 ) return false;


	bool oilInj = injectant.isOil();

	if(oilInj)     // Layer refcreation
	{

	  if(isInOilFloodVec())  {
		for(int i = 0; i < shyp->numCorners(); ++i)
			if(!shyp->oilLayerConst()[i].isInOilFloodVec() &&
				 ! shyp->oilLayerConst()[i].exists()
				 && entryPc() > shyp->oilLayerConst()[i].entryPc()
			  )
				cornToAdd.push_back(i);
	  }
	  else
	  {
		for(int i = 0; i < shyp->numCorners(); ++i)
			if(!shyp->oilLayerConst()[i].isInOilFloodVec() &&
				 ! shyp->oilLayerConst()[i].exists()
			  )
				cornToAdd.push_back(i);
		}

		if(!cornToAdd.empty())
			for(size_t j = 0; j < cnctions_.size(); ++j)
				if(cnctions_[j]->model()->conductsAnyOil())
					return true;

		return false;
	}
	else        // Layer collapse
	{
		for(int i = 0; i < shyp->numCorners(); ++i)
			if(!shyp->oilLayerConst()[i].isInWatFloodVec() &&
				shyp->oilLayerConst()[i].exists() &&
				trapIndexOil_.first < 0 &&
				(trapIndexWatFilm_.first < 0 || trapIndexWatBulk_.first < 0)
			 )
				cornToAdd.push_back(i);

		if(!cornToAdd.empty())
		   return true;
		else
			return false;
	}
}




void Elem::calcCentreEntryPrsOilInj()  {
	entryPc_=model_->centreEntryPrsOilInj();
}

void Elem::calcCentreEntryPrsWatInj()  {
	entryPc_=model_->centreEntryPrsWatInj();
}




//#include "ElementConstSolver.cpp"
//#include "ElementPore.cpp"
//#include "ElementThroat.cpp"
