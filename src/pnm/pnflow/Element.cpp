#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include <vector>
#include <stack>
#include <cassert>

#include "fluid.h"
#include "Element.h"
#include "inputData.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"


using namespace std;

bool Element::USE_GRAV_IN_KR = false;
double Element::COND_CUT_OFF = 1.0e-29;
const double Element::PI = acos(-1.0);

int Element::nErrs = 0;

double Element::RRR() const {return model_->RRR();}

Element::Element(const CommonData& comn, int indx, dbl3 nod, double radius, double volume,
				 double volClay, double shapeFact, int connNum, bool amPore, int typ) 
	 :  Apex(), iAmAPore_(amPore), rockIndex_(typ),comn_(comn), flowVolume_(volume), clayVolume_(volClay), connectionNum_(connNum), 
	trapIndexOil_(-1, 0.0),
	trapIndexWatBulk_(-1, 0.0),
	trapIndexWatFilm_(-1, 0.0),
	numOilCentreFeederNeis_(0),
	numWatCentreFeederNeis_(connNum), index_(indx), node_(nod)
{
	waterSaturation_ = 1.0;
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



	if(shapeFact <=  sqrt(3.0)/36.0+0.00001)                   // Triangular:  0 >= G >= sqrt(3.0)/36.0
	{
		model_ = new Triangle(*this, comn_, radius, min(shapeFact, sqrt(3.0)/36.0-0.00005), connectionNum_, typ);
		 comn.countTriangle();

	}
	else if(shapeFact < 0.07)                         // Square:      G == 1.0/16/0
	{
		model_ = new Square(*this, comn_, radius, connectionNum_, typ);
		comn.countSquare();
	}
	else                                                // Circular:    G == 1.0/4.0*PI
	{
		model_ = new Circle(*this, comn_, radius, connectionNum_, typ);
		comn.countCircle();
	}





	///.       set parent
	setConnections(model_,-1);


}


bool Element::convertToMicroPorosityForSven(bool entry1Exit0)
{
		return false;
}

int Element::ffaz() const {return int(model_->containCOil())+1;} ///. Viz only
//double Element::saturation() const 
//{ return
	//(flowVolume_*waterSaturation()+
	//connections_[0]->waterSaturation()*connections_[0]->flowVolume()/connections_[0]->connectionNum()+ 
	//connections_[1]->waterSaturation()*connections_[1]->flowVolume()/connections_[1]->connectionNum() ) /
	//(flowVolume_+
	 //connections_[0]->flowVolume()/connections_[0]->connectionNum()+
	 //connections_[1]->flowVolume()/connections_[1]->connectionNum()
	//);
//}

Element::~Element(){delete model_;}


void Element::checkConnections() const
{
	if(connectionNum_ != int(connections_.size()))
	{
		cout<< "========================================" << endl
			<< "Error: The connection number is incorrect" << endl
			<< "=========================================" << endl;				exit(-1);//ERROR_STATE = true;
	}

	for(int ij = 0; ij < connectionNum_; ++ij)
	{
		if(connections_[ij] == NULL)
		{ 
			cout<< "==================================="         << endl
				<< "Error: missing network connecetions"          << endl
				<< "Interanl pointer remains NULL      "          << endl
				<< "Element index: " << indexOren()               << endl
				<< "Element type: " << typeid(*this).name()       << endl
				<< "Connection: " << ij                         << endl
				<< "Total connection number: " << connectionNum_ << endl
				<< "==================================="          << endl;	//exit(-1);//ERROR_STATE = true;
		}

		Element* nextElem = connections_[ij];

		bool connExist(false);
		for(int ii = 0; ii < nextElem->connectionNum(); ++ii)
		{
			connExist = (nextElem->connection(ii) == this);
			if (connExist) break;
		}

		if(!connExist)
		{
			cout<< "==================================="                               << endl
				<< "Error: mismatch in network connecetions"                                << endl
				<< "Interanal links do not match up    "                                << endl
				<< typeid(*this).name() <<" index: " << indexOren()                     << endl
				<< "Connection: " << ij                                               << endl
				<< "Neighbouring " <<typeid(*nextElem).name() << " index: " << nextElem->indexOren()    << endl
				<< "Neighbouring " <<typeid(*nextElem).name() << "  number of connections: " << nextElem->connectionNum()    << endl
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
void Element::findMarkTrappedOilGanglia(double prs, vector<Element*>& trappingStorageOil,
								   double& elapsed, TrappingCriteria criteria)
{
	if(!isEntryOrExitRes() && model_->conductsAnyOil())
	{
		clock_t startTrapRoutine(clock());
		if(criteria == escapeToBoth)
		{
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
inline void Element::trapOil(double prs)
{
	trapIndexOil_.first = comn_.newOilTrappingIndex();
	trapIndexOil_.second = prs;
}




	/**
	* During depth first graph traversal we need to retrive the next sucsessor
	* after an elemnt. This will be an element that is not marked and contains
	* the fluid in question (oil in this case)
	*/
	inline Element* Element::nextUntrappedOil(TrappingCriteria criteria)
	{
		if(criteria == escapeToOutlet)  ///.  escapeToInlet does not make any difference
		{
			for(short i = 0; i < connectionNum_; ++i)
			{
				if(!connections_[i]->isEntryRes() &&
					!connections_[i]->isTrappedOil() &&
					connections_[i]->model()->conductsAnyOil())
					return connections_[i];
			}
		}    
		else if(criteria == escapeToEither)
		{
			for(short i = 0; i < connectionNum_; ++i)
			{
				if(// !connections_[i]->isEntryOrExitRes() &&
					!connections_[i]->isTrappedOil() &&
					connections_[i]->model()->conductsAnyOil())
					return connections_[i];
			}
		}
		else if(criteria == escapeToInlet)
		{
			for(short i = connectionNum_-1; i >= 0; --i)
			{
				if(!connections_[i]->isExitRes() &&
					!connections_[i]->isTrappedOil() &&
					connections_[i]->model()->conductsAnyOil())///. assumes elem oil is connected to all connections
					return connections_[i];
			}
		}
		else
		{
			for(short i = 0; i < connectionNum_; ++i)
			{
				if(//!connections_[i]->isEntryOrExitRes() &&
					!connections_[i]->isTrappedOil() &&
					connections_[i]->model()->conductsAnyOil())
					return connections_[i];
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
bool Element::foundEscapePathOil_trapOtherwise(double pc, vector<Element*>& trappingStorage, TrappingCriteria criteria)
{

	Element* elemPtr = this;
	stack<Element*> elemStack;                      // Simulate recursive behavior using a stack
	//elemStack.push(this);
	double datumPc(pc+gravityCorrection());

	while(elemPtr)
	{
		if(  (elemPtr->isConnectedToEntry() && (criteria == escapeToInlet || criteria == escapeToEither))
		  || (elemPtr->isConnectedToExit() && (criteria == escapeToOutlet || criteria == escapeToEither))  )
		{
			for(size_t i = 0; i < trappingStorage.size(); ++i)
			{
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




void Element::findMarkTrappedWaterGanglia(double prs, FluidBlob startPt, vector< pair<Element*,FluidBlob> >& trappingStorageWat,
								   double& elapsed, TrappingCriteria criteria)
{
	if(!isEntryOrExitRes() && model_->conductAnyWaterBlob(startPt))
	{
		clock_t startTrapRoutine(clock());
		if(criteria == escapeToBoth)
		{
			if  (foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, escapeToOutlet)) ;
			else foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, escapeToInlet );            
		}
		else     foundEscapePathWat_trapOtherwise(prs, startPt, trappingStorageWat, criteria);

		elapsed += (double)(clock() - startTrapRoutine) / CLOCKS_PER_SEC;
	}
}




	inline void Element::trapWat(double prs, FluidBlob blob)
	{
		if(blob == bulkBlob)
		{
			trapIndexWatBulk_.first = comn_.newWatTrappingIndex();
			trapIndexWatBulk_.second = prs;
		}
		else
		{
			trapIndexWatFilm_.first = comn_.newWatTrappingIndex();
			trapIndexWatFilm_.second = prs;
		}
	}





	inline Element* Element::nextSuccessorWat(TrappingCriteria criteria, FluidBlob& blob)
	{
		if(criteria == escapeToInlet)
		{
			for(short i = connectionNum_; i > 0; --i)
			{
				if(!connections_[i-1]->isEntryOrExitRes() &&
					!connections_[i-1]->isTrappedWat(blob) &&
					connections_[i-1]->model()->conductAnyWaterBlob(blob))
					//(connections_[i-1]->model()->conductsAnyWaterBlob(blob) || connections_[i-1]->connectedToEntryOrExit()))
					return connections_[i-1];
			}
		}
		else
		{
			for(short i = 0; i < connectionNum_; ++i)
			{
				if(!connections_[i]->isEntryOrExitRes() &&
					!connections_[i]->isTrappedWat(blob) &&
					connections_[i]->model()->conductAnyWaterBlob(blob))
					//(connections_[i]->model()->conductsAnyWaterBlob(blob) || connections_[i]->connectedToEntryOrExit()))
					return connections_[i];
			}
		}

		FluidBlob blobbo(blob == filmBlob ? bulkBlob: filmBlob);
		if(!model_->disConectedCentreWCornerW() && model_->conductAnyWaterBlob(blobbo) && !isTrappedWat(blobbo))
		{
			blob = blobbo;
			return this;
		}

		return NULL;
	}

inline bool Element::foundEscapePathWat_trapOtherwise(double pc, FluidBlob startPt, vector< pair<Element*,FluidBlob> >& trappingStorage,
										TrappingCriteria criteria)
{
	pair<Element*, FluidBlob> elem(this, startPt);
	stack< pair<Element*, FluidBlob> > elemStack;
	elemStack.push(elem);
	double datumPc(pc+gravityCorrection());

	while(elem.first)
	{
		if((elem.first->isConnectedToEntry() && (criteria == escapeToInlet || criteria == escapeToEither)) ||
			(elem.first->isConnectedToExit() && (criteria == escapeToOutlet || criteria == escapeToEither)))
		{
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

		while(!elem.first)
		{
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
void Element::identifyConnectedPoreElems()
{
	setConnectedToNetwork(true);      // Do allow escapes through entry reservoir

	for(int throat = 0; throat < connectionNum_; ++throat)
	{
		if(!connections_[throat]->connectedToNetwork())                // Might have been marked previously
		{
			set< Element * > frontier, oldFrontier;
			frontier.insert(connections_[throat]);

			while(!frontier.empty())
			{
				oldFrontier = frontier;
				frontier.clear();

				for(auto elm:oldFrontier)
				{
					if(!elm->connectedToNetwork())                   // Flag is not set
					{
						elm->setConnectedToNetwork(true);               // Set flag
						for(int i = 0; i < elm->connectionNum(); ++i)
							frontier.insert(elm->connection(i));
					}
				}
			}
		}
	}
}



/// We have the option to not drain pores that are connected to the network through a
/// single throat. These pores (and throats) are removed completley from the network
/// initially. This is found to greatly improve solver convergence
int Element::removeFromNetwork()
{
	int numRemoved(1);
	for(int i = 0; i < connectionNum_; ++i)
	{
		if(connections_[i]->connectedToNetwork())
		{
			if(connections_[i]->connectionNum() > 2)  connections_[i]->severConnection(this);
			else                                      numRemoved += connections_[i]->removeFromNetwork();
		}
	}
	connections_.clear(); connectionNum_=0; connectedToNetwork_=false;
	return numRemoved;
}

/// Physically remove the single connetion to the network. This cluster is now isolated
void Element::severConnection(Element* elm)
{
	ensure(connectionNum_ > 0);
	vector<Element*>::iterator delm = find(connections_.begin(), connections_.end(), elm);
	if(*delm == elm)  connections_.erase(delm);
	else              cerr << endl << "oooops..." << endl;
	--connectionNum_;
}


















/**
// The corner conductancees are computed using the suggested procedure by Patzek SPE59312
// Elements that are trapped have their capillary pressure frozen at the time they're
// identified as being trapped.
*/
double Element::updateSat_calcR(double cappPrs)
{
	if(isEntryOrExitRes()) return 0.0;

	waterSaturation_ = model_->calcR(cappPrs);

	return isInCalcBox_ ? waterSaturation_ * flowVolume_ + clayVolume_ : 0.0;
}


/**
// As an elemnt is drained, its neighbours have an increase in number of elems
// containing oil.
*/
void Element::fillElemCentreWithOilRemoveLayers()
{
	//ensure((comn_.water()));
	//double warninnng;
	if(model_->bulkFluid()->isOil()) return;


	if (!isInOilFloodVec_ && debugLevel>0) cout<<"Qo";
	isInOilFloodVec_ = false;
	model_->fillCentreWithOilRemoveLayers();

	for(int i = 0; i < connectionNum_; ++i)
	{
		connections_[i]->IncreaseNumOilCentreFeederNeis();
		if (!rockIndex())	connections_[i]->ReduceNumWatCentreFeederNeis();
	}

	///. rubbish 
	if(model_->displacementType() == 'S')
	{
		//model_->addFillingEventToHistory(2);
		eventI_ = -1;
	}
	else
	{
		//model_->addFillingEventToHistory(1);
		eventI_ = connectionNum_ - numOilCentreFeederNeis();
	}
}

/**
// As an elemt is imbided we need to update the imbibition entry pressures of
// neighbouring elements.
*/
void Element::fillElemCentreWithWaterCreateLayers(bool snapOffOverRide)
{
	if (debugLevel>0)
	{
		if (!snapOffOverRide && !isInWatFloodVec_) cout<<"Qw";
		else if (snapOffOverRide && isInWatFloodVec_) cout<<"Qi";
	}

	isInWatFloodVec_ = false;

	ensure(model_->bulkFluid()->isOil(), "RefillingWithWater ");

	model_->fillCentreWithWaterCreateLayers(snapOffOverRide);

	for(int i = 0; i < connectionNum_; ++i)
	{
		connections_[i]->ReduceNumOilCentreFeederNeis(); ///. 
		if (!rockIndex())	connections_[i]->IncreaseNumWatCentreFeederNeis(); 
	}

	if(model_->displacementType() == 'S')
	{
		//model_->addFillingEventToHistory(4);
		eventI_ = -1;
	}
	else
	{
		//model_->addFillingEventToHistory(3);
		eventI_ = numOilCentreFeederNeis();
	}
}



bool Element::canBeAddedToEventVec(const Fluid& injectant) const
{
	if(model_->canNOTReconfigure(injectant) || !connectedToNetwork_)
		return false;

   if(injectant.isOil())
   {
		if (model_->containCOil() && rockIndex() == 0) 
		{
			if (debugLevel>0) cout<<" sgsobqbd ";
			return false;//         || nonTrappedOilCentreNeighbour(injectant));
		}

		if( isInOilFloodVec_ || (trapIndexWatBulk_.first > -1) ) return false;

		if (model_->conductsAnyOil()) return true;//         || nonTrappedOilCentreNeighbour(injectant));


		for(int i = 0; i < connectionNum_; ++i)
		{
			if(connections_[i]->model()->conductCOil() && !connections_[i]->isTrappedOil()) /// second is extra to delete
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


		for(int i = 0; i < connectionNum_; ++i)
			if(connections_[i]->model()->conductCWater())
				return true;


		//if(numWatCentreFeederNeis_ == 0) return false;

		return false;
	}
 
}





///.  rare funcs
bool Element::addToLayerVec(const Fluid& injectant, Polygon* shyp, std::vector< int >& cornToAdd) const
{
	if( !shyp->hasOilLayer() || shyp->containCOil() || trapIndexOil_.first != -1 ) return false;


	bool oilInj = injectant.isOil();

	if(oilInj)     // Layer refcreation  
	{

	  if(isInOilFloodVec())
	  {
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
			for(size_t j = 0; j < connections_.size(); ++j)
				if(connections_[j]->model()->conductsAnyOil()) 
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




void Element::calcCentreEntryPrsOilInj()
{
	entryPc_=model_->centreEntryPrsOilInj();
}

void Element::calcCentreEntryPrsWatInj()
{
	entryPc_=model_->centreEntryPrsWatInj();
}




//#include "ElementConstSolver.cpp"
//#include "ElementPore.cpp"
//#include "ElementThroat.cpp"

