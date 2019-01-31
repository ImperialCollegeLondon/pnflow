#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include <vector>
#include <stack>
#include <cassert>
using namespace std;

//#include "node.h"
#include "fluid.h"
//#include "compareFuncs.h"
#include "Element.h"
#include "inputData.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"


//bool Element::ERROR_STATE = false;
bool Element::USE_GRAV_IN_KR = false;
double Element::COND_CUT_OFF = 0.0;
const double Element::PI = acos(-1.0);

int Element::nErrs = 0;

double Element::RRR() const {return m_elemModel->radius();}

Element::Element(CommonData& common, const Oil& oil, const Water& water, double radius, double volume,
                 double volClay, double shapeFact, int connNum, bool amPore, std::string elemData) 
	 :  Apex(), m_randomNum(rand()), m_iAmAPore(amPore), m_comn(common), m_flowVolume(volume), m_clayVolume(volClay), m_connectionNum(connNum), 
    m_trappingIndexOil(-1, 0.0),
    m_trappingIndexWatBulk(-1, 0.0),
    m_trappingIndexWatFilm(-1, 0.0),
    m_numOilCentreFeederNeis(0),
    m_numWatCentreFeederNeis(connNum)
{
    m_waterSaturation = 1.0;
    m_isInsideSolverBox = false;
    m_isInsideSatBox = false;
    m_isInWatFloodVec = false;
    m_isInOilFloodVec = false;
    m_isOnInletSlvrBdr = false;
    m_isOnOutletSlvrBdr = false;
    m_isExitRes = false;
    m_isEntryRes = false;
    m_isInWatFloodVec = false;
    m_isInOilFloodVec = false;
    m_connectedToNetwork = false;
    m_isConnectedToExit = false;
    m_isConnectedToEntry = false;
    m_connectedToEntryOrExit = false;
    m_canWFilmBePassedToSolver = false;
    m_canWBulkBePassedToSolver = false;
    m_canOBulkBePassedToSolver = false;
    m_isWFilmTouchedInSearch = false;
    m_isWBulkTouchedInSearch = false;
    m_isOBulkTouchedInSearch = false;


    m_poreToPoreCond = 0.0;
    m_fillingEventRecord = -2;


	istringstream inputData;
	inputData.str( elemData );
	inputData>>m_iRockType;

	char hetroData(';');
	inputData>>hetroData;

	if(shapeFact <=  sqrt(3.0)/36.0+0.00001)                   // Triangular:  0 >= G >= sqrt(3.0)/36.0
    {
        m_elemModel = new Triangle(*this, m_comn, radius, min(shapeFact, sqrt(3.0)/36.0-0.00005),
            m_connectionNum, inputData);
         common.countTriangle();

    }
    else if(shapeFact < 0.07)                         // Square:      G == 1.0/16/0
    {
        m_elemModel = new Square(*this, m_comn, radius, m_connectionNum, inputData);
		common.countSquare();
    }
    else                                                // Circular:    G == 1.0/4.0*PI
    {
        m_elemModel = new Circle(*this, m_comn, radius, m_connectionNum, inputData);
		common.countCircle();
    }





	///.       set parent
	setConnections(m_elemModel,-1);


}


bool Element::convertToMicroPorosityForSven(bool entry1Exit0)
{
		return false;
}

int Element::ffaz() const {return int(m_elemModel->containCOil())+1;} ///. Viz only
double Element::saturation() const 
{ return
	(m_flowVolume*waterSaturation()+
	m_connections[0]->waterSaturation()*m_connections[0]->flowVolume()/m_connections[0]->connectionNum()+ 
	m_connections[1]->waterSaturation()*m_connections[1]->flowVolume()/m_connections[1]->connectionNum() ) /
	(m_flowVolume+
	 m_connections[0]->flowVolume()/m_connections[0]->connectionNum()+
	 m_connections[1]->flowVolume()/m_connections[1]->connectionNum()
	);
}

Element::~Element(){delete m_elemModel;}


void Element::checkConnections() const
{
    if(m_connectionNum != int(m_connections.size()))
    {
        cout<< "========================================" << endl
            << "Error: The connection number is incorrect" << endl
            << "=========================================" << endl;				exit(-1);//ERROR_STATE = true;
    }

    for(int conn = 0; conn < m_connectionNum; ++conn)
    {
        if(m_connections[conn] == NULL)
        { 
            cout<< "==================================="         << endl
                << "Error: Missing network connecetions"          << endl
                << "Interanl pointer remains NULL      "          << endl
                << "Element index: " << orenIndex()               << endl
                << "Element type: " << typeid(*this).name()       << endl
                << "Connection: " << conn                         << endl
                << "Total connection number: " << m_connectionNum << endl
                << "==================================="          << endl;	exit(-1);//ERROR_STATE = true;
        }

        Element* nextElem = m_connections[conn];

        bool connExist(false);
        for(int nxtConn = 0; nxtConn < nextElem->connectionNum(); ++nxtConn)
        {
            connExist = (nextElem->connection(nxtConn) == this);
            if (connExist) break;
        }

        if(!connExist)
        {
            cout<< "==================================="                               << endl
                << "Error: Missing network connecetions"                                << endl
                << "Interanal links do not match up    "                                << endl
                << typeid(*this).name() <<" index: " << orenIndex()                     << endl
                << "Connection: " << conn                                               << endl
                //<< "Element type: " << typeid(*this).name()                             << endl
                << "Neighbouring " <<typeid(*nextElem).name() << " index: " << nextElem->orenIndex()    << endl
                //<< "Opposing type: " << typeid(*nextElem).name()                        << endl
                //<< "Total connection number: " <<                        << endl
                << "Neighbouring " <<typeid(*nextElem).name() << "  number of connections: " << nextElem->connectionNum()    << endl
                << "===================================" << endl;			exit(-1);//ERROR_STATE = true;
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
    if(!m_isExitRes && !m_isEntryRes && m_elemModel->conductsAnyOil())
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
    m_trappingIndexOil.first = m_comn.newOilTrappingIndex();
    m_trappingIndexOil.second = prs;
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
			for(short i = 0; i < m_connectionNum; ++i)
			{
				if(!m_connections[i]->isEntryRes() &&
					!m_connections[i]->isTrappedOil() &&
					m_connections[i]->model()->conductsAnyOil())
					return m_connections[i];
			}
		}    
		else if(criteria == escapeToEither)
		{
			for(short i = 0; i < m_connectionNum; ++i)
			{
				if(// !m_connections[i]->isEntryOrExitRes() &&
					!m_connections[i]->isTrappedOil() &&
					m_connections[i]->model()->conductsAnyOil())
					return m_connections[i];
			}
		}
		else if(criteria == escapeToInlet)
		{
			for(short i = m_connectionNum-1; i >= 0; --i)
			{
				if(!m_connections[i]->isExitRes() &&
					!m_connections[i]->isTrappedOil() &&
					m_connections[i]->model()->conductsAnyOil())///. assumes elem oil is connected to all connections
					return m_connections[i];
			}
		}
		else
		{
			for(short i = 0; i < m_connectionNum; ++i)
			{
				if(//!m_connections[i]->isEntryOrExitRes() &&
					!m_connections[i]->isTrappedOil() &&
					m_connections[i]->model()->conductsAnyOil())
					return m_connections[i];
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
// This consern was stong enough to warrant seperate implementations for each fluid
*/
bool Element::foundEscapePathOil_trapOtherwise(double pc, vector<Element*>& trappingStorage, TrappingCriteria criteria)
{

    Element* elemPtr = this;
    stack<Element*> elemStack;                      // Simulate recursive behavior using a stack
    //elemStack.push(this);
    double datumPc(pc+m_elemModel->rhogh());

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

        double localPc(datumPc - elemPtr->model()->rhogh());
        softAssert(!elemPtr->isEntryOrExitRes());
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
    if(!m_isExitRes && !m_isEntryRes && m_elemModel->conductAnyWaterBlob(startPt))
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
			m_trappingIndexWatBulk.first = m_comn.newWatTrappingIndex();
			m_trappingIndexWatBulk.second = prs;
		}
		else
		{
			m_trappingIndexWatFilm.first = m_comn.newWatTrappingIndex();
			m_trappingIndexWatFilm.second = prs;
		}
	}





	inline Element* Element::nextSuccessorWat(TrappingCriteria criteria, FluidBlob& blob)
	{
		if(criteria == escapeToInlet)
		{
			for(short i = m_connectionNum; i > 0; --i)
			{
				if(!m_connections[i-1]->isEntryOrExitRes() &&
					!m_connections[i-1]->isTrappedWat(blob) &&
					m_connections[i-1]->model()->conductAnyWaterBlob(blob))
					//(m_connections[i-1]->model()->conductsAnyWaterBlob(blob) || m_connections[i-1]->connectedToEntryOrExit()))
					return m_connections[i-1];
			}
		}
		else
		{
			for(short i = 0; i < m_connectionNum; ++i)
			{
				if(!m_connections[i]->isEntryOrExitRes() &&
					!m_connections[i]->isTrappedWat(blob) &&
					m_connections[i]->model()->conductAnyWaterBlob(blob))
					//(m_connections[i]->model()->conductsAnyWaterBlob(blob) || m_connections[i]->connectedToEntryOrExit()))
					return m_connections[i];
			}
		}

		FluidBlob blobbo(blob == filmBlob ? bulkBlob: filmBlob);
		if(!m_elemModel->disConectedCentreWCornerW() && m_elemModel->conductAnyWaterBlob(blobbo) && !isTrappedWat(blobbo))
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
    double datumPc(pc+m_elemModel->rhogh());

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

        softAssert(!elem.first->isEntryOrExitRes());
        double localPc(datumPc-elem.first->model()->rhogh());
        elem.first->trapWat(localPc, elem.second);
        trappingStorage.push_back(elem);
        softAssert(elem.first);
        elem.first = elem.first->nextSuccessorWat(criteria, elem.second);

        while(!elem.first)
        {
            elemStack.pop();
            if(elemStack.empty()) break;
            elem = elemStack.top();
            softAssert(elem.first);
            elem.first = elem.first->nextSuccessorWat(criteria, elem.second);
        }

        if(elem.first) elemStack.push(elem);
    }
    return true;
}


/**
// It is possible that some parts of the network in fact is physically isolated from the
// network due to diagenesis processes. These parts will remain water filled but are not
// available for displacement. Identify these elements during initialization by doing a
// graph traversal from the inlet.
*/
void Element::identifyConnectedPoreElems()
{
    m_connectedToNetwork = true;                                            // Do allow escapes through entry reservoir
    set<Element*> investigatedThroats;

    for(int throat = 0; throat < m_connectionNum; ++throat)
    {
        if(investigatedThroats.count(m_connections[throat]) == 0 &&     // Already failed to find outlet for this
           !m_connections[throat]->connectedToNetwork())                // Might have been marked previously
        {
            set< Element * > frontier, oldFrontier;
            frontier.insert(m_connections[throat]);

            while(!frontier.empty())
            {
                oldFrontier = frontier;
                frontier.clear();

                for(ItrSet itrElem = oldFrontier.begin(); itrElem != oldFrontier.end(); ++itrElem)
                {
                    if(!(*itrElem)->connectedToNetwork())                   // Flag is not set
                    {
                        (*itrElem)->setConnectedToNetwork(true);               // Set flag
                        for(int i = 0; i < (*itrElem)->connectionNum(); ++i)
                            frontier.insert((*itrElem)->connection(i));
                    }
                }
            }
            for(int j = 0; j < m_connectionNum; ++j)
            {
                if(m_connections[j]->connectedToNetwork())              // Make sure we don't traverse same region twice. If the
                    investigatedThroats.insert(m_connections[j]);       // elem was not connected it would have been unmarked
            }
        }
    }
}



/**
// We have the option to not drain pores that are connected to the network through a
// single throat. These pores (and throats) are removed completley from the network
// initially. This is found to greatly improve solver convergence
*/
int Element::removeFromNetwork()
{
    int numRemoved(1);
    m_connectedToNetwork = false;
    for(int i = 0; i < m_connectionNum; ++i)
    {
        if(m_connections[i]->connectedToNetwork())
        {
            if(m_connections[i]->connectionNum() > 2)
            {
                m_connections[i]->severConnection(this);
                severConnection(m_connections[i]);
            }
            else
                numRemoved += m_connections[i]->removeFromNetwork();
        }
    }
    return numRemoved;
}

/**
// Physically remove the single connetion to the network. This cluster is now isolated
*/
void Element::severConnection(Element* connection)
{
    softAssert(m_connectionNum > 0);
    vector<Element*>::iterator delCandidate;
    delCandidate = find(m_connections.begin(), m_connections.end(), connection);
    if(*delCandidate == connection)
        m_connections.erase(delCandidate);
    else
    {
        cerr << endl << "oooops..." << endl;
        //exit(-1);
    }
    --m_connectionNum;
    //m_elemModel->reduceNeighbourCount();
}


















/**
// The corner conductancees are computed using the suggested procedure by Patzek SPE59312
// Elements that are trapped have their capillary pressure frozen at the time they're
// identified as being trapped.
*/
double Element::updateSat_calcR(double cappPrs)
{
    if(m_isEntryRes || m_isExitRes) return 0.0;

    m_waterSaturation = m_elemModel->calcR(cappPrs);

    return m_isInsideSatBox ? m_waterSaturation * m_flowVolume + m_clayVolume : 0.0;
}


/**
// As an elemnt is drained, its neighbours have an increase in number of elems
// containing oil.
*/
void Element::fillElemCentreWithOilRemoveLayers()
{
    //softAssert((m_comn.water()));
    //double warninnng;
    if(m_elemModel->bulkFluid()->isOil()) return;
    
    
    if (!m_isInOilFloodVec && m_comn.debugMode>0) cout<<"Qo";
	m_isInOilFloodVec = false;
    m_elemModel->fillCentreWithOilRemoveLayers();

    for(int i = 0; i < m_connectionNum; ++i)
    {
        m_connections[i]->IncreaseNumOilCentreFeederNeis();
        if (!iRockType())	m_connections[i]->ReduceNumWatCentreFeederNeis();
	}

	///. rubbish 
    if(m_elemModel->displacementType() == 'S')
    {
        //m_elemModel->addFillingEventToHistory(2);
        m_fillingEventRecord = -1;
    }
    else
    {
        //m_elemModel->addFillingEventToHistory(1);
        m_fillingEventRecord = m_connectionNum - numOilCentreFeederNeis();
    }
}

/**
// As an elemt is imbided we need to update the imbibition entry pressures of
// neighbouring elements.
*/
void Element::fillElemCentreWithWaterCreateLayers(bool snapOffOverRide)
{
	if (m_comn.debugMode>0)
	{
		if (!snapOffOverRide && !m_isInWatFloodVec) cout<<"Qw";
		else if (snapOffOverRide && m_isInWatFloodVec) cout<<"Qi";
	}
    
	m_isInWatFloodVec = false;
	
    if(!m_elemModel->bulkFluid()->isOil()) 
    {
		//bool accurat= (m_elemModel[1000000].shapeFactor()>0);
		//cout<<accurat;
		 cout<<" Qisadasklaaqwghhwetj ";
	}

    m_elemModel->fillCentreWithWaterCreateLayers(snapOffOverRide);

    for(int i = 0; i < m_connectionNum; ++i)
    {
        m_connections[i]->ReduceNumOilCentreFeederNeis(); ///. 
        if (!iRockType())	m_connections[i]->IncreaseNumWatCentreFeederNeis(); 
	}

    if(m_elemModel->displacementType() == 'S')
    {
        //m_elemModel->addFillingEventToHistory(4);
        m_fillingEventRecord = -1;
    }
    else
    {
        //m_elemModel->addFillingEventToHistory(3);
        m_fillingEventRecord = numOilCentreFeederNeis();
    }
}



bool Element::canBeAddedToEventVec(const Fluid* injectant) const
{
    if(m_elemModel->canNOTReconfigure(injectant) || !m_connectedToNetwork)
        return false;

   if(injectant->isOil())
   {
        if (m_elemModel->containCOil() && iRockType() == 0) 
        {
			if (m_comn.debugMode>0) cout<<" sgsobqbd ";
			return false;//         || nonTrappedOilCentreNeighbour(injectant));
		}

        if( m_isInOilFloodVec || (m_trappingIndexWatBulk.first > -1) ) return false;

        if (m_elemModel->conductsAnyOil()) return true;//         || nonTrappedOilCentreNeighbour(injectant));


		for(int i = 0; i < m_connectionNum; ++i)
		{
			if(m_connections[i]->model()->conductCOil() && !m_connections[i]->isTrappedOil()) /// second is extra to delete
				return true;
		}
		
		if(m_numOilCentreFeederNeis == 0) return false;
		
		return false;
	}
	else
    {
        if (!m_elemModel->conductsAnyOil() && iRockType() == 0) 
        {
			if (m_comn.debugMode>0) cout<<" sjswbqld ";
			return false;//         || nonTrappedOilCentreNeighbour(injectant));
		}
		
        if(m_isInWatFloodVec || (m_trappingIndexOil.first > -1) ) return false;
        
        if (m_elemModel->conductAnyWater() ) return true;//         || nonTrappedWatNeighbour(injectant));
        

		for(int i = 0; i < m_connectionNum; ++i)
			if(m_connections[i]->model()->conductCWater() && !m_connections[i]->isTrappedWat(bulkBlob)) /// second is extra to delete
				return true;


		//if(m_numWatCentreFeederNeis == 0) return false;

		return false;
	}
 
}





///.  rare funcs
bool Element::addToLayerVec(const Fluid* injectant, Polygon* polyShape, std::vector< int >& cornToAdd) const
{
    if( !polyShape->hasOilLayer() || polyShape->containCOil() || m_trappingIndexOil.first != -1 ) return false;


    bool oilInj = injectant->isOil();

    if(oilInj)     // Layer refcreation  
    {

      if(isInOilFloodVec())    
        for(int i = 0; i < polyShape->numCorners(); ++i)
            if(!polyShape->oilLayerConst()[i].isInOilFloodVec() &&
                 ! polyShape->oilLayerConst()[i].exists() 
                 && entryPc() > polyShape->oilLayerConst()[i].entryPc()
              )
                cornToAdd.push_back(i);
      else
        for(int i = 0; i < polyShape->numCorners(); ++i)
            if(!polyShape->oilLayerConst()[i].isInOilFloodVec() &&
                 ! polyShape->oilLayerConst()[i].exists() 
              )
                cornToAdd.push_back(i);

        if(!cornToAdd.empty())
            for(size_t j = 0; j < m_connections.size(); ++j)
                if(m_connections[j]->model()->conductsAnyOil()) 
                    return true;

        return false;
    }
    else        // Layer collapse
    {
        for(int i = 0; i < polyShape->numCorners(); ++i)
            if(!polyShape->oilLayerConst()[i].isInWatFloodVec() && 
                polyShape->oilLayerConst()[i].exists() &&
                m_trappingIndexOil.first < 0 &&
                (m_trappingIndexWatFilm.first < 0 || m_trappingIndexWatBulk.first < 0)
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
	m_entryPc=m_elemModel->centreEntryPrsOilInj();
}

void Element::calcCentreEntryPrsWatInj()
{
	m_entryPc=m_elemModel->centreEntryPrsWatInj();
}




//#include "ElementConstSolver.cpp"
//#include "ElementPore.cpp"
//#include "ElementThroat.cpp"

