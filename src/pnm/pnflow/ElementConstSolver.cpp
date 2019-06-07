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
//#include "elem_hetroPorous.h"
//#include "elem_porous.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"











/**
* Returns the neighbouring pore, going through the throat. Returns a NULL pointer
* if there is no fluid conductance between the pores or if the next pore does not
* exist (ie current pore is on a boundary). The resulting conductance between the
* pores is also calculated. Conductance is passed in as a reference.
*/
Element* Pore::getConnectionProp(int conn, double& conductance, double& deltaGrav, const Fluid *fluid, bool resistivitySolve) const
{
	if (conn >= m_connectionNum)
	{
		cerr << "=========================================" << endl
			<< "Error: Requesting network connection that" << endl
			<< "does not exist.                          " << endl
			<< "=========================================" << endl;        exit(-1);
	}

	Throat *throat = dynamic_cast< Throat* >(m_connections[conn]);
	ensure(throat);

	Element *nextPore;
	double nextLength, thisLength, throatLength(throat->length());
	if(throat->connection(0) == this)
	{
		nextPore = throat->connection(1);
		nextLength = throat->poreLength(1);
		thisLength = throat->poreLength(0);
	}
	else
	{
		nextPore = throat->connection(0);
		nextLength = throat->poreLength(0);
		thisLength = throat->poreLength(1);
	}

	bool filmConnection(
		(!fluid->isOil() && 
		 canWFilmBePassedToSolver() && throat->canWFilmBePassedToSolver() && (   nextPore->canWFilmBePassedToSolver()
	 || ( nextPore->isOnSlvrBdr() && (nextPore->isEntryOrExitRes() || nextPore->model()->conductsAny(fluid)) )     )
	)
	);  ///.wrong condtainsCFluid
 
    bool bulkConnection( 
		 fluid->isOil() ?
		 (  canOBePassedToSolver(/*bulkBlob*/) && throat->canOBePassedToSolver(/*bulkBlob*/) && (  nextPore->canOBePassedToSolver(/*bulkBlob*/) || 
           (nextPore->isOnSlvrBdr() && (nextPore->model()->conductsAny(fluid) || nextPore->isEntryOrExitRes())) ) 
         )
         :
		 (  
		   canWBulkBePassedToSolver() && throat->canWBulkBePassedToSolver() && (  nextPore->canWBulkBePassedToSolver() || 
           (nextPore->isOnSlvrBdr() && (nextPore->model()->conductsAny(fluid) || nextPore->isEntryOrExitRes())) ) 
         )
       );

	bool paralell( !fluid->isOil() && 
		throat->model()->disConectedCentreWCornerW() && filmConnection && bulkConnection &&
		m_elemModel->disConectedCentreWCornerW() && nextPore->model()->disConectedCentreWCornerW());
 
	if (filmConnection || bulkConnection)
	{
		double flowResistance(0.0);
		bool neighbourToExit(nextPore->isEntryOrExitRes());

		if (resistivitySolve)
		{
			double thisCond, throatCond, nextCond;
			//double resistivity = fluid->resistivity();
			thisCond = m_elemModel->electricalConductance(fluid);
			throatCond = throat->model()->electricalConductance(fluid);
			nextCond = nextPore->model()->electricalConductance(fluid, neighbourToExit);

			//cout<<thisCond <<" "<< throatCond <<" "<< nextCond <<" "<<nextPore->isEntryOrExitRes()<< endl;
			ensure(thisCond > 0.0 &&  throatCond > 0.0 && nextCond > 0.0);

			flowResistance = (throatLength/throatCond + thisLength/thisCond + nextLength/nextCond);
			double A_clay = (throat->clayVolume() + clayVolume() + nextPore->clayVolume()) / (nextLength + thisLength + throatLength);
			double condd = 1.0 / flowResistance + A_clay / (nextLength + thisLength + throatLength) / (m_comn.clay().resistivity());
			flowResistance = 1.0 / condd;
 
		}

		else if(paralell)
		{
			double thisCond, throatCond, nextCond;
			ensure(dynamic_cast< const Water* >(fluid));
			thisCond = m_elemModel->getWaterConductance(filmBlob);
			throatCond = throat->model()->getWaterConductance(filmBlob);
			nextCond = nextPore->model()->getWaterConductance(filmBlob, neighbourToExit);
			ensure(thisCond*throatCond*nextCond != 0.0);

			if(thisCond*throatCond*nextCond > 0.0)
				flowResistance = (throatLength/throatCond + thisLength/thisCond + nextLength/nextCond);
			else
				flowResistance = 1.0e32;

		        double thisCondBulk = m_elemModel->getWaterConductance(bulkBlob);
		        double throatCondBulk = throat->model()->getWaterConductance(bulkBlob);
		        double nextCondBulk = nextPore->model()->getWaterConductance(bulkBlob, neighbourToExit);
		        ensure(thisCondBulk*throatCondBulk*nextCondBulk != 0.0);

		        double flowResistBulk = 1.0e32;

				if(thisCondBulk*throatCondBulk*nextCondBulk > 0.0)
					flowResistBulk = (throatLength/throatCondBulk+thisLength/thisCondBulk+nextLength/nextCondBulk);

		        flowResistance = 1.0/(1.0/flowResistance + 1.0/flowResistBulk);
		}
		else
		{
			double thisCond, throatCond, nextCond;
			thisCond = m_elemModel->getConductance(fluid);
			throatCond = throat->model()->getConductance(fluid);
			nextCond = nextPore->model()->getConductance(fluid, neighbourToExit);



			if(thisCond*throatCond*nextCond <=  1.0e-150 )
			{   if(m_comn.debugMode>10)
				{
					if  ( fluid ==  &m_comn.water())
					cout
					<<" filmSolverFlag  " << canWFilmBePassedToSolver() <<" "<<(throat->canWFilmBePassedToSolver()) <<" "<< (nextPore->canWFilmBePassedToSolver()) << "\n"
					<<" bulkSolverFlag  " << canWBulkBePassedToSolver()<<" "<< (throat->canWBulkBePassedToSolver())<<" "<< (nextPore->canWBulkBePassedToSolver()) << "\n"
					<<" rockType        "<<iRockType()<<" "<<throat->iRockType()<<" "<<nextPore->iRockType()<<"\n"; 
					else
					cout
					//<<" filmSolverFlag  " << canOBePassedToSolver(filmBlob) <<" "<<(throat->canOBePassedToSolver(filmBlob)) <<" "<< (nextPore->canOBePassedToSolver(filmBlob)) << "\n"
					<<" bulkSolverFlag  " << canOBePassedToSolver(/*bulkBlob*/)<<" "<< (throat->canOBePassedToSolver(/*bulkBlob*/))<<" "<< (nextPore->canOBePassedToSolver(/*bulkBlob*/)) << "\n"
					<<" rockType        "<<iRockType()<<" "<<throat->iRockType()<<" "<<nextPore->iRockType()<<"\n"; 

					cout<< "cond: this = " << thisCond <<" throat = "<< throatCond <<" next = "<< nextCond <<",  area = "<<throat->model()->area()<< "\n\n";
					//cout<<nextPore->isTrappedWat(filmBlob)<<nextPore->isTrappedWat(bulkBlob);
				}
				flowResistance = 5.0e49;
			}
			else
		        flowResistance = (throatLength/throatCond + thisLength/thisCond + nextLength/nextCond);
		}

		ensure(flowResistance > 0.0);

		conductance = 1.0 / flowResistance;

        if(nextPore->isEntryOrExitRes() || resistivitySolve)
            deltaGrav = 0.0;
        else
        {
            deltaGrav = rhogh(fluid->density(), m_node->xPos()-nextPore->node()->xPos(),
                m_node->yPos()-nextPore->node()->yPos(), m_node->zPos()-nextPore->node()->zPos());
        }

        throat->setPoreToPoreCond(conductance, fluid, resistivitySolve);

	}
	else
	{
		conductance = 0.0;
	}


	return nextPore;
}
























//inline bool Element::isTouchedInSearch(FluidBlob blob) const
//{
	//if(blob == bulkBlob)
        //return m_isTouchedInSearch.second;
	//else
        //return m_isTouchedInSearch.first;
//}

//inline bool Element::canWBePassedToSolver() const
//{
		//return  m_canWBulkBePassedToSolver || m_canWFilmBePassedToSolver;
//}


//bool Element::canOBePassedToSolver() const
//{
	//return m_canOBulkBePassedToSolver;
//}




/**
// Passes back the prevous pressure determined when solving the laplace eq.
*/
bool Pore::prevSolvrRes(const Fluid* fluid, int resistSolve, double loc, double& res, double& flowRate) const
{
	
	
    if(fluid->isOil())
    {
        if(resistSolve)
            res = m_oilSolverVolt;
        else
            res = m_oilSolverPrs;
    }    
    else
    {
        if(resistSolve)
            res = m_watSolverVolt;
        else
            res = m_watSolverPrs;
    }
    
	if ( res<-1.0e+12 )
	{
		res = 0.0;
	}


    return true;
}


void Pore::setSolverResults(const Fluid* fluid, int resistSolve, double res) const
{
    //ensure(///. we add -1.5e24 for those with dummy value
    //((canBePassedToSolver() || m_isOnInletSlvrBdr || m_isOnOutletSlvrBdr) && res>-1.0) || res<-1.0);

		
    if(fluid->isOil())
    {
		if ( res<-1.0e-12 && (canOBePassedToSolver() || m_isOnInletSlvrBdr || m_isOnOutletSlvrBdr) )
		{
			if ( res>-1.4e24 || res<-1.6e24 )
			{
				if (m_comn.debugMode>100) 
				{ 
					std::cout<<"E";
				}
				++nErrs;
			}
			//else
				//res = 0;
		}

        if(resistSolve)
            m_oilSolverVolt = res;
        else
            m_oilSolverPrs = res;
    }    
    else
    {
		if ( res<-1.0e-12 && (canWBePassedToSolver() || m_isOnInletSlvrBdr || m_isOnOutletSlvrBdr) )
		{
			if ( res>-1.4e24 || res<-1.6e24 )
			{
				//if (m_comn.debugMode>100) std::cout<<"E";
				++nErrs;
			}
		}

		if(resistSolve)
            m_watSolverVolt = res;
        else
            m_watSolverPrs = res;
    }

}


/**
// Passes back the prevous pressure determined when solving the laplace eq.
// The routine for throats interpolates between the connecting pore pressures and
// also calculates the flowrate between the pores
*/
bool Throat::prevSolvrRes(const Fluid* fluid, int resistSolve, double loc, double& res, double& flowRate) const
{
    double resZero(0.0), resOne(0.0), tmp(0.0), gravCorrZero(0.0), gravCorrOne(0.0);

    m_connections[0]->prevSolvrRes(fluid, resistSolve, loc, resZero, tmp);
    m_connections[1]->prevSolvrRes(fluid, resistSolve, loc, resOne, tmp);

    if(USE_GRAV_IN_KR && !resistSolve)
    {
        gravCorrZero = rhogh(fluid->density(), m_connections[0]->node()->xPos(),
            m_connections[0]->node()->yPos(), m_connections[0]->node()->zPos());
        gravCorrOne = rhogh(fluid->density(), m_connections[1]->node()->xPos(),
            m_connections[1]->node()->yPos(), m_connections[1]->node()->zPos());
    }

    if(m_connections[0]->isEntryOrExitRes())                        // Datum correct in/outlet reservoirs. These are always assumed
    {                                                               // to be at the same lavel as connecting pore.
        resZero += gravCorrOne-gravCorrZero;
        gravCorrZero = gravCorrOne;
        if (resistSolve==2) resZero = resOne;
    }

    if(m_connections[1]->isEntryOrExitRes())
    {
        resOne += gravCorrZero-gravCorrOne;
        gravCorrOne = gravCorrZero;
        if (resistSolve==2) resOne = resZero;
    }

    double locOne(m_connections[0]->node()->xPos());
    double locTwo(m_connections[1]->node()->xPos());

    ensure(m_poreToPoreCond>=0.0);

    if(m_poreToPoreCond > 0.0)
    {
        flowRate = locOne < locTwo ?
            (resZero-gravCorrZero-resOne+gravCorrOne)*m_poreToPoreCond:
            (resOne-gravCorrOne-resZero+gravCorrZero)*m_poreToPoreCond;
    }
    else
        flowRate = 0.0;

    resZero -= gravCorrZero;                          // Pressures passed back are datum levelled
    resOne -= gravCorrOne;
    if (abs(locTwo-locOne) <= 1.0e-32)
		res = resZero;
    else
		res = resZero + (resOne-resZero) * (loc-locOne) / (locTwo-locOne);

    ensure((loc >= locOne && loc <=  locTwo) || (loc >= locTwo && loc <=  locOne));
    ensure(res >= min(resZero, resOne) && res <=  max(resZero, resOne));

    return m_poreToPoreCond > 0.0;
}


inline void Element::setWSolverFlag(FluidBlob blob) const
{
    if(blob == bulkBlob)
    {
		//m_canBePassedToSolver.second = true;
        //m_isTouchedInSearch.second = true;
		m_canWBulkBePassedToSolver = true;
        m_isWBulkTouchedInSearch = true;
    }
	else
    {
		//m_canBePassedToSolver.first = true;
        //m_isTouchedInSearch.first = true;
		m_canWFilmBePassedToSolver = true;
        m_isWFilmTouchedInSearch = true;
    }
}

inline void Element::clearWSolverFlag(FluidBlob blob) const
{
	//if(m_iRockType>0) return;
	int Warning;
    if(blob == bulkBlob)
		m_canWBulkBePassedToSolver = false;
	else
		m_isWFilmTouchedInSearch = false;

}

void Element::clearAllWSolverFlag() const
{
    //m_canBePassedToSolver.first = false;
    //m_isTouchedInSearch.first = false;
    //m_canBePassedToSolver.second = false;
    //m_isTouchedInSearch.second = false;
		m_canWBulkBePassedToSolver = false;
        m_isWBulkTouchedInSearch = false;
		m_canWFilmBePassedToSolver = false;
        m_isWFilmTouchedInSearch = false;
         
    m_poreToPoreCond = 0.0;
}

inline void Element::setOSolverFlag(FluidBlob blob) const
{
		m_canOBulkBePassedToSolver = true;
        m_isOBulkTouchedInSearch = true;

}

inline void Element::clearOSolverFlag(FluidBlob blob) const
{
		m_canOBulkBePassedToSolver = false;
}

void Element::clearAllOSolverFlag() const
{
		m_canOBulkBePassedToSolver = false;
        m_isOBulkTouchedInSearch = false;
        
    m_poreToPoreCond = 0.0;
}








	inline const Element* Element::nextSuccSolvrOil(bool& outletFound) const
	{
		for(short i = 0; i < m_connectionNum; ++i)
		{
			if(!m_connections[i]->isOTouchedInSearch(/*bulkBlob*/) &&
				m_connections[i]->isInsideSolverBox() &&
				((m_connections[i]->model()->conductsAnyOil() && m_connections[i]->model()->oilConductance() > COND_CUT_OFF) || m_connections[i]->connectedToEntryOrExit()) 
				//// && (m_connections[i]->model()->oilConductance() > COND_CUT_OFF || m_connections[i]->connectedToEntryOrExit())
				)
			{
				return m_connections[i];
			}
			else if(m_connections[i]->isOnOutletSlvrBdr() &&
				(m_connections[i]->model()->conductsAnyOil() || m_connections[i]->isEntryOrExitRes()))
			{
				outletFound = true;
			}
		}

		return NULL;
	}


/**
// Before the solver can solve the pressure field we need to ensure that there is pressure
// communication. Previous trapping routines are not sufficient for this as region connected
// through oil layers might have become trapped without having been identified. This is done
// through a depth first traversal from the elements at the inlet boundary. All connected
// elements (through the defined fluid) are marked. These marked elements can then be passed to
// the solver. The routine is somewhat more involved than the trapping routine since we do not
// want connections through in/outlet to be marked. What we want is only the region of nodes
// that are connected both to the in and outlet.
*/
bool Element::connectedToOutlet(const Fluid* fluid) const
{
    if(!m_elemModel->conductsAny(fluid) && !isEntryOrExitRes()) return false;


	if(dynamic_cast<const Water*>(fluid))
		for(int throat = 0; throat < m_connectionNum; ++throat)
		{            
			pair<const Element *, FluidBlob> elem(m_connections[throat], filmBlob);

            if(markWaterElemForSolver(elem)) return true;

            elem.second = bulkBlob;
            if(markWaterElemForSolver(elem)) return true;
		}
	else
		for(int throat = 0; throat < m_connectionNum; ++throat)
		{  
            //if(markOilElemForSolver(m_connections[throat])) return true;
			const Element* elem = m_connections[throat];
			
			if( (!elem->model()->conductsAnyOil() &&	!elem->connectedToEntryOrExit()) 
				||	elem->isOTouchedInSearch() || !elem->isInsideSolverBox() )
			{
				continue;
			}

			vector<const Element*> trappingStorage;
			stack<const Element*> elemStack;
			elemStack.push(elem);
			bool outletFound(false);

			while(elem)
			{
				elem->setOSolverFlag();
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
					trappingStorage[k]->clearOSolverFlag(bulkBlob);
			else
				return true;
            
            
        }

    return false;
}




		inline const Element* Element::nextSuccSolvrWat(bool& outletFound, FluidBlob& blob) const
		{
			for(short i = 0; i < m_connectionNum; ++i)
			{
				if(!m_connections[i]->isWTouchedInSearch(blob) && m_connections[i]->isInsideSolverBox() &&
					(m_connections[i]->model()->conductAnyWaterBlob(blob) || m_connections[i]->connectedToEntryOrExit()))
				{
					return m_connections[i];
				}
				else if(m_connections[i]->isOnOutletSlvrBdr() &&
					(m_connections[i]->model()->conductAnyWater() || m_connections[i]->isEntryOrExitRes()))
				{
					outletFound = true;
				}
			}

			FluidBlob blobbo(blob == filmBlob ? bulkBlob: filmBlob);
			if(!m_elemModel->disConectedCentreWCornerW() && m_elemModel->conductAnyWaterBlob(blobbo) && !isWTouchedInSearch(blobbo))
			{
				blob = blobbo;
				return this;
			}

			return NULL;
		}
		
	bool Element::markWaterElemForSolver(pair<const Element*, FluidBlob> elem) const
	{
		if( (!elem.first->model()->conductAnyWaterBlob(elem.second) && !elem.first->connectedToEntryOrExit())
			 || elem.first->isWTouchedInSearch(elem.second) ||	!elem.first->isInsideSolverBox())
		{
			return false;
		}

		vector< pair<const Element*,FluidBlob> > trappingStorage;
		stack< pair<const Element*, FluidBlob> > elemStack;
		elemStack.push(elem);
		bool outletFound(false);

		while(elem.first)
		{
			elem.first->setWSolverFlag(elem.second);
			trappingStorage.push_back(elem);
			elem.first = elem.first->nextSuccSolvrWat(outletFound, elem.second);

			while(!elem.first)
			{
				elemStack.pop();
				if(elemStack.empty()) break;
				elem = elemStack.top();
				elem.first = elem.first->nextSuccSolvrWat(outletFound, elem.second);
			}

			if(elem.first) elemStack.push(elem);
		}

		if(!outletFound)
		{
			for(size_t k = 0; k < trappingStorage.size(); ++k)         // If region was isolated it should not be passed
			{                                                       // to solver
				trappingStorage[k].first->clearWSolverFlag(trappingStorage[k].second);
			}
		}

		return outletFound;
	}





