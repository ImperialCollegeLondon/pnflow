#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
using namespace std;

#include "inputData.h"
#include "Element.h"
#include "solver.h"
#include "netStatistics.h"


#include "netsim.h"




inline void Netsim::insertReCalcImbibeEntryPrs(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, Element* elem, double cappPrs)
{
	if(elem->isInWatFloodVec())
	{
		if (m_eventsCh.remove(elem))
		{   ///. order is important
				elem->calcCentreEntryPrsWatInj();
				m_eventsCh.insert(elem);
                if (!elem->model()->bulkFluid()->isOil()) cout<< " deaw " ;
		}
		else
		{
			elem->setInWatFloodVec(false);
			elem->calcCentreEntryPrsWatInj();
			if (elem->canBeAddedToEventVec(&m_water))
			{

				m_eventsCh.insert(elem);		elem->setInWatFloodVec(true);
				//cout<<"|";
			}
			else
			{
					cout<<"~";
			}
		}
	}
	else
	{
		if (elem->canBeAddedToEventVec(&m_water))
		{
			elem->ChModel()->initWaterInjection(cappPrs);
			elem->calcCentreEntryPrsWatInj();
			m_eventsCh.insert(elem);
			elem->setInWatFloodVec(true);
		}
	}

	//elem->ChModel()->calcR(cappPrs);///. candidate to be removed

} 



inline void Netsim::addElemTo_layerImbibeEvents(SortedEvents<Apex*,PceImbCmp>& m_layerEventsCh, Element* elem)
{
    Polygon* polyShape = dynamic_cast< Polygon* >(elem->ChModel());
    if(polyShape)
    {
		vector<int> addCrns;
		//double rhogh(elem->model()->rhogh());
		if(elem->addToLayerVec(invadingFluid, polyShape, addCrns))
		{
			softAssert(!addCrns.empty());
			for(size_t i = 0; i < addCrns.size(); ++i)
			{
				softAssert(!polyShape->oilLayerConst()[addCrns[i]].isInWatFloodVec());
				int dd45454jnn;
				m_layerEventsCh.insert(&polyShape->oilLayerCh()[addCrns[i]]);
				polyShape->oilLayerCh()[addCrns[i]].setInWatFloodVec(true);
			}
		}
	}
}


inline void Netsim::clearTrappedOilFromEvents(SortedEvents<Apex*,PceImbCmp>&    m_eventsCh)
{
    while(!m_eventsCh.empty() && m_eventsCh.peek()->parentModel()->eleman()->isTrappedOil())
    {
        m_eventsCh.pop()->setInWatFloodVec(false);
    }

    while(!m_layerEventsCh.empty() && m_layerEventsCh.peek()->trappingCL().first>-1)
    {
        m_layerEventsCh.pop()->setInWatFloodVec(false);
    }
}





/**
* Inject water (imbibition) until target water saturation is reached or alternatively that the queue
* containing elements to be popped is empty. After drainage throats connected to the outlets have a
* maximum of one oil neighbour. If all elements were drained water injection will then occur at both
* faces. We can remove an injection face by increasing the oil neighbour count in throats connected to
* that that face. If we remove both faces and all elements were previously drained, the firts imbibition
* event will be a snap off.
 */
void Netsim::Imbibition(InputData& input, double& Sw, double& Pc, double requestedFinalSw, double requestedFinalPc,
              double deltaSw, double deltaPc, double deltaPcIncFactor, bool calcKr, bool calcI, 
			  bool entreL, bool entreR, bool exitL, bool exitR, bool swOut)
{

	m_out << "\n********************* Water-injection --  cycle "<<m_comn.floodingCycle()+1<<" *********************"<<endl;

    SortedEvents<Apex*,PceImbCmp>    m_eventsCh;

	//if (useHypre) 
     m_solver = new hypreSolver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_comn.debugMode, "solverImbibe", m_writeSlvMatrixAsMatlab);
    //else 
     //m_solver = new amg_solver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_maxNonZeros, m_comn.debugMode,
		//m_matrixFileName + "_imbcycle_"+to_string(m_comn.floodingCycle()), m_writeSlvMatrixAsMatlab);



    Netsim::initializeImbibition(m_eventsCh, calcKr, calcI, entreL, entreR, exitL, exitR, input);


	if (m_comn.debugMode>100) 
	{	cout<<m_rockLattice[0]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<m_rockLattice[0]->model()->Pc_pistonTypeRec()<<endl;
		cout<<m_rockLattice[m_numPores+1]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<m_rockLattice[m_numPores+1]->model()->Pc_pistonTypeRec()<<endl;
	}


    bool satCompress(false), compressWat(false), compressOil(false);
    double krThreshold, newDeltaSw, dSw(deltaSw);
    input.relPermCompression(satCompress, krThreshold, newDeltaSw, compressWat, compressOil);
    if(satCompress && calcKr && compressWat) dSw = newDeltaSw;
    double SwTarget = min(requestedFinalSw, Sw + dSw*0.5);
    double PcTarget = max(requestedFinalPc, Pc - (deltaPc+abs(Pc)*deltaPcIncFactor)*0.1);
    bool residualSat(false);


    while(/*!residualSat &&*/ Sw <=  requestedFinalSw && Pc > requestedFinalPc)
    {

		Netsim::singleImbibeStep(m_eventsCh, SwTarget, PcTarget, residualSat);

		double krw = m_watFlowRate / (m_singlePhaseWaterQ+1.0e-200);
		double kro = m_oilFlowRate / (m_singlePhaseOilQ+1.0e-200);;
		//double resIdx = m_resistivityIdx;

        Sw = m_satWater;
        Pc = m_cappPress;
        if(satCompress && calcKr && compressWat && krw > krThreshold) dSw = deltaSw;
        if(satCompress && calcKr && compressOil && kro < krThreshold) dSw = newDeltaSw;
        SwTarget = min(requestedFinalSw+1.0e-15, round((m_satWater + 0.75*dSw)/dSw)*dSw);
        PcTarget = max(requestedFinalPc-1.0e-7, Pc - (deltaPc+abs(Pc)*deltaPcIncFactor+1.0e-16));
        //if(swOut) Netsim::recordWaterSatMap();

    }

    Netsim::writeResultData(m_wantRelPerm, m_wantResIdx);

	m_solver=NULL;
    Netsim::finaliseImbibition();
}




/**
 * At the end of draiange imbibition displacement is initialized as max Pc is set. This
 * involves determening entry pressures for all elements.
 */
void Netsim::initializeImbibition(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR, InputData& input)
{
	m_dd = -1;

    invadingFluid = &m_water;
    retreatingFluid = &m_oil;

    m_comn.incrementFloodCycle();
    m_comn.injectant(invadingFluid);
    m_comn.setDrainageCycle(false);

	dbgFile<<std::endl<<std::endl<<"initializeImbibition cycle "<<m_comn.floodingCycle()<<" "<<endl; dbgFile.flush();

    m_cpuTimeTotal = 0.0;
    m_cpuTimeCoal = 0.0;
    m_cpuTimeKrw = 0.0;
    m_cpuTimeKro = 0.0;
    m_cpuTimeTrapping = 0.0;
    m_cpuTimeResIdx = 0.0;
    m_totNumFillings = 0;
    m_minCycleCappPress = m_cappPress;
    m_comn.setMaxPcLastDrainCycle(m_maxCycleCappPress);

    m_maxOilFlowErr = 0.0;
    m_maxWatFlowErr = 0.0;
    m_maxResIdxErr = 0.0;


	m_injAtLeftRes = entreL;
	m_injAtRightRes = entreR;
	
    if(m_sourceNode>0)
    {
        if(m_sourceNode > m_numPores)
        {
            cerr<<"\nError: Source pore ("<< m_sourceNode<<") out of range [1,"<<m_numPores<<"]."<<endl;            exit(-1);
        }
        m_trappingCriteria = escapeToBoth;
        //m_injAtLeftRes = false;
        //m_injAtRightRes = false;
        if(m_injAtLeftRes || m_injAtRightRes)	m_out << "Warning injecting from both inlet/outlet boundary and source node"<<endl;
    }
    else if(!m_injAtLeftRes && !m_injAtRightRes) m_out << "Error: injection boundary is not specified"<<endl;

    if(exitL && exitR)    m_trappingCriteria = escapeToEither;
    else if(exitL)        m_trappingCriteria = escapeToInlet;
    else if(exitR)        m_trappingCriteria = escapeToOutlet;
	else    			  m_out << "Error: exit boundary is not specified "<<exitL << exitR<<endl;
    

    if(m_comn.floodingCycle() == 2) 
    {
		Netsim::applyFWettabilityChange(input);
	}


    if(!m_wantRelPerm) solve_forRelPermResIndex(wantRelPerm, false);     // Create starting points if not previously
    if(!m_wantResIdx) solve_forRelPermResIndex(false, wantResIdx);       // available

    Netsim::recordRes(wantRelPerm, wantResIdx);
    m_wantRelPerm = wantRelPerm;
    m_wantResIdx = wantResIdx;

	///. initialize inlet and outlets
    if(m_injAtRightRes)       ///. Right
    {
        ((InOutBoundaryPore*)m_rockLattice[m_numPores+1])->fillElemCentreWithWaterCreateLayersIO(m_cappPress+0.1);      ///.  inlet BC
	}
    else if(!m_injAtRightRes)
    {
      ((InOutBoundaryPore*)m_rockLattice[m_numPores+1])->fillElemCentreWithOilRemoveLayersIO(m_cappPress-1000000000.1);     ///.  outlet BC
	}

    if(m_injAtLeftRes)    ///. Left
    {
        ((InOutBoundaryPore*)m_rockLattice[0])->fillElemCentreWithWaterCreateLayersIO(m_cappPress+0.1);                 ///.  inlet BC
	}
    else if(!m_injAtLeftRes)
        ((InOutBoundaryPore*)m_rockLattice[0])->fillElemCentreWithOilRemoveLayersIO(m_cappPress-1000000000.1);            ///.  outlet BC

    if(m_sourceNode != 0 && m_rockLattice[m_sourceNode]->model()->bulkFluid() != invadingFluid)
        m_rockLattice[m_sourceNode]->fillElemCentreWithWaterCreateLayers();     ///.  source BC

    //m_eventsCh.clear();
    //m_layerEventsCh.clear();

	dbgFile<<" 2"; dbgFile.flush();

	///. initialize all elements, except inlet and outlet BCs
    for(size_t i = m_numPores+2; i < m_rockLattice.size(); ++i)
    {
        if(m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->initWaterInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());
        }
    }
    for(int i = 1; i <=  m_numPores; ++i)
    {
        if(m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->initWaterInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());
        }
    }
    

	int nInWaterFlood0 = 0;
    for(int i = 1; i < int(m_rockLattice.size()); ++i)
    {
        if(i != m_numPores + 1 && m_rockLattice[i]->connectedToNetwork())
        {
            if(m_rockLattice[i]->canBeAddedToEventVec(&m_water))
            {
                if (!m_rockLattice[i]->model()->bulkFluid()->isOil())  		cout<< " depn" ; 

				m_rockLattice[i]->calcCentreEntryPrsWatInj();
                m_eventsCh.quickInsert(m_rockLattice[i]);
                m_rockLattice[i]->setInWatFloodVec(true);
				if( m_rockLattice[i]->iRockType() == 0 && m_rockLattice[i]->isInWatFloodVec() )		++nInWaterFlood0;
            }
			addElemTo_layerImbibeEvents(m_layerEventsCh, m_rockLattice[i]);
        }
    }

	dbgFile<<"3"; dbgFile.flush();

	///.  sort events to find the highest priority element to start injecting from
    m_eventsCh.sortEvents();
    m_layerEventsCh.sortEvents();



	int nInWaterFlood = 0, count1WF = 0, nNotInWaterFlood = 0, count1NIWF = 0;    
	for(int i = 1; i < int(m_rockLattice.size()); ++i)
    {
        if(i != m_numPores + 1 && m_rockLattice[i]->connectedToNetwork())
        {   register Element* elem = m_rockLattice[i];
            if(elem->iRockType() == 0)
            {
				if(elem->isInWatFloodVec())			nInWaterFlood++;
				else 								nNotInWaterFlood++;
			}
			else  if(elem->isInWatFloodVec())		count1WF++;
				else 								count1NIWF++;
		}
	}
	cout<<"\n nInWaterFlood0 "<<nInWaterFlood0<<endl;
	cout<<" nNotInWaterFlood "<<nNotInWaterFlood<<endl;
	cout<<" nInWaterFlood "<<nInWaterFlood<<endl;
	cout<<" count1WF "<<count1WF<<endl;
	cout<<" count1NIWF "<<count1NIWF<<endl;
	cout<<" n_EventsCh "<<m_eventsCh.size()<<endl;
	cout<<" n_rockLattice "<<m_rockLattice.size()<<endl;



	int WARNING_TOCKECK;
	/*///. connect oil to all elements connected to inlet, ERROR TODO DELETE TOTEST
    //if(m_injAtLeftRes)
    //{
        //for(int inT = 0; inT < m_rockLattice[0]->connectionNum(); ++inT)
        //{
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[0]->connection(inT), filmBlob);
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[0]->connection(inT), bulkBlob);
        //}
    //}
//
    //if(m_injAtRightRes)
    //{
        //for(int outT = 0; outT < m_rockLattice[m_numPores+1]->connectionNum(); ++outT)
        //{
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[m_numPores+1]->connection(outT), filmBlob);
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[m_numPores+1]->connection(outT), bulkBlob);
        //}
    //}
//
    //if(m_sourceNode != 0)
    //{
        //for(int sourceT = 0; sourceT < m_rockLattice[m_sourceNode]->connectionNum(); ++sourceT)
        //{
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[m_sourceNode]->connection(sourceT), filmBlob);
            //Netsim::untrap_WaterGanglia(m_eventsCh, m_rockLattice[m_sourceNode]->connection(sourceT), bulkBlob);
        //}
    //}*/

    m_amottDataImbibition[0] = m_satWater;
    m_amottDataImbibition[1] = -1.0;

    m_usbmDataImbibition.clear();
    recordUSBMData(false);

    if(m_writeImbList)
    {
        ostringstream fileName;
        fileName << "fill_imbcycle_" << m_comn.floodingCycle() << ".m";
        m_imbListOut.open(fileName.str().c_str());
        m_imbListOut << "% The backbone identifies which pores/throats are water filled at the start of water flooding." << endl
            << "% The first row is pore/throat index, followed by 1 for pores and 0 for thoats." << endl;
        m_imbListOut << "backbone = [";
        for(size_t i = 0; i < m_rockLattice.size(); ++i)
        {
            if(!m_rockLattice[i]->isEntryOrExitRes() && !m_rockLattice[i]->model()->containCOil())
            {
                bool isAPore(dynamic_cast< Pore* >(m_rockLattice[i]) != 0);
                m_imbListOut << m_rockLattice[i]->orenIndex() << ", ";
                m_imbListOut << isAPore << "; ..." << endl;
            }
        }
        m_imbListOut << "];" << endl;
        m_imbListOut << endl << "% The filling list identifies the order through which pores/throats get filled by water" << endl;
        m_imbListOut << "fill = [";
    }
    

        Netsim::checkUntrapOilIfUnstableConfigsImb(m_eventsCh);



	dbgFile<<" . ";
	cout<<endl;

}




/**
 * Do a single displacement with water and recalculate saturation and relative permeability.
 * We do not allow for empty filling Vec as this just might create a mess in imbibition. The
 * reason is that as Pc increases (negatively) more and more regions might get trapped as oil
 * layers collapse, and we don't want to spend huge amounts of time checking on this.
 */
void Netsim::singleImbibeStep(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double SwTarget, double& PcTarget, bool& residualSat)
{	

	clock_t startClock(clock());


	if(SwTarget < m_satWater || PcTarget > m_cappPress)
	{
		m_out << "================================================="              << endl
			<< "Nothing to be done:"                                            << endl;
		if(SwTarget < m_satWater)			m_out << "Target water saturation (" << SwTarget << ")   is lower than current saturation (" << m_satWater << ")"   << endl;
		else			m_out << "Target capillary pressure (" << PcTarget << ")    is higher than current pressure (" << m_cappPress << ")"  << endl;
		m_out << "=================================================="             << endl;

		return;
	}


	double lastStepWaterSat = m_satWater;
	
	//m_out << "Sw = "<< setw(6) << m_satWater << " -> ";

    int numSteps(0), totNumFill(0);
    int fillTarget = max(m_minNumFillings,
						int(m_initStepSize*(m_numPores+m_numThroats)*(SwTarget-m_satWater)));
	double oldCappPress(m_cappPress);

    while(m_satWater <= SwTarget && m_cappPress > PcTarget-1.0e-32 /*&& !m_eventsCh.empty()*/)
    {
        double oldWaterSat(m_satWater);
        int numInv(0), invInsideBox(0);
        bool insideBox(false);

		dbgFile<<"\n Sw = "<<m_satWater<<":"<<"  Pc = "<<m_cappPress<<": ";dbgFile.flush();

        while(invInsideBox < fillTarget && !m_eventsCh.empty() && nextCentrInjPc(m_eventsCh) >= PcTarget)
        {
            Netsim::popUpdateWaterInj(m_eventsCh,insideBox, m_cappPressCh, PcTarget);
            m_comn.GuessCappPress(m_cappPress);                   // Use local maxima for Pc
            ++numInv;
            if(insideBox) ++invInsideBox;
            if(m_cappPress==0.0 || m_cappPress*oldCappPress < 0.0)
            {

                Netsim::checkUntrapOilIfUnstableConfigsImb(m_eventsCh);
				oldCappPress=m_cappPress-1.0e-12;
				updateSatAndConductances(m_cappPress);
                m_out << "Pc cross over at Sw = " << setw(6) << m_satWater << "; ";
                m_amottDataImbibition[1] = m_satWater;//Netsim::recordAmottData(false);
                Netsim::solve_forRelPermResIndex(m_wantRelPerm, m_wantResIdx);
                Netsim::recordRes(m_wantRelPerm, m_wantResIdx);
               	m_out << endl;
            }

///. Third inner loop, only when option m_StableFilling, until no layer ready for pop(?)
            if(m_StableFilling)   
             while(nextCentrInjPc(m_eventsCh) >= m_cappPress-1.0e-32)
                {	//cout<<" StableFilling ";
                    Netsim::popUpdateWaterInj(m_eventsCh, insideBox, m_cappPressCh, PcTarget);
                    ++numInv;
                    if(insideBox) ++invInsideBox;
                 }
        }

        if (nextCentrInjPc(m_eventsCh) < PcTarget && (m_cappPress > PcTarget) )	
        {
			m_cappPressCh = PcTarget;
			
			Netsim::updateSatAndConductances(m_cappPress);
			//if (m_satWater < SwTarget)
				//m_cappPressCh = satTrack->predictPc( m_cappPress,m_satWater,PcTarget,SwTarget, m_comn.floodingCycle()-1); ///. Warning hardcopied iCycle

			if(m_cappPress==0.0 || m_cappPress*oldCappPress < 0.0)
			{
				m_out << "Pc cross over at Sw = " << setw(6) << m_satWater << "\n";
				m_amottDataImbibition[1] = m_satWater;//Netsim::recordAmottData(false);
				oldCappPress=m_cappPress-1.0e-12;
			}

		}

        Netsim::checkUntrapOilIfUnstableConfigsImb(m_eventsCh);
		Netsim::updateSatAndConductances(m_cappPress);

		Netsim::recordUSBMData(false);
        fillTarget = max(m_minNumFillings, (int)min((fillTarget*m_maxFillIncrease),
            (m_extrapCutBack*(invInsideBox/(m_satWater-oldWaterSat))*(SwTarget-m_satWater)) )  );
        totNumFill += numInv;
        ++numSteps;
    }

    residualSat = m_eventsCh.empty();
    m_minCycleCappPress = min(m_cappPress, m_minCycleCappPress);

	cout.precision(3);
	m_out << "Sw: " << setw(8) << std::left << m_satWater << " Pc: " << setw(10) << round(m_cappPress)  << " "; 
    m_out << setw(2) << std::right << numSteps << " steps " << setw(5) << totNumFill << " invasions " ;

    Netsim::solve_forRelPermResIndex(m_wantRelPerm, m_wantResIdx);
    m_totNumFillings += totNumFill;

	Netsim::recordRes(m_wantRelPerm, m_wantResIdx);


	m_out << endl;
	m_cpuTimeTotal += cpuTimeElapsed(startClock);
}



/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains water. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void Netsim::popUpdateWaterInj(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, bool& insideBox, double& localPc, double localPcTarget)
{

	while (m_eventsCh.peek()->subIndex()>=0 && !m_eventsCh.empty() && m_eventsCh.peek()->gravCorrectedEntryPress() > localPcTarget)
    {
		Netsim::popUpdate_layerInj_Water(m_eventsCh.pop(), m_eventsCh, localPc);

		clearTrappedOilFromEvents(m_eventsCh); ///. skip trapped regions
		if( m_eventsCh.empty() ) return;
	}

    if(m_eventsCh.empty() ||  m_eventsCh.peek()->gravCorrectedEntryPress() < localPcTarget) return;

	{
		Element *currElemCh = m_eventsCh.pop()->parentModel()->ChParent();
		Netsim::popUpdateCentreInj_Water(currElemCh, m_eventsCh, localPc);
		
		if(m_writeImbList)
		{
			bool isAPore(dynamic_cast< Pore* >(currElemCh) != 0);
			m_imbListOut << currElemCh->orenIndex() << ", ";
			m_imbListOut << isAPore << "; ..." << endl;
		}
		
		clearTrappedOilFromEvents(m_eventsCh);
		insideBox = currElemCh->isInsideSatBox();	
    }
}



/**
 * Since a single element have many possible imbibition entry pressures, a priority Vec
 * is unsuitable for keeping track of the filling sequence. Hence we use a sorted list.
 * For each imbibition event the entry pressures for the neighbours might change since
 * the number of oil filled neighbours have decreased. We need to remove these elements
 * from the list and then reinsert them with the new entry pressure, while keeping the
 * computational cost as low as possible. Just resorting the list would cost O(N) whreas
 * extracting the elements, O(logN), and then reinserting them, O(logN), have a total
 * cost of 2*connectionNum*O(logN). For large networks this should be cheaper.
 */
void Netsim::popUpdateCentreInj_Water(Element* currElemCh, SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double& localPc)
{

    const Element *currElem = currElemCh;
    double currentPc = currElem->model()->gravCorrectedEntryPress();
     if (m_eventsCh.remove(currElemCh)) cout<<" dblInsWInj " ;
	if(!currElemCh->model()->conductCOil())
	{
		currElemCh->setInWatFloodVec(false);
		cout<<"dsfgfg "<<currElemCh->canBeAddedToEventVec(&m_water)<<endl;
		if (m_eventsCh.remove(currElemCh)) cout<<" dblInsElm" ;
		return ;
	} 
	
	if(!m_eventsCh.empty() && currentPc < m_eventsCh.peek()->entryPc()) cout<<"\n\nError: currentPc < m_eventsCh.peek()->entryPc() \n"<<endl;

    localPc = min(localPc, currentPc);
    //m_boundPress = max(m_boundPress, currentPc);
    dbgFile<<currElem->model()->displacementType()<<currElem->iRockType();dbgFile.flush();
 
    bool HadNoInvading(!currElem->model()->conductAnyWater());

    currElemCh->fillElemCentreWithWaterCreateLayers();         ///. popUpdateCentreInj_Water    inject through centre



    if(currElem->isTrappedWat(filmBlob) && !currElem->model()->disConectedCentreWCornerW())
       untrap_WaterGanglia(m_eventsCh, currElemCh, filmBlob);

    for(int j = 0; j < currElem->connectionNum(); ++j)
    {
        if(HadNoInvading) Netsim::untrap_WaterGanglia(m_eventsCh, currElem->connection(j), filmBlob);///. invading fluid coalescence
        untrap_WaterGanglia(m_eventsCh, currElem->connection(j), bulkBlob);
    }

    Netsim::addElemTo_layerImbibeEvents(m_layerEventsCh, currElemCh);///. for film and centre coalescence

    bool disconnectRetreading(!currElem->model()->conductsAnyOil());

    for(int i = 0; i < currElem->connectionNum(); ++i)
    {
        Element *connection = currElem->connection(i);
        if(!connection->isEntryOrExitRes())///. ! in out, check for water trapping after oil injection
        {///. THE RULES TO BE CHECKED
            if(disconnectRetreading)
            {
                 Netsim::findMarkStoreTrappedOil(m_layerEventsCh,connection,localPc);
			}

				insertReCalcImbibeEntryPrs(m_eventsCh, connection,localPc);
				addElemTo_layerImbibeEvents(m_layerEventsCh,connection);
        }
    }

    clearTrappedOilFromEvents(m_eventsCh);

}


/** trapping routine
 */
void Netsim::findMarkStoreTrappedOil(SortedEvents<Apex*,PceImbCmp>& m_layerEventsCh, Element* elem, double localPc)
{
    vector<Element*> trappingStorage;
    double rhogh(elem->model()->rhogh());

    elem->findMarkTrappedOilGanglia(localPc-rhogh, trappingStorage, m_cpuTimeTrapping, m_trappingCriteria);

    if(!trappingStorage.empty())
    {
		dbgFile<<' '<<trappingStorage.size()<<'t'<<'o'<<elem->iRockType()<<" ";dbgFile.flush();
        for(size_t elm = 0; elm < trappingStorage.size(); ++elm)
        {
			if(trappingStorage[elm]->isInWatFloodVec())
			{
				m_layerEventsCh.remove(trappingStorage[elm]);
				trappingStorage[elm]->setInWatFloodVec(false);
			}
			Polygon* polyShape = dynamic_cast< Polygon* >(trappingStorage[elm]->ChModel());
			if(polyShape)
			{
				for(int cor = 0; cor < polyShape->numCorners(); ++cor)
				{
					if(polyShape->oilLayerConst()[cor].isInWatFloodVec())
					{
						m_layerEventsCh.remove(&polyShape->oilLayerCh()[cor]);
						polyShape->oilLayerCh()[cor].setInWatFloodVec(false);
					}
				}

				polyShape->calcOilLayerPc_syncTrappings(localPc-polyShape->rhogh());
				//addElemTo_layerImbibeEvents(m_layerEventsCh,connection);
			}
				//insertReCalcImbibeEntryPrs(m_eventsCh, connection,m_cappPress);
				//trappingStorage[elm]->parentModel()->ChParent()->ChModel()->calcR(localPc);///. candidate to be removed


        }

		//sort(trappingStorage.begin(), trappingStorage.end(), TrappingOilStorageCmp());  // Having it sorted hels us when we want to coalesce the blob
        m_comn.addTrappedRegionOil(trappingStorage);
    }
    ///. trapped elements are left in the events storage to be cleared later just before they pop up
}






/**
 * Done with imbibition  => clean up
 */
void Netsim::finaliseImbibition()
{

    if(m_maxOilFlowErr > 0.1 || m_maxWatFlowErr > 0.1 || m_maxResIdxErr > 0.1)
    {
        m_out << endl
            << "==================================================== "   << endl
            << "Warning: For some rel perm calculations there were"     << endl
            << "more than 10% difference between flowrates computed"    << endl
            << "at the inlet and outlet. Perhaps try to reduce"         << endl
            << "solver tolerance"                                       << endl
            << "Max water flowrate error:    "  << m_maxWatFlowErr      << endl
            << "Max oil flowrate error:      "  << m_maxOilFlowErr      << endl
            << "Max resistivity index error: "  << m_maxResIdxErr       << endl
            << "==================================================== "   << endl
            << endl;
    }

    vector< int > sumfillingEventsP(6);
    vector< int > imbEventT(3);

    for(int elm = 1; elm < static_cast< int >(m_rockLattice.size()); ++elm)
    {
        int fillingEventRecord = m_rockLattice[elm]->fillingEventRecord();

        if(elm <=  m_numPores)
        {
            if(fillingEventRecord == -1) ++sumfillingEventsP[0];
            else if(fillingEventRecord == 0 || fillingEventRecord == 1) ++sumfillingEventsP[1];
            else if(fillingEventRecord == 2) ++sumfillingEventsP[2];
            else if(fillingEventRecord == 3) ++sumfillingEventsP[3];
            else if(fillingEventRecord > 3) ++sumfillingEventsP[4];
            else ++sumfillingEventsP[5];
        }
        else if (elm > m_numPores+1)
        {
            if(fillingEventRecord == -1) ++imbEventT[0];
            else if(fillingEventRecord == 0 || fillingEventRecord == 1) ++imbEventT[1];
            else ++imbEventT[2];
        }
    }

m_out << endl
        << "===================== Imbibition  ===================== "                                   << endl
        << "Total elapsed time for imbibition:       "  << m_cpuTimeTotal                               << endl
        << "Solving for water relative permeability: "  << m_cpuTimeKrw                                 << endl
        << "Solving for oil relative permeability:   "  << m_cpuTimeKro                                 << endl
        << "Solving for Resistivity Index:           "  << m_cpuTimeResIdx                              << endl
        << "Identifying trapped elements:            "  << m_cpuTimeTrapping                            << endl
        << "Coalesceing trapped water:               "  << m_cpuTimeCoal                                << endl
        << "Max water flowrate error:                "  << m_maxWatFlowErr*100.0 << " %"                << endl
        << "Max oil flowrate error:                  "  << m_maxOilFlowErr*100.0 << " %"                << endl
        << "Max resistivity index error:             "  << m_maxResIdxErr*100.0 << " %"                 << endl
        << endl
        << "===================== Network State  ===================== "                                << endl
        << "Minimum capillary pressure reached (Pa): "  << m_cappPress                                  << endl
        << "Water saturation:                        "  << m_satWater                                   << endl
        << "Number of elements invaded:              "  << m_totNumFillings                             << endl
        << "Number of trapped water regions:         "  << (int)m_comn.numTrappedWatRegions()    << endl
        << "Total number of trapped water elements:  "  << (int)m_comn.numTrappedWatElems()      << endl
        << "Number of trapped oil regions:           "  << (int)m_comn.numTrappedOilRegions()    << endl
        << "Total number of trapped oil elements:    "  << (int)m_comn.numTrappedOilElems()      << endl
        << endl
        << "==================   Pore Filling Process  ================="                                << endl
        << "Uninvaded:                               "  << sumfillingEventsP[5]                                 << endl
        << "Snap off:                                "  << sumfillingEventsP[0]                                 << endl
        << "Piston type displacement:                "  << sumfillingEventsP[1]                                 << endl
        << "Pore body filling, I2:                   "  << sumfillingEventsP[2]                                 << endl
        << "Pore body filling, I3:                   "  << sumfillingEventsP[3]                                 << endl
        << "Pore body filling, I4+:                  "  << sumfillingEventsP[4]                                 << endl
        << endl
        << "=================Throat Filling Process  ================="                               << endl
        << "Uninvaded:                               "  << imbEventT[2]                                 << endl
        << "Snap off:                                "  << imbEventT[0]                                 << endl
        << "Piston type displacement:                "  << imbEventT[1]                                 << endl
        << endl;

    m_amottDataImbibition[2] = m_satWater;
    if(m_cappPress > 0.0)
        m_amottWaterIdx = 1.0;
    else if(m_cappPress < 0.0 && m_amottDataImbibition[1] < 0.0)
        m_amottWaterIdx = 0.0;
    else
    {
        m_amottWaterIdx = (m_amottDataImbibition[1]-m_amottDataImbibition[0])/
            (m_amottDataImbibition[2]-m_amottDataImbibition[0]);
    }

    if(m_comn.floodingCycle() > 1)
    {
        m_out << endl
            << "=================Wettability State  ================ "                     << endl
            << "Amott water index, Iw:                  " << m_amottWaterIdx                << endl
            << "Amott oil index, Io:                    " << m_amottOilIdx                  << endl
            << "Amott wettability index, I = Iw - Io:   " << m_amottWaterIdx-m_amottOilIdx  << endl
            << "USBM wettability index:                 " << calcUSBMindex()                << endl
            << endl;
    }


    if(m_writeImbList)
    {
        m_imbListOut << "];" << endl;
        m_imbListOut.close();
    }

    //if(m_prtPressureProfile && m_wantRelPerm)
    //{
        //writePrsProfileData(invadingFluid, m_floodCycleResultsOut); 
        //writePrsProfileData(retreatingFluid, m_floodCycleResultsOut);
    //}


    for(int i = 1; i < int(m_rockLattice.size()); ++i)
    {
        if(i != m_numPores + 1 && m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->finitWaterInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());        
        }
    }


	m_out<<"\n:/"<<endl<<endl;

}











void Netsim::popUpdate_layerInj_Water(Apex*apex,  SortedEvents<Apex*,PceImbCmp>& m_eventsCh, double & localPc)
{///. oil layer collapse

		dbgFile<<'L'<<'w';dbgFile.flush();

		Polygon* polyCh = (Polygon*)apex->parentModel();
        const int & cornerIndex = apex->subIndex();
		const Polygon* poly = polyCh;

		softAssert(!poly->eleman()->isTrappedOil());


	   double currentPc = polyCh->Polygon::Pc_pin_disconnectOilLayer(cornerIndex)+poly->rhogh();
	   localPc = min(currentPc, localPc);
		softAssert(m_eventsCh.empty() || currentPc >= m_eventsCh.peek()->entryPc());
 
		if(poly->numLayers() == poly->numCorners()-1)
		{  // Oblique corner  => water in center and corner is joined, check for coalescing
            Netsim::untrap_WaterGanglia(m_eventsCh, poly->ChParent(), filmBlob);
            Netsim::untrap_WaterGanglia(m_eventsCh, poly->ChParent(), bulkBlob);
        }
		else if(poly->numLayers()==0)//cornerIndex == 0)   ///. last element in W-W, first element in O-W OInj 
		{           //Sharpest corner  => no more oil, check for oil trapping
			if (poly->numLayers()==0) cout<<"dsknb";
			for(int i = 0; i < poly->eleman()->connectionNum(); ++i)
			{
				Netsim::findMarkStoreTrappedOil(m_layerEventsCh, poly->eleman()->connection(i),localPc);
			}
		}
}



void Netsim::untrap_WaterGanglia(SortedEvents<Apex*,PceImbCmp>& m_eventsCh, Element* elem, FluidBlob blob)
{
    pair<int, double> trapWatInfo = elem->trappingWat(blob);

    if(trapWatInfo.first > -1)                                                      // We need to untrap this water blob
    {
        clock_t startCoalesce(clock());
        const vector< pair<Element*,FluidBlob> >& newElems = m_comn.trappedRegionsWat(trapWatInfo.first);
        trapWatInfo.second += elem->model()->rhogh();
		double localPc = trapWatInfo.second;

		dbgFile<<'['<<newElems.size()<<'u'<<'w';dbgFile.flush();

			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				if(newElems[elm].first->isInWatFloodVec())
				{
					if (!m_eventsCh.remove(newElems[elm].first)) cout<<" sykr1 ";cout.flush();
					newElems[elm].first->setInWatFloodVec(false);
				}
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[elm].first->ChModel());
				if(polyShape)
				{
					for(int i = 0; i < polyShape->numCorners(); ++i)
					{
						if(polyShape->oilLayerConst()[i].isInWatFloodVec())
						{
							if (!m_layerEventsCh.remove(&polyShape->oilLayerCh()[i])) cout<<" gonw ";
							polyShape->oilLayerCh()[i].setInWatFloodVec(false);
						}
					}
					softAssert(!newElems[elm].first->isInWatFloodVec());
					softAssert(!polyShape->oilLayerConst()[0].isInWatFloodVec());
				}
				///. start from the local ganglia Pc and gradually equlibriate
				//newElems[elm].first->ChModel()->finitWaterInjection(max(localPc, mcappPress)-newElems[elm].first->ChModel()->rhogh());
				newElems[elm].first->unTrapWat(newElems[elm].second);
				if(polyShape)
				{
					polyShape->calcOilLayerPc_markUntrappedFilms(localPc-polyShape->rhogh());
				}
			}


		if(localPc < m_cappPress) ///. lower curvature, we need to inject oil from centre to equilibrate
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 
			
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm].first->ChModel()->finitWaterInjection(localPc-newElems[elm].first->ChModel()->rhogh());
				newElems[elm].first->ChModel()->initOilInjection(localPc-newElems[elm].first->ChModel()->rhogh());
			}

			SortedEvents< Apex*, PceDrainCmp > oilFillingEvents;

			for(size_t i = 0; i < newElems.size(); ++i)               // First stage: check for unstable configurations
			{                                                       // when the pressure in the coalesced blob is
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[i].first->ChModel());  // equilibrated with the rest
				if(polyShape)
				{
					//Netsim::reCalcWaterFilmPc_markUnrappedFilms(oilFillingEvents, newElems[i].first, m_cappPress);
					//polyShape->calcOilLayerPc _markUntrappedFilms(localPc-polyShape->rhogh());
					polyShape->insertOilSnapEvent_IfSnapPcLgPc(oilFillingEvents, m_cappPress-polyShape->rhogh());       //  => insert oil snap off
					if(polyShape->hasOilLayer_TrappedOutside_PcHsnapPc(m_cappPress))
					{
						dbgFile<<'D';dbgFile.flush();
						cout<<"D";cout.flush();
					}
				}
			}

			if (!oilFillingEvents.empty())   
			{   //do the events 
				oilFillingEvents.sortEvents();
				//Netsim::increasePressureInCoalescedBlobTrue(oilFillingEvents,localPc, m_cappPress);  
				while(!oilFillingEvents.empty() &&  oilFillingEvents.peek()->gravCorrectedEntryPress() < m_cappPress)
				{
					dbgFile<<'^';dbgFile.flush();

					bool tmp;
					popUpdateOilInj(oilFillingEvents, tmp, localPc, m_cappPress-0.000001);
				}

			}
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
					newElems[elm].first->ChModel()->finitOilInjection(m_cappPress-newElems[elm].first->ChModel()->rhogh());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 

			//double highestPc=max(trapWatInfo.second,m_cappPress);
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm].first->ChModel()->initWaterInjection(localPc-newElems[elm].first->ChModel()->rhogh());
			}

		}




        for(size_t elm = 0; elm < newElems.size(); ++elm)                         // Final stage: Now that all is peach, we can
        {                                                                   // possibly start adding new elements to the
			//newElems[elm].first->calcCentreEntryPrsWatInj();
            // if(!newElems[elm]->isTrappedOil())	                           // displacement vectors
            {
                insertReCalcImbibeEntryPrs(m_eventsCh, newElems[elm].first,localPc);
                addElemTo_layerImbibeEvents(m_layerEventsCh, newElems[elm].first);
                for(int k = 0; k < newElems[elm].first->connectionNum(); ++k)
                {
					insertReCalcImbibeEntryPrs(m_eventsCh, newElems[elm].first->connection(k),localPc);
					addElemTo_layerImbibeEvents(m_layerEventsCh, newElems[elm].first->connection(k));
				}
			}
		}

		m_comn.removeTrappedRegionWat(trapWatInfo.first);
		m_cpuTimeCoal += cpuTimeElapsed(startCoalesce);
		
		
		while(!m_eventsCh.empty() && m_eventsCh.peek()->gravCorrectedEntryPress() >= m_cappPress)
		{
			bool tmp;
			popUpdateWaterInj(m_eventsCh, tmp, localPc, m_cappPress);
		}
	dbgFile<<']';dbgFile.flush();
    }
}




///. never does anything TO DELETE
void Netsim::checkUntrapOilIfUnstableConfigsImb(SortedEvents<Apex*,PceImbCmp>& m_eventsCh)
{
    dbgFile<<"[";dbgFile.flush();
    for(size_t i = 0; i < m_rockLattice.size(); ++i)
    {
        Element* elem = m_rockLattice[i];
        double rhogh(elem->model()->rhogh());
        if(elem->model()->waterLayer_UntrappedCorner_PcLsnapPc(m_cappPress-rhogh))
        {
			dbgFile<<'I';dbgFile.flush();
			cout<<'I';cout.flush();

           if(elem->isInWatFloodVec())
            {
                if (m_eventsCh.remove(elem)) elem->setInWatFloodVec(false);
                else  cout<<" clei "; cout.flush();
            }
			softAssert(elem->model()->bulkFluid()->isOil());
			
            double capillaryPressure=1000000000.0;
		   Netsim::popUpdateCentreInj_Water(elem, m_eventsCh, capillaryPressure);

		}
    }
    dbgFile<<"]"<<"  ";
}

 



