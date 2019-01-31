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


#include "netsim.h"






/**
 * Continue injecting oil until a target saturation is reached, draiange queue is empty
 * or capillary pressure target is reached. Writing results to output file. Oil injection
 * needs to occur at least at one face.
 */
void Netsim::Drainage(InputData& input, double& Sw, double& Pc, double requestedFinalSw, double requestedFinalPc,
              double deltaSw, double deltaPc, double deltaPcIncFactor, bool calcKr, bool calcI, 
			  bool entreL, bool entreR, bool exitL, bool exitR, bool swOut)
{

	m_out << "\n********************* Oil-injection -- cycle "<<m_comn.floodingCycle()+1<<" *********************"<<endl;

    SortedEvents<Apex*,PceDrainCmp>    m_eventsCh;
    //SortedEvents<Apex*, PceDrainCmp>&   m_layerEventsCh(m_eventsCh);


    //if (useHypre) 
     m_solver = new hypreSolver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_comn.debugMode, "solverDrain", m_writeSlvMatrixAsMatlab);
    //else 
     //m_solver = new amg_solver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_maxNonZeros, m_comn.debugMode,
		//m_matrixFileName + "_draincycle_"+to_string(m_comn.floodingCycle()), m_writeSlvMatrixAsMatlab);


    Netsim::initializeDrainage(m_eventsCh,calcKr, calcI, entreL, entreR, exitL, exitR);


	if (m_comn.debugMode>100) 
	{	cout<<m_rockLattice[0]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<m_rockLattice[0]->model()->Pc_pistonTypeRec()<<endl;
		cout<<m_rockLattice[m_numPores+1]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<m_rockLattice[m_numPores+1]->model()->Pc_pistonTypeRec()<<endl;
	}


    bool satCompress(false), compressWat(false), compressOil(false);
    double krThreshold, newDeltaSw, dSw(deltaSw);
    input.relPermCompression(satCompress, krThreshold, newDeltaSw, compressWat, compressOil);
    if(satCompress && calcKr && compressOil) dSw = newDeltaSw;
    double SwTarget = max(requestedFinalSw, Sw - dSw*0.5);
    double PcTarget = min(requestedFinalPc, Pc + (deltaPc+abs(Pc)*deltaPcIncFactor)*0.1);
    bool residualSat(false);


    while(/*!residualSat &&*/ Sw >= requestedFinalSw && Pc < requestedFinalPc)
    {

		Netsim::singleDrainStep(m_eventsCh, SwTarget, PcTarget, residualSat);

		double krw = m_watFlowRate / (m_singlePhaseWaterQ+1.0e-200);
		double kro = m_oilFlowRate / (m_singlePhaseOilQ+1.0e-200);;
		//double resIdx = m_resistivityIdx;
        Sw = m_satWater;
        Pc = m_cappPress;
        if(satCompress && calcKr && compressOil && kro > krThreshold) dSw = deltaSw;
        if(satCompress && calcKr && compressWat && krw < krThreshold) dSw = newDeltaSw;
        SwTarget = max(requestedFinalSw-1.0e-15, round((m_satWater - 0.75*dSw)/dSw)*dSw);
        PcTarget = min(requestedFinalPc+1.0e-7, m_cappPress + (deltaPc+abs(Pc)*deltaPcIncFactor));
        //if(swOut) Netsim::recordWaterSatMap();

    }

    Netsim::writeResultData(m_wantRelPerm, m_wantResIdx);

	delete m_solver;
	m_solver=NULL;
    Netsim::finaliseDrainage();
}



/**
 * Do a single displacement with oil and recalculate saturation and relative permeability
 */
void Netsim::singleDrainStep(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double SwTarget, double PcTarget, bool& residualSat)
{	

	clock_t startClock(clock());


	if(SwTarget >= m_satWater || PcTarget <= m_cappPress)
	{
		m_out << "================================================="              << endl
			<< "Nothing to be done:"                                            << endl;
		if(SwTarget > m_satWater)			m_out << "Target water saturation (" << SwTarget << ")   is higher than current saturation (" << m_satWater << ")"   << endl;
		else			m_out << "Target capillary pressure (" << PcTarget << ")   is lower than current pressure (" << m_cappPress << ")"  << endl;
		m_out << "=================================================="             << endl;

		return;
	}


	///HASAN  add for loop for time stepping
	
    int numSteps(0), totNumFill(0);
    int fillTarget = max(m_minNumFillings,
						int(m_initStepSize*(m_numPores+m_numThroats)*(SwTarget-m_satWater)));
    double oldCappPress(m_cappPress);

    while(m_satWater > SwTarget && m_cappPress < PcTarget+1.0e-32 /*&& !m_eventsCh.empty()*/) ///. outer filling loop
    {
        double oldWaterSat(m_satWater);
        int numInv(0), invInsideBox(0);
        bool insideBox(false);

		dbgFile<<"\n Sw = "<<m_satWater<<":"<<"  Pc = "<<m_cappPress<<"   "<<"  "<<nextCentrInjPc(m_eventsCh)<<": ";dbgFile.flush();

        while((invInsideBox < fillTarget && !m_eventsCh.empty() && nextCentrInjPc(m_eventsCh) <= PcTarget) )  ///. inner filling loop
        {
            Netsim::popUpdateOilInj(m_eventsCh,insideBox, m_cappPressCh, PcTarget); ///.  popUpdateCentreInj _Oil
            m_comn.GuessCappPress(m_cappPress); /// used as initialGuess
            ++numInv;
            if(insideBox) ++invInsideBox;
            if(m_cappPress==0.0 || m_cappPress*oldCappPress < 0.0)
            {

                Netsim::checkUntrapWaterIfUnstableConfigsDrain(m_eventsCh);

				updateSatAndConductances(m_cappPress);
                m_out << "Pc cross over at Sw = " << setw(6) << m_satWater << "; ";
                m_amottDataDrainage[1] = m_satWater;//Netsim::recordAmottData(true);
                Netsim::solve_forRelPermResIndex(m_wantRelPerm, m_wantResIdx);
                Netsim::recordRes(m_wantRelPerm, m_wantResIdx);
		        oldCappPress = m_cappPress+1.0e-12;
               	m_out << endl;
            }


///. Third inner loop, only when option m_StableFilling, until no layer ready for pop(?)
            if(m_StableFilling)   
             while(nextCentrInjPc(m_eventsCh) <= m_cappPress+1.0e-32)
                {	//cout<<" StableFilling ";
                    Netsim::popUpdateOilInj(m_eventsCh, insideBox, m_cappPressCh, PcTarget);
                    ++numInv;
                    if(insideBox) ++invInsideBox;
                 }
        }

        if (nextCentrInjPc(m_eventsCh) > PcTarget && (m_cappPress < PcTarget) )	
        {
			//cout<<"  * "<<m_cappPress<<":"<<satTrack->predictPc( m_cappPress,m_satWater,PcTarget,SwTarget, 0)<<":"<<PcTarget<<" * ";
			m_cappPressCh = PcTarget;
			//m_cappPressCh = min(PcTarget,m_cappPress+((1.0+double(numSteps))/50.0)*(PcTarget-m_cappPress));/// with additional underrelaxation to  check for saturation
			
			Netsim::updateSatAndConductances(m_cappPress);
			//if (m_satWater > SwTarget)
				//m_cappPressCh = satTrack->predictPc( m_cappPress,m_satWater,PcTarget,SwTarget, m_comn.floodingCycle()-1);
				//m_cappPressCh = min(PcTarget,m_cappPress+
				//abs(oldCappPress-m_cappPress)/(abs(m_satWater-oldWaterSat)+1.0e-6)*(oldWaterSat-SwTarget+0.1*(lastStepWaterSat-SwTarget)+1.1e-6)*0.1
				//);
			if(m_cappPress==0.0 || m_cappPress*oldCappPress < 0.0)
			{				
				m_out << "Pc cross over at Sw = " << setw(6) << m_satWater << "\n";
				m_amottDataDrainage[1] = m_satWater;//Netsim::recordAmottData(true);
				oldCappPress = m_cappPress+1.0e-12;
			}


		}


        Netsim::checkUntrapWaterIfUnstableConfigsDrain(m_eventsCh);
		Netsim::updateSatAndConductances(m_cappPress);

		Netsim::recordUSBMData(true);
        fillTarget = max(m_minNumFillings, (int)min((fillTarget*m_maxFillIncrease),
            (m_extrapCutBack*(invInsideBox/(m_satWater-oldWaterSat))*(SwTarget-m_satWater)) )  );
        totNumFill += numInv;
        ++numSteps;
    }

    residualSat = m_eventsCh.empty();
    m_maxCycleCappPress = max(m_cappPress, m_maxCycleCappPress);

	cout.precision(3);
	m_out << "Sw: " << setw(8) << std::left << m_satWater << " Pc: " << setw(10) << round(m_cappPress)  << " "; 
    m_out << setw(2) << std::right << numSteps << ":" << setw(5) << totNumFill << " invasions " ;

    Netsim::solve_forRelPermResIndex(m_wantRelPerm, m_wantResIdx);
    m_totNumFillings += totNumFill;

	Netsim::recordRes(m_wantRelPerm, m_wantResIdx);



	if (m_comn.debugMode>0 && Element::nErrs>0) 
		cout<<"   nErrs:"<<Element::nErrs<<"  ";
	Element::nErrs = 0;

	m_out << endl;
	m_cpuTimeTotal += cpuTimeElapsed(startClock);
}



inline bool Netsim::insertReCalcDrainEntryPrs(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem, double cappPrs)
{
	if(elem->isInOilFloodVec())
	{
		if (m_eventsCh.remove(elem))
		{   ///. order is important
		
				elem->calcCentreEntryPrsOilInj();///. order is important
				m_eventsCh.insert(elem);
                if (elem->model()->bulkFluid()->isOil()) cout<< " deab " ;
		}
		else
		{
			elem->setInOilFloodVec(false);
			elem->calcCentreEntryPrsOilInj();///. order is important
			if (elem->canBeAddedToEventVec(&m_oil))
			{
				m_eventsCh.insert(elem);	elem->setInOilFloodVec(true);
				if (elem->model()->bulkFluid()->isOil()) cout<< " dwpb " ;
				//cout<<"|";
			}
			else
			{
					cout<<"~";
			}
			//elem->Ch S hape()->calc R(cappPrs);///. candidate to be removed
			return false;
		}
	}
	else
	{
		if (elem->canBeAddedToEventVec(&m_oil))
		{
			elem->ChModel()->initOilInjection(cappPrs);
			elem->calcCentreEntryPrsOilInj();
			m_eventsCh.insert(elem);
			elem->setInOilFloodVec(true);
		}
	}
	//elem->ChS hape()->ca lcR(cappPrs);///. candidate to be removed
	return true;
} 

template<class compType>
inline void Netsim::reCalcOilLayerPc_markTrappedFilms(SortedEvents<Apex*,compType>& m_layerEventsCh, Element* elem, double cappPrs)
{///. pc is just used for error handling
    Polygon* polyShape = dynamic_cast< Polygon* >(elem->ChModel());
    if(polyShape)
    {
        for(int i = 0; i < polyShape->numCorners(); ++i)
        {
            if(polyShape->oilLayerConst()[i].isInOilFloodVec())
            {
                if (!m_layerEventsCh.remove(&polyShape->oilLayerCh()[i])) cout<<" 454545886 ";
                //m_layerEventsCh->setInOilFloodVec(false);
            }
        }
 
        polyShape->calcOilLayerPc_syncTrappings(cappPrs-polyShape->rhogh());
 
        for(int j = 0; j < polyShape->numCorners(); ++j)
        {
            if(polyShape->oilLayerConst()[j].isInOilFloodVec())
			{
				if (polyShape->oilLayerConst()[j].trappedOLayer().first<0)
					m_layerEventsCh.insert(&polyShape->oilLayerCh()[j]);
				else
				{
					polyShape->oilLayerCh()[j].setInOilFloodVec(false);
					cout<<" nlsua ";
				}
					
			}
        }
    }
}


inline void Netsim::addElemTo_layerDrainEvents(SortedEvents<Apex*,PceDrainCmp>& m_layerEventsCh, Element* elem)
{
    Polygon* polyShape = dynamic_cast< Polygon* >(elem->ChModel());
    if(polyShape)
    {
		vector<int> addCrns;
		double rhogh(elem->model()->rhogh());
		if(elem->addToLayerVec(invadingFluid, polyShape, addCrns))
		{
			for(size_t i = 0; i < addCrns.size(); ++i)
			{
				softAssert(!polyShape->oilLayerConst()[addCrns[i]].isInOilFloodVec());
				m_layerEventsCh.insert(&polyShape->oilLayerCh()[addCrns[i]]);
				polyShape->oilLayerCh()[addCrns[i]].setInOilFloodVec(true);
			}
		}
	}
}



inline void Netsim::clearTrappedWatFromEvents(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh)
{
    while(!m_eventsCh.empty() && m_eventsCh.peek()->parentModel()->eleman()->isTrappedWat(bulkBlob) )
    {
        m_eventsCh.pop()->setInOilFloodVec(false);
    }

    while(!m_layerEventsCh.empty() && m_layerEventsCh.peek()->trappingCL().first>-1)
    {
        m_layerEventsCh.pop()->setInOilFloodVec(false);
    }
}





/**
 * Sets up the solver, fills the sorted lists etc
 */
void Netsim::initializeDrainage(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, bool wantRelPerm, bool wantResIdx, bool entreL, bool entreR, bool exitL, bool exitR)
{
	m_dd = 1;

    invadingFluid = &m_oil;
    retreatingFluid = &m_water;

    m_comn.incrementFloodCycle();
    m_comn.injectant(invadingFluid);
    m_comn.setDrainageCycle(true);

	dbgFile<<std::endl<<std::endl<<"initializeDrainage cycle "<<m_comn.floodingCycle()<<" "<<endl;

    m_cpuTimeTotal = 0.0;
    m_cpuTimeCoal = 0.0;
    m_cpuTimeKrw = 0.0;
    m_cpuTimeKro = 0.0;
    m_cpuTimeTrapping = 0.0;
    m_cpuTimeResIdx = 0.0;
    m_totNumFillings = 0;
    m_maxCycleCappPress = m_cappPress;
    m_comn.setMinPcLastCycle(m_minCycleCappPress);

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
        if(m_injAtLeftRes || m_injAtRightRes)	m_out << "Injecting from inlet/outlet and source node "<<m_sourceNode<<endl;
        else									m_out << "Injecting from node "<<m_sourceNode<<endl;
    }
    else if(!m_injAtLeftRes && !m_injAtRightRes) m_out << "Error: injection boundary is not specified"<<endl;

    if(exitL && exitR)    m_trappingCriteria = escapeToEither;
    else if(exitL)        m_trappingCriteria = escapeToInlet;
    else if(exitR)        m_trappingCriteria = escapeToOutlet;
	else    			  m_out << "Error: exit boundary is not specified "<<exitL << exitR<<endl;
    





    if(!m_wantRelPerm) solve_forRelPermResIndex(wantRelPerm, false);     // Create starting points if not previously
    if(!m_wantResIdx) solve_forRelPermResIndex(false, wantResIdx);       // available

    Netsim::recordRes(wantRelPerm, wantResIdx);
    m_wantRelPerm = wantRelPerm;
    m_wantResIdx = wantResIdx;

	///. initialize inlet and outlets
    if(m_injAtRightRes)       ///. Right
    {
        ((InOutBoundaryPore*)m_rockLattice[m_numPores+1])->fillElemCentreWithOilRemoveLayersIO(m_cappPress-0.1);        ///.  inlet BC
	}
    else if(!m_injAtRightRes)
    {
		((InOutBoundaryPore*)m_rockLattice[m_numPores+1])->fillElemCentreWithWaterCreateLayersIO(m_cappPress+1000000000.1);             ///.  outlet BC
	}

    if(m_injAtLeftRes)    ///. Left
    {
        ((InOutBoundaryPore*)m_rockLattice[0])->fillElemCentreWithOilRemoveLayersIO(m_cappPress-0.1);            ///.  inlet BC
	}
    else if(!m_injAtLeftRes)
        ((InOutBoundaryPore*)m_rockLattice[0])->fillElemCentreWithWaterCreateLayersIO(m_cappPress+1000000000.1);                ///.  outlet BC

    if(m_sourceNode != 0 && m_rockLattice[m_sourceNode]->model()->bulkFluid() != invadingFluid)
        m_rockLattice[m_sourceNode]->fillElemCentreWithOilRemoveLayers();   ///.  source node

    //m_eventsCh.clear();
    //m_layerEventsCh.clear();

	dbgFile<<" 2"; dbgFile.flush();

	///. initialize all elements, except inlet and outlet BCs
    for(size_t i = m_numPores+2; i < m_rockLattice.size(); ++i)
    {
        if(m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->initOilInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());
        }
    }
    for(int i = 1; i <=  m_numPores; ++i)
    {
        if(m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->initOilInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());
        }
    }
    

    for(int i = 1; i < int(m_rockLattice.size()); ++i)
    {
        if(i != m_numPores + 1 && m_rockLattice[i]->connectedToNetwork())
        {
            if(m_rockLattice[i]->canBeAddedToEventVec(&m_oil))
            {
                if (m_rockLattice[i]->model()->bulkFluid()->isOil())  		cout<< " depb " ; 

				m_rockLattice[i]->calcCentreEntryPrsOilInj();
               m_eventsCh.quickInsert(m_rockLattice[i]);
                m_rockLattice[i]->setInOilFloodVec(true);

            }
			addElemTo_layerDrainEvents(m_layerEventsCh, m_rockLattice[i]);
        }
    }

	dbgFile<<"3"; dbgFile.flush();

	///.  sort events to find the highest priority element to start injecting from
    m_eventsCh.sortEvents();
    m_layerEventsCh.sortEvents();


	/*/////. connect oil to all elements connected to inlet, ERROR TODO DELETE TOTEST
    //if(m_injAtLeftRes)
    //{
        //for(int inT = 0; inT < m_rockLattice[0]->connectionNum(); ++inT)
        //{
            //Netsim::untrap_OilGanglia(m_eventsCh, m_rockLattice[0]->connection(inT));
        //}
    //}
//
    //if(m_injAtRightRes)
    //{
        //for(int outT = 0; outT < m_rockLattice[m_numPores+1]->connectionNum(); ++outT)
        //{
            //Netsim::untrap_OilGanglia(m_eventsCh, m_rockLattice[m_numPores+1]->connection(outT));
        //}
    //}
//
    //if(m_sourceNode != 0)
    //{
        //for(int sourceT = 0; sourceT < m_rockLattice[m_sourceNode]->connectionNum(); ++sourceT)
        //{
            //Netsim::untrap_OilGanglia(m_eventsCh, m_rockLattice[m_sourceNode]->connection(sourceT));
        //}
    //}*/

    m_amottDataDrainage[0] = m_satWater;
    m_amottDataDrainage[1] = -1.0;

    m_usbmDataDrainage.clear();
    recordUSBMData(true);

    if(m_writeDrainList)
    {
        ostringstream fileName;
        fileName << m_baseFileName << "fill_cycle" << m_comn.floodingCycle() << "_drain.m";
        m_drainListOut.open(fileName.str().c_str());
        m_drainListOut << "% The backbone identifies which pores/throats are oil filled at the start of drainage." << endl
            << "% The first row is pore/throat index, followed by 1 for pores and 0 for thoats." << endl;
        m_drainListOut << "backbone = [";
        for(size_t i = 0; i < m_rockLattice.size(); ++i)
        {
            if(!m_rockLattice[i]->isEntryOrExitRes() && m_rockLattice[i]->model()->containCOil())
            {
                bool isAPore(dynamic_cast< Pore* >(m_rockLattice[i]) != 0);
                m_drainListOut << m_rockLattice[i]->orenIndex() << ", ";
                m_drainListOut << isAPore << "; ..." << endl;
            }
        }
        m_drainListOut << "];" << endl;
        m_drainListOut << endl << "% The filling list identifies the order through which pores/throats get filled by oil" << endl;
        m_drainListOut << "fill = [";
    }
    

        Netsim::checkUntrapWaterIfUnstableConfigsDrain(m_eventsCh);



	dbgFile<<" . ";
	cout<<endl;

}




/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains oil. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void Netsim::popUpdateOilInj(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, bool& insideBox, double& localPc, double localPcTarget)
{

	while (m_eventsCh.peek()->subIndex()>=0 && !m_eventsCh.empty() && m_eventsCh.peek()->gravCorrectedEntryPress() < localPcTarget)
    {
		Netsim::popUpdate_layerInjs_Oil(m_eventsCh.pop(), m_eventsCh, localPc);

		clearTrappedWatFromEvents(m_eventsCh); ///. skip trapped regions
		if( m_eventsCh.empty() ) return;
	}

    if(m_eventsCh.empty() ||  m_eventsCh.peek()->gravCorrectedEntryPress() > localPcTarget) return;

	{
		Element *currElemCh = m_eventsCh.pop()->parentModel()->ChParent();
		Netsim::popUpdateCentreInj_Oil(currElemCh, m_eventsCh, localPc);
		if(m_writeDrainList)
		{
			bool isAPore(dynamic_cast< Pore* >(currElemCh) != 0);
			m_drainListOut << currElemCh->orenIndex() << ", ";
			m_drainListOut << isAPore << "; ..." << endl;
		}
		
		clearTrappedWatFromEvents(m_eventsCh);
		insideBox = currElemCh->isInsideSatBox();	
    }
}



/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains oil. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void Netsim::popUpdateCentreInj_Oil(Element *currElemCh, SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double& localPc)
{

    const Element *currElem = currElemCh;
    double currentPc = currElem->model()->gravCorrectedEntryPress();
     if(currElem->model()->bulkFluid()->isOil()) 
     {		 cout<<"\n bypd "<<currElem->latticeIndex()<<"  "<<currElem->isInOilFloodVec()<<"  "<<currElem->model()->displacementType()<<endl;		 return; //exit(-1);
	 }

    localPc = max(localPc, currentPc);
    //m_boundPress = min(m_boundPress, currentPc);
    dbgFile<<currElem->model()->displacementType()<<currElem->iRockType();dbgFile.flush();
 
    bool HadNoOil(!currElem->model()->conductsAnyOil());//any oil

    currElemCh->fillElemCentreWithOilRemoveLayers(); 
    if (m_eventsCh.remove(currElemCh)) cout<<" dblInsOInj ";     ///. popUpdateCentreInj _Oil  inject through centre
    if (m_eventsCh.checkIfThere(currElemCh)) cout<<" dblInsOIndsdaadsj ";     ///. popUpdateCentreInj _Oil  inject through centre
	if (currElem->isInOilFloodVec())  cout<<" CWVN "; 


    if(HadNoOil) 
    {
		for(int j = 0; j < currElem->connectionNum(); ++j)
		{
			Netsim::untrap_OilGanglia(m_eventsCh, currElem->connection(j)); ///. invading fluid coalescence
		}
	}

    bool disconnectRetreading(!currElem->model()->conductAnyWater());//any water

    for(int i = 0; i < currElem->connectionNum(); ++i)
    {
        Element *connection = currElem->connection(i);
        if(!connection->isEntryOrExitRes())///. ! in out, check for water trapping after oil injection
        {   if(connection->model()->conductCWater())
			{///. THE RULES TO BE CHECKED
                if(disconnectRetreading || connection->model()->disConectedCentreWCornerW() )
                {
                    Netsim::findMarkStoreTrappedWater_reCalcOlPc(m_layerEventsCh,connection, bulkBlob, localPc); 
					if(disconnectRetreading && connection->model()->disConectedCentreWCornerW())
						Netsim::findMarkStoreTrappedWater_reCalcOlPc(m_layerEventsCh,connection, filmBlob, localPc);
				}
            }
            else if(disconnectRetreading && connection->model()->conductAnyWater())
            {//&& connection->model()->bulk Fluid() == oil -> no need to re-add to drain events
                Netsim::findMarkStoreTrappedWater_reCalcOlPc(m_layerEventsCh,connection, filmBlob, localPc);
            }

				insertReCalcDrainEntryPrs(m_eventsCh, connection, localPc);
				addElemTo_layerDrainEvents(m_layerEventsCh,connection);
        }
    }

    clearTrappedWatFromEvents(m_eventsCh);
}


/** trapping routine
 */
void Netsim::findMarkStoreTrappedWater_reCalcOlPc(SortedEvents<Apex*,PceDrainCmp>& m_layerEventsCh, Element* elem, FluidBlob startPt, double localPc) 
{
    vector< pair<Element*, FluidBlob> > trappingStorage;
    double rhogh(elem->model()->rhogh());

    elem->findMarkTrappedWaterGanglia(localPc-rhogh, startPt, trappingStorage, m_cpuTimeTrapping, m_trappingCriteria);///. marks trapOilIndx ///unTrapWat(

    if(!trappingStorage.empty())
    {
		int WarningImproveEfficiency;
		dbgFile<<' '<<trappingStorage.size()<<'t'<<'w'<<elem->iRockType()<<" ";dbgFile.flush();
        for(size_t elm = 0; elm < trappingStorage.size(); ++elm)
        {
			if(trappingStorage[elm].first->isInOilFloodVec())
			{
				m_layerEventsCh.remove(trappingStorage[elm].first); ///. Candidate to be removed to improve efficiency
				trappingStorage[elm].first->setInOilFloodVec(false);
			}
			Polygon* polyShape = dynamic_cast< Polygon* >(trappingStorage[elm].first->ChModel());
			if(polyShape)
			{
				for(int cor = 0; cor < polyShape->numCorners(); ++cor)
				{
					if(polyShape->oilLayerConst()[cor].isInOilFloodVec())
					{
						m_layerEventsCh.remove(&polyShape->oilLayerCh()[cor]); ///. Candidate to be removed to improve efficiency
						polyShape->oilLayerCh()[cor].setInOilFloodVec(false);
					}
				}
				Netsim::reCalcOilLayerPc_markTrappedFilms(m_layerEventsCh, trappingStorage[elm].first, localPc);///. markstrapIndside/Outside from trapOilIndx
				Netsim::addElemTo_layerDrainEvents(m_layerEventsCh, trappingStorage[elm].first);
			}
        }

		sort(trappingStorage.begin(), trappingStorage.end(), TrappingWatStorageCmp());  // Having it sorted hels us when we want to coalesce the blob
        m_comn.addTrappedRegionWat(trappingStorage);
    }
    ///. trapped elements are left in the events storage to be cleared later just before they pop up
}






/**
 * Done with drainage  => clean up
 */
void Netsim::finaliseDrainage()
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

    m_out << endl
        << "===================== Drainage   ===================== "                                   << endl
        << "Total elapsed time for  drainage:        "  << m_cpuTimeTotal                               << endl
        << "Solving for water relative permeability: "  << m_cpuTimeKrw                                 << endl
        << "Solving for oil relative permeability:   "  << m_cpuTimeKro                                 << endl
        << "Solving for Resistivity Index:           "  << m_cpuTimeResIdx                              << endl
        << "Identifying trapped elements:            "  << m_cpuTimeTrapping                            << endl
        << "Coalesceing trapped oil:                 "  << m_cpuTimeCoal                                << endl
        << "Max water flowrate error:                "  << m_maxWatFlowErr*100.0 << " %"                << endl
        << "Max oil flowrate error:                  "  << m_maxOilFlowErr*100.0 << " %"                << endl
        << "Max resistivity index error:             "  << m_maxResIdxErr*100.0 << " %"                 << endl
        << endl
        << "===================== Network State  ===================== "                                << endl
        << "Maximum capillary pressure reached (Pa): "  << m_maxCycleCappPress                               << endl
        << "Water saturation:                        "  << m_satWater                                   << endl
        << "Number of elements invaded:              "  << m_totNumFillings                             << endl
        << "Remaining uninvaded elements:            "  << (int)m_rockLattice.size()-2-m_totNumFillings << endl
        << "Number of trapped water regions:         "  << (int)m_comn.numTrappedWatRegions()    << endl
        << "Total number of trapped water elements:  "  << (int)m_comn.numTrappedWatElems()      << endl
        << "Number of trapped oil regions:           "  << (int)m_comn.numTrappedOilRegions()    << endl
        << "Total number of trapped oil elements:    "  << (int)m_comn.numTrappedOilElems()      << endl
        << endl;

    m_amottDataDrainage[2] = m_satWater;
    if(m_cappPress < 0.0 && m_amottDataDrainage[1] < 0.0)
        m_amottOilIdx = 1.0;
    else if(m_amottDataDrainage[1] < 0.0)
        m_amottOilIdx = 0.0;
    else
    {
        m_amottOilIdx = (m_amottDataDrainage[0]-m_amottDataDrainage[1])/
            (m_amottDataDrainage[0]-m_amottDataDrainage[2]);
    }
    //if(m_cappPress < 0.0 ) 		m_amottOilIdx = 1.0;
    //else 						m_amottOilIdx = 0.0;


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


    if(m_writeDrainList)
    {
        m_drainListOut << "];" << endl;
        m_drainListOut.close();
    }

    //if(m_prtPressureProfile && m_wantRelPerm)
    //{
        //writePrsProfileData(retreatingFluid, m_floodCycleResultsOut); 
        //writePrsProfileData(invadingFluid, m_floodCycleResultsOut);
    //}


    for(int i = 1; i < int(m_rockLattice.size()); ++i)
    {
        if(i != m_numPores + 1 && m_rockLattice[i]->connectedToNetwork())
        {
            m_rockLattice[i]->ChModel()->finitOilInjection(m_cappPress-m_rockLattice[i]->model()->rhogh());        
        }
    }

	m_out<<"\n:/"<<endl<<endl;

}











void Netsim::popUpdate_layerInjs_Oil(Apex*apex,  SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, double & localPc)
{///. oil layer growth

		dbgFile<<'l'<<'o';dbgFile.flush();

		Polygon* polyCh = (Polygon*)apex->parentModel();
        const int & cornerIndex = apex->subIndex();
		const Polygon* poly = polyCh;

		softAssert(!poly->eleman()->isTrappedOil());
		//softAssert(!poly->conductCOil());
     if(polyCh->bulkFluid()->isOil() )
     {

		if(m_comn.debugMode > 0) 
		{
			cout<<" bydsdspd ";			cout<<apex->isInOilFloodVec()<<"  ";			cout<< apex->subIndex()<<endl;
			apex->setInOilFloodVec(false);
			vector<int> addCrns;
			cout<<(apex->parentModel()->eleman()->addToLayerVec(invadingFluid, polyCh, addCrns));

			//cout<<currElem->latticeIndex()<<"  ";
			cout<<apex->parentModel()->eleman()->canBeAddedToEventVec(&m_oil)<<"  ";
			//cout<<currElem->model()->displacementType()<<endl;
			//exit(-1);
		}
		return;
	 }
		
	   double currentPc = apex->entryPc();
	   if ( !polyCh->Pc_growStableOilLayerDrain_UseLess(localPc,cornerIndex)+poly->rhogh() ) return;
	   localPc = max(currentPc, localPc);
		//softAssert(currentPc < nextEventPc);
		softAssert(m_eventsCh.empty() || currentPc <= m_eventsCh.peek()->entryPc());

		if(poly->numLayers() == poly->numCorners()-1)
		{  // Oblique corner  => water in center and corner is separated, check for trapping
            Netsim::findMarkStoreTrappedWater_reCalcOlPc(m_layerEventsCh, polyCh->ChParent(), filmBlob, localPc);
            Netsim::findMarkStoreTrappedWater_reCalcOlPc(m_layerEventsCh, polyCh->ChParent(), bulkBlob, localPc);
        }
		else if(poly->numLayers() == 0)   ///. last element in W-W, first element in O-W OInj 
		{  if (poly->numLayers()==0) cout<<"dsssdknb";
          // Sharpest corner  => new oil connection need to check for oil coalescence
			Netsim::insertReCalcDrainEntryPrs(m_layerEventsCh,polyCh->ChParent(),localPc);
			//Netsim::addElemToDrainEvents(m_layerEventsCh,polyCh->ChParent());
			for(int i = 0; i < poly->eleman()->connectionNum(); ++i)
			{
				Netsim::untrap_OilGanglia(m_eventsCh, polyCh->eleman()->connection(i) );
				Netsim::addElemTo_layerDrainEvents(m_layerEventsCh, polyCh->eleman()->connection(i));
			}
		}
} 




void Netsim::untrap_WaterGanglia(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem, FluidBlob blob)
{
    pair<int, double> trapWatInfo = elem->trappingWat(blob);

    if(trapWatInfo.first > -1)                                                      // We need to untrap this water blob
    {
        clock_t startCoalesce(clock());
        const vector< pair<Element*,FluidBlob> >& newElems = m_comn.trappedRegionsWat(trapWatInfo.first);
        trapWatInfo.second += elem->model()->rhogh();
		double localPc = trapWatInfo.second;

		dbgFile<<newElems.size()<<'u'<<'W';dbgFile.flush();

			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				if(newElems[elm].first->isInOilFloodVec())
				{
					if (!m_eventsCh.remove(newElems[elm].first)) cout<<" sykr1 ";cout.flush();
					newElems[elm].first->setInOilFloodVec(false);
				}
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[elm].first->ChModel());
				if(polyShape)
				{
					for(int i = 0; i < polyShape->numCorners(); ++i)
					{
						if(polyShape->oilLayerConst()[i].isInOilFloodVec())
						{
							if (!m_layerEventsCh.remove(&polyShape->oilLayerCh()[i])) cout<<" gonw ";
							polyShape->oilLayerCh()[i].setInOilFloodVec(false);
						}
					}
					softAssert(!newElems[elm].first->isInOilFloodVec());
					softAssert(!polyShape->oilLayerConst()[0].isInOilFloodVec());
				}
				///. start from the local ganglia Pc and gradually equlibriate
				//newElems[elm].first->ChModel()->finitOilInjection(max(localPc, mcappPress)-newElems[elm].first->ChModel()->rhogh());
				newElems[elm].first->unTrapWat(newElems[elm].second);
				if(polyShape)
				{
					polyShape->calcOilLayerPc_markUntrappedFilms(localPc-polyShape->rhogh());
				}
			}


		if(localPc > m_cappPress) 
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm].first->ChModel()->finitOilInjection(localPc-newElems[elm].first->ChModel()->rhogh());
				newElems[elm].first->ChModel()->initWaterInjection(localPc-newElems[elm].first->ChModel()->rhogh());
			}

			SortedEvents< Apex*, PceImbCmp > waterFillingEvents;

			for(size_t i = 0; i < newElems.size(); ++i)               // First stage: check for unstable configurations
			{                                                       // when the pressure in the coalesced blob is
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[i].first->ChModel());  // equilibrated with the rest
				if(polyShape)
				{
					//Netsim::reCalcWaterFilmPc_markUnrappedFilms(waterFillingEvents, newElems[i].first, m_cappPress);
					//polyShape->calcOilLayerPc _markUntrappedFilms(localPc-polyShape->rhogh());
					polyShape->insertWatSnapEvent_IfSnapPcHgPc(waterFillingEvents, m_cappPress-polyShape->rhogh());   //  => insert water snap off + oil layer collapse
				}
			}


			if (!waterFillingEvents.empty())                                  // have been identified it's time to increase/reduce
			{   //do the events                                               // pressure in the blob. Dropping pressure will in fact be
				waterFillingEvents.sortEvents();                              // water injection process (water snap off, layer collapse).
				while(!waterFillingEvents.empty() && waterFillingEvents.peek()->gravCorrectedEntryPress() > m_cappPress) 
				{
					dbgFile<<'V'<<'o'<<'!';dbgFile.flush();
					bool tmp;
					popUpdateWaterInj(waterFillingEvents, tmp, localPc, m_cappPress);
				}

			}
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
					newElems[elm].first->ChModel()->finitWaterInjection(localPc-newElems[elm].first->ChModel()->rhogh());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 

			//double lowestPc=max(trapWatInfo.second,m_cappPress);
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm].first->ChModel()->initOilInjection(localPc);
			}

		}




        for(size_t elm = 0; elm < newElems.size(); ++elm)                         // Final stage: Now that all is peach, we can
        {                                                                   // possibly start adding new elements to the
			newElems[elm].first->calcCentreEntryPrsOilInj();
            // if(!newElems[elm]->isTrappedOil())	                           // displacement vectors
            {
                insertReCalcDrainEntryPrs(m_eventsCh, newElems[elm].first,localPc);
                addElemTo_layerDrainEvents(m_layerEventsCh, newElems[elm].first);
                for(int k = 0; k < newElems[elm].first->connectionNum(); ++k)
                {
					insertReCalcDrainEntryPrs(m_eventsCh, newElems[elm].first->connection(k),localPc);
					addElemTo_layerDrainEvents(m_layerEventsCh, newElems[elm].first->connection(k));
				}
			}
		}

		m_comn.removeTrappedRegionWat(trapWatInfo.first);
		m_cpuTimeCoal += cpuTimeElapsed(startCoalesce);
		
		
		while(!m_eventsCh.empty() && m_eventsCh.peek()->gravCorrectedEntryPress() <= m_cappPress)
		{
			bool tmp;
			popUpdateOilInj(m_eventsCh, tmp, localPc, m_cappPress);
		}

    }
}




///. never does anything TO DELETE
void Netsim::checkUntrapWaterIfUnstableConfigsDrain(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh)
{
    dbgFile<<"[";dbgFile.flush();
    for(size_t i = 0; i < m_rockLattice.size(); ++i)
    {
        Element* elem = m_rockLattice[i];
        double rhogh(elem->model()->rhogh());
        if(elem->model()->hasOilLayer_TrappedOutside_PcHsnapPc(m_cappPress-rhogh))
        {
			dbgFile<<'D';dbgFile.flush();
			cout<<"D";cout.flush();

           if(elem->isInOilFloodVec())
            {
                if (m_eventsCh.remove(elem)) elem->setInOilFloodVec(false);
                else  cout<<" cled "; cout.flush();
            }

			softAssert(!elem->model()->bulkFluid()->isOil());
			
            double capillaryPressure=-1000000000.0;
		   Netsim::popUpdateCentreInj_Oil(elem, m_eventsCh, capillaryPressure);

		}
    }
    dbgFile<<"]"<<"  ";
}

 



void Netsim::untrap_OilGanglia(SortedEvents<Apex*,PceDrainCmp>& m_eventsCh, Element* elem)
{
    pair<int, double> trapOilInfo = elem->trappingOil();

    if(trapOilInfo.first > -1)                                                      // We need to untrap this oil blob
    {
        clock_t startCoalesce(clock());
        vector<Element*> newElems = m_comn.trappedRegionsOil(trapOilInfo.first);
		trapOilInfo.second += elem->model()->rhogh();
		double localPc = trapOilInfo.second;

		dbgFile<<newElems.size()<<'u'<<'o';dbgFile.flush();

			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				if(newElems[elm]->isInOilFloodVec())///. useless
				{
					if (!m_eventsCh.remove(newElems[elm])) cout<<" sykr ";cout.flush();
					newElems[elm]->setInOilFloodVec(false);
				}
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[elm]->ChModel());
				if(polyShape)
				{
					for(int i = 0; i < polyShape->numCorners(); ++i)
					{
						if(polyShape->oilLayerConst()[i].isInOilFloodVec())
						{
							if (!m_layerEventsCh.remove(&polyShape->oilLayerCh()[i])) cout<<" gonw ";
							polyShape->oilLayerCh()[i].setInOilFloodVec(false);
						}
					}
					softAssert(!newElems[elm]->isInOilFloodVec());
					softAssert(!polyShape->oilLayerConst()[0].isInOilFloodVec());
				}

				//newElems[elm]->ChModel()->finitOilInjection(localPc-newElems[elm]->ChModel()->rhogh());
				newElems[elm]->unTrapOil();
				if(polyShape)
				{
					polyShape->calcOilLayerPc_markUntrappedFilms(m_cappPress-polyShape->rhogh());
				}
			}
 

		if(localPc > m_cappPress) ///. higher curvature, we need to inject water from films to equilibrate
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 
			
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm]->ChModel()->finitOilInjection(localPc-newElems[elm]->ChModel()->rhogh());
				newElems[elm]->ChModel()->initWaterInjection(localPc-newElems[elm]->ChModel()->rhogh());
			}

			SortedEvents< Apex*, PceImbCmp > waterFillingEvents;///. < element , corner indx or -1forCentre, localPc >

			for(size_t i = 0; i < newElems.size(); ++i)   //store the events            // First stage: check for unstable configurations
			{
				Polygon* polyShape = dynamic_cast< Polygon* >(newElems[i]->ChModel());
				if(polyShape)
				{
					//Netsim::reCalcOilLayerPc_markUntrappedFilms(waterFillingEvents, newElems[i], m_cappPress);
					//polyShape->calcOilLayerPc _markUntrappedFilms(localPc-polyShape->rhogh());
					polyShape->insertWatSnapEvent_IfSnapPcHgPc(waterFillingEvents, m_cappPress-polyShape->rhogh());   //  => insert water snap off + oil layer collapse
        if(polyShape->waterLayer_UntrappedCorner_PcLsnapPc(m_cappPress))
        {
			dbgFile<<'I';dbgFile.flush();
			cout<<'I';cout.flush();
		}
				}
			}

			if (!waterFillingEvents.empty())                                                          	// Second stage: Once all possible unstable configurations
			{   //do the events                                              // have been identified it's time to increase/reduce
				waterFillingEvents.sortEvents();                            // pressure in the blob. Dropping pressure will in fact be
				while(!waterFillingEvents.empty() && waterFillingEvents.peek()->gravCorrectedEntryPress() > m_cappPress)
				{
					dbgFile<<'V'<<'o';dbgFile.flush();

					bool tmp;
					popUpdateWaterInj(waterFillingEvents, tmp, localPc, m_cappPress);
				}
				
			}
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
					newElems[elm]->ChModel()->finitWaterInjection(localPc-newElems[elm]->ChModel()->rhogh());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			m_comn.injectant(invadingFluid); //================================================= 

			//double lowestPc=min(trapOilInfo.second,m_cappPress);
			for(size_t elm = 0; elm < newElems.size(); ++elm)
			{
				newElems[elm]->ChModel()->initOilInjection(localPc-newElems[elm]->ChModel()->rhogh());
			}

		}




        for(size_t elm = 0; elm < newElems.size(); ++elm)                         // Final stage: Now that all is peach, we can
        {                                                                   // possibly start adding new elements to the
			//newElems[elm]->calcCentreEntryPrsOilInj();
            // if(!newElems[elm]->isTrappedOil())	                           // displacement vectors
            {
                insertReCalcDrainEntryPrs(m_eventsCh, newElems[elm],localPc);
                addElemTo_layerDrainEvents(m_layerEventsCh, newElems[elm]);
                for(int k = 0; k < newElems[elm]->connectionNum(); ++k)
                {
					insertReCalcDrainEntryPrs(m_eventsCh, newElems[elm]->connection(k),localPc);
					addElemTo_layerDrainEvents(m_layerEventsCh, newElems[elm]->connection(k));
				}
			}
		}

		m_comn.removeTrappedRegionOil(trapOilInfo.first);
		m_cpuTimeCoal += cpuTimeElapsed(startCoalesce);
		
		
		while(!m_eventsCh.empty() && m_eventsCh.peek()->gravCorrectedEntryPress() <= m_cappPress)
		{
			bool tmp;
			popUpdateOilInj(m_eventsCh, tmp, localPc, m_cappPress);
		}

    }
}



