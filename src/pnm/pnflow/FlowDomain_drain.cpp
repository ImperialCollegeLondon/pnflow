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
#include "hypreSolver.h"


#include "FlowDomain.h"




/**
 * Continue injecting oil until a target saturation is reached, draiange queue is empty
 * or capillary pressure target is reached. Writing results to output file. Oil injection
 * needs to occur at least at one face.
 */
void FlowDomain::Drainage(double requestedFinalSw, double requestedFinalPc,  double deltaSw, bool wantRelPerm, bool wantResIdx,    bool entreL, bool entreR, bool exitL, bool exitR)
{

	out_ << "\n********************  Oil-injection -- cycle "<<comn_.dispCycle()+1<<"  ********************"<<endl;

	SortedEvents<Apex*,PceDrainCmp>    evnts;


	solver_ = new hypreSolver(elemans_, krInletBoundary_, krOutletBoundary_, nBSs_,nBpPors_, debugLevel, "solverDrain", writeSlvMatrixAsMatlab_);
	{
		dd_ = 1;

		invadingFluid = &oil_;
		retreatingFluid = &water_;

		comn_.incrementFloodCycle();
		comn_.injectant(invadingFluid);

		outD<<std::endl<<std::endl<<"initializeDrainage cycle "<<comn_.dispCycle()<<" "<<endl;

		cpuTimeTotal_ = 0.0;
		cpuTimeCoal_ = 0.0;
		cpuTimeKrw_ = 0.0;
		cpuTimeKro_ = 0.0;
		cpuTimeTrapping_ = 0.0;
		cpuTimeResIdx_ = 0.0;
		totNumFillings_ = 0;
		maxCyclePc_ = Pc_;
		comn_.setMinPcLastCycle(minCyclePc_);

		maxOilFlowErr_ = 0.0;		maxWatFlowErr_ = 0.0;		maxResIdxErr_ = 0.0;



		if(sourceNode_>0)
		{
			if(sourceNode_ > nPors_)
			{
				cerr<<"\nError: Source pore ("<< sourceNode_<<") out of range [1,"<<nPors_<<"]."<<endl;            exit(-1);
			}
			trappingCriteria_ = escapeToBoth;
			if(entreL || entreR) out_ << "Injecting from inlet/outlet and source node "<<sourceNode_<<endl;
			else                                out_ << "Injecting from node "<<sourceNode_<<endl;
		}
		else if(!entreL && !entreR) out_ << "Error: injection boundary is not specified"<<endl;

		if(exitL && exitR)  trappingCriteria_ = escapeToEither;
		else if(exitL)      trappingCriteria_ = escapeToInlet;
		else if(exitR)      trappingCriteria_ = escapeToOutlet;
		else                out_ << "Error: exit boundary is not specified "<<exitL << exitR<<endl;




		setElemProps(input_, elemans_, nBpPors_, out_, comn_);



		solve_forRelPermResIndex(wantRelPerm, wantResIdx);   // Create starting points if not previously available


		wantRelPerm_ = wantRelPerm;
		wantResIdx_ = wantResIdx;

		///. initialize inlet and outlets
		if(entreR)       ///. Right
		{
			((InOutBoundary*)elemans_[OutI])->fillElemCentreWithOilRemoveLayersIO(Pc_-0.1);        ///.  inlet BC
		}
		else if(!entreR)
		{
			((InOutBoundary*)elemans_[OutI])->fillElemCentreWithWaterCreateLayersIO(Pc_+1000000000.1);             ///.  outlet BC
		}

		if(entreL)    ///. Left
		{
			((InOutBoundary*)elemans_[0])->fillElemCentreWithOilRemoveLayersIO(Pc_-0.1);            ///.  inlet BC
		}
		else if(!entreL)
			((InOutBoundary*)elemans_[0])->fillElemCentreWithWaterCreateLayersIO(Pc_+1000000000.1);                ///.  outlet BC

		if(sourceNode_ != 0 && elemans_[sourceNode_]->model()->bulkFluid() != invadingFluid)
			elemans_[sourceNode_]->fillElemCentreWithOilRemoveLayers();   ///.  source node

		//evnts.clear();
		//evnts.clear();

		outD<<" 2"; outD.flush();

		///. initialize all elements, except inlet and outlet BCs
		for(size_t i = nBpPors_; i < elemans_.size(); ++i)
		{
			if(elemans_[i]->connectedToNetwork())
			{
				elemans_[i]->ChModel()->initOilInjection(Pc_-elemans_[i]->gravityCorrection());
			}
		}
		for(int i = nBSs_; i <  nBpPors_; ++i)
		{
			if(elemans_[i]->connectedToNetwork())
			{
				elemans_[i]->ChModel()->initOilInjection(Pc_-elemans_[i]->gravityCorrection());
			}
		}


		for(int i = nBSs_; i < int(elemans_.size()); ++i)
		{
			if(elemans_[i]->connectedToNetwork())
			{
				if(elemans_[i]->canBeAddedToEventVec(oil_))
				{
					if (elemans_[i]->model()->bulkFluid()->isOil())  		cout<< " depb " ;

					elemans_[i]->calcCentreEntryPrsOilInj();
				   evnts.quickInsert(elemans_[i]);
					elemans_[i]->setInOilFloodVec(true);

				}
				addElemTo_layerDrainEvents(evnts, elemans_[i]);
			}
		}

		outD<<"3"; outD.flush();

		///.  sort events to find the highest priority element to start injecting from
		evnts.sortEvents();
		evnts.sortEvents();


		/*/////. connect oil to all elements connected to inlet, ERROR TODO DELETE TOTEST
		//if(entreL)
		//{
			//for(int inT = 0; inT < elemans_[0]->connectionNum(); ++inT)
			//{
				//FlowDomain::untrap_OilGanglia(evnts, elemans_[0]->connection(inT));
			//}
		//}
	//
		//if(entreR)
		//{
			//for(int outT = 0; outT < elemans_[OutI]->connectionNum(); ++outT)
			//{
				//FlowDomain::untrap_OilGanglia(evnts, elemans_[OutI]->connection(outT));
			//}
		//}
	//
		//if(sourceNode_ != 0)
		//{
			//for(int sourceT = 0; sourceT < elemans_[sourceNode_]->connectionNum(); ++sourceT)
			//{
				//FlowDomain::untrap_OilGanglia(evnts, elemans_[sourceNode_]->connection(sourceT));
			//}
		//}*/

		amottDataDrainage_[0] = Sw_;
		amottDataDrainage_[1] = -1.0;

		usbmDataDrainage_.clear();
		recordUSBMData(true);

		if(writeDrainList_)
		{
			ostringstream fileName;
			fileName << title_ << "fill_cycle" << comn_.dispCycle() << "_drain.m";
			drainListOut_.open(fileName.str().c_str());
			drainListOut_ << "% The backbone identifies which pores/throats are oil filled at the start of drainage." << endl
				<< "% The first row is pore/throat index, followed by 1 for pores and 0 for thoats." << endl;
			drainListOut_ << "backbone = [";
			for(size_t i = nBSs_; i < elemans_.size(); ++i)
			{
				if(elemans_[i]->model()->containCOil())
				{
					bool isAPore(dynamic_cast< Pore* >(elemans_[i]) != 0);
					drainListOut_ << elemans_[i]->indexOren() << ", ";
					drainListOut_ << isAPore << "; ..." << endl;
				}
			}
			drainListOut_ << "];" << endl;
			drainListOut_ << endl << "% The filling list identifies the order through which pores/throats get filled by oil" << endl;
			drainListOut_ << "fill = [";
		}


		FlowDomain::checkUntrapWaterIfUnstableConfigsDrain(evnts);



		outD<<" . ";
		cout<<endl;

	}



	if (debugLevel>100)
	{	cout<<elemans_[0]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<elemans_[0]->model()->Pc_pistonTypeRec()<<endl;
		cout<<elemans_[OutI]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<elemans_[OutI]->model()->Pc_pistonTypeRec()<<endl;
	}


	bool satCompress(false), compressWat(false), compressOil(false);
	double krThreshold, newDeltaSw, dSw(deltaSw);
	input_.relPermCompression(satCompress, krThreshold, newDeltaSw, compressWat, compressOil);
	if(satCompress && wantRelPerm && compressOil) dSw = newDeltaSw;
	double SwTarget = max(requestedFinalSw, Sw_ - dSw*0.5);
	double PcTarget = min(requestedFinalPc, Pc_ + (KcIncr_+abs(Pc_)*KcIncrFactor_)*0.1);
	bool residualSat(false);


	while(Sw_ >= requestedFinalSw && Pc_ < requestedFinalPc)
	{

		FlowDomain::singleDrainStep(evnts, SwTarget, PcTarget, residualSat);

		double krw = watFlowRate_ /singlePhaseWaterQ_;
		double kro = oilFlowRate_ /singlePhaseOilQ_;

		if(satCompress && wantRelPerm && compressOil && kro > krThreshold) dSw = deltaSw;
		if(satCompress && wantRelPerm && compressWat && krw < krThreshold) dSw = newDeltaSw;
		SwTarget = max(requestedFinalSw-1.0e-15, round((Sw_ - 0.75*dSw)/dSw)*dSw);
		PcTarget = min(requestedFinalPc+1.0e-7, Pc_ + (KcIncr_+abs(Pc_)*KcIncrFactor_));
	}

	FlowDomain::writeResultData(wantRelPerm, wantResIdx);

	delete solver_;
	solver_=NULL;
	//FlowDomain::finaliseDrainage();

	/// Done with drainage  => clean up
	//void FlowDomain::finaliseDrainage()
	{

		if(maxOilFlowErr_ > 0.1 || maxWatFlowErr_ > 0.1 || maxResIdxErr_ > 0.1)
		{
			out_ << endl
				<< "==================================================== "  << endl
				<< "Warning: For some rel perm calculations there were"     << endl
				<< "more than 10% difference between flowrates computed"    << endl
				<< "at the inlet and outlet. "                              << endl
				<< "Max water flowrate error:    "  << maxWatFlowErr_       << endl
				<< "Max oil flowrate error:"  << maxOilFlowErr_       << endl
				<< "Max resistivity index error: "  << maxResIdxErr_        << endl
				<< "==================================================== "  << endl
				<< endl;
		}

		out_ << endl
			<< "===================== Drainage   ===================== "                   << endl
			<< "Total elapsed time for  drainage:        " << cpuTimeTotal_                << endl
			<< "Solving for water relative permeability: " << cpuTimeKrw_                  << endl
			<< "Solving for oil relative permeability:   " << cpuTimeKro_                  << endl
			<< "Solving for resistivity index:           " << cpuTimeResIdx_               << endl
			<< "Identifying trapped elements:            " << cpuTimeTrapping_             << endl
			<< "Coalesceing trapped oil:                 " << cpuTimeCoal_                 << endl
			<< "Max water flow solver residual:          " << maxWatFlowErr_*100.0 << " %" << endl
			<< "Max oil flow solver residual:            " << maxOilFlowErr_*100.0 << " %" << endl
			<< "Max resistivity index residual:          " << maxResIdxErr_*100.0 << " %"  << endl
			<< endl
			<< "===================== Network State  ===================== "                         << endl
			<< "Maximum capillary pressure reached (Pa): " << maxCyclePc_                            << endl
			<< "Water saturation:                        " << Sw_                                    << endl
			<< "Number of elements invaded:              " << totNumFillings_                        << endl
			<< "Remaining uninvaded elements:            " << int(elemans_.size())-2-totNumFillings_ << endl
			<< "Number of trapped water regions:         " << comn_.numTrappedWatRegions()      << endl
			<< "Total number of trapped water elements:  " << comn_.numTrappedWatElems()        << endl
			<< "Number of trapped oil regions:           " << comn_.numTrappedOilRegions()      << endl
			<< "Total number of trapped oil elements:    " << comn_.numTrappedOilElems()        << endl
			<< endl;

		amottDataDrainage_[2] = Sw_;

		if(Pc_ < 0.0 && amottDataDrainage_[1] < 0.0)   amottOilIdx_ = 1.0;
		else if(amottDataDrainage_[1] < 0.0)           amottOilIdx_ = 0.0;
		else  amottOilIdx_ = (amottDataDrainage_[0]-amottDataDrainage_[1])/(amottDataDrainage_[0]-amottDataDrainage_[2]);

		if(comn_.dispCycle() > 1)
		{
			out_ << endl
				<< "=================Wettability State  ================ "                     << endl
				<< "Amott water index, Iw:                  " << amottWaterIdx_                << endl
				<< "Amott oil index, Io:                    " << amottOilIdx_                  << endl
				<< "Amott wettability index, I = Iw - Io:   " << amottWaterIdx_-amottOilIdx_  << endl
				<< "USBM wettability index:                 " << calcUSBMindex()                << endl
				<< endl;
		}


		if(writeDrainList_)
		{
			drainListOut_ << "];" << endl;
			drainListOut_.close();
		}

		//if(prtPressureProfile_ && wantRelPerm_)
		//{
			//writePrsProfileData(retreatingFluid, floodCycleResultsOut_);
			//writePrsProfileData(invadingFluid, floodCycleResultsOut_);
		//}


		for(int i = nBSs_; i < int(elemans_.size()); ++i)
			if(elemans_[i]->connectedToNetwork())
				elemans_[i]->ChModel()->finitOilInjection(Pc_-elemans_[i]->gravityCorrection());       

		out_<<"\n:/\n\n"<<endl;

	}


}



/**
 * Do a single displacement with oil and recalculate saturation and relative permeability
 */
void FlowDomain::singleDrainStep(SortedEvents<Apex*,PceDrainCmp>& evnts, double SwTarget, double PcTarget, bool& residualSat)
{

	clock_t startClock(clock());


	if(SwTarget >= Sw_ || PcTarget <= Pc_)
	{
		out_ << "================================================="<<endl << "Nothing to be done:"  <<endl;
		if(SwTarget>Sw_)  out_<<"Target water saturation ("<< SwTarget <<")   is higher than current Sw ("<<Sw_<<")" <<endl;
		else              out_<<"Target capillary pressure ("<< PcTarget <<")   is lower than current Pc ("<<Pc_<<")" <<endl;
		out_ << "=================================================="  <<endl;
		return;
	}


	///HASAN  add for loop for time stepping

	int numSteps(0), totNumFill(0);
	int fillTarget = max(minNumFillings_, int(initStepSize_*(nPors_+nTrots_)*(SwTarget-Sw_)));
	double oldCappPress(Pc_);

	while(Sw_ > SwTarget && Pc_ < PcTarget+1.0e-32 /*&& !evnts.empty()*/) ///. outer filling loop
	{
		double oldWaterSat(Sw_);
		int numInv(0), invInsideBox(0);
		bool insideBox(false);

		outD<<"\n Sw = "<<Sw_<<":"<<"  Pc = "<<Pc_<<"   "<<"  "<<nextCentrInjPc(evnts)<<": ";outD.flush();

		while((invInsideBox < fillTarget && !evnts.empty() && nextCentrInjPc(evnts) <= PcTarget) )  ///. inner filling loop
		{
			FlowDomain::popUpdateOilInj(evnts,insideBox, cappPressCh_, PcTarget); ///.  popUpdateCentreInj _Oil
			comn_.GuessCappPress(Pc_); /// used as initialGuess
			++numInv;
			if(insideBox) ++invInsideBox;
			if(Pc_==0.0 || Pc_*oldCappPress < 0.0)
			{

				FlowDomain::checkUntrapWaterIfUnstableConfigsDrain(evnts);

				updateSatAndConductances(Pc_);
				out_ << "Pc cross over at Sw = " << setw(6) << Sw_ << "; ";
				amottDataDrainage_[1] = Sw_;//FlowDomain::recordAmottData(true);
				FlowDomain::solve_forRelPermResIndex(wantRelPerm_, wantResIdx_);
				oldCappPress = Pc_+1.0e-12;
			   	out_ << endl;
			}


///. Third inner loop, only when option stableFilling_, until no layer ready for pop(?)
			if(stableFilling_)  
			 while(nextCentrInjPc(evnts) <= Pc_+1.0e-32)
				{	//cout<<" StableFilling ";
					FlowDomain::popUpdateOilInj(evnts, insideBox, cappPressCh_, PcTarget);
					++numInv;
					if(insideBox) ++invInsideBox;
				 }
		}

		if (nextCentrInjPc(evnts) > PcTarget && (Pc_ < PcTarget) )
		{
			//cout<<"  * "<<Pc_<<":"<<satTrack_->predictPc( Pc_,Sw_,PcTarget,SwTarget, 0)<<":"<<PcTarget<<" * ";
			cappPressCh_ =  PcTarget;
			//cappPressCh_ = min(PcTarget,Pc_+((1.0+double(numSteps))/50.0)*(PcTarget-Pc_));/// with additional underrelaxation to  check for saturation

			FlowDomain::updateSatAndConductances(Pc_);

			if(Pc_==0.0 || Pc_*oldCappPress < 0.0)
			{
				out_ << "Pc cross over at Sw = " << setw(6) << Sw_ << "\n";
				amottDataDrainage_[1] = Sw_;//FlowDomain::recordAmottData(true);
				oldCappPress = Pc_+1.0e-12;
			}


		}


		FlowDomain::checkUntrapWaterIfUnstableConfigsDrain(evnts);
		FlowDomain::updateSatAndConductances(Pc_);

		FlowDomain::recordUSBMData(true);
		fillTarget = max(minNumFillings_, (int)min((fillTarget*maxFillIncrease_),
			(extrapCutBack_*(invInsideBox/(Sw_-oldWaterSat))*(SwTarget-Sw_)) )  );
		totNumFill += numInv;
		++numSteps;
	}

	residualSat = evnts.empty();
	maxCyclePc_ = max(Pc_, maxCyclePc_);

	cout.precision(3);
	out_ << "Sw: " << setw(8) << std::left << Sw_ << " Pc: " << setw(10) << round(Pc_)  << " ";
	out_ << setw(2) << std::right << numSteps << ":" << setw(5) << totNumFill << " invasions " ;

	totNumFillings_ += totNumFill;
	FlowDomain::solve_forRelPermResIndex(wantRelPerm_, wantResIdx_);



	if (debugLevel>0 && Element::nErrs>0)
		cout<<"   nErrs:"<<Element::nErrs<<"  ";
	Element::nErrs = 0;

	out_ << endl;
	cpuTimeTotal_ += cpuTimeElapsed(startClock);
}



inline bool FlowDomain::insertReCalcDrainEntryPrs(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, double cappPrs)
{
	if(elem->isInOilFloodVec())
	{
		if (evnts.remove(elem))
		{   ///. order is important

				elem->calcCentreEntryPrsOilInj();///. order is important
				evnts.insert(elem);
				if (elem->model()->bulkFluid()->isOil()) cout<< " deab " ;
		}
		else
		{
			elem->setInOilFloodVec(false);
			elem->calcCentreEntryPrsOilInj();///. order is important
			if (elem->canBeAddedToEventVec(oil_))
			{
				evnts.insert(elem);	elem->setInOilFloodVec(true);
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
		if (elem->canBeAddedToEventVec(oil_))
		{
			elem->ChModel()->initOilInjection(cappPrs);
			elem->calcCentreEntryPrsOilInj();
			evnts.insert(elem);
			elem->setInOilFloodVec(true);
		}
	}
	//elem->ChS hape()->ca lcR(cappPrs);///. candidate to be removed
	return true;
}

template<class compType>
inline void FlowDomain::reCalcOilLayerPc_markTrappedFilms(SortedEvents<Apex*,compType>& evnts, Element* elem, double cappPrs)
{///. pc is just used for error handling
	Polygon* shyp = dynamic_cast< Polygon* >(elem->ChModel());
	if(shyp)
	{
		for(int i = 0; i < shyp->numCorners(); ++i)
		{
			if(shyp->oilLayerConst()[i].isInOilFloodVec())
			{
				if (!evnts.remove(&shyp->oilLayerCh()[i])) cout<<" 454545886 ";
				//evnts->setInOilFloodVec(false);
			}
		}

		shyp->calcOilLayerPc_syncTrappings(cappPrs-elem->gravityCorrection());

		for(int j = 0; j < shyp->numCorners(); ++j)
		{
			if(shyp->oilLayerConst()[j].isInOilFloodVec())
			{
				if (shyp->oilLayerConst()[j].trappedOLayer().first<0)
					evnts.insert(&shyp->oilLayerCh()[j]);
				else
				{
					shyp->oilLayerCh()[j].setInOilFloodVec(false);
					cout<<" nlsua ";
				}

			}
		}
	}
}


inline void FlowDomain::addElemTo_layerDrainEvents(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem)
{
	Polygon* shyp = dynamic_cast< Polygon* >(elem->ChModel());
	if(shyp)
	{
		vector<int> addCrns;
		if(elem->addToLayerVec(oil_, shyp, addCrns))
		{
			for(size_t i = 0; i < addCrns.size(); ++i)
			{
				ensure(!shyp->oilLayerConst()[addCrns[i]].isInOilFloodVec());
				evnts.insert(&shyp->oilLayerCh()[addCrns[i]]);
				shyp->oilLayerCh()[addCrns[i]].setInOilFloodVec(true);
			}
		}
	}
}



inline void FlowDomain::clearTrappedWatFromEvents(SortedEvents<Apex*,PceDrainCmp>& evnts)
{
	while(!evnts.empty() && evnts.peek()->parentModel()->eleman()->isTrappedWat(bulkBlob) )
	{
		evnts.pop()->setInOilFloodVec(false);
	}

	while(!evnts.empty() && evnts.peek()->trappingCL().first>-1)
	{
		evnts.pop()->setInOilFloodVec(false);
	}
}







/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains oil. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void FlowDomain::popUpdateOilInj(SortedEvents<Apex*,PceDrainCmp>& evnts, bool& insideBox, double& localPc, double localPcTarget)
{

	while (evnts.peek()->subIndex()>=0 && !evnts.empty() && evnts.peek()->gravCorrectedEntryPress() < localPcTarget)
	{
		FlowDomain::popUpdate_layerInjs_Oil(evnts.pop(), evnts, localPc);

		clearTrappedWatFromEvents(evnts); ///. skip trapped regions
		if( evnts.empty() ) return;
	}

	if(evnts.empty() ||  evnts.peek()->gravCorrectedEntryPress() > localPcTarget) return;

	{
		Element *currElemCh = evnts.pop()->parentModel()->ChParent();
		FlowDomain::popUpdateCentreInj_Oil(currElemCh, evnts, localPc);
		if(writeDrainList_)
		{
			bool isAPore(dynamic_cast< Pore* >(currElemCh) != 0);
			drainListOut_ << currElemCh->indexOren() << ", ";
			drainListOut_ << isAPore << "; ..." << endl;
		}

		clearTrappedWatFromEvents(evnts);
		insideBox = currElemCh->isInCalcBox();
	}
}



/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains oil. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void FlowDomain::popUpdateCentreInj_Oil(Element *currElemCh, SortedEvents<Apex*,PceDrainCmp>& evnts, double& localPc)
{

	const Element *currElem = currElemCh;
	double currentPc = currElem->gravCorrectedEntryPress();
	 if(currElem->model()->bulkFluid()->isOil())
	 {		 cout<<"\n bypd "<<currElem->index()<<"  "<<currElem->isInOilFloodVec()<<"  "<<currElem->model()->displacementType()<<endl;		 return; //exit(-1);
	 }

	localPc = max(localPc, currentPc);
	//boundPress_ = min(boundPress_, currentPc);
	outD<<currElem->model()->displacementType()<<currElem->rockIndex();outD.flush();

	bool HadNoOil(!currElem->model()->conductsAnyOil());//any oil

	currElemCh->fillElemCentreWithOilRemoveLayers();
	if (evnts.remove(currElemCh)) cout<<" dblInsOInj ";     ///. popUpdateCentreInj _Oil  inject through centre
	if (evnts.checkIfThere(currElemCh)) cout<<" dblInsOIndsdaadsj ";     ///. popUpdateCentreInj _Oil  inject through centre
	if (currElem->isInOilFloodVec())  cout<<" CWVN ";


	if(HadNoOil)
	{
		for(int j = 0; j < currElem->connectionNum(); ++j)
		{
			FlowDomain::untrap_OilGanglia(evnts, currElem->connection(j)); ///. invading fluid coalescence
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
					FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(evnts,connection, bulkBlob, localPc);
					if(disconnectRetreading && connection->model()->disConectedCentreWCornerW())
						FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(evnts,connection, filmBlob, localPc);
				}
			}
			else if(disconnectRetreading && connection->model()->conductAnyWater())
			{//&& connection->model()->bulk Fluid() == oil -> no need to re-add to drain events
				FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(evnts,connection, filmBlob, localPc);
			}

				insertReCalcDrainEntryPrs(evnts, connection, localPc);
				addElemTo_layerDrainEvents(evnts,connection);
		}
	}

	clearTrappedWatFromEvents(evnts);
}


/** trapping routine
 */
void FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, FluidBlob startPt, double localPc)
{
	vector< pair<Element*, FluidBlob> > trappingStorage;

	elem->findMarkTrappedWaterGanglia(localPc-elem->gravityCorrection(), startPt, trappingStorage, cpuTimeTrapping_, trappingCriteria_);///. marks trapOilIndx ///unTrapWat(

	if(!trappingStorage.empty())
	{
		//int WarningImproveEfficiency;
		outD<<' '<<trappingStorage.size()<<'t'<<'w'<<elem->rockIndex()<<" ";outD.flush();
		for(auto& elm:trappingStorage)
		{
			if(elm.first->isInOilFloodVec())
			{
				evnts.remove(elm.first); ///. Candidate to be removed to improve efficiency
				elm.first->setInOilFloodVec(false);
			}
			Polygon* shyp = dynamic_cast< Polygon* >(elm.first->ChModel());
			if(shyp)
			{
				for(int cor = 0; cor < shyp->numCorners(); ++cor)
				{
					if(shyp->oilLayerConst()[cor].isInOilFloodVec())
					{
						evnts.remove(&shyp->oilLayerCh()[cor]); ///. Candidate to be removed to improve efficiency
						shyp->oilLayerCh()[cor].setInOilFloodVec(false);
					}
				}
				FlowDomain::reCalcOilLayerPc_markTrappedFilms(evnts, elm.first, localPc);///. markstrapIndside/Outside from trapOilIndx
				FlowDomain::addElemTo_layerDrainEvents(evnts, elm.first);
			}
		}

		sort(trappingStorage.begin(), trappingStorage.end(), TrappingWatStorageCmp());  // Having it sorted hels us when we want to coalesce the blob
		comn_.addTrappedRegionWat(trappingStorage);
	}
	///. trapped elements are left in the events storage to be cleared later just before they pop up
}















void FlowDomain::popUpdate_layerInjs_Oil(Apex*apex,  SortedEvents<Apex*,PceDrainCmp>& evnts, double & localPc)
{///. oil layer growth

		outD<<'l'<<'o';outD.flush();

		Polygon* polyCh = (Polygon*)apex->parentModel();
		const int & cornerIndex = apex->subIndex();
		const Polygon* poly = polyCh;

		ensure(!poly->eleman()->isTrappedOil());
		//ensure(!poly->conductCOil());
	 if(polyCh->bulkFluid()->isOil() )
	 {

		if(debugLevel > 0)
		{
			cout<<" bydsdspd ";			cout<<apex->isInOilFloodVec()<<"  ";			cout<< apex->subIndex()<<endl;
			apex->setInOilFloodVec(false);
			vector<int> addCrns;
			cout<<(apex->parentModel()->eleman()->addToLayerVec(oil_, polyCh, addCrns));

			//cout<<currElem->index()<<"  ";
			cout<<apex->parentModel()->eleman()->canBeAddedToEventVec(oil_)<<"  ";
			//cout<<currElem->model()->displacementType()<<endl;
			//exit(-1);
		}
		return;
	 }

	   double currentPc = apex->entryPc();
	   if ( !polyCh->Pc_growStableOilLayerDrain_UseLess(localPc,cornerIndex)+poly->eleman()->gravityCorrection() ) return;
	   localPc = max(currentPc, localPc);
		//ensure(currentPc < nextEventPc);
		ensure(evnts.empty() || currentPc <= evnts.peek()->entryPc());

		if(poly->numLayers() == poly->numCorners()-1)
		{  // Oblique corner  => water in center and corner is separated, check for trapping
			FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(evnts, polyCh->ChParent(), filmBlob, localPc);
			FlowDomain::findMarkStoreTrappedWater_reCalcOlPc(evnts, polyCh->ChParent(), bulkBlob, localPc);
		}
		else if(poly->numLayers() == 0)   ///. last element in W-W, first element in O-W OInj
		{  if (poly->numLayers()==0) cout<<"dsssdknb";
		  // Sharpest corner  => new oil connection need to check for oil coalescence
			FlowDomain::insertReCalcDrainEntryPrs(evnts,polyCh->ChParent(),localPc);
			//FlowDomain::addElemToDrainEvents(evnts,polyCh->ChParent());
			for(int i = 0; i < poly->eleman()->connectionNum(); ++i)
			{
				FlowDomain::untrap_OilGanglia(evnts, polyCh->eleman()->connection(i) );
				FlowDomain::addElemTo_layerDrainEvents(evnts, polyCh->eleman()->connection(i));
			}
		}
}




void FlowDomain::untrap_WaterGanglia(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem, FluidBlob blob)
{
	pair<int, double> trapWatInfo = elem->trappingWat(blob);

	if(trapWatInfo.first > -1)                                                      // We need to untrap this water blob
	{
		clock_t startCoalesce(clock());
		const vector< pair<Element*,FluidBlob> >& newElems = comn_.trappedRegionsWat(trapWatInfo.first);
		trapWatInfo.second += elem->gravityCorrection();
		double localPc = trapWatInfo.second;

		outD<<newElems.size()<<'u'<<'W';outD.flush();

			for(auto& elm:newElems)
			{
				if(elm.first->isInOilFloodVec())
				{
					if (!evnts.remove(elm.first)) (cout<<" sykr1 ").flush();
					elm.first->setInOilFloodVec(false);
				}
				Polygon* shyp = dynamic_cast< Polygon* >(elm.first->ChModel());
				if(shyp)
				{
					for(int i = 0; i < shyp->numCorners(); ++i)
					{
						if(shyp->oilLayerConst()[i].isInOilFloodVec())
						{
							if (!evnts.remove(&shyp->oilLayerCh()[i])) cout<<" gonw1 ";
							shyp->oilLayerCh()[i].setInOilFloodVec(false);
						}
					}
					ensure(!elm.first->isInOilFloodVec());
					ensure(!shyp->oilLayerConst()[0].isInOilFloodVec());
				}
				///. start from the local ganglia Pc and gradually equlibriate
				//elm.first->ChModel()->finitOilInjection(max(localPc, mcappPress)-elm.first->gravityCorrection());
				elm.first->unTrapWat(elm.second);
				if(shyp)
				{
					shyp->calcOilLayerPc_markUntrappedFilms(localPc-elem->gravityCorrection());
				}
			}


		if(localPc > Pc_)
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			for(auto& elm:newElems)
			{
				elm.first->ChModel()->finitOilInjection(localPc-elm.first->gravityCorrection());
				elm.first->ChModel()->initWaterInjection(localPc-elm.first->gravityCorrection());
			}

			SortedEvents< Apex*, PceImbCmp > waterFillingEvents;

			for(size_t i = 0; i < newElems.size(); ++i)              // First stage: check for unstable configurations
			{                                                       // when the pressure in the coalesced blob is
				Polygon* shyp = dynamic_cast< Polygon* >(newElems[i].first->ChModel());  // equilibrated with the rest
				if(shyp)
				{
					shyp->insertWatSnapEvent_IfSnapPcHgPc(waterFillingEvents, Pc_-newElems[i].first->gravityCorrection());   //  => insert water snap off + oil layer collapse
				}
			}

			if (!waterFillingEvents.empty())                                  // have been identified it's time to increase/reduce
			{   //do the events                                               // pressure in the blob. Dropping pressure will in fact be
				waterFillingEvents.sortEvents();                              // water injection process (water snap off, layer collapse).
				while(!waterFillingEvents.empty() && waterFillingEvents.peek()->gravCorrectedEntryPress() > Pc_)
				{
					outD<<'V'<<'o'<<'!';outD.flush();
					bool tmp;
					popUpdateWaterInj(waterFillingEvents, tmp, localPc, Pc_);
				}

			}
			for(auto& elm:newElems)
			{
					elm.first->ChModel()->finitWaterInjection(localPc-elm.first->gravityCorrection());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			//double lowestPc=max(trapWatInfo.second,Pc_);
			for(auto& elm:newElems)
			{
				elm.first->ChModel()->initOilInjection(localPc);
			}

		}




		for(auto& elm:newElems)                         // Final stage: Now that all is peach, we can
		{                                                                   // possibly start adding new elements to the
			elm.first->calcCentreEntryPrsOilInj();
			// if(!elm->isTrappedOil())							   // displacement vectors
			{
				insertReCalcDrainEntryPrs(evnts, elm.first,localPc);
				addElemTo_layerDrainEvents(evnts, elm.first);
				for(int k = 0; k < elm.first->connectionNum(); ++k)
				{
					insertReCalcDrainEntryPrs(evnts, elm.first->connection(k),localPc);
					addElemTo_layerDrainEvents(evnts, elm.first->connection(k));
				}
			}
		}

		comn_.removeTrappedRegionWat(trapWatInfo.first);
		cpuTimeCoal_ += cpuTimeElapsed(startCoalesce);


		while(!evnts.empty() && evnts.peek()->gravCorrectedEntryPress() <= Pc_)
		{
			bool tmp;
			popUpdateOilInj(evnts, tmp, localPc, Pc_);
		}

	}
}




///. never does anything TO DELETE
void FlowDomain::checkUntrapWaterIfUnstableConfigsDrain(SortedEvents<Apex*,PceDrainCmp>& evnts)
{
	outD<<"[";outD.flush();
	for(size_t i = 0; i < elemans_.size(); ++i)
	{
		Element* elem = elemans_[i];
		if(elem->model()->hasOilLayer_TrappedOutside_PcHsnapPc(Pc_-elem->gravityCorrection()))
		{
			outD<<'D';outD.flush();
			cout<<"D";cout.flush();

		   if(elem->isInOilFloodVec())
			{
				if (evnts.remove(elem)) elem->setInOilFloodVec(false);
				else  (cout<<" cled ").flush();
			}

			ensure(!elem->model()->bulkFluid()->isOil());

			double capillaryPressure=-1000000000.0;
		   FlowDomain::popUpdateCentreInj_Oil(elem, evnts, capillaryPressure);

		}
	}
	outD<<"]"<<"  ";
}





void FlowDomain::untrap_OilGanglia(SortedEvents<Apex*,PceDrainCmp>& evnts, Element* elem)
{
	pair<int, double> trapOilInfo = elem->trappingOil();

	if(trapOilInfo.first > -1)                                                      // We need to untrap this oil blob
	{
		clock_t startCoalesce(clock());
		vector<Element*> newElems = comn_.trappedRegionsOil(trapOilInfo.first);
		trapOilInfo.second += elem->gravityCorrection();
		double localPc = trapOilInfo.second;

		outD<<newElems.size()<<'u'<<'o';outD.flush();

		for(auto& elm:newElems)
		{
			if(elm->isInOilFloodVec())///. useless
			{
				if (!evnts.remove(elm)) (cout<<" sykr ").flush();
				elm->setInOilFloodVec(false);
			}
			Polygon* shyp = dynamic_cast< Polygon* >(elm->ChModel());
			if(shyp)
			{
				for(int i = 0; i < shyp->numCorners(); ++i)
				{
					if(shyp->oilLayerConst()[i].isInOilFloodVec())
					{
						if (!evnts.remove(&shyp->oilLayerCh()[i])) cout<<" evol2 ";
						shyp->oilLayerCh()[i].setInOilFloodVec(false);
					}
				}
				ensure(!elm->isInOilFloodVec());
				ensure(!shyp->oilLayerConst()[0].isInOilFloodVec());
			}

			//elm->ChModel()->finitOilInjection(localPc-elm->gravityCorrection());
			elm->unTrapOil();
			if(shyp)
			{
				shyp->calcOilLayerPc_markUntrappedFilms(Pc_-elm->gravityCorrection());
			}
		}


		if(localPc > Pc_) ///. higher curvature, we need to inject water from films to equilibrate
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			for(auto& elm:newElems)
			{
				elm->ChModel()->finitOilInjection(localPc-elm->gravityCorrection());
				elm->ChModel()->initWaterInjection(localPc-elm->gravityCorrection());
			}

			SortedEvents< Apex*, PceImbCmp > waterFillingEvents;///. < element , corner indx or -1forCentre, localPc >

			for(size_t i = 0; i < newElems.size(); ++i)   //store the events            // First stage: check for unstable configurations
			{
				Polygon* shyp = dynamic_cast< Polygon* >(newElems[i]->ChModel());
				if(shyp)
				{
					shyp->insertWatSnapEvent_IfSnapPcHgPc(waterFillingEvents, Pc_-newElems[i]->gravityCorrection());   //  => insert water snap off + oil layer collapse
					if(shyp->waterLayer_UntrappedCorner_PcLsnapPc(Pc_))
					{
						outD<<'I';outD.flush();
						cout<<'I';cout.flush();
					}
				}
			}

			if (!waterFillingEvents.empty())                                                          	// Second stage: Once all possible unstable configurations
			{   //do the events                                              // have been identified it's time to increase/reduce
				waterFillingEvents.sortEvents();                            // pressure in the blob. Dropping pressure will in fact be
				while(!waterFillingEvents.empty() && waterFillingEvents.peek()->gravCorrectedEntryPress() > Pc_)
				{
					outD<<'V'<<'o';outD.flush();

					bool tmp;
					popUpdateWaterInj(waterFillingEvents, tmp, localPc, Pc_);
				}

			}
			for(auto& elm:newElems)
			{
					elm->ChModel()->finitWaterInjection(localPc-elm->gravityCorrection());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			//double lowestPc=min(trapOilInfo.second,Pc_);
			for(auto& elm:newElems)
			{
				elm->ChModel()->initOilInjection(localPc-elm->gravityCorrection());
			}

		}




		for(auto& elm:newElems)                         // Final stage: Now that all is peach, we can
		{                                                                   // possibly start adding new elements to the
			//elm->calcCentreEntryPrsOilInj();
			// if(!elm->isTrappedOil())							   // displacement vectors
			{
				insertReCalcDrainEntryPrs(evnts, elm,localPc);
				addElemTo_layerDrainEvents(evnts, elm);
				for(int k = 0; k < elm->connectionNum(); ++k)
				{
					insertReCalcDrainEntryPrs(evnts, elm->connection(k),localPc);
					addElemTo_layerDrainEvents(evnts, elm->connection(k));
				}
			}
		}

		comn_.removeTrappedRegionOil(trapOilInfo.first);
		cpuTimeCoal_ += cpuTimeElapsed(startCoalesce);


		while(!evnts.empty() && evnts.peek()->gravCorrectedEntryPress() <= Pc_)
		{
			bool tmp;
			popUpdateOilInj(evnts, tmp, localPc, Pc_);
		}

	}
}



