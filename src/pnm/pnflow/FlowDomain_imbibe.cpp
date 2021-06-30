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
#include "netStatistics.h"


#include "FlowDomain.h"




/**
* Inject water (imbibition) until target water saturation is reached or alternatively that the queue
* containing elements to be popped is empty. After drainage throats connected to the outlets have a
* maximum of one oil neighbour. If all elements were drained water injection will then occur at both
* faces. We can remove an injection face by increasing the oil neighbour count in throats connected to
* that that face. If we remove both faces and all elements were previously drained, the firts imbibition
* event will be a snap off.
 */
void FlowDomain::Imbibition(double requestedFinalSw, double requestedFinalPc,  double deltaSw, bool wantRelPerm, bool wantResIdx,  bool entreL, bool entreR, bool exitL, bool exitR)  {

	out_ << "\n********************  Water-injection -- cycle "<<comn_.dispCycle()+1<<"  ********************"<<endl;


	Events<Apex*,PceImbCmp>    evnts;

	 solver_ = new hypreSolver(elemans_, krInletBoundary_, krOutletBoundary_, nBSs_,nBpPors_, debugLevel, title_+"_solverImb", writeSlvMatrixAsMatlab_);

	
	///At the end of draiange imbibition displacement is initialized as max Pc is set. This
	///involves determening entry pressures for all elements.
	{
		dd_ = -1;

		invadingFluid = &water_;
		retreatingFluid = &oil_;

		comn_.incrementFloodCycle();
		comn_.injectant(invadingFluid);

		outD<<std::endl<<std::endl<<"initializeImbibition cycle "<<comn_.dispCycle()<<" "<<endl;

		cpuTimeTotal_ = 0.;
		cpuTimeCoal_ = 0.;
		cpuTimeKrw_ = 0.;
		cpuTimeKro_ = 0.;
		cpuTimeTrapping_ = 0.;
		cpuTimeResIdx_ = 0.;
		totNumFillings_ = 0;
		minCyclePc_ = Pc_;
		comn_.setMaxPcLastDrainCycle(maxCyclePc_);

		maxOilFlowErr_ = 0.;		maxWatFlowErr_ = 0.;		maxResIdxErr_ = 0.;




		if(sourceNode_>0)  {
			if(sourceNode_ > nPors_)  {
				cerr<<"\nError: Source pore ("<< sourceNode_<<") out of range [1,"<<nPors_<<"]."<<endl;            exit(-1);
			}
			trappingCriteria_ = escapeToBoth;
			if(entreL || entreR)	out_ << "Warning injecting from both inlet/outlet boundary and source node"<<endl;
		}
		else if(!entreL && !entreR) out_ << "Error: injection boundary is not specified"<<endl;

		if(exitL && exitR)    trappingCriteria_ = escapeToEither;
		else if(exitL)        trappingCriteria_ = escapeToInlet;
		else if(exitR)        trappingCriteria_ = escapeToOutlet;
		else    			  out_ << "Error: exit boundary is not specified "<<exitL << exitR<<endl;
	   

		setElemProps(input_, elemans_, nBpPors_, out_, comn_);



		solve_forRelPermResIndex(wantRelPerm, wantResIdx);   // Create starting points if not previously available

		wantRelPerm_ = wantRelPerm;
		wantResIdx_ = wantResIdx;

		///. initialize inlet and outlets
		if(entreR)       ///. Right
			((InOutBoundary*)elemans_[OutI])->fillElemCentreWithWaterCreateLayersIO(Pc_+0.1);      ///.  inlet BC
		else if(!entreR)
			((InOutBoundary*)elemans_[OutI])->fillElemCentreWithOilRemoveLayersIO(Pc_-1000000000.1);     ///.  outlet BC

		if(entreL)    ///. Left
			((InOutBoundary*)elemans_[0])->fillElemCentreWithWaterCreateLayersIO(Pc_+0.1);                 ///.  inlet BC
		else if(!entreL)
			((InOutBoundary*)elemans_[0])->fillElemCentreWithOilRemoveLayersIO(Pc_-1000000000.1);            ///.  outlet BC

		if(sourceNode_ != 0 && elemans_[sourceNode_]->model()->bulkFluid() != invadingFluid)
			elemans_[sourceNode_]->fillElemCentreWithWaterCreateLayers();     ///.  source BC


		outD<<" 2"; outD.flush();

		///. initialize all elements, except inlet and outlet BCs
		for(size_t i = nBpPors_; i < elemans_.size(); ++i)
			if(elemans_[i]->connectedToNetwork())
				elemans_[i]->ChModel()->initWaterInjection(Pc_-elemans_[i]->gravityCorrection());

		for(int i = nBSs_; i <  nBpPors_; ++i)
			if(elemans_[i]->connectedToNetwork())
				elemans_[i]->ChModel()->initWaterInjection(Pc_-elemans_[i]->gravityCorrection());
	   

		int nInWaterFlood0 = 0;
		for(int i = nBSs_; i < int(elemans_.size()); ++i)  {
			if(elemans_[i]->connectedToNetwork())  {
				if(elemans_[i]->canBeAddedToEventVec(water_))  {
					if (!elemans_[i]->model()->bulkFluid()->isOil())  		cout<< " depn" ;

					elemans_[i]->calcCentreEntryPrsWatInj();
					evnts.quickInsert(elemans_[i]);
					elemans_[i]->setInWatFloodVec(true);
					if( elemans_[i]->rockIndex() == 0 && elemans_[i]->isInWatFloodVec() )		++nInWaterFlood0;
				}
				addElemTo_layerImbibeEvents(evnts, elemans_[i]);
			}
		}

		outD<<"3"; outD.flush();

		///.  sort events to find the highest priority element to start injecting from
		evnts.sortEvents();



		int nInWaterFlood = 0, count1WF = 0, nNotInWaterFlood = 0, count1NIWF = 0;   
		for(int i = nBSs_; i < int(elemans_.size()); ++i)  {
			if(elemans_[i]->connectedToNetwork())  {   Elem* elem = elemans_[i];
				if(elem->rockIndex() == 0)  {
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
		cout<<" n_EventsCh "<<evnts.size()<<endl;
		(cout<<" n_rockLattice "<<elemans_.size()).flush();



		//int WARNING_TOCKECK;
		/*///. connect oil to all elements connected to inlet, ERROR TODO DELETE TOTEST
		//if(injAtLeftRes_)
		//{
			//for(int inT = 0; inT < elemans_[0]->nCncts(); ++inT)
			//{
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[0]->neib(inT), filmBlob);
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[0]->neib(inT), bulkBlob);
			//}
		//}
	//
		//if(entreR)
		//{
			//for(int outT = 0; outT < elemans_[OutI]->nCncts(); ++outT)
			//{
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[OutI]->neib(outT), filmBlob);
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[OutI]->neib(outT), bulkBlob);
			//}
		//}
	//
		//if(sourceNode_ != 0)
		//{
			//for(int sourceT = 0; sourceT < elemans_[sourceNode_]->nCncts(); ++sourceT)
			//{
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[sourceNode_]->neib(sourceT), filmBlob);
				//FlowDomain::untrap_WaterGanglia(evnts, elemans_[sourceNode_]->neib(sourceT), bulkBlob);
			//}
		//}*/

		amottDataImbibition_[0] = Sw_;
		amottDataImbibition_[1] = -1.;

		usbmDataImbibition_.clear();
		recordUSBMData(false);

		if(input_.fillingList(1))  {
			imbListOut_.open("fill_imbcycle_"+_s(comn_.dispCycle())+".m");
			imbListOut_ << "% The backbone identifies which pores/throats are water filled at the start of water flooding.\n"
			            << "% The first row is pore/throat index, followed by 1 for pores and 0 for thoats.\n";
			imbListOut_ << "backbone = [";
			for(size_t i = nBSs_; i < elemans_.size(); ++i)
				if(!elemans_[i]->model()->containCOil())
					imbListOut_ << elemans_[i]->indexOren() << ", "<< (dynamic_cast<Pore*>(elemans_[i]) != 0) << "; ..." << endl;
			imbListOut_ << "];\n\n" 
			            << "% The filling list identifies the order through which pores/throats get filled by water\n"
			            << "fill = [";
		}


		FlowDomain::checkUntrapOilIfUnstableConfigsImb(evnts);



		outD<<" . ";
		cout<<endl;

	}



	if (debugLevel>100)  {	cout<<elemans_[0]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<elemans_[0]->model()->Pc_pistonTypeRec()<<endl;
		cout<<elemans_[OutI]->model()->Pc_pistonTypeAdv()<<"  ";
		cout<<elemans_[OutI]->model()->Pc_pistonTypeRec()<<endl;
	}


	bool satCompress(false), compressWat(false), compressOil(false);
	double krThreshold, newDeltaSw, dSw(deltaSw);
	input_.relPermCompression(satCompress, krThreshold, newDeltaSw, compressWat, compressOil);
	if(satCompress && wantRelPerm && compressWat) dSw = newDeltaSw;
	double SwTarget = min(requestedFinalSw, Sw_+ dSw*0.5);
	double PcTarget = max(requestedFinalPc, Pc_ - (KcIncr_+abs(Pc_)*KcIncrFactor_)*0.1);
	bool residualSat(false);


	while(/*!residualSat &&*/ Sw_ <=  requestedFinalSw && Pc_ > requestedFinalPc)  {

		FlowDomain::singleImbibeStep(evnts, SwTarget, PcTarget, residualSat);

		double krw = watFlowRate_ / (singlePhaseWaterQ_+1e-200);
		double kro = oilFlowRate_ / (singlePhaseOilQ_+1e-200);;
		//double resIdx = resistivityIdx_;


		if(satCompress && wantRelPerm && compressWat && krw > krThreshold) dSw = deltaSw;
		if(satCompress && wantRelPerm && compressOil && kro < krThreshold) dSw = newDeltaSw;
		SwTarget = min(requestedFinalSw+1e-15, round((Sw_ + 0.75*dSw)/dSw)*dSw);
		PcTarget = max(requestedFinalPc-1e-7, Pc_ - (KcIncr_+abs(Pc_)*KcIncrFactor_+1e-16));
	}

	FlowDomain::writeResultData(wantRelPerm, wantResIdx);

	solver_=NULL;
	//FlowDomain::finaliseImbibition();

	/// - Done with imbibition  => clean up
	//void FlowDomain::finaliseImbibition()
	{

		if(maxOilFlowErr_ > 0.1 || maxWatFlowErr_ > 0.1 || maxResIdxErr_ > 0.1)  {
			out_ << endl
				<< "==================================================== "   << endl
				<< "Warning: For some rel perm calculations there were"     << endl
				<< "more than 10% difference between flowrates computed"    << endl
				<< "at the inlet and outlet. "                                       << endl
				<< "Max water flowrate error:    "  << maxWatFlowErr_      << endl
				<< "Max oil flowrate error:"  << maxOilFlowErr_      << endl
				<< "Max resistivity index error: "  << maxResIdxErr_       << endl
				<< "==================================================== "   << endl
				<< endl;
		}

		vector< int > sumfillingEventsP(6);
		vector< int > imbEventT(3);

		for(int ii = nBSs_; ii < static_cast< int >(elemans_.size()); ++ii)  {
			int evntI = elemans_[ii]->eventI();

			if(ii <  nBpPors_)  {
				if(evntI == -1) ++sumfillingEventsP[0];
				else if(evntI == 0 || evntI == 1) ++sumfillingEventsP[1];
				else if(evntI == 2) ++sumfillingEventsP[2];
				else if(evntI == 3) ++sumfillingEventsP[3];
				else if(evntI > 3)  ++sumfillingEventsP[4];
				else ++sumfillingEventsP[5];
			}
			else if (ii >= nBpPors_)  {
				if(evntI == -1) ++imbEventT[0];
				else if(evntI == 0 || evntI == 1) ++imbEventT[1];
				else ++imbEventT[2];
			}
		}

		out_<< endl
			<< "===================== Imbibition  ===================== "                                   << endl
			<< "Total elapsed time for imbibition:       "  << cpuTimeTotal_                               << endl
			<< "Solving for water relative permeability: "  << cpuTimeKrw_                                 << endl
			<< "Solving for oil relative permeability:   "  << cpuTimeKro_                                 << endl
			<< "Solving for resistivity index:           "  << cpuTimeResIdx_                              << endl
			<< "Identifying trapped elements:            "  << cpuTimeTrapping_                            << endl
			<< "Coalesceing trapped water:               "  << cpuTimeCoal_                                << endl
			<< "Max water flow solver residual:          "  << maxWatFlowErr_*100. << " %"                << endl
			<< "Max oil flow solver residual:            "  << maxOilFlowErr_*100. << " %"                << endl
			<< "Max resistivity index residual:          "  << maxResIdxErr_*100. << " %"                 << endl
			<< endl
			<< "===================== Network State  ===================== "                                << endl
			<< "Minimum capillary pressure reached (Pa): "  << Pc_                                  << endl
			<< "Water saturation:                        "  << Sw_                                   << endl
			<< "Number of elements invaded:              "  << totNumFillings_                             << endl
			<< "Number of trapped water regions:         "  << (int)comn_.numTrappedWatRegions()    << endl
			<< "Total number of trapped water elements:  "  << (int)comn_.numTrappedWatElems()      << endl
			<< "Number of trapped oil regions:           "  << (int)comn_.numTrappedOilRegions()    << endl
			<< "Total number of trapped oil elements:    "  << (int)comn_.numTrappedOilElems()      << endl
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

		amottDataImbibition_[2] = Sw_;
		if(Pc_ > 0.)
			amottWaterIdx_ = 1.;
		else if(Pc_ < 0. && amottDataImbibition_[1] < 0.)
			amottWaterIdx_ = 0.;
		else  {
			amottWaterIdx_ = (amottDataImbibition_[1]-amottDataImbibition_[0])/
				(amottDataImbibition_[2]-amottDataImbibition_[0]);
		}

		if(comn_.dispCycle() > 1)  {
			out_ << endl
				<< "=================Wettability State  ================ "                     << endl
				<< "Amott water index, Iw:                  " << amottWaterIdx_                << endl
				<< "Amott oil index, Io:                    " << amottOilIdx_                  << endl
				<< "Amott wettability index, I = Iw - Io:   " << amottWaterIdx_-amottOilIdx_  << endl
				<< "USBM wettability index:                 " << calcUSBMindex()                << endl
				<< endl;
		}


		if(imbListOut_)  {
			imbListOut_ << "];" << endl;
			imbListOut_.close();
		}

		//if(prtPressureProfile_ && wantRelPerm_)
		//{
			//writePrsProfileData(invadingFluid, floodCycleResultsOut_);
			//writePrsProfileData(retreatingFluid, floodCycleResultsOut_);
		//}


		for(int i = 2; i < int(elemans_.size()); ++i)
			if(elemans_[i]->connectedToNetwork())
				elemans_[i]->ChModel()->finitWaterInjection(Pc_-elemans_[i]->gravityCorrection());       

		out_<<"\n:/\n\n"<<endl;

	}


}



/**
 * Do a single displacement with water and recalculate saturation and relative permeability.
 * We do not allow for empty filling Vec as this just might create a mess in imbibition. The
 * reason is that as Pc increases (negatively) more and more regions might get trapped as oil
 * layers collapse, and we don't want to spend huge amounts of time checking on this.
 */
void FlowDomain::singleImbibeStep(Events<Apex*,PceImbCmp>& evnts, double SwTarget, double& PcTarget, bool& residualSat)  {

	clock_t startClock(clock());


	if(SwTarget < Sw_ || PcTarget > Pc_)  {
		out_ << "================================================="              << endl
			<< "Nothing to be done:"                                            << endl;
		if(SwTarget < Sw_)			out_ << "Target water saturation (" << SwTarget << ")   is lower than current saturation (" << Sw_ << ")"   << endl;
		else			out_ << "Target capillary pressure (" << PcTarget << ")    is higher than current pressure (" << Pc_ << ")"  << endl;
		out_ << "=================================================="             << endl;

		return;
	}


	//double lastStepWaterSat = Sw_;

	//out_ << "Sw = "<< setw(6) << Sw_ << " -> ";

	int numSteps(0), totNumFill(0);
	int fillTarget = max(minNumFillings_,
						int(initStepSize_*(nPors_+nTrots_)*(SwTarget-Sw_)));
	double oldCappPress(Pc_);

	while(Sw_ <= SwTarget && Pc_ > PcTarget-1e-32 /*&& !evnts.empty()*/)  {
		double oldWaterSat(Sw_);
		int numInv(0), invInsideBox(0);
		bool insideBox(false);

		outD<<"\n Sw = "<<Sw_<<":"<<"  Pc = "<<Pc_<<": ";outD.flush();

		while(invInsideBox < fillTarget && !evnts.empty() && nextCentrInjPc(evnts) >= PcTarget)  {
			FlowDomain::popUpdateWaterInj(evnts,insideBox, cappPressCh_, PcTarget);
			comn_.GuessCappPress(Pc_);                   // Use local maxima for Pc
			++numInv;
			if(insideBox) ++invInsideBox;
			if(Pc_==0. || Pc_*oldCappPress < 0.)  {

				FlowDomain::checkUntrapOilIfUnstableConfigsImb(evnts);
				oldCappPress=Pc_-1e-12;
				updateSatAndConductances(Pc_);
				out_ << "Pc cross over at Sw = " << setw(6) << Sw_ << "; ";
				amottDataImbibition_[1] = Sw_;//FlowDomain::recordAmottData(false);
				FlowDomain::solve_forRelPermResIndex(wantRelPerm_, wantResIdx_);
			   	out_ << endl;
			}

		///. Third inner loop, only when option stableFilling_, until no layer ready for pop(?)
			if(stableFilling_)  
			 while(nextCentrInjPc(evnts) >= Pc_-1e-32)  {	//cout<<" StableFilling ";
					FlowDomain::popUpdateWaterInj(evnts, insideBox, cappPressCh_, PcTarget);
					++numInv;
					if(insideBox) ++invInsideBox;
				 }
		}

		if (nextCentrInjPc(evnts) < PcTarget && (Pc_ > PcTarget) )  {
			cappPressCh_ = PcTarget;

			FlowDomain::updateSatAndConductances(Pc_);

			if(Pc_==0. || Pc_*oldCappPress < 0.)  {
				out_ << "Pc cross over at Sw = " << setw(6) << Sw_ << "\n";
				amottDataImbibition_[1] = Sw_;//FlowDomain::recordAmottData(false);
				oldCappPress=Pc_-1e-12;
			}

		}

		FlowDomain::checkUntrapOilIfUnstableConfigsImb(evnts);
		FlowDomain::updateSatAndConductances(Pc_);

		FlowDomain::recordUSBMData(false);
		fillTarget = max(minNumFillings_, (int)min((fillTarget*maxFillIncrease_),
			(extrapCutBack_*(invInsideBox/(Sw_-oldWaterSat))*(SwTarget-Sw_)) )  );
		totNumFill += numInv;
		++numSteps;
	}

	residualSat = evnts.empty();
	minCyclePc_ = min(Pc_, minCyclePc_);

	cout.precision(3);
	out_ << "Sw: " << setw(8) << std::left << Sw_ << " Pc: " << setw(10) << round(Pc_)  << " ";
	out_ << setw(2) << std::right << numSteps << " steps " << setw(5) << totNumFill << " invasions " ;

	totNumFillings_ += totNumFill;
	FlowDomain::solve_forRelPermResIndex(wantRelPerm_, wantResIdx_);


	out_ << endl;
	cpuTimeTotal_ += cpuTimeElapsed(startClock);
}



inline void FlowDomain::insertReCalcImbibeEntryPrs(Events<Apex*,PceImbCmp>& evnts, Elem* elem, double cappPrs)  {
	if(elem->isInWatFloodVec())  {
		if (evnts.remove(elem))  {   ///. order is important
				elem->calcCentreEntryPrsWatInj();
				evnts.insert(elem);
				if (!elem->model()->bulkFluid()->isOil()) cout<< " deaw " ;
		}
		else
		{
			elem->setInWatFloodVec(false);
			elem->calcCentreEntryPrsWatInj();
			if (elem->canBeAddedToEventVec(water_))  {

				evnts.insert(elem);		elem->setInWatFloodVec(true);
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
		if (elem->canBeAddedToEventVec(water_))  {
			elem->ChModel()->initWaterInjection(cappPrs);
			elem->calcCentreEntryPrsWatInj();
			evnts.insert(elem);
			elem->setInWatFloodVec(true);
		}
	}

	//elem->ChModel()->calcR(cappPrs);///. candidate to be removed

}



inline void FlowDomain::addElemTo_layerImbibeEvents(Events<Apex*,PceImbCmp>& evnts, Elem* elem)  {
	Polygon* shyp = dynamic_cast< Polygon* >(elem->ChModel());
	if(shyp)  {
		vector<int> addCrns;
		if(elem->addToLayerVec(water_, shyp, addCrns))  {
			ensure(!addCrns.empty());
			for(size_t i = 0; i < addCrns.size(); ++i)  {
				ensure(!shyp->oilLayerConst()[addCrns[i]].isInWatFloodVec());
				evnts.insert(&shyp->oilLayerCh()[addCrns[i]]);
				shyp->oilLayerCh()[addCrns[i]].setInWatFloodVec(true);
			}
		}
	}
}


inline void FlowDomain::clearTrappedOilFromEvents(Events<Apex*,PceImbCmp>&    evnts)  {
	while(!evnts.empty() && evnts.peek()->parentModel()->eleman()->isTrappedOil())  {
		evnts.pop()->setInWatFloodVec(false);
	}

	while(!evnts.empty() && evnts.peek()->trappingCL().first>-1)  {
		evnts.pop()->setInWatFloodVec(false);
	}
}







/**
 * Pops the elements with the highest priority acording to the compare function and updates
 * that element which now contains water. New elements now available for injection are inserted
 * into the queue. The entry pressure rquired to do that invasion is returned.
 */
void FlowDomain::popUpdateWaterInj(Events<Apex*,PceImbCmp>& evnts, bool& insideBox, double& localPc, double localPcTarget)  {

	while (evnts.peek()->subIndex()>=0 && !evnts.empty() && evnts.peek()->gravCorrectedEntryPress() > localPcTarget)  {
		FlowDomain::popUpdate_layerInj_Water(evnts.pop(), evnts, localPc);

		clearTrappedOilFromEvents(evnts); ///. skip trapped regions
		if( evnts.empty() ) return;
	}

	if(evnts.empty() ||  evnts.peek()->gravCorrectedEntryPress() < localPcTarget) return;

	{
		Elem *currElemCh = evnts.pop()->parentModel()->ChParent();
		FlowDomain::popUpdateCentreInj_Water(currElemCh, evnts, localPc);

		if(imbListOut_)  {
			bool isAPore(dynamic_cast< Pore* >(currElemCh) != 0);
			imbListOut_ << currElemCh->indexOren() << ", ";
			imbListOut_ << isAPore << "; ..." << endl;
		}

		clearTrappedOilFromEvents(evnts);
		insideBox = currElemCh->isInCalcBox();
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
 * cost of 2*nCncts*O(logN). For large networks this should be cheaper.
 */
void FlowDomain::popUpdateCentreInj_Water(Elem* currElemCh, Events<Apex*,PceImbCmp>& evnts, double& localPc)  {

	const Elem *currElem = currElemCh;
	double currentPc = currElem->gravCorrectedEntryPress();
	 if (evnts.remove(currElemCh)) cout<<" dblInsWInj " ;
	if(!currElemCh->model()->conductCOil())  {
		currElemCh->setInWatFloodVec(false);
		cout<<"dsfgfg "<<currElemCh->canBeAddedToEventVec(water_)<<endl;
		if (evnts.remove(currElemCh)) cout<<" dblInsElm" ;
		return ;
	}

	if(!evnts.empty() && currentPc < evnts.peek()->entryPc()) cout<<"\n\nError: currentPc < evnts.peek()->entryPc() \n"<<endl;

	localPc = min(localPc, currentPc);
	//boundPress_ = max(boundPress_, currentPc);
	outD<<currElem->model()->displacementType()<<currElem->rockIndex();outD.flush();

	bool HadNoInvading(!currElem->model()->conductAnyWater());

	currElemCh->fillElemCentreWithWaterCreateLayers();         ///. popUpdateCentreInj_Water    inject through centre



	if(currElem->isTrappedWat(filmBlob) && !currElem->model()->disConectedCentreWCornerW())
	   untrap_WaterGanglia(evnts, currElemCh, filmBlob);

	for(int j = 0; j < currElem->nCncts(); ++j)  {
		if(HadNoInvading) FlowDomain::untrap_WaterGanglia(evnts, currElem->neib(j), filmBlob);///. invading fluid coalescence
		untrap_WaterGanglia(evnts, currElem->neib(j), bulkBlob);
	}

	FlowDomain::addElemTo_layerImbibeEvents(evnts, currElemCh);///. for film and centre coalescence

	bool disconnectRetreading(!currElem->model()->conductsAnyOil());

	for(int i = 0; i < currElem->nCncts(); ++i)  {
		Elem *connection = currElem->neib(i);
		if(!connection->isEntryOrExitRes())///. ! in out, check for water trapping after oil injection
		{///. THE RULES TO BE CHECKED
			if(disconnectRetreading)  {
				 FlowDomain::findMarkStoreTrappedOil(evnts,connection,localPc);
			}

				insertReCalcImbibeEntryPrs(evnts, connection,localPc);
				addElemTo_layerImbibeEvents(evnts,connection);
		}
	}

	clearTrappedOilFromEvents(evnts);

}


/** trapping routine
 */
void FlowDomain::findMarkStoreTrappedOil(Events<Apex*,PceImbCmp>& evnts, Elem* elem, double localPc)  {
	vector<Elem*> trappingStorage;

	elem->findMarkTrappedOilGanglia(localPc-elem->gravityCorrection(), trappingStorage, cpuTimeTrapping_, trappingCriteria_);

	if(!trappingStorage.empty())  {
		outD<<' '<<trappingStorage.size()<<'t'<<'o'<<elem->rockIndex()<<" ";outD.flush();
		for(auto& elm:trappingStorage)  {
			if(elm->isInWatFloodVec())  {
				evnts.remove(elm);
				elm->setInWatFloodVec(false);
			}
			Polygon* shyp = dynamic_cast< Polygon* >(elm->ChModel());
			if(shyp)  {
				for(int cor = 0; cor < shyp->numCorners(); ++cor)  {
					if(shyp->oilLayerConst()[cor].isInWatFloodVec())  {
						evnts.remove(&shyp->oilLayerCh()[cor]);
						shyp->oilLayerCh()[cor].setInWatFloodVec(false);
					}
				}

				shyp->calcOilLayerPc_syncTrappings(localPc-elm->gravityCorrection());
				//addElemTo_layerImbibeEvents(evnts,connection);
			}
				//insertReCalcImbibeEntryPrs(evnts, connection,Pc_);
				//elm->parentModel()->ChParent()->ChModel()->calcR(localPc);///. candidate to be removed


		}

		//sort(trappingStorage.begin(), trappingStorage.end(), TrappingOilStorageCmp());  // Having it sorted hels us when we want to coalesce the blob
		comn_.addTrappedRegionOil(trappingStorage);
	}
	///. trapped elements are left in the events storage to be cleared later just before they pop up
}















void FlowDomain::popUpdate_layerInj_Water(Apex*apex,  Events<Apex*,PceImbCmp>& evnts, double & localPc)  {///. oil layer collapse

		outD<<'L'<<'w';outD.flush();

		Polygon* polyCh = (Polygon*)apex->parentModel();
		const int & cornerIndex = apex->subIndex();
		const Polygon* poly = polyCh;

		ensure(!poly->eleman()->isTrappedOil());


	   double currentPc = polyCh->Polygon::Pc_pin_disconnectOilLayer(cornerIndex)+poly->eleman()->gravityCorrection();
	   localPc = min(currentPc, localPc);
		ensure(evnts.empty() || currentPc >= evnts.peek()->entryPc());

		if(poly->numLayers() == poly->numCorners()-1)  {  // Oblique corner  => water in center and corner is joined, check for coalescing
			FlowDomain::untrap_WaterGanglia(evnts, poly->ChParent(), filmBlob);
			FlowDomain::untrap_WaterGanglia(evnts, poly->ChParent(), bulkBlob);
		}
		else if(poly->numLayers()==0)//cornerIndex == 0)   ///. last element in W-W, first element in O-W OInj
		{           //Sharpest corner  => no more oil, check for oil trapping
			if (poly->numLayers()==0) cout<<"dsknb";
			for(int i = 0; i < poly->eleman()->nCncts(); ++i)  {
				FlowDomain::findMarkStoreTrappedOil(evnts, poly->eleman()->neib(i),localPc);
			}
		}
}



void FlowDomain::untrap_WaterGanglia(Events<Apex*,PceImbCmp>& evnts, Elem* elem, FluidBlob blob)  {
	pair<int, double> trapWatInfo = elem->trappingWat(blob);

	if(trapWatInfo.first > -1)                                                      // We need to untrap this water blob
	{
		clock_t startCoalesce(clock());
		const vector< pair<Elem*,FluidBlob> >& newElems = comn_.trappedRegionsWat(trapWatInfo.first);
		trapWatInfo.second += elem->gravityCorrection();
		double localPc = trapWatInfo.second;

		outD<<'['<<newElems.size()<<'u'<<'w';outD.flush();

		for(auto& elm:newElems)  {
			if(elm.first->isInWatFloodVec())  {
				if (!evnts.remove(elm.first)) (cout<<" sykr1 ").flush();
				elm.first->setInWatFloodVec(false);
			}
			Polygon* shyp = dynamic_cast< Polygon* >(elm.first->ChModel());
			if(shyp)  {
				for(int i = 0; i < shyp->numCorners(); ++i)  {
					if(shyp->oilLayerConst()[i].isInWatFloodVec())  {
						if (!evnts.remove(&shyp->oilLayerCh()[i])) cout<<" gonwi ";
						shyp->oilLayerCh()[i].setInWatFloodVec(false);
					}
				}
				ensure(!elm.first->isInWatFloodVec());
				ensure(!shyp->oilLayerConst()[0].isInWatFloodVec());
			}
			///. start from the local ganglia Pc and gradually equlibriate
			//elm.first->ChModel()->finitWaterInjection(max(localPc, mcappPress)-elm.first->gravityCorrection());
			elm.first->unTrapWat(elm.second);
			if(shyp)  {
				shyp->calcOilLayerPc_markUntrappedFilms(localPc-elm.first->gravityCorrection());
			}
		}


		if(localPc < Pc_) ///. lower curvature, we need to inject oil from centre to equilibrate
		{

			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			for(auto& elm:newElems)  {
				elm.first->ChModel()->finitWaterInjection(localPc-elm.first->gravityCorrection());
				elm.first->ChModel()->initOilInjection(localPc-elm.first->gravityCorrection());
			}

			Events< Apex*, PceDrainCmp > oilFillingEvents;

			for(size_t i = 0; i < newElems.size(); ++i)               // First stage: check for unstable configurations
			{                                                       // when the pressure in the coalesced blob is
				Polygon* shyp = dynamic_cast< Polygon* >(newElems[i].first->ChModel());  // equilibrated with the rest
				if(shyp)  {
					shyp->insertOilSnapEvent_IfSnapPcLgPc(oilFillingEvents, Pc_-newElems[i].first->gravityCorrection());       //  => insert oil snap off
					if(shyp->hasOilLayer_TrappedOutside_PcHsnapPc(Pc_))  {
						outD<<'D';outD.flush();
						cout<<"D";cout.flush();
					}
				}
			}

			if (!oilFillingEvents.empty())  
			{   //do the events
				oilFillingEvents.sortEvents();
				//FlowDomain::increasePressureInCoalescedBlobTrue(oilFillingEvents,localPc, Pc_); 
				while(!oilFillingEvents.empty() &&  oilFillingEvents.peek()->gravCorrectedEntryPress() < Pc_)  {
					outD<<'^';outD.flush();

					bool tmp;
					popUpdateOilInj(oilFillingEvents, tmp, localPc, Pc_-0.000001);
				}

			}
			for(auto& elm:newElems)  {
					elm.first->ChModel()->finitOilInjection(Pc_-elm.first->gravityCorrection());
			}
			{const Fluid* OldInv=invadingFluid; invadingFluid=retreatingFluid; retreatingFluid=OldInv;}
			comn_.injectant(invadingFluid); //=================================================

			//double highestPc=max(trapWatInfo.second,Pc_);
			for(auto& elm:newElems)  {
				elm.first->ChModel()->initWaterInjection(localPc-elm.first->gravityCorrection());
			}

		}




		for(auto& elm:newElems)                         // Final stage: Now that all is peach, we can
		{                                                                   // possibly start adding new elements to the
			//elm.first->calcCentreEntryPrsWatInj();
			// if(!elm->isTrappedOil())							   // displacement vectors
			{
				insertReCalcImbibeEntryPrs(evnts, elm.first,localPc);
				addElemTo_layerImbibeEvents(evnts, elm.first);
				for(int k = 0; k < elm.first->nCncts(); ++k)  {
					insertReCalcImbibeEntryPrs(evnts, elm.first->neib(k),localPc);
					addElemTo_layerImbibeEvents(evnts, elm.first->neib(k));
				}
			}
		}

		comn_.removeTrappedRegionWat(trapWatInfo.first);
		cpuTimeCoal_ += cpuTimeElapsed(startCoalesce);


		while(!evnts.empty() && evnts.peek()->gravCorrectedEntryPress() >= Pc_)  {
			bool tmp;
			popUpdateWaterInj(evnts, tmp, localPc, Pc_);
		}
	outD<<']';outD.flush();
	}
}




///. never does anything TO DELETE
void FlowDomain::checkUntrapOilIfUnstableConfigsImb(Events<Apex*,PceImbCmp>& evnts)  {
	outD<<"[";outD.flush();
	for(Elem* elem:elemans_)  {
		if(elem->model()->waterLayer_UntrappedCorner_PcLsnapPc(Pc_-elem->gravityCorrection()))  {
			outD<<'I';outD.flush();
			cout<<'I';cout.flush();

			if(elem->isInWatFloodVec())  {
				if (evnts.remove(elem)) elem->setInWatFloodVec(false);
				else  (cout<<" clei ").flush();
			}
			ensure(elem->model()->bulkFluid()->isOil());

			double capillaryPressure=1000000000.;
		   FlowDomain::popUpdateCentreInj_Water(elem, evnts, capillaryPressure);

		}
	}
	outD<<"]"<<"  ";
}





