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
#include <array>
#include <algorithm>
#include <functional>
#include <map>

#include "FlowDomain.h"

using namespace std;
#include "hypreSolver.h"
#include "NetworkTransform.h"
#include "netStatistics.h"


#include "InOutCNM.cpp"
#include "setGVoidProps.h"





string title2(const InputFile& input, const string& title)
{
	string titl2 = input.getOr(title+"_res/","outputDir");
	if (titl2.back()=='/' || titl2.back()=='\\')
	{
		if (input.keyvals("writeMatrix").size() || input.keyvals("visuaLight").size() || input.keyvals("visualize").size() )
			cout<<"creating folder: "<<titl2
			#if defined(_WIN32)
				<< mkdir(titl2.c_str()) ///. check also _mkdir
			#else
				<< mkdir(titl2.c_str(), 0733) // note that 777 is different from 0777
			#endif
			<<endl;
		else if (input.keyvals("writeXmf").size()) titl2=title;
	}
	return titl2;
}

/// FlowDomain constructor
FlowDomain::FlowDomain(InputFile & input, unsigned int randSeed)
:	input_(input), comn_(input,randSeed), oil_(comn_.oil()), water_(comn_.water()), Pc_(cappPressCh_),
	title_(input.outputName()),  out_(title_ + "_pnflow.prt"),
	results3D_(input,&comn_,title2(input,title_), reinterpret_cast<const std::vector<Element const *>*>(&elemans_)) // need to be init()ed later
{
	srand(input_.randSeed());
	inletSolverPrs_= 1.0;   outletSolverPrs_= 0.0;
	minCyclePc_= 0.0;       maxCyclePc_= 0.0;
	Sw_= 1.0;               cappPressCh_= 0.0;
	flowVolume_= 0.0;       clayVolume_= 0.0;
	maxOilFlowErr_= 0.0;    maxWatFlowErr_= 0.0;  maxResIdxErr_= 0.0;
	//StableFilling_= true;
	solver_= NULL;
	amottOilIdx_= 0.0;       amottWaterIdx_= 1.0;
	singlePhaseDprs_= 1.0;   singlePhaseDvolt_= 1.0;
	oilFlowRate_= 0.0;       watFlowRate_= 0.0;  current_= 0.0;
	writeWatMatrix_= false;  writeOilMatrix_= false;    writeResMatrix_= false;
	writeWatVelocity_= false;  writeOilVelocity_= false;  writeResVelocity_= false;
	writeSlvMatrixAsMatlab_= false;
	amottDataDrainage_.resize(3);   amottDataImbibition_.resize(3);
	deltaPo_= 0.0;  deltaPw_= 0.0;

	if(debugLevel) outD.open(title_ + "_pnflow.dbg");


	input.echoKeywords(out_.fileStream());

	cout<<  "Initializing network:"  << endl;



	bool createLocData(false), writeMatInitOnly(false), verboseSlvr(false);

	input_.calcBox(solverBoxStart_, solverBoxEnd_);//input.calcBox(satBoxStart_, satBoxEnd_);


	//input_.satConvergence(minNumFillings_, initStepSize_, extrapCutBack_, maxFillIncrease_, StableFilling_);
	double  eps(0.0);	double scaleFact(1.0e18), condCutOff(0.0); int slvrOutput(0);
	input_.solverTune(eps, scaleFact, slvrOutput, verboseSlvr, condCutOff);

	input_.prsBdrs(useAvrPrsAsBdr_, prtPressureProfile_, numPressurePlanes_);
	reportMaterialBal_= input.lookupOr("MAT_BAL", false);

	sourceNode_= input.lookupOr("POINT_SOURCE", 0);

	input_.solverDebug(writeWatMatrix_, writeOilMatrix_, writeResMatrix_, writeWatVelocity_, writeOilVelocity_, writeResVelocity_, writeSlvMatrixAsMatlab_, writeMatInitOnly);

	input_.fillingList(writeDrainList_, writeImbList_, createLocData);
	input_.prsDiff(inletSolverPrs_, outletSolverPrs_, includeGravityInRelPerm_);
	deltaPw_= deltaPo_= inletSolverPrs_-outletSolverPrs_;

	Element::set_useGravInKr(includeGravityInRelPerm_);
	Element::conductanceCutOff(condCutOff);




	nBSs_=2;// not properly implemented yet
	coreFrac_=1.0;
	X0_=dbl3(0.0, 0.0, 0.0);
	string fileName;
	initNetworkOld(input_);


	double medianLenToRadRatio;
	{
		vector< double > LToRRatio(nTrots_);
		for(int ti = 0; ti < nTrots_; ++ti)
		{	Throat& tr = *static_cast<Throat*>(elemans_[ti+nBpPors_]);  LToRRatio[ti] = (tr.poreLength(0)+tr.poreLength(1)+tr.throatLength())/tr.RRR();  }
		sort(LToRRatio.begin(), LToRRatio.end());
		medianLenToRadRatio = LToRRatio[nTrots_/2];
	}

	//for(size_t ei = 0; ei < nBSs_; ++ei) dynamic_cast<InOutBoundary*>(elemans_[ei])->prepare2();

	satBoxVolume_= box_.x * box_.y * box_.z * (solverBoxEnd_ - solverBoxStart_);

	modifyConnNum_removeConnectionsFromNetwork(input_);  // Reduce connection number

	//elemans_[0]->identifyConnectedPoreElems();    // Possibly remove singlets
	for(auto elm:elemans_)  elm->prepare2();
	//for(auto elm:elemans_)  elm->setConnectedToNetwork(true); // to allow boundary element removal
	//for(int ie=2; ie<nBSs_; ++ie)  elemans_[ie]->removeFromNetwork();
	for(auto elm:elemans_)  elm->setConnectedToNetwork(false); // to recompute connectivity
	elemans_[0]->identifyConnectedPoreElems();    // Possibly remove singlets

	//for(int ie=2; ie<nBSs_; ++ie) elemans_[ie]->setConnectedToNetwork(false);
	//nBSs_=2; 

	bool inOrOutThroats(false);
	if(input.getVar(inOrOutThroats, "convertInOutThroatsToMicroPorosity"))
	{	int converted = 0;
		for(auto& el:elemans_)
			if(el->connectedToNetwork() && el->convertToMicroPorosityForSven(inOrOutThroats)) ++converted;
		out_ << endl << converted << " elements are converted to micro-porosity type 1\n"<<endl;
	}


	int numSingletsRemoved(0);
	if( ! input.lookupOr("DRAIN_SINGLETS",true) )
	{
		for(int pr = nBSs_; pr < nBpPors_; ++pr)
			if(elemans_[pr]->connectedToNetwork() && elemans_[pr]->connectionNum() == 1)
				numSingletsRemoved += elemans_[pr]->removeFromNetwork();
	}


	int numIsolatedElems=0,  sumNConnection=0;
	for(int pr = nBSs_; pr < nBpPors_; ++pr) sumNConnection+=elemans_[pr]->connectionNum();
	for(auto elm:elemans_)
	{
		if(elm->connectedToNetwork())
		{
			if(elm->isOnInletSlvrBdr())  krInletBoundary_.push_back(elm);
			if(elm->isOnOutletSlvrBdr()) krOutletBoundary_.push_back(elm);
			elm->sortConnectingElems_DistToExit();
			if(elm->isInCalcBox()) {
				flowVolume_ += elm->flowVolume();
				clayVolume_ += elm->clayVolume();
			}
		} 
		else ++numIsolatedElems;
	}



	if(useAvrPrsAsBdr_ && numPressurePlanes_ < 2) numPressurePlanes_= 2;
	pressurePlanes_.resize(numPressurePlanes_);

	for(int p = 0; p < numPressurePlanes_; ++p)
	{
		double loc = solverBoxStart_ + p * (solverBoxEnd_ - solverBoxStart_) / (numPressurePlanes_ - 1);
		pressurePlanesLoc_.push_back(loc);
		for(size_t t = nBpPors_; t < elemans_.size(); ++t)
		{
			if(elemans_[t]->crossesPlaneAt(loc*box_.x) && elemans_[t]->connectedToNetwork())
				pressurePlanes_[p].push_back(elemans_[t]);
		}
	}
	comn_.setNumElem(nPors_, nTrots_);

	modifyNetwork(input_);
	results3D_.init(nBSs_,nBpPors_);

	NetworkTransform(input_.randSeed()+rand()/RAND_MAX*1000, out_).modify(input_, elemans_,nPors_, box_.x*box_.y*box_.z);



	flowVolume_= 0.0;
	clayVolume_= 0.0;
	for(size_t i = nBSs_; i < elemans_.size(); ++i)
		if(elemans_[i]->isInCalcBox())
		{
			flowVolume_ += elemans_[i]->flowVolume();
			clayVolume_ += elemans_[i]->clayVolume();
		}


	//set contact angles
	//vector<Element*> pores2Set;pores2Set.reserve(nBpPors_);//(elemans_.begin(), elemans_.begin()+nBpPors_);
	//vector<Element*> throats2Set;throats2Set.reserve(elemans_.size()-nBpPors_);//(elemans_.begin()+nBpPors_, elemans_.end());

	//for(int el = nBSs_; el < nBpPors_; ++el)
		//if(elemans_[el]->rockIndex() == 0)
		//{
			//pores2Set.push_back(elemans_[el]);			//absVolume += elemans_[el]->flowVolume();
		//}
	//for(size_t el = nBpPors_; el < elemans_.size(); ++el)
		//if(elemans_[el]->rockIndex() == 0)
		//{
			//throats2Set.push_back(elemans_[el]);			//absVolume += elemans_[el]->flowVolume();
		//}
	//{
		//int	   CAMdl = 1;
		//string	angDistScheme;
		//double minAng(0.0), maxAng(0.0),  deltaExp, etaExp, modelTwoSepAng(0.0);
		//setContactAngles(pores2Set, throats2Set, minAng, maxAng, deltaExp, etaExp, CAMdl, modelTwoSepAng,  angDistScheme);
	//}

	SolveSinglePhaseFlow(input_); /// A single phase solve is conducted to have flowrates to scale relperms against.



	writeNetworkToFile(input);

	results3D_.write3D(Pc_, comn_.sigmaOW());
	if(input.lookupOr("writeStatistics", false))
	{
		const std::vector<Element const *>* rockLatticeConst = reinterpret_cast<const std::vector<Element const *>*>(&elemans_);
		printRadiusStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printShapeFactorStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printCornerNumStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printCornerAngStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printAspectRatioStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printCoordinaNumStatistics(*rockLatticeConst, nBSs_,nBpPors_);
		printDistanceMapStatistics(*rockLatticeConst, nBSs_,nBpPors_);
	}



	if(createLocData) createMatlabLocationData();

	double permeability = (singlePhaseWaterQ_ * water_.viscosity() * box_.x * (solverBoxEnd_-solverBoxStart_))
		/ (box_.y * box_.z * (inletSolverPrs_ - outletSolverPrs_));
	double formationfactor = ((box_.y * box_.z) * (inletSolverPrs_ - outletSolverPrs_))
		/ (water_.resistivity() * (singlePhaseCurrent_+1.0e-200) * box_.x * (solverBoxEnd_-solverBoxStart_));

	out_<< "\n" << endl
		<< "\n Number of pores:                " << nPors_
		<< "\n Number of throats:              " << nTrots_ << endl
		<<   " Average connection number:      " << (double)(sumNConnection-nPors_)/nPors_
		<< "\n Number of inlet connections:    " << elemans_[0]->connectionNum()
		<< "\n Number of outlet connections:   " << elemans_[OutI]->connectionNum()
		<< "\n Number of isolated elements:    " << numIsolatedElems
		<< "\n Number of singlets removed:     " << numSingletsRemoved
		<< "\n Number of triangular elements:  " << comn_.numTriangles()
		<< "\n Number of square     elements:  " << comn_.numSquares()
		<< "\n Number of circular   elements:  " << comn_.numCircles()
		<< "\n Median of throat length/radius: " << medianLenToRadRatio
		<< "\n Clay bound porosity:            " << clayVolume_ / satBoxVolume_
		<< "\n Total porosity:        " << (flowVolume_+clayVolume_) / satBoxVolume_
		<< "\n Absolute permeability: " << permeability<<" (m2), "<<permeability/9.869233E-16<<" (mD) "
		<< "\n Formation factor:      " << formationfactor     << endl
		<< endl;

	 //coreFrac_=1.0; // TODO implement
	comn_.KcSwKrsQsss_.reserve(10);
	comn_.KcSwKrsQsss_.push_back(vector<vector<double> >(1,vector<double>(8,0.0))); //-Note to sync: n_KcSwKrsQsss=8
	comn_.KcSwKrsQsss_[0][0][0]=1.0;//comn_.water().interfacialTen(); ///. needed for post-processing, valid cycle numbers start from 1 //	vector<double>KcSwKrsQs(8,0.0)
	comn_.KcSwKrsQsss_[0][0][1]=(flowVolume_+clayVolume_) / satBoxVolume_ ;
	comn_.KcSwKrsQsss_[0][0][2]=permeability/coreFrac_ ;
	comn_.KcSwKrsQsss_[0][0][3]=permeability/coreFrac_ ;
	comn_.KcSwKrsQsss_[0][0][4]=formationfactor*coreFrac_ ;

	singlePhaseWaterQ_=max(singlePhaseWaterQ_,1.0e-64);
	singlePhaseOilQ_= singlePhaseWaterQ_ * (water_.viscosity() / oil_.viscosity());

	writeResultData(true,true);




		KcIncr_=1000.0; KcIncrFactor_=.2; minNumFillings_=nBpPors_/200+1, initStepSize_=0.05; extrapCutBack_=0.8; maxFillIncrease_=5.0;
		istringstream data;
		if(input.getData(data,"SAT_CONVERGENCE"))
		{
			data >> minPcStep_ >> KcIncr_ >> KcIncrFactor_ >> minNumFillings_ >> initStepSize_ >> extrapCutBack_ >> maxFillIncrease_;
			//KcIncr_ /= comn_.watoil().sigma();
			stableFilling_ = readBoolOr("T",data);
			cout<<"SAT_CONVERGENCE:  minPcStep PcIncr  PcIncrFactor  minVolFillings  initStepSize  extrapCutBack  maxFillIncrease  StableFilling"<<endl;
			cout<<"SAT_CONVERGENCE: "<< minPcStep_ <<" "<< KcIncr_ <<" "<< KcIncrFactor_ <<" "<< minNumFillings_ <<" "<< initStepSize_ <<" "<< extrapCutBack_ <<" "<< maxFillIncrease_<<" "<<(stableFilling_?'T':'F')<<endl;
			//input.checkEndOfData(data,"SAT_CONVERGENCE");
		}





	Apex::nErrors = 0;
}



/// Destructor
FlowDomain::~FlowDomain()
{
	for(size_t i = 0; i < elemans_.size(); ++i)		delete elemans_[i];
}




void FlowDomain::modifyConnNum_removeConnectionsFromNetwork(InputData& input)
{
	string modelToChangeCN;
	double targetConnNum(-1.0);
	input.modifyConnNum(targetConnNum, modelToChangeCN);

	if(targetConnNum <= 0.0)   return;



	int totNumConn(0);
	for(int i = 2; i < nBpPors_; ++i)  totNumConn += elemans_[i]->connectionNum();

	if(targetConnNum > double(totNumConn)/nPors_)
	{
		out_ << endl
			<< "==========================================================="	<< endl
			<< "Warning: The requested connection number is higher than the"	<< endl
			<< "original. Connections can only be removed, not added."		  << endl
			<< "==========================================================="	<< endl
			<< endl;
	}

	vector<Element*> throats(elemans_.begin()+nBpPors_, elemans_.end());

	if(modelToChangeCN[1] == 'o' || modelToChangeCN[1] == 'O')	  // Net Volume 'volume'
		sort(throats.begin(), throats.end(), ElemVolCmpRed());
	else if(modelToChangeCN[1] == 'h' || modelToChangeCN[1] == 'H') // ElemModel Factor 'shape'
		sort(throats.begin(), throats.end(), ElemGCmpRed());
	else if(modelToChangeCN[2] == 'n' || modelToChangeCN[2] == 'N') // Random 'rand'
		sort(throats.begin(), throats.end(), ElemGCmpRed());
	else															// Radius 'radius'
		sort(throats.begin(), throats.end(), ElemRadCmpRed());

	int numToRemove((totNumConn-targetConnNum*nPors_)/2);
	while(numToRemove > 0 && !throats.empty())
	{
		Element *damned = throats.back();
		throats.pop_back();

		if(!damned->connectedToEntryOrExit())
		{
			Element *poreOne = damned->connection(0);
			Element *poreTwo = damned->connection(1);

			poreOne->severConnection(damned);
			poreTwo->severConnection(damned);

			ensure(elemans_[damned->index()] == damned);
			elemans_[damned->index()] = 0;
			delete damned;
			--numToRemove;
			--nTrots_;
		}
	}
	Element *ghost = NULL;
	elemans_.erase(remove(elemans_.begin()+nBpPors_, elemans_.end(), ghost),
		elemans_.end());
	ensure(static_cast< int >(elemans_.size()) == nBpPors_+nTrots_);

	for(int neoTIdx = nBpPors_; neoTIdx < static_cast< int >(elemans_.size()); ++neoTIdx)
		elemans_[neoTIdx]->updateLatticeIndex(neoTIdx);
}

void FlowDomain::modifyNetwork(InputData& input)
{

	vector<Element*> pores(elemans_.begin()+nBSs_, elemans_.begin()+nBpPors_);
	vector<Element*> throats(elemans_.begin()+nBpPors_, elemans_.end());

	///. modify shape factors
	{
		int throatGModel(0), poreGModel(0);
		int numPtsGDist(0);
		string throatGOptions, poreGOptions;
		bool  writeGDistToFile(false);
		input.getModifyGDist(throatGModel, poreGModel, throatGOptions, poreGOptions, writeGDistToFile,
			numPtsGDist);
		if(poreGModel == -1 || throatGModel == -1)
		{
			string fileName("ShapeFactDist.csv");
			vector<Element*> elems(elemans_.begin()+nBSs_, elemans_.begin()+nBpPors_);
			elems.insert(elems.end(), elemans_.begin()+nBpPors_, elemans_.end());
			if(throatGModel == -1)
				modifyShapeFactor(poreGModel, poreGOptions, elems, fileName, writeGDistToFile, numPtsGDist);
			else
				modifyShapeFactor(throatGModel, throatGOptions, elems, fileName, writeGDistToFile, numPtsGDist);
		}
		else
		{
			string fileNamePores("ShapeFactDist_pores.csv"), fileNameThroats("ShapeFactDist_throats.csv");
			modifyShapeFactor(throatGModel, throatGOptions, throats, fileNameThroats, writeGDistToFile, numPtsGDist);
			modifyShapeFactor(poreGModel, poreGOptions, pores, fileNamePores, writeGDistToFile, numPtsGDist);
		}
	}

  ///. scale
	{
		double scaleFactor(-1.0);
		input.getModifyModelSize(scaleFactor);
		if(scaleFactor > 0.0)
		{
			for(size_t elm = 0; elm < elemans_.size(); ++elm)
			{
				double oldRad(elemans_[elm]->model()->RRR());
				elemans_[elm]->ChModel()->setRadius(oldRad*scaleFactor);
			}
		}
	}


	double netPoroTrgt(-1.0), clayPoroTrgt(-1.0);


	///. MODIFY_RAD_DIST: modify radius and areas, skip changing volumes for now, they
	///. will be updated below, or overwritten based on keyword MODIFY_PORO
	{
		double oldMedianLenToRadRatio;
		{
			vector< double > LToRRatio(nTrots_);
			for(int ti = 0; ti < nTrots_; ++ti)
			{	Throat& tr = *static_cast<Throat*>(elemans_[ti+nBpPors_]);  LToRRatio[ti] = (tr.poreLength(0)+tr.poreLength(1)+tr.throatLength())/tr.RRR();  }
			sort(LToRRatio.begin(), LToRRatio.end());
			oldMedianLenToRadRatio = LToRRatio[nTrots_/2];
		}

		string throatRadOptions, poreRadOptions;
		bool writeRDistToFile(false), maintainLtoR(false);
		int  throatRadModel(0), poreRadModel(0), numPtsRDist(0);
		input.getModifyRadDistOptions(throatRadModel, poreRadModel, throatRadOptions, poreRadOptions,
			maintainLtoR, writeRDistToFile, numPtsRDist);
		if(poreRadModel == -1 || throatRadModel == -1)
		{
			string fileName("RadDist.csv");
			vector<Element*> elems(elemans_.begin()+nBSs_, elemans_.begin()+nBpPors_);
			elems.insert(elems.end(), elemans_.begin()+nBpPors_, elemans_.end());
			if(throatRadModel == -1)
				modifyInscribedRadii(poreRadModel, poreRadOptions, elems, fileName, writeRDistToFile, numPtsRDist);
			else
				modifyInscribedRadii(throatRadModel, throatRadOptions, elems, fileName, writeRDistToFile, numPtsRDist);
		}
		else if(throatRadModel != 2)
		{
			string fileNamePores("RadDist_pores.csv"), fileNameThroats("RadDist_throats.csv");
			modifyInscribedRadii(throatRadModel, throatRadOptions, throats, fileNameThroats, writeRDistToFile, numPtsRDist);
			modifyInscribedRadii(poreRadModel, poreRadOptions, pores, fileNamePores, writeRDistToFile, numPtsRDist, true);
		}
		else //if(throatRadModel == 2)
		{
			string fileNamePores("RadDist_pores.csv"), fileNameThroats("RadDist_throats.csv");
			modifyInscribedRadii(poreRadModel, poreRadOptions, pores, fileNamePores, writeRDistToFile, numPtsRDist, true);
			modifyInscribedRadii(throatRadModel, throatRadOptions, throats, fileNameThroats, writeRDistToFile, numPtsRDist);
		}

		double medianLenToRadRatio;
		{
			vector< double > LToRRatio(nTrots_);
			for(int ti = 0; ti < nTrots_; ++ti)
			{	Throat& tr = *static_cast<Throat*>(elemans_[ti+nBpPors_]);  LToRRatio[ti] = (tr.poreLength(0)+tr.poreLength(1)+tr.throatLength())/tr.RRR();  }
			sort(LToRRatio.begin(), LToRRatio.end());
			medianLenToRadRatio = LToRRatio[nTrots_/2];
		}

		if(poreRadModel && maintainLtoR)
		//if(Rad_scaleFactor > 0.0)
		{
			double Rad_scaleFactor = oldMedianLenToRadRatio/medianLenToRadRatio;
			for(int i = 0; i < nBpPors_; ++i)
				elemans_[i]->node()*=Rad_scaleFactor;

			for(size_t j = nBpPors_; j < elemans_.size(); ++j)
				elemans_[j]->modifyLength(Rad_scaleFactor);

			if(netPoroTrgt < 0) netPoroTrgt = flowVolume_/satBoxVolume_;
			if(clayPoroTrgt < 0) clayPoroTrgt = clayVolume_/satBoxVolume_;
			box_.x *= Rad_scaleFactor;
			box_.y *= Rad_scaleFactor;
			box_.z *= Rad_scaleFactor;
			satBoxVolume_= box_.x*box_.y*box_.z*(solverBoxEnd_-solverBoxStart_);
		}
	}



	///. modify volumes and porosity; based on MODIFY_PORO keyword, or based on MODIFY_RAD_DIST if MODIFY_PORO not set,
	{

		input.getModifyPoro(netPoroTrgt, clayPoroTrgt);///. will not change anything if keyword MODIFY_PORO not found

		double targetNetVol = netPoroTrgt >= 0.0 ? satBoxVolume_*netPoroTrgt: flowVolume_;
		double targetClayVol = clayPoroTrgt >= 0.0 ? satBoxVolume_*clayPoroTrgt: clayVolume_;

		if(netPoroTrgt >= 0.0 || clayPoroTrgt >= 0.0)
		{
			double changeNetFrac((targetNetVol-flowVolume_)/flowVolume_);
			double changeClayFrac((targetClayVol-clayVolume_)/clayVolume_);
			double changeClayFromNet((targetClayVol-flowVolume_)/flowVolume_);

			for(size_t i = nBSs_; i < elemans_.size(); ++i)
			{
				double neoNetVol = elemans_[i]->flowVolume()*(1.0+changeNetFrac);

				double neoClayVol(0.0);
				if(clayPoroTrgt >= 0.0 && clayVolume_ > 0.0)
					neoClayVol = elemans_[i]->clayVolume()*(1.0+changeClayFrac);
				else if(clayPoroTrgt >= 0.0 && clayVolume_ < 1.0e-12)
					neoClayVol = elemans_[i]->flowVolume()*(1.0+changeClayFromNet);

				elemans_[i]->adjustVolume(neoNetVol, neoClayVol);
			}
		}
	} // volumes are out of sync here

}



/// Does the first flow through the network
void FlowDomain::SolveSinglePhaseFlow(InputData& input)
{

	out_ << "solving for single-phase flow";cout.flush();
	wantRelPerm_=true; wantResIdx_=true;
	updateSatAndConductances(Pc_); ///.  updateSatAndConductances

	double flwErr(0.0), cpuTmp(0.0), crntErr(0.0);

	if((writeWatVelocity_ || writeOilVelocity_) && writeSlvMatrixAsMatlab_) createMatlabLocationData();

	for(size_t it = nBpPors_; it< elemans_.size(); ++it) { auto tr=static_cast<Throat*>(elemans_[it]);  tr->calcR2(comn_.water());  }


	hypreSolver slvr(elemans_, krInletBoundary_, krOutletBoundary_, nBSs_, nBpPors_, debugLevel, title_ + "_mat_init", writeSlvMatrixAsMatlab_);

	out_<< " "; cout.flush();

	singlePhaseWaterQ_= slvr.flowrate(inletSolverPrs_, outletSolverPrs_, water_, flwErr, cpuTmp, 1.0, writeWatVelocity_, writeWatMatrix_);

	if(useAvrPrsAsBdr_)
	{
		FlowDomain::prsOrVoltDrop(water_, singlePhaseDprs_);
		deltaPw_= singlePhaseDprs_;
	}
	//if(prtPressureProfile_ && wantRelPerm_) recordPrsProfiles(water_);

	watFlowRate_= singlePhaseWaterQ_;
	oilFlowRate_= 0.0;



	out_ << " "; cout.flush();
	singlePhaseCurrent_= slvr.flowrate(inletSolverPrs_, outletSolverPrs_, comn_.elec(), crntErr, cpuTmp, 1.0, writeResVelocity_, writeResMatrix_);
	current_= singlePhaseCurrent_;
	if(useAvrPrsAsBdr_) FlowDomain::prsOrVoltDrop(comn_.elec(), singlePhaseDvolt_);
	out_ << " \n";cout.flush();

	if(fabs(flwErr) > 0.05 || fabs(crntErr) > 0.05)
	{
		out_ << endl
			<< "====================================================="  << endl
			<< "Warning. Large errors were detected when solving for"   << endl
			<< "single phase conditions."							   << endl
			<< "Flow Error: " << flwErr						 << endl
			<< "Resistivity Error: " << crntErr					  << endl
			<< "====================================================="  << endl;
	}





}


/**
 * Gets the pressure drop between two planes in the network, used for kr
 */
bool FlowDomain::prsOrVoltDrop(const Fluid& fluid, double& prsDrop) const
{
	double prsOut(0.0), stdOut(0.0);
	int numOut(0);

	bool btOut = FlowDomain::avrPrsOrVolt(fluid, numPressurePlanes_-1, prsOut, stdOut, numOut);

	double prsIn(0.0), stdIn(0.0);
	int numIn(0);

	bool btIn = FlowDomain::avrPrsOrVolt(fluid, 0, prsIn, stdIn, numIn);

	prsDrop = prsIn - prsOut;
	return btOut && btIn;
}

/**
 * Calculates the average pressure in throats crossing a plane. Assumes the pore
 * pressures to be point distributed. Will also calculate the number of pressure
 * points along with standard deviation
 */
bool FlowDomain::avrPrsOrVolt(const Fluid& fluid, int prsPlane, double& res,  double& stdDev, int& numVal) const
{
	double loc(pressurePlanesLoc_[prsPlane]*box_.x), resSum(0.0), resSumSq(0.0), flowSum(0.0), flowTarget;
	vector< pair< double, double > > pqArray;
	string phase("water");

	if(fluid.ff()==ELEC)
	{
		flowTarget = current_;
		phase = "current";
	}
	else if(fluid.isOil())
	{
		flowTarget = oilFlowRate_;
		phase = "oil";
	}
	else
		flowTarget = watFlowRate_;


	flowSum = 0.0;
	for(size_t i = 0; i < pressurePlanes_[prsPlane].size(); ++i)
	{
		double val(0.0), flowRate(0.0);
		if(pressurePlanes_[prsPlane][i]->prevSolvrRes(fluid, loc, val, flowRate))
		{
			pair< double, double > valPair(val, flowRate);
			pqArray.push_back(valPair);
			flowSum += flowRate;
		}
	}

	double flowError(fabs((flowTarget-flowSum)/flowSum));

	if(fluid.ff()==ELEC && flowTarget > 0.0)          maxResIdxErr_  = max(maxResIdxErr_, flowError);
	else if(fluid.isOil() && flowTarget > 0.0)  maxOilFlowErr_= max(maxOilFlowErr_, flowError);
	else if(flowTarget > 0.0)                    maxWatFlowErr_= max(maxWatFlowErr_, flowError);


	const double  MAX_FLOW_ERR = 0.02;
	if(flowError > MAX_FLOW_ERR && flowTarget > 0.0) cout<< "Large flow error (" << flowError*100.0 << "%) for " << phase << " at " << loc << endl;

	int numPts(static_cast< int >(pqArray.size()));
	numVal = 0;
	for(int j = 0; j < numPts; ++j)
	{
		resSum += pqArray[j].first;
		resSumSq += pqArray[j].first*pqArray[j].first;
		++numVal;
	}

	if(numVal > 0) res = resSum / numVal;
	if(numVal > 1) stdDev = sqrt((numVal*resSumSq - resSum*resSum)/(numVal*(numVal-1)));

	return numVal > 0;
}








/// Do a single calculation of relative permeability
void FlowDomain::solve_forRelPermResIndex(bool wantRelPerm, bool wantResIdx)
{


	vector<double>KcSwKrsQs(8,0.0);
	KcSwKrsQs[0]= Pc_;
	KcSwKrsQs[1]= Sw_;

	if(wantRelPerm)
	{
		double oilErr(0.0), watErr(0.0);

		watFlowRate_= solver_->flowrate(inletSolverPrs_, outletSolverPrs_, water_, watErr, cpuTimeKrw_, Sw_, writeWatVelocity_, writeWatMatrix_);
		double relPermWater = watFlowRate_ /singlePhaseWaterQ_;
		if(watFlowRate_ != 0.0 && useAvrPrsAsBdr_)
		{
			bool chk = FlowDomain::prsOrVoltDrop(water_, deltaPw_);
			ensure(chk);
			relPermWater = (watFlowRate_ * singlePhaseDprs_) / (singlePhaseWaterQ_ * deltaPw_);
		}
		KcSwKrsQs[2]=relPermWater;

		oilFlowRate_= solver_->flowrate(inletSolverPrs_, outletSolverPrs_, oil_, oilErr, cpuTimeKro_, Sw_, writeOilVelocity_, writeOilMatrix_);
		double relPermOil = oilFlowRate_ /singlePhaseOilQ_;
		if(oilFlowRate_ != 0.0 && useAvrPrsAsBdr_)
		{
				bool chk = FlowDomain::prsOrVoltDrop(oil_, deltaPo_);
				ensure(chk);
				relPermOil = (oilFlowRate_ * singlePhaseDprs_) / (singlePhaseOilQ_ * deltaPo_);
		}
		KcSwKrsQs[3]=relPermOil;

		out_ << " krw: " << setw(9) << std::left << relPermWater << " kro: " << setw(9) << std::left << relPermOil<< " ";

		maxOilFlowErr_= max(maxOilFlowErr_, oilErr);
		maxWatFlowErr_= max(maxWatFlowErr_, watErr);

		resultWaterFlowRate_.push_back({watFlowRate_, deltaPw_});
		resultOilFlowRate_.push_back({oilFlowRate_, deltaPo_});
	}

	if(wantResIdx)
	{
		double resErr(0.0);


		current_= solver_->flowrate(inletSolverPrs_, outletSolverPrs_, comn_.elec(), resErr, cpuTimeResIdx_, Sw_, writeResVelocity_, writeResMatrix_);
		double resIndex = singlePhaseCurrent_ / (abs(current_)+1.0e-32);
		if(current_ > 0.0 && useAvrPrsAsBdr_)
		{
				bool chk = FlowDomain::prsOrVoltDrop(comn_.elec(), deltaV_);
				resIndex = singlePhaseCurrent_ / (abs(current_)+1.0e-32) * deltaV_/singlePhaseDvolt_;
				ensure(chk);
		}
		KcSwKrsQs[4]=resIndex;

		resultResistivityIdx_.push_back(resIndex);
		out_ << "RI: " << resIndex;;

		maxResIdxErr_= max(maxResIdxErr_, resErr);
	}



	resultWaterSat_.push_back(Sw_);
	resultCappPress_.push_back(Pc_);

	if(reportMaterialBal_)
	{
		resultWaterMass_.push_back((flowVolume_ + clayVolume_)*Sw_*water_.density());
		resultOilMass_.push_back((flowVolume_ + clayVolume_)*(1 - Sw_)*oil_.density());
	}


	comn_.addKcSwKrsQs(KcSwKrsQs);
	results3D_.write3D(Pc_, comn_.sigmaOW());

}




/// Iterates across all network elements (not in/outlets) and computes the water saturation
void FlowDomain::updateSatAndConductances(double Pc)
{
	double volWater(0.0);
	for(Element* el:elemans_)  volWater += el->updateSat_calcR(Pc);
	if(wantRelPerm_)
	 for(size_t it = nBpPors_; it< elemans_.size(); ++it) { auto tr=static_cast<Throat*>(elemans_[it]);  tr->calcR2(comn_.water());   tr->calcR2(comn_.oil());  }
	if(wantResIdx_)
	 for(size_t it = nBpPors_; it< elemans_.size(); ++it) { auto tr=static_cast<Throat*>(elemans_[it]);  tr->calcR2(comn_.elec());  }
	Sw_= volWater / (flowVolume_ + clayVolume_);
}





string FlowDomain::calcUSBMindex() const
{
	ostringstream out;
	out.flags(ios::showpoint);
	out.flags(ios::fixed);
	out.precision(3);

	size_t numDrainPts(usbmDataDrainage_.size()), numImbPts(usbmDataImbibition_.size());
	double areaDrain(0.0), areaImb(0.0);

	if(numDrainPts > 1)
	{
		for(size_t i = 1; i < numDrainPts; ++i)
		{
			pair< double, double > ptOne(usbmDataDrainage_[i]), ptTwo(usbmDataDrainage_[i-1]); ///. Sw - Pc
			areaDrain += (ptTwo.second-ptOne.second) * ((ptTwo.first+ptOne.first)*0.5);
		}
	}

	if(numImbPts > 1)
	{
		for(size_t i = 1; i < numImbPts; ++i)
		{
			pair< double, double > ptOne(usbmDataImbibition_[i-1]), ptTwo(usbmDataImbibition_[i]);  ///. Sw - Pc
			areaImb += (ptTwo.second-ptOne.second)*(-(ptTwo.first+ptOne.first)*0.5);
		}
	}

	if(areaDrain == 0.0)
		out << "-INF";
	else if(areaImb == 0.0)
		out << "+INF";
	else
		out << log10(areaDrain/areaImb);

	return out.str();
}



