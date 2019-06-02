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
using namespace std;

#include "inputData.h"
#include "Element.h"
#include "solver.h"
#include "netStatistics.h"
#include "netsim.h"
#include "NetworkTransform.h"



const double	Netsim::MAX_FLOW_ERR = 0.02;
const int	   Netsim::DUMMY_IDX = -99;
OnDemandStream    outD;  ///. alias to OnDemandStream::dbgFile  //thread_local

/**
 * Netsim constructor
 */
Netsim::Netsim(InputData & input)
: 	m_input(input),
	m_comn(input),
	m_oil(m_comn.oil()),	
	m_water(m_comn.water()),
	m_cappPress(m_cappPressCh),	
	m_baseFileName(input.title()),
	m_out(m_baseFileName + "_pnflow.prt"),
	m_vtkWriter(input.keywordData("visualize"),this,m_baseFileName+"_res/")

{
	srand(input.randSeed());
	m_inletSolverPrs = 1.0;
	m_outletSolverPrs = 0.0;
	m_minCycleCappPress = 0.0;
	m_maxCycleCappPress = 0.0;
	m_satWater = 1.0;
	m_cappPressCh = 0.0;
	//m_clayAdjust = 0.0;
	m_totalFlowVolume = 0.0;
	m_totalClayVolume = 0.0;
	m_maxNonZeros = 0;
	m_maxOilFlowErr = 0.0;
	m_maxWatFlowErr = 0.0;
	m_maxResIdxErr = 0.0;
	m_numIsolatedElems = 0;
	m_totNumFillings = 0;
	m_keepFraction = 1.0;
	m_wantRelPerm = true;
	m_wantResIdx = true;
	m_StableFilling = true;
	m_injAtLeftRes = false;
	m_injAtRightRes = false;
	m_reportMaterialBal = false;
	m_includeGravityInRelPerm = false;
	m_solver = NULL;
	m_solverBoxStart = 0.0;
	m_amottOilIdx = 0.0;
	m_amottWaterIdx = 1.0;
	m_solverBoxEnd = 1.0;
	m_singlePhaseDprs = 1.0;
	m_singlePhaseDvolt = 1.0;
	m_numPressurePlanes = 0;
	m_oilFlowRate = 0.0;
	m_watFlowRate = 0.0;
	m_current = 0.0;
	m_writeWatMatrix = false;
	m_writeOilMatrix = false;
	m_writeResMatrix = false;
	m_writeWatVelocity = false;
	m_apexPrsReported = false;
	m_writeOilVelocity = false;
	m_writeResVelocity = false;
	m_writeSlvMatrixAsMatlab = false;
	m_amottDataDrainage.resize(3);
	m_amottDataImbibition.resize(3);
	m_deltaPo = 0.0;
	m_deltaPw = 0.0;
	if(m_comn.debugMode) outD.open(m_baseFileName + "_pnflow.dbg");

	init(input);


	Apex::debugMode = m_comn.debugMode;
	Apex::nErrors = 0;
}

/**
 * Copy constructor takes a copy of another netsim instance. This might be useful in
 * say simulations where one would use multiple instances of the network model in
 * different blocks.
 */

/**
 * Destructor
 */
Netsim::~Netsim()
{
	for(size_t i = 0; i < m_rockLattice.size(); ++i)		delete m_rockLattice[i];
}

/**
 * The network is created and initialized by reading in the connection files that also contains
 * all other relevant information about the pores/throats. Pointers to all elements are contained
 * in a vector. Element 0 and (n_pores + 1) are inlet/outlet. Throats follows after the pores.
 * Since all rock elements contain pointers to connecting elements rather than vector indicies,
 * initialization has to be in correct order: throats, pores, in/outlet and finally finishing the
 * throats. A single phase solve is also conducted to have flowrates to scale relperms against.
 */
void Netsim::init(InputData& input)
{
	input.echoKeywords(m_out.fileStream());
	
	cout<<  "Initializing network:"  << endl;



	double /*gravX(0.0), gravY(0.0), gravZ(-9.81),*/ eps(0.0);
	double scaleFact(1.0e18); int slvrOutput(0);
	vector< string > solverOptions;
	string modPoroOptions;
	bool createLocData(false), writeMatInitOnly(false), verboseSlvr(false);
	double condCutOff(0.0);

	input.calcBox(m_satBoxStart, m_satBoxEnd);
	input.network(m_numPores, m_numThroats, m_xSize, m_ySize, m_zSize);
	

	input.satConvergence(m_minNumFillings, m_initStepSize, m_extrapCutBack, m_maxFillIncrease, m_StableFilling);
	input.solverTune(eps, scaleFact, slvrOutput, verboseSlvr, condCutOff);
	//m_clayAdjust = input.clayEdit();
	//input.resFormat(m_matlabFormat, m_excelFormat);
	input.prsBdrs(m_useAvrPrsAsBdr, m_prtPressureProfile, m_numPressurePlanes);
	input.matBal(m_reportMaterialBal);
	input.apexPrs(m_apexPrsReported);
	input.aCloseShave(m_keepFraction);

	//input.trapping(drainSinglets, KrwatcornAtSw0);
	input.sourceNode(m_sourceNode);

		//TrappingCriteria trpCrit(escapeToEither);	   // Anything connected to the entry reservoir is
		//if(m_injAtLeftRes && !m_injAtRightRes)		  // by default connected. ie we need less restrictive
			//trpCrit = escapeToInlet;					// trapping criteria.
		//else if(!m_injAtLeftRes && m_injAtRightRes)
			//trpCrit = escapeToOutlet;

	input.solverDebug(m_writeWatMatrix, m_writeOilMatrix, m_writeResMatrix, m_writeWatVelocity, m_writeOilVelocity,
		m_writeResVelocity, m_matrixFileName, m_writeSlvMatrixAsMatlab, writeMatInitOnly);
	if(writeMatInitOnly)
	{
		m_writeWatMatrix = false;
		m_writeOilMatrix = false;
		m_writeResMatrix = false;
		m_writeWatVelocity = false;
		m_writeOilVelocity = false;
		m_writeResVelocity = false;
	}

	input.fillingList(m_writeDrainList, m_writeImbList, createLocData);
	input.prsDiff(m_inletSolverPrs, m_outletSolverPrs, m_includeGravityInRelPerm);
	m_deltaPw = m_deltaPo = m_inletSolverPrs-m_outletSolverPrs;

	Element::set_useGravInKr(m_includeGravityInRelPerm);
	Element::conductanceCutOff(condCutOff);


	if(!m_useAvrPrsAsBdr)						   // If we don't average pressures to obtain
	{											   // rel perm (ie we're moving the boundary)
		m_solverBoxStart = m_satBoxStart;		   // then the lattice we pass to the solver
		m_solverBoxEnd = m_satBoxEnd;			   // will be the same as that we compute sat across
	}

	initNetwork(input);

	m_comn.setNumElem(m_numPores, m_numThroats);
	modifyNetwork(input);

	NetworkTransform(input.randSeed()+rand()/RAND_MAX*1000, m_out).modify(input, m_rockLattice,m_numPores, m_xSize*m_ySize*m_zSize);


	{  // m_elemans is not used much yet, it will  replace m_rockLattice
		m_elemans.resize(m_rockLattice.size());
		                              m_elemans[0]=m_rockLattice[0];    m_elemans[1]=m_rockLattice[1+m_numPores];
		for(int el = 1; el < 1+m_numPores; ++el)		                m_elemans[el+1]=m_rockLattice[el];
		for(size_t el = 2+m_numPores; el < m_rockLattice.size(); ++el)	m_elemans[el]=m_rockLattice[el];
	}

	m_totalFlowVolume = 0.0;
	m_totalClayVolume = 0.0;
	for(size_t i = 0; i < m_rockLattice.size(); ++i)
	{
		if(m_rockLattice[i]->isInsideSatBox())
		{
			m_totalFlowVolume += m_rockLattice[i]->flowVolume();
			m_totalClayVolume += m_rockLattice[i]->clayVolume();
		}
	}

	updateSatAndConductances(m_cappPress);

	//set contact angles
	vector<Element*> pores2Set;pores2Set.reserve(m_numPores+2);//(m_rockLattice.begin(), m_rockLattice.begin()+2+m_numPores);
	vector<Element*> throats2Set;throats2Set.reserve(m_rockLattice.size()-m_numPores+2);//(m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());
	
	for(int el = 0; el < 2+m_numPores; ++el)
		if(m_rockLattice[el]->iRockType() == 0)
		{
			pores2Set.push_back(m_rockLattice[el]);
			//absVolume += m_rockLattice[el]->flowVolume();
		}

	for(size_t el = 2+m_numPores; el < m_rockLattice.size(); ++el)
		if(m_rockLattice[el]->iRockType() == 0)
		{
			throats2Set.push_back(m_rockLattice[el]);
			//absVolume += m_rockLattice[el]->flowVolume();
		}

	int	   wettingClass = 1;
	string	angDistScheme;
	double minEqConAng(0.0), maxEqConAng(0.0),  wettDelta, wettEta, modelTwoSepAng(0.0);
	input.initConAng(wettingClass, minEqConAng, maxEqConAng, wettDelta, wettEta, angDistScheme, modelTwoSepAng);
	setContactAngles(pores2Set, throats2Set, minEqConAng, maxEqConAng, wettDelta, wettEta, wettingClass, modelTwoSepAng,  angDistScheme);


	useHypre=true;
	input.getVar(useHypre,"useHypre");
	if(useHypre) m_out << "using hypre solver"<<endl;
	else 	m_out << "using netlib amg-1.5 (1995) solver"<<endl;
	//if(useHypre)
		//amg_solver::initSolver(eps, scaleFact, slvrOutput, verboseSlvr, m_includeGravityInRelPerm);

	SolveSinglePhaseFlow(input);



	writeNetworkToFile(input);
	
//	m_vtkWriter.vtuWrite(&m_elemans, m_numPores, m_cappPress, m_oil.interfacialTen());
	m_vtkWriter.vtuWrite(reinterpret_cast<const std::vector<Element const *>*>(&m_rockLattice), m_numPores, m_cappPress, m_oil.interfacialTen());
	const std::vector<Element const *>* rockLatticeConst = reinterpret_cast<const std::vector<Element const *>*>(&m_rockLattice);
	printCornerNumStatistics(rockLatticeConst, m_numPores);
	printCornerAngStatistics(rockLatticeConst, m_numPores);
	printShapeFactorStatistics(rockLatticeConst, m_numPores);
	printRadiusStatistics(rockLatticeConst, m_numPores);
	printAspectRatioStatistics(rockLatticeConst, m_numPores);
	printCoordinaNumStatistics(rockLatticeConst, m_numPores);
	printDistanceMapStatistics(*rockLatticeConst, m_numPores);



	if(createLocData) createMatlabLocationData();
	
	double absPermeability = (m_singlePhaseWaterQ * m_water.viscosity() * m_xSize * (m_solverBoxEnd-m_solverBoxStart))
		/ (m_ySize * m_zSize * 9.869233E-16 * (m_inletSolverPrs - m_outletSolverPrs));
	double formationFactor = ((m_ySize * m_zSize) * (m_inletSolverPrs - m_outletSolverPrs))
		/ (m_water.resistivity() * (m_singlePhaseCurrent+1.0e-200) * m_xSize * (m_solverBoxEnd-m_solverBoxStart));
		
	m_out << endl << endl
		<<  " Number of pores:						   " << m_numPores									 << endl
		<<  " Number of throats:						 " << m_numThroats								   << endl
		<<  " Average connection number:				 " << (double)(m_maxNonZeros-m_numPores)/m_numPores  << endl
		<<  " Number of connections to inlet:			" << m_rockLattice[0]->connectionNum()			  << endl
		<<  " Number of connections to outlet:		   " << m_rockLattice[m_numPores+1]->connectionNum()   << endl
		<<  " Number of physically isolated elements:	" << m_numIsolatedElems							 << endl
		//<<  " Number of singlets removed:				" << numSingletsRemoved							 << endl
		<<  " Number of triangular shaped elements:	  " << m_comn.numTriangles()			   << endl
		<<  " Number of square shaped elements:		  " << m_comn.numSquares()				 << endl
		<<  " Number of circular shaped elements:		" << m_comn.numCircles()				 << endl
		<<  " Median throat length to radius ratio:	  " << m_rockLattice[m_numPores+m_numThroats/2]->lenToRadRatio()						  << endl
		<<  " Net porosity:							  " << m_totalFlowVolume / m_satBoxVolume					 << endl
		<<  " Clay bound porosity:					   " << m_totalClayVolume / m_satBoxVolume					<< endl
		<<  " Absolute permeability (mD):				" << absPermeability							  << endl
		<<  " Absolute permeability (m2):				" << absPermeability * 9.869233E-16			   << endl
		<<  " Formation factor:						  " << formationFactor							  << endl
		<< endl;
		
}

/**
 * Creates all the pores and throats
 */
void Netsim::initNetwork(InputData& input)
{


	vector< pair< int, double > > insidePoreHashs(m_numPores);
	int newNumPores = setupInsidePoreHashings(input, insidePoreHashs);
	m_rockLattice.resize(newNumPores + 2);
	
	vector< int > insideThroatHashs(m_numThroats, DUMMY_IDX);
	vector< pair<int, int> > newTrotNeibors(m_numThroats);
	vector<Element*> throatsToInlet;								 // All throats connected to the in/outlet are
	vector<Element*> throatsToOutlet;								// recorded while crating the throats

	int newNumThroats =
   readAndCreateThroats(input, newTrotNeibors, throatsToInlet, throatsToOutlet, insidePoreHashs, insideThroatHashs, newNumPores);
   readAndCreatePores(input, /*newTrotNeibors, insidePoreHashs,*/ insideThroatHashs, newNumPores);
	input.clearNetworkData(); //functions above create a lot of additional working data not needed anymore, lets clean them to free the memory

	m_xSize *= m_keepFraction;
	m_numPores = newNumPores;
	m_satBoxVolume = m_xSize * m_ySize * m_zSize * (m_satBoxEnd - m_satBoxStart);

	double yMid(m_ySize/2.0), zMid(m_zSize/2.0);
	createInAndOutletPore(0, -1.0E-15, yMid, zMid, throatsToInlet);
	createInAndOutletPore(m_numPores + 1, m_xSize+1.0E-15, yMid, zMid, throatsToOutlet);

	vector< double > lenToRadRatio(newNumThroats);
	int runIdx(0);
	for(int i = 0; i < m_numThroats; ++i)						// Adding the pore pointers to the throats had
	{															// to be delayed until after having created
		int poreIndex1 = reIndex(newTrotNeibors[i].first);	  // the pores. The network should now be
		int poreIndex2 = reIndex(newTrotNeibors[i].second);	 // properly initialized.

		if(poreIndex1 >= 0 && poreIndex2 >= 0)
		{
			ensure(poreIndex1 < m_numPores+2 && poreIndex2 < m_numPores+2);
			Element* pore1 = m_rockLattice[poreIndex1];
			Element* pore2 = m_rockLattice[poreIndex2];
			ensure(pore1 != NULL && pore2 != NULL);

			m_rockLattice[m_numPores + 2 + runIdx]->addConnections(pore1, pore2, m_xSize*m_solverBoxStart,
				m_xSize*m_solverBoxEnd, !m_useAvrPrsAsBdr);
			lenToRadRatio[runIdx] = m_rockLattice[m_numPores + 2 + runIdx]->lenToRadRatio();
			++runIdx;
		}
	}

	m_numThroats = newNumThroats;
	sort(lenToRadRatio.begin(), lenToRadRatio.end());
	//m_medianLenToRadRatio = lenToRadRatio[m_numThroats/2];





	modifyConnNum_removeConnectionsFromNetwork(input);								// Reduce connection number

	m_rockLattice[0]->identifyConnectedPoreElems();					 // Possibly remove singlets


	bool inOrOutThroats(false);
	if(input.getVar(inOrOutThroats, "convertInOutThroatsToMicroPorosity"))
	{
		int converted = 0;
		for(size_t elem = 0; elem < m_rockLattice.size(); ++elem)
		{
			if(m_rockLattice[elem]->connectedToNetwork())
			{
				if(m_rockLattice[elem]->convertToMicroPorosityForSven(inOrOutThroats))
				++converted;
			}
		}
		m_out << endl << converted << " elements are converted to micro-porosity type 1"<<endl<<endl;
	}


	int numSingletsRemoved(0);
	bool drainSinglets(true);
	input.getVar(drainSinglets, "DRAIN_SINGLETS");
	if(!drainSinglets)
	{
		for(int pr = 1; pr <= m_numPores; ++pr)
		{
			if(m_rockLattice[pr]->connectedToNetwork() && m_rockLattice[pr]->connectionNum() == 1)
				numSingletsRemoved += m_rockLattice[pr]->removeFromNetwork();
		}
	}
	
	m_out <<  " Number of singlets removed:				" << numSingletsRemoved							 << endl;

	for(size_t elem = 0; elem < m_rockLattice.size(); ++elem)
	{
		ensure(m_rockLattice[elem] != NULL);
		m_rockLattice[elem]->calcVolume_CheckIntegrity(m_totalFlowVolume, m_totalClayVolume, m_maxNonZeros, m_numIsolatedElems);
		if(m_rockLattice[elem]->connectedToNetwork())
		{
			if(m_rockLattice[elem]->isOnInletSlvrBdr()) m_krInletBoundary.push_back(m_rockLattice[elem]);
			if(m_rockLattice[elem]->isOnOutletSlvrBdr()) m_krOutletBoundary.push_back(m_rockLattice[elem]);
			m_rockLattice[elem]->sortConnectingElems_DistToExit();
		}
	}



	if(m_useAvrPrsAsBdr && m_numPressurePlanes < 2) m_numPressurePlanes = 2;
	m_pressurePlanes.resize(m_numPressurePlanes);

	for(int p = 0; p < m_numPressurePlanes; ++p)
	{
		double loc = m_satBoxStart + p * (m_satBoxEnd - m_satBoxStart) / (m_numPressurePlanes - 1);
		m_pressurePlanesLoc.push_back(loc);
		for(size_t t = m_numPores + 2; t < m_rockLattice.size(); ++t)
		{
			if(m_rockLattice[t]->crossesPlaneAt(loc*m_xSize) && m_rockLattice[t]->connectedToNetwork())
				m_pressurePlanes[p].push_back(m_rockLattice[t]);
		}
	}
}


void Netsim::modifyConnNum_removeConnectionsFromNetwork(InputData& input)
{
	string modelToChangeCN;
	double targetConnNum(-1.0);
	int totNumConn(0);

	input.modifyConnNum(targetConnNum, modelToChangeCN);

	for(int i = 1; i <= m_numPores; ++i)
		totNumConn += m_rockLattice[i]->connectionNum();

	if(targetConnNum <= 0.0)
		return;
	else if(targetConnNum > double(totNumConn)/m_numPores)
	{
		ostringstream out;
		out << endl
			<< "==========================================================="	<< endl
			<< "Warning: The requested connection number is higher than the"	<< endl
			<< "original. Connections can only be removed, not added."		  << endl
			<< "==========================================================="	<< endl
			<< endl;
		writePrtData(out);
	}

	vector<Element*> throats(m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());

	if(modelToChangeCN[1] == 'o' || modelToChangeCN[1] == 'O')	  // Net Volume 'volume'
		sort(throats.begin(), throats.end(), ElemVolCmpRed());
	else if(modelToChangeCN[1] == 'h' || modelToChangeCN[1] == 'H') // ElemModel Factor 'shape'
		sort(throats.begin(), throats.end(), ElemGCmpRed());
	else if(modelToChangeCN[2] == 'n' || modelToChangeCN[2] == 'N') // Random 'rand'
		sort(throats.begin(), throats.end(), ElemGCmpRed());
	else															// Radius 'radius'
		sort(throats.begin(), throats.end(), ElemRadCmpRed());

	int numToRemove((totNumConn-targetConnNum*m_numPores)/2);
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

			ensure(m_rockLattice[damned->latticeIndex()] == damned);
			m_rockLattice[damned->latticeIndex()] = 0;
			delete damned;
			--numToRemove;
			--m_numThroats;
		}
	}
	Element *ghost = NULL;
	m_rockLattice.erase(remove(m_rockLattice.begin()+m_numPores+2, m_rockLattice.end(), ghost),
		m_rockLattice.end());
	ensure(static_cast< int >(m_rockLattice.size()) == m_numPores+m_numThroats+2);

	for(int newTIdx = m_numPores+2; newTIdx < static_cast< int >(m_rockLattice.size()); ++newTIdx)
		m_rockLattice[newTIdx]->updateLatticeIndex(newTIdx);
}

void Netsim::modifyNetwork(InputData& input)
{

	vector<Element*> pores(m_rockLattice.begin()+1, m_rockLattice.begin()+1+m_numPores);
	vector<Element*> throats(m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());

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
			vector<Element*> elems(m_rockLattice.begin()+1, m_rockLattice.begin()+1+m_numPores);
			elems.insert(elems.end(), m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());
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
			for(size_t elm = 0; elm < m_rockLattice.size(); ++elm)
			{
				double oldRad(m_rockLattice[elm]->model()->radius());
				m_rockLattice[elm]->ChModel()->setRadius(oldRad*scaleFactor);
			}
		}
	}


	double netPoroTrgt(-1.0), clayPoroTrgt(-1.0);


	///. MODIFY_RAD_DIST: modify radius and areas, skip changing volumes for now, they
	///. will be updated below, or overwritten based on keyword MODIFY_PORO
	{
		sort(throats.begin(), throats.end(), lenToRadRatioCmp());
		double oldmMedianLenToRadRatio(throats[throats.size()/2]->lenToRadRatio());
		
		string throatRadOptions, poreRadOptions;
		bool writeRDistToFile(false), maintainLtoR(false);
		int  throatRadModel(0), poreRadModel(0), numPtsRDist(0);
		input.getModifyRadDistOptions(throatRadModel, poreRadModel, throatRadOptions, poreRadOptions,
			maintainLtoR, writeRDistToFile, numPtsRDist);
		if(poreRadModel == -1 || throatRadModel == -1)
		{
			string fileName("RadDist.csv");
			vector<Element*> elems(m_rockLattice.begin()+1, m_rockLattice.begin()+1+m_numPores);
			elems.insert(elems.end(), m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());
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
		sort(throats.begin(), throats.end(), lenToRadRatioCmp());
		//double oldmMedianLenToRadRatio(m_medianLenToRadRatio);
		double MedianLenToRadRatio = throats[throats.size()/2]->lenToRadRatio();
		
		if(poreRadModel && maintainLtoR)
		//if(Rad_scaleFactor > 0.0)
		{
			double Rad_scaleFactor = oldmMedianLenToRadRatio/MedianLenToRadRatio;
			for(int i = 0; i < m_numPores + 2; ++i)
				m_rockLattice[i]->node()->rePosition(Rad_scaleFactor);

			for(size_t j = m_numPores + 2; j < m_rockLattice.size(); ++j)
				m_rockLattice[j]->modifyLength(Rad_scaleFactor);

			if(netPoroTrgt < 0) netPoroTrgt = m_totalFlowVolume/m_satBoxVolume;
			if(clayPoroTrgt < 0) clayPoroTrgt = m_totalClayVolume/m_satBoxVolume;
			m_xSize *= Rad_scaleFactor;
			m_ySize *= Rad_scaleFactor;
			m_zSize *= Rad_scaleFactor;
			m_satBoxVolume = m_xSize*m_ySize*m_zSize*(m_satBoxEnd-m_satBoxStart);
		}
	}



	///. modify volumes and porosity; based on MODIFY_PORO keyword, or based on MODIFY_RAD_DIST if MODIFY_PORO not set,
	{

		input.getModifyPoro(netPoroTrgt, clayPoroTrgt);///. will not change anything if keyword MODIFY_PORO not found

		double targetNetVol = netPoroTrgt >= 0.0 ? m_satBoxVolume*netPoroTrgt: m_totalFlowVolume;
		double targetClayVol = clayPoroTrgt >= 0.0 ? m_satBoxVolume*clayPoroTrgt: m_totalClayVolume;

		if(netPoroTrgt >= 0.0 || clayPoroTrgt >= 0.0)
		{
			double changeNetFrac((targetNetVol-m_totalFlowVolume)/m_totalFlowVolume);
			double changeClayFrac((targetClayVol-m_totalClayVolume)/m_totalClayVolume);
			double changeClayFromNet((targetClayVol-m_totalFlowVolume)/m_totalFlowVolume);
			
			double netVolSum(0.0), clayVolSum(0.0);
			for(size_t i = 0; i < m_rockLattice.size(); ++i)
			{
				double newNetVol = m_rockLattice[i]->flowVolume()*(1.0+changeNetFrac);
				
				double newClayVol(0.0);
				if(clayPoroTrgt >= 0.0 && m_totalClayVolume > 0.0)
					newClayVol = m_rockLattice[i]->clayVolume()*(1.0+changeClayFrac);
				else if(clayPoroTrgt >= 0.0 && m_totalClayVolume < 1.0e-12)
					newClayVol = m_rockLattice[i]->flowVolume()*(1.0+changeClayFromNet);

				m_rockLattice[i]->adjustVolume(newNetVol, newClayVol, netVolSum, clayVolSum);
			}
			
			m_totalFlowVolume = netVolSum;
			m_totalClayVolume = clayVolSum;
		}
	}

}


/**
 * Change the advancing contact angle depending on several criteria set out in
 * input file. Wettability alterations will only occur in elems that have been
 * drained
 */
void Netsim::applyFWettabilityChange(InputData& input)
{



	//set equil contact angles
	vector<Element*> pores2Set;pores2Set.reserve(m_numPores+2);//(m_rockLattice.begin(), m_rockLattice.begin()+2+m_numPores);
	vector<Element*> throats2Set;throats2Set.reserve(m_rockLattice.size()-m_numPores+2);//(m_rockLattice.begin()+2+m_numPores, m_rockLattice.end());
	
	for(int el = 0; el < 2+m_numPores; ++el)
		if(m_rockLattice[el]->iRockType() == 0)
		{
			pores2Set.push_back(m_rockLattice[el]);
			// TotalVol += m_rockLattice[el]->flowVolume();
		}

	for(size_t el = 2+m_numPores; el < m_rockLattice.size(); ++el)
		if(m_rockLattice[el]->iRockType() == 0)
		{
			throats2Set.push_back(m_rockLattice[el]);
			// TotalVol += m_rockLattice[el]->flowVolume();
		}
	
	
	
	
	string  angDistScheme;
	double 	modelTwoSepAng(0.0);
	int	   wettingClass = 1;

	//input.equilConAng(wettingClass, minEqConAng, maxEqConAng, wettDelta, wettEta, angDistScheme, modelTwoSepAng);
	const string keywordEqui("EQUIL_CON_ANG");
	const string& equiStr=m_input.keywordData(keywordEqui); ;
	if(!equiStr.empty())
	{
		double minEqConAng(0.0), maxEqConAng(0.0),  wettDelta, wettEta, modelTwoSepAng(0.0);		

		istringstream data;
		data.str(equiStr);
		//cout<< "Reading " << keywordEqui << endl;
		data >> wettingClass >> minEqConAng >> maxEqConAng >> wettDelta >> wettEta;
		if(!data) m_input.errorMsg(keywordEqui);

		data >> angDistScheme;				  // Maintain backward compatibility
		if(!data) angDistScheme = "rand";	   // by allowing both model and
											// separation angle not to be present
		data >> modelTwoSepAng;
		if(!data) modelTwoSepAng = 25.2;

		m_input.errorInDataCheck(data, keywordEqui);

		minEqConAng *= acos(-1.0) / 180.0;
		maxEqConAng *= acos(-1.0) / 180.0;
		modelTwoSepAng *= acos(-1.0) / 180.0;
		
		
		setContactAngles(pores2Set, throats2Set, minEqConAng, maxEqConAng, wettDelta, wettEta, wettingClass, modelTwoSepAng,  angDistScheme);
	}
	else
	{
		m_input.missingDataErr(keywordEqui+" not found" );
	}




	//input.fracWetting(applyFracWetting, fraction, volBased, minAng, maxAng, deltaExp, etaExp,
		//fracWetModel, clustDiam, oilClustInWat);
	const string frac_con_ang("FRAC_CON_ANG");
	const string& dataStr=m_input.keywordData(frac_con_ang);
	if(!dataStr.empty())
	{
		bool volBased(false), oilClustInWat(true);
		double fraction, minAng, maxAng, deltaExp, etaExp, oilInvadedVol(0.0), TotalVol(0.0), fracWetted(0.0);
		int totElem(0), clustDiam(0);
		string fracWetModel;

		{
			istringstream data;
			data.str(dataStr);

			char volModel;
			//cout<< "Reading " << frac_con_ang << endl;
			data >> fraction >> volModel  >> minAng >> maxAng >> deltaExp >> etaExp >> fracWetModel;
			minAng *= acos(-1.0) / 180.0;
			maxAng *= acos(-1.0) / 180.0;
			volBased = (volModel == 'T' || volModel == 't');
			if(!(volModel == 'T' || volModel == 't' || volModel == 'F' || volModel == 'f'))
				cout<<"	  Error Wrong Model, expected T/t or F/f"<<endl;
			if(!data) m_input.errorMsg(frac_con_ang);

			char oInW('T');
			data >> clustDiam >> oInW;
			oilClustInWat = (oInW == 'T' || oInW == 't');
			if(clustDiam<1 && (fracWetModel[1] == 'o' || fracWetModel[1] == 'O')) cout<<" please provide a value larger than 1 for cluster diameter in keyword "<< frac_con_ang <<endl<<dataStr<<endl;

			m_input.errorInDataCheck(data, frac_con_ang);
			//cout<<"FRAC_CON_ANG read"<<endl;
			//cout<< " " << fraction << " " << volModel  << " " << minAng << " " << maxAng << " " << deltaExp << " " << etaExp << " " << fracWetModel<<endl;
			//cout<< " " << clustDiam << " " << oilClustInWat<<endl;		
		}






		vector<Element*> pores2BAltered, throats2BAltered;



		for(size_t el = 0; el < m_rockLattice.size(); ++el)
		{
			if(m_rockLattice[el]->iRockType() == 0)
			{
				if(m_rockLattice[el]->model()->conductsAnyOil())
				{
					if(el < static_cast< size_t >(m_numPores+2)) ++totElem;						// Only count pores
					oilInvadedVol += m_rockLattice[el]->flowVolume();
				}
				TotalVol += m_rockLattice[el]->flowVolume();
			}
		}
		
		if(oilInvadedVol/TotalVol>1.0e-12)
		{
		  if(fracWetModel[1] == 'o' || fracWetModel[1] == 'O')	// Spatial correlation approach
		  {
			int numFracWetted(0), clusterIdx(1);
			bool allSystemsAreGo(true);
			double targetFrac = oilClustInWat ? fraction: 1.0-fraction;
			long long counter(0);
			while(allSystemsAreGo)
			{
				if(++counter>m_numPores*100) {m_out<<"Warning: could not apply the requested wetting fraction"<<endl; break;}
				double clusterVol(0.0);
				int clusterElem(0);
				set< Element * > frontier, oldFrontier;
				double randNum = double(rand()) / RAND_MAX;
				//double randNum = 0.5;   // delete me
				int randPoreIdx = randNum*m_numPores;
				randPoreIdx = max(1, randPoreIdx);
				if(m_rockLattice[randPoreIdx]->iRockType() == 0)
				{
				  frontier.insert(m_rockLattice[randPoreIdx]);
				  typedef set<Element*>::iterator ItrSet;
				  int frontExpansions(0);
				  while(!frontier.empty() && allSystemsAreGo && frontExpansions < clustDiam)  // Front expansions will go over thraots as well as pores, wheras
				  {																		   // cluster diam is wrt pores  => Don't divide by 2
					oldFrontier = frontier;
					frontier.clear();
					++frontExpansions;
					for(ItrSet itrElem = oldFrontier.begin(); itrElem != oldFrontier.end(); ++itrElem)
					{

						if( (*itrElem)->model()->clusterIndex() != clusterIdx && allSystemsAreGo)
						{
							if((*itrElem)->model()->clusterIndex() == 0 && (*itrElem)->model()->conductsAnyOil())
							{
								if(oilClustInWat && dynamic_cast< Pore* >(*itrElem) != 0)
								{
									ensure((*itrElem)->iAmAPore());
									pores2BAltered.push_back(*itrElem);
								}
								else if(oilClustInWat)
								{
									ensure(!(*itrElem)->iAmAPore());
									throats2BAltered.push_back(*itrElem);
								}
								(*itrElem)->ChModel()->setClusterIndex(clusterIdx);					// Set flag
								clusterVol += (*itrElem)->flowVolume();
								fracWetted += (*itrElem)->flowVolume();
								++clusterElem;
								if(dynamic_cast< Pore* >(*itrElem) != 0) ++numFracWetted;   // Only count pores

								//m_wettingFraction = volBased ? fracWetted/oilInvadedVol: numFracWetted/totElem;
								allSystemsAreGo = ((volBased && fracWetted/oilInvadedVol < targetFrac) ||
									(!volBased && static_cast<double>(numFracWetted)/totElem < targetFrac));
								if(!allSystemsAreGo)
									break;
							}

							for(int i = 0; i < (*itrElem)->connectionNum(); ++i)
							{
								if(!(*itrElem)->connection(i)->isEntryOrExitRes() &&
									(*itrElem)->connection(i)->model()->clusterIndex() != (clusterIdx) &&
									 ((*itrElem)->connection(i)->iRockType() == 0)
									)
								{
									frontier.insert((*itrElem)->connection(i));
								}
							}			
						}
					}
				  }
				  ++clusterIdx;
			   }
			}
			if(!oilClustInWat)
			{
				numFracWetted = 0;
				fracWetted = oilInvadedVol-fracWetted;
				for(size_t elm = 0; elm < m_rockLattice.size(); ++elm)
				{
					if(m_rockLattice[elm]->model()->conductsAnyOil() && m_rockLattice[elm]->model()->clusterIndex() == 0 &&
									 (m_rockLattice[elm]->iRockType() == 0) )
					{
						++numFracWetted;
						if(static_cast< int >(elm) < m_numPores+2)
						{
							ensure(m_rockLattice[elm]->iAmAPore());
							pores2BAltered.push_back(m_rockLattice[elm]);
						}
						else
						{
							ensure(!m_rockLattice[elm]->iAmAPore());
							throats2BAltered.push_back(m_rockLattice[elm]);
						}
					}
				}
			}

			m_out << "================================================================ "	<< endl
				<< "Fractional wetting was applied"										 << endl
				<< "Number of correlated regions: " << clusterIdx-1						 << endl
				<< "Altered volume (of oil invaded volume): " << fracWetted/oilInvadedVol	   << endl
				<< "Altered volume (of total net pore volume): " << fracWetted/TotalVol	<< endl
				<< "Number of pores altered: " << numFracWetted							 << endl
				<< "================================================================="	  << endl;
		 }
		  else
		  {
			oilClustInWat = true;   // Ensure this is set for non-spatial approach
			vector< pair<double, Element*> > toBeAltered;

			for(int i = 1; i <= m_numPores; ++i)
			{
				if(m_rockLattice[i]->model()->conductsAnyOil() && (m_rockLattice[i]->iRockType() == 0))
				{
					pair<double, Element*> rockEntry;
					rockEntry.second = m_rockLattice[i];
					if(fracWetModel[1] == 'a' || fracWetModel[1] == 'A')		///. rand // Random
					{
						rockEntry.first = double(rand()) / RAND_MAX;
						//rockEntry.first = 0.5;  // delete me
					}
					else if(fracWetModel[0] == 'r' || fracWetModel[0] == 'R')   ///. rMax/rMin // Radius
						rockEntry.first = m_rockLattice[i]->model()->radius();
					else
					{
						cerr << endl
							<< "=============================================================== "   << endl
							<< "Error: Did not recognize fractional wetting model: " << fracWetModel  << endl
							<< "=============================================================== "   << endl << endl;   exit(-1);
					}
					toBeAltered.push_back(rockEntry);
				}
			}
			if(fracWetModel[2] == 'a' || fracWetModel[2] == 'A') ///. rM'a'x
				sort(toBeAltered.begin(), toBeAltered.end(), FracWettInc());
			else
				sort(toBeAltered.begin(), toBeAltered.end(), FracWettDec()); ///. rMin/rand

			int numElem(0), targetNum(fraction*toBeAltered.size());
			double poreVol(0.0);

			while(!toBeAltered.empty() &&
				((volBased && poreVol/oilInvadedVol < fraction) || (!volBased && numElem < targetNum)))
			{
				++numElem;
				Element* elem = toBeAltered.back().second;
				poreVol += elem->flowVolume();
				pores2BAltered.push_back(elem);
				elem->ChModel()->setClusterIndex(1);
				toBeAltered.pop_back();
				for(int t = 0; t < elem->connectionNum(); ++t)
				{
					Element* throat = elem->connection(t);
					if(throat->model()->clusterIndex() == 0 && throat->model()->conductsAnyOil())
					{
						double randNum = double(rand()) / double(RAND_MAX);
						//double randNum = 0.5;   // delete me
						if(randNum < fraction)
						{
							throats2BAltered.push_back(throat);
							throat->ChModel()->setClusterIndex(1);
							poreVol += throat->flowVolume();
						}
					}
				}
			}
			m_out << "================================================================ "   << endl
				<< "Fractional wetting was applied"									 << endl
				<< "Altered volume (of oil invaded volume): " << poreVol/oilInvadedVol	  << endl
				<< "Altered volume (of total net pore volume): " << poreVol/TotalVol   << endl
				<< "Number of pores altered: " << numElem							   << endl
				<< "================================================================="  << endl;
	  	  }
		}
		else
			m_out << "================================================================ "	<< endl
			<< "Warning:  fractional wetting was NOT applied, "										 << endl
			<< "too low oil invaded volume/total volume: " << oilInvadedVol/TotalVol	<< endl
			<< "================================================================="	  << endl;
	
		setContactAngles(pores2BAltered, throats2BAltered, minAng, maxAng, deltaExp, etaExp, wettingClass, modelTwoSepAng,  angDistScheme, oilClustInWat);
	}

}

/**
 * Does the first flow through of the network
 */
void Netsim::SolveSinglePhaseFlow(InputData& input)
{

	m_out << "solving for single-phase flow";cout.flush();
	updateSatAndConductances(0.0); ///.  updateSatAndConductances

	double singlePhaseErr(0.0), cpuTmp(0.0), currentErr(0.0);

	if((m_writeWatVelocity || m_writeOilVelocity) && m_writeSlvMatrixAsMatlab) createMatlabLocationData();

	string matFile = m_matrixFileName + "_init";

	//if(useHypre)
	 m_solver = new hypreSolver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_comn.debugMode, matFile, m_writeSlvMatrixAsMatlab);
	//else
	 //m_solver = new amg_solver(m_rockLattice, m_krInletBoundary, m_krOutletBoundary, m_numPores+1, m_maxNonZeros, m_comn.debugMode, matFile, m_writeSlvMatrixAsMatlab);
	m_out << " ";cout.flush();

	m_singlePhaseWaterQ = m_solver->flowrate(m_inletSolverPrs, m_outletSolverPrs, m_water, singlePhaseErr,
		cpuTmp, 1.0, m_writeWatVelocity, m_writeWatMatrix);

	m_singlePhaseOilQ = m_singlePhaseWaterQ * (m_water.viscosity() / m_oil.viscosity());
	if(m_useAvrPrsAsBdr)
	{
		Netsim::prsOrVoltDrop(&m_water, 0, m_singlePhaseDprs);
		m_deltaPw = m_singlePhaseDprs;
	}
	//if(m_prtPressureProfile && m_wantRelPerm) recordPrsProfiles(m_water);

	m_watFlowRate = m_singlePhaseWaterQ;
	m_oilFlowRate = 0.0;



	m_out << " ";cout.flush();
	m_singlePhaseCurrent = m_solver->flowrate(m_inletSolverPrs, m_outletSolverPrs, m_water, currentErr, cpuTmp, 1.0,
		m_writeResVelocity, m_writeResMatrix, true);
	m_current = m_singlePhaseCurrent;
	if(m_useAvrPrsAsBdr) Netsim::prsOrVoltDrop(&m_water, 1, m_singlePhaseDvolt);
	m_out << " \n";cout.flush();

	if(fabs(singlePhaseErr) > 0.05 || fabs(currentErr) > 0.05)
	{
		m_out << endl
			<< "====================================================="  << endl
			<< "Warning. Large errors were detected when solving for"   << endl
			<< "single phase conditions."							   << endl
			<< "Flow Error: " << singlePhaseErr						 << endl
			<< "Resistivity Error: " << currentErr					  << endl
			<< "====================================================="  << endl;
	}





	delete m_solver;
	m_solver = NULL;
}


void Netsim::createMatlabLocationData() const
{
	ofstream outp("poreLocation.m");
	outp.flags(ios::scientific);
	outp.precision(3);
	outp << "% The physical dimensions (m) of the network is: " << endl;
	outp << "modSize(1) = " << m_xSize << ";" << endl;
	outp << "modSize(2) = " << m_ySize << ";" << endl;
	outp << "modSize(3) = " << m_zSize << ";" << endl;

	outp << endl << "% The location of individual pores (m). The array is 1-based" << endl;
	outp << "poreLoc = [";
	for(int i = 1; i <= m_numPores; ++i)
	{
		outp << m_rockLattice[i]->node()->xPos() << ", "
			<< m_rockLattice[i]->node()->yPos() << ", "
			<< m_rockLattice[i]->node()->zPos() << "; ..."
			<< endl;
	}
	outp << "];" << endl;
	outp.close();

	ofstream outt("throatConnection.m");
	outt << "% Pores to which the throats are connected. The indexing is" << endl
		<< "% 1-based, which is the same as that used by the Oren format" << endl;

	outt << "throatConn = [";
	for(size_t j = m_numPores+2; j < m_rockLattice.size(); ++j)
	{
		outt << m_rockLattice[j]->connection(0)->orenIndex() << ", ";
		if(m_rockLattice[j]->connectionNum() == 2)
			outt << m_rockLattice[j]->connection(1)->orenIndex() << "; " << endl;
		else
			outt << m_rockLattice[j]->connection(0)->orenIndex() << "; " << endl;
	}
	outt << "];" << endl;
	outt.close();
}


/**
 * Gets the pressure drop between two planes in the network, used for kr
 */
bool Netsim::prsOrVoltDrop(const Fluid* fluid, int resistSolve, double& prsDrop) const
{
	double prsOut(0.0), stdOut(0.0);
	int numOut(0);

	bool btOut = Netsim::avrPrsOrVolt(fluid, resistSolve, m_numPressurePlanes-1, prsOut, stdOut, numOut);

	double prsIn(0.0), stdIn(0.0);
	int numIn(0);

	bool btIn = Netsim::avrPrsOrVolt(fluid, resistSolve, 0, prsIn, stdIn, numIn);

	prsDrop = prsIn - prsOut;
	return btOut && btIn;
}

/**
 * Calculates the average pressure in throats crossing a plane. Assumes the pore
 * pressures to be point distributed. Will also calculate the number of pressure
 * points along with standard deviation
 */
bool Netsim::avrPrsOrVolt(const Fluid* fluid, int resistSolve, int prsPlane, double& res,
						  double& stdDev, int& numVal) const
{
	double loc(m_pressurePlanesLoc[prsPlane]*m_xSize), resSum(0.0), resSumSq(0.0), flowSum(0.0), flowTarget;
	vector< pair< double, double > > valArray;
	string phase("water");

	if(resistSolve)
	{
		flowTarget = m_current;
		phase = "current";
	}
	else if(dynamic_cast< const Oil* >(fluid) != NULL)
	{
		flowTarget = m_oilFlowRate;
		phase = "oil";
	}
	else
		flowTarget = m_watFlowRate;


	flowSum = 0.0;
	for(size_t i = 0; i < m_pressurePlanes[prsPlane].size(); ++i)
	{
		double val(0.0), flowRate(0.0);
		if(m_pressurePlanes[prsPlane][i]->prevSolvrRes(fluid, resistSolve, loc, val, flowRate))
		{
			pair< double, double > valPair(val, flowRate);
			valArray.push_back(valPair);
			flowSum += flowRate;
		}
	}

	double flowError(fabs((flowTarget-flowSum)/flowSum));

	if(resistSolve && flowTarget > 0.0)
		m_maxResIdxErr = max(m_maxResIdxErr, flowError);
	else if(fluid->isOil() && flowTarget > 0.0)
		m_maxOilFlowErr = max(m_maxOilFlowErr, flowError);
	else if(flowTarget > 0.0)
		m_maxWatFlowErr = max(m_maxWatFlowErr, flowError);

	if(flowError > MAX_FLOW_ERR && flowTarget > 0.0)
		cout<< "Large flow error (" << flowError*100.0 << "%) for " << phase << " at " << loc << endl;

	int numPts(static_cast< int >(valArray.size()));
	numVal = 0;
	for(int j = 0; j < numPts; ++j)
	{
		resSum += valArray[j].first;
		resSumSq += valArray[j].first*valArray[j].first;
		++numVal;
	}

	if(numVal > 0) res = resSum / numVal;
	if(numVal > 1) stdDev = sqrt((numVal*resSumSq - resSum*resSum)/(numVal*(numVal-1)));

	return numVal > 0;
}


/**
 * The data for the throats are read from the link files. Since the pores are not yet created their
 * indicies are stored in a temporary vector and the pointers will be initialized later. Since
 * in/outlet does not have separete entries in the data files, we keep track of connecting throats.
 * The strucure of the link files are as follows:
 *
 * *_link1.dat:
 * index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
 *
 * *_link2.dat:
 * index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
 */
int Netsim::readAndCreateThroats(InputData& input, vector< pair<int, int> >& newTrotNeibors,
								  vector<Element*>& throatsToInlet, vector<Element*>& throatsToOutlet,
								  const vector< pair< int, double> >& insidePoreHashs, vector< int >& insideThroatHashs, int newNumPores)
{
	cout<<"Reading throats"<<endl;

	vector<std::string> throatTypes(m_numThroats,"0");
	if(!input.keywordData("microporosity").empty())
	   loadRockTypeData(input.netFileBaseName()+ "_link3.dat",throatTypes);
	double clayPorosityToAdd = input.clayEdit();
	int numLengthErrors(0);
	int newNumThroats(0);
	for(int iTrot = 0; iTrot < m_numThroats; ++iTrot)
	{

		int pore1Idx, pore2Idx;
		double radius, shapeFactor, lenTot, lenPore1, lenPore2, lenThroat, volume, clayVolume;

		input.throatData(iTrot+1, pore1Idx, pore2Idx, volume, clayVolume, radius, shapeFactor,
			lenPore1, lenPore2, lenThroat, lenTot);


		///. allow initializing pores/throats with zero-length/volume
		radius = max(radius,1.0e-32); volume = max(volume,1.0e-96); shapeFactor = max(shapeFactor,1.0e-32);
		lenThroat = max(lenThroat,1.0e-32); lenPore1 = max(lenPore1,1.0e-32); lenPore2 = max(lenPore2,1.0e-32);
		lenTot = max(lenTot,3.0e-32);


		if(pore1Idx > 0)  newTrotNeibors[iTrot].first = insidePoreHashs[pore1Idx-1].first;
		else				newTrotNeibors[iTrot].first = pore1Idx;

		if(pore2Idx > 0)  newTrotNeibors[iTrot].second = insidePoreHashs[pore2Idx-1].first;
		else				newTrotNeibors[iTrot].second = pore2Idx;

		if(newTrotNeibors[iTrot].first > 0 || newTrotNeibors[iTrot].second > 0)
		{
			++newNumThroats;
			insideThroatHashs[iTrot] = newNumThroats;



			//double initConAng = weibull(minConAng, maxConAng, delta, eta);
			//double adjustingVol(input.clayEdit()*(volume+clayVolume));
			//adjustingVol = min(adjustingVol, volume);
			//adjustingVol = -min(-adjustingVol, clayVolume);

			//volume -= adjustingVol;
			//clayVolume += adjustingVol;
			clayVolume += clayPorosityToAdd*m_xSize*m_ySize*m_zSize/(m_numThroats+m_numPores);

			if(fabs(lenPore1+lenPore2+lenThroat-lenTot)/(lenTot+lenPore1+lenPore2+lenThroat) > 0.01)	++numLengthErrors;



			Element *throat = new Throat(m_comn, m_oil, m_water, radius, volume, clayVolume, shapeFactor,
							  lenThroat, lenPore1, lenPore2, newNumPores+1+newNumThroats,throatTypes[iTrot]);


			m_rockLattice.push_back(throat);

			if(newTrotNeibors[iTrot].first == DUMMY_IDX)
			{
				if(insidePoreHashs[pore1Idx-1].second < m_xSize/2.0)
					  newTrotNeibors[iTrot].first = -1;
				else  newTrotNeibors[iTrot].first = 0;
			}
			else if(newTrotNeibors[iTrot].second == DUMMY_IDX)
			{
				if(insidePoreHashs[pore2Idx-1].second < m_xSize/2.0)
					  newTrotNeibors[iTrot].second = -1;
				else  newTrotNeibors[iTrot].second = 0;
			}			
			if(newTrotNeibors[iTrot].first == -1 || newTrotNeibors[iTrot].second == -1)
				throatsToInlet.push_back(throat);
			else if(newTrotNeibors[iTrot].first == 0 || newTrotNeibors[iTrot].second == 0)
				throatsToOutlet.push_back(throat);
		}
	}

	if(numLengthErrors > 0)
	{


		m_out << endl
			<< "=================================================== "				<< endl
			<< "Warning: For " << numLengthErrors << " throats the lengths of the"  << endl
			<< "pore-throat-pore did not match the total length."				   << endl
			<< "This is generally only an artifact of the network"				  << endl
			<< "reconstruction process, and is not serious."						<< endl
			<< "=================================================== "				<< endl
			<< endl;


	}

	return newNumThroats;
}

/**
 * We're using different indicies for out and inlet. Paal-Eric uses -1 and 0 wheras we use 0 and (numPores+1), hence
 * we need to renumber these. The reason for this is that -1 is not a good index when storing the element
 * pointers in a vector.
 */
int Netsim::reIndex(int index) const
{
	if(index == -1)
		return 0;
	else if(index == 0)
		return m_numPores + 1;
	else
		return index;
}

int Netsim::setupInsidePoreHashings(InputData& input, vector< pair< int, double > >& insidePoreHashs) const
{
	int numPores(0);
	for(int index = 1; index <= m_numPores; ++index)
	{
		pair< int, double > entry(DUMMY_IDX, 0.0);
		input.poreLocation(index, entry.second);
		if(entry.second >= m_xSize*(1.0-m_keepFraction)/2.0 && entry.second <= m_xSize-m_xSize*(1.0-m_keepFraction)/2.0)
		{
			entry.first = ++numPores;
		}
		insidePoreHashs[index-1] = entry;
	}
	return numPores;
}


/**
 * The pore data is read from the node files. At this point the throats are already created and the pointers
 * can be set. The strucure of the node files are as follows:
 *
 * *_node1.dat:
 * index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
 *
 * *_node2.dat:
 * index, volume, radius, shape factor, clay volume
 */
void Netsim::readAndCreatePores(InputData& input, /*vector< pair<int, int> >& newTrotNeibors,const vector< pair< int, double> >& insidePoreHashs,*/
								 const vector< int >& insideThroatHashs, int newNumPores)
{
	cout<<"Reading pores\n";cout.flush();
	double  shaveOff(m_xSize*(1.0-m_keepFraction)/2.0);

	vector<std::string> poreTypes(m_numPores,"0");
	if(!input.keywordData("microporosity").empty())
	   loadRockTypeData(input.netFileBaseName()+ "_node3.dat",poreTypes);
	double clayPorosityToAdd = input.clayEdit();

	int newIndex(0);
	for(int iPore = 1; iPore <= m_numPores; ++iPore)
	{
		int connNumber;
		double xPos, yPos, zPos, volume, radius, shapeFactor, clayVolume;
		vector< int > connThroats, connPores;
		vector<Element*> connectingThroats;

		input.poreData(iPore, xPos, yPos, zPos, connNumber, connThroats, connPores, volume, clayVolume,
			radius, shapeFactor);

		if(xPos >= shaveOff && xPos <= m_xSize-shaveOff)
		{
			++newIndex;
			
			//double adjustingVol(input.clayEdit()*(volume+clayVolume));
			//adjustingVol = min(adjustingVol, volume);
			//adjustingVol = -min(-adjustingVol, clayVolume);
			//volume -= adjustingVol;
			//clayVolume += adjustingVol;
			clayVolume += clayPorosityToAdd*m_xSize*m_ySize*m_zSize/(m_numThroats+m_numPores);

			connectingThroats.resize(connNumber);
			for(int j = 0; j < connNumber; ++j)
			{
				//int hashedThroatIdx = insideThroatHashs[connThroats[j]-1];

				/////.  not used
				//int hashedPoreIdx = connPores[j];
				//if(hashedPoreIdx > 0 && insidePoreHashs[connPores[j]-1].first != DUMMY_IDX)
					//hashedPoreIdx = insidePoreHashs[connPores[j]-1].first;
				//else if(hashedPoreIdx > 0 && insidePoreHashs[connPores[j]-1].second < m_xSize/2.0)
					//hashedPoreIdx = -1;
				//else if(hashedPoreIdx > 0 && insidePoreHashs[connPores[j]-1].second > m_xSize/2.0)
					//hashedPoreIdx = 0;
				//ensure(hashedPoreIdx == newTrotNeibors[connThroats[j]-1].first
					//|| hashedPoreIdx == newTrotNeibors[connThroats[j]-1].second);


				//ensure(newIndex == newTrotNeibors[connThroats[j]-1].first || newIndex == newTrotNeibors[connThroats[j]-1].second);

				connectingThroats[j] = m_rockLattice[ newNumPores+1+insideThroatHashs[connThroats[j]-1] ];
			}

			double initSolvPrs = (m_outletSolverPrs + m_inletSolverPrs)/2.0;

			Node *currNode = new Node(newIndex, newNumPores, xPos-shaveOff, yPos, zPos, m_xSize*m_keepFraction); ///. shifted to left by half of curtailed fraction
			bool insideSlvrBox(currNode->isInsideBox(m_solverBoxStart, m_solverBoxEnd));
			bool insideSatBox(currNode->isInsideBox(m_satBoxStart, m_satBoxEnd));

			m_rockLattice[newIndex] = new Pore(m_comn, currNode, m_oil, m_water, radius, volume, clayVolume,
							 shapeFactor, insideSlvrBox, insideSatBox, initSolvPrs, connectingThroats,poreTypes[iPore-1]);
		}
	}
	cout<<" "<<endl;

}

//+++++++++ rubbish ++++++++++++++



/**
* Writes the optimized network to file
*/
void Netsim::writeNetworkToFile(const InputData& input) const
{

	bool writeNetToFile(false), writeInBinary(false);
	string fileNameBase;

	input.writeNetwork(writeNetToFile, writeInBinary, fileNameBase);

	if(!writeNetToFile) return;


	//if(!drainSinglets)
	//{
	cout<< "==================================================" << endl
		<< "CATION: The option not to drain dangling ends	   " << endl
		<< "DRAIN_SINGLETS not be used together with the option" << endl
		<< "to write the optimized network to file." << endl
		<< "==================================================" << endl;
	//exit(-1);
	//}

	if(writeInBinary)
	{
		string poreFileName(fileNameBase + "_node.bin"), throatFileName(fileNameBase + "_link.bin");
		ofstream pOut(poreFileName.c_str(), ios::binary), tOut(throatFileName.c_str(), ios::binary);
		if(!pOut || !tOut)
		{
			cout<< endl
				<< "=====================================================" << endl
				<< "Warning: Could not open " << poreFileName << endl
				<< "for writing. Network is not written to file." << endl
				<< "=====================================================" << endl
				<< endl;
		}

		pOut.write((char *)(&m_numPores), sizeof(int));
		pOut.write((char *)(&m_xSize), sizeof(double));
		pOut.write((char *)(&m_ySize), sizeof(double));
		pOut.write((char *)(&m_zSize), sizeof(double));
		tOut.write((char *)(&m_numThroats), sizeof(int));

		for (int i = 1; i <= m_numPores; ++i)
			m_rockLattice[i]->writeNetworkDataBinary(pOut);

		for (size_t j = m_numPores + 2; j < m_rockLattice.size(); ++j)
			m_rockLattice[j]->writeNetworkDataBinary(tOut);

		pOut.close();
		tOut.close();
	}
	else
	{
		string pOut1FileName(fileNameBase + "_node1.dat"), pOut2FileName(fileNameBase + "_node2.dat");

		ofstream pOut1, pOut2;
		pOut1.open(pOut1FileName.c_str());
		pOut2.open(pOut2FileName.c_str());
		pOut1.flags(ios::showpoint);
		pOut1.flags(ios::scientific);
		if(!pOut1 || !pOut2)
		{
			cout<< endl
				<< "=====================================================" << endl
				<< "Warning: Could not open " << pOut1FileName << endl
				<< "for writing. Network is not written to file." << endl
				<< "=====================================================" << endl
				<< endl;
		}

		pOut1 << m_numPores << "   " << m_xSize << "   " << m_ySize << "   " << m_zSize << '\n';

		for (int i = 1; i <= m_numPores; ++i)
			m_rockLattice[i]->writeNetworkData(pOut1, pOut2);

		pOut1.close();
		pOut2.close();

		string tOut1FileName(fileNameBase + "_link1.dat"), tOut2FileName(fileNameBase + "_link2.dat");
		ofstream tOut1, tOut2;
		tOut1.open(tOut1FileName.c_str());
		tOut2.open(tOut2FileName.c_str());

		tOut1 << m_numThroats << "   " << '\n';

		for (size_t j = m_numPores + 2; j < m_rockLattice.size(); ++j)
			m_rockLattice[j]->writeNetworkData(tOut1, tOut2);

		tOut1.close();
		tOut2.close();
	}
}


/**
 * In and outlet only need to know what throats are connected to it. Their indicies are 0 and (n_pores + 1)
 */
void Netsim::createInAndOutletPore(int index, double xPos, double yPos, double zPos, vector<Element*>& connThroats)
{
	if(xPos > 0.0 && xPos < m_xSize)
	{   cerr << "\n\nError: Entry and exit reservoirs cannot be within the network model area.\n" << endl;		exit(-1);
	}

	Node *currNode = new Node(index, m_numPores, xPos, yPos, zPos, m_xSize);
	Element *pore = new InOutBoundaryPore(m_comn, currNode, m_oil, m_water, connThroats);
	m_rockLattice[currNode->index()] = pore;
}

/*
//double Netsim::trappedWWoil() const
//{
	//double volOil(0.0);
	//for(size_t i = 0; i < m_rockLattice.size(); ++i)
	//{
		//Element* el = m_rockLattice[i];
		//if(el->isInsideSatBox() &&
			//el->model()->waterWet() &&
			//el->model()->bulkOil() &&
			//el->isTrappedOil())
		//{
			//volOil += el->oilVolume();
		//}
	//}
	//return volOil / (m_totalFlowVolume + m_totalClayVolume);
//}
*/


/**
 * Do a single calculation of relative permeability
 */
void Netsim::solve_forRelPermResIndex(bool wantRelPerm, bool wantResIdx)
{
	if(wantRelPerm)
	{
		Netsim::solveforRelPerm();
	}

	if(wantResIdx)
	{
		Netsim::solveforResIndex();
	}


}

void Netsim::solveforRelPerm()
{
	double oilErr(0.0), watErr(0.0);

	m_watFlowRate = m_solver->flowrate(m_inletSolverPrs, m_outletSolverPrs, m_water, watErr,
						   m_cpuTimeKrw, m_satWater, m_writeWatVelocity, m_writeWatMatrix);
	double relPermWater = m_watFlowRate / (m_singlePhaseWaterQ+1.0e-100);
	if(m_watFlowRate != 0.0)
	{
		if(m_useAvrPrsAsBdr)
		{
			bool chk = Netsim::prsOrVoltDrop(&m_water, 0, m_deltaPw);
			ensure(chk);
			relPermWater = (m_watFlowRate * m_singlePhaseDprs) / ((m_singlePhaseWaterQ+1.0e-100) * m_deltaPw);
		}
	}

	m_oilFlowRate = m_solver->flowrate(m_inletSolverPrs, m_outletSolverPrs, m_oil, oilErr, m_cpuTimeKro,
		m_satWater, m_writeOilVelocity, m_writeOilMatrix);

	double relPermOil = m_oilFlowRate / (m_singlePhaseOilQ+1.0e-100);
	if(m_oilFlowRate != 0.0)
	{
		if(m_useAvrPrsAsBdr)
		{
			bool chk = Netsim::prsOrVoltDrop(&m_oil, 0, m_deltaPo);
			ensure(chk);
			relPermOil = (m_oilFlowRate * m_singlePhaseDprs) / ((m_singlePhaseOilQ+1.0e-100) * m_deltaPo);
		}
	}

	m_out << " krw: " << setw(9) << std::left << relPermWater << " kro: " << setw(9) << std::left << relPermOil<< " ";

	m_maxOilFlowErr = max(m_maxOilFlowErr, oilErr);
	m_maxWatFlowErr = max(m_maxWatFlowErr, watErr);


}

void Netsim::solveforResIndex()
{
	double resErr(0.0);


	m_current = m_solver->flowrate(m_inletSolverPrs, m_outletSolverPrs, m_water, resErr, m_cpuTimeResIdx, m_satWater,
		m_writeResVelocity, m_writeResMatrix, true);
	double resIndex = m_singlePhaseCurrent / (abs(m_current)+1.0e-100);
	if(m_current > 0.0)
	{
		if(m_useAvrPrsAsBdr)
		{
			bool chk = Netsim::prsOrVoltDrop(&m_water, 1, m_deltaV);
			resIndex = m_singlePhaseCurrent / (abs(m_current)+1.0e-100) * m_deltaV/m_singlePhaseDvolt;

			ensure(chk);
		}
	}

	m_out << "RI: " << resIndex;;

	m_maxResIdxErr = max(m_maxResIdxErr, resErr);

}






/**
 * Iterates across all network elements (not in/outlets) and computes the water saturation
 */
void Netsim::updateSatAndConductances(double Pc)
{
	double volWater(0.0);
	for(size_t elem = 0; elem < m_rockLattice.size(); ++elem)
	{
		Element& poreI = *(m_rockLattice[elem]);//for the sake of doxygen only
		volWater += poreI.updateSat_calcR(Pc);
	}

	m_satWater = volWater / (m_totalFlowVolume + m_totalClayVolume);
}



void Netsim::writeResultData(bool relPermIncluded, bool resIdxIncluded)
{

	ostringstream resFileName;
	resFileName << m_baseFileName;
	
	bool matlabFormat, excelFormat, MCPMode(false);
	m_input.resFormat(matlabFormat, excelFormat, MCPMode);
	formatResults(matlabFormat, excelFormat, MCPMode);
	ostringstream out;
	string legend;
	int cycle;

	if(m_comn.injectant() == &m_oil)
	{
		legend = matlabFormat ? "draincycle_" : " drainage, cycle  ";
		cycle = m_comn.floodingCycle();
		resFileName << "_cycle"<<floodingCycle()<<"_drain";
	}
	else
	{
		legend = matlabFormat ? "imb_": " imbcycle ";
		cycle = m_comn.floodingCycle();
		resFileName << "_cycle"<<floodingCycle()<<"_imb";
	}

	if(!MCPMode)
	{
		char seperator = excelFormat ? ',': ' ';
		if(!matlabFormat)
		{
			out << m_input.netFileBaseName()<<"  "<<legend<<"\n";
			out << (m_results.size()) << seperator << "data points" << endl;
			out << "Sw,			Pc(Pa)			";

			if(m_apexPrsReported)
			{
				out << seperator << "Pc_b(Pa)	";
			}


			if(relPermIncluded)
			{
				out << "   ,Krw,		  	   ,Kro			   ";
			}

			if(resIdxIncluded)
			{
				out << "   ,RI	       	  CC ";
			}

			if(m_reportMaterialBal)
			{
				out << seperator << "Mass_w (kg) " << legend << cycle << ")"
					<< seperator << "Mass_o (kg) " << legend << cycle << ")";
			}
			out << endl;
		}
		else
		{
			out << "DataLegend = {'Sw " << legend << cycle << "' 'Pc (Pa) " << legend << cycle << "'";
			if(relPermIncluded) out << "'Krw " << legend << cycle << "' 'Kro " << legend << cycle << "'";
			if(resIdxIncluded)  out << "'I " << legend << cycle << "'";
			if(m_reportMaterialBal) out << "'Mass_w (kg) " << legend << cycle << "' 'Mass_o (kg) " << legend << cycle << "'";
			out << "};" << endl << "Res_" << legend << cycle << " = [";
		}
	}

	
	if(matlabFormat || excelFormat || !MCPMode)
	{
		string extension(".out");
		if(matlabFormat) extension = ".m";
		else if(excelFormat) extension = ".csv";

		ofstream of((resFileName.str()+extension).c_str());

		if(!of)
		{
			cerr << "\n\n *** Error: Could not open " << (resFileName.str()+extension) << " for writing ***\n\n"  << endl;		exit(-1);
		}

		of << out.str();
		for(size_t j = 0; j < m_results.size(); ++j)
			of << m_results[j] << endl;
		if(matlabFormat) of << "];" << endl;
	}





		if(MCPMode)
		{
			string mcpFilename=m_baseFileName+"_upscaled.dat";
			ofstream of;
			if(cycle==1)
			{
				of.open(mcpFilename.c_str());

				of<<endl<<"\nhomogeneous: \t"<<m_baseFileName<<";"<<endl<<endl;

				double absPermeability = (m_singlePhaseWaterQ * m_water.viscosity() * m_xSize * (m_solverBoxEnd-m_solverBoxStart))
					/ (m_ySize * m_zSize * (m_inletSolverPrs - m_outletSolverPrs));
				double formationFactor = ((m_ySize * m_zSize) * (m_inletSolverPrs - m_outletSolverPrs))
					/ (m_water.resistivity() * (m_singlePhaseCurrent+1.0e-200) * m_xSize * (m_solverBoxEnd-m_solverBoxStart));
					

				of<<endl<<m_baseFileName<<"_permeability: \t" <<absPermeability<<";"<<endl;

				of<<endl<<m_baseFileName<<"_porosity: \t"   <<m_totalFlowVolume / m_satBoxVolume<<";"<<endl;

				of<<endl<<m_baseFileName<<"_formationfactor: \t" <<formationFactor<<";"<<endl<<endl;
			}
			else
				of.open(mcpFilename.c_str(), std::ios_base::app);

			of<<"\n\n\n"<<m_baseFileName<<"_SwPcKrwKroRI_cycle"<<cycle<<endl;
			of << "%"<< legend  << endl;
			of << "%Sw	  \t   Pc(Pa)   \t   Krw	   \t   Kro	   \t   RI"<<endl;

			for(size_t j = 0; j < m_results.size(); ++j)
				of << m_results[j] << endl;
			of << "\n\n" << endl;
			of.close();

		}

	m_results.clear();
	m_resultWaterFlowRate.clear();
	m_resultOilFlowRate.clear();
	m_resultWaterSat.clear();
	m_resultCappPress.clear();
	//m_resultBoundingPc.clear();
	m_resultResistivityIdx.clear();
	m_resultWaterMass.clear();
	m_resultOilMass.clear();
	
}





string Netsim::calcUSBMindex() const
{
	ostringstream out;
	out.flags(ios::showpoint);
	out.flags(ios::fixed);
	out.precision(3);

	size_t numDrainPts(m_usbmDataDrainage.size()), numImbPts(m_usbmDataImbibition.size());
	double areaDrain(0.0), areaImb(0.0);

	if(numDrainPts > 1)
	{
		for(size_t i = 1; i < numDrainPts; ++i)
		{
			pair< double, double > ptOne(m_usbmDataDrainage[i]), ptTwo(m_usbmDataDrainage[i-1]); ///. Sw - Pc
			areaDrain += (ptTwo.second-ptOne.second) * ((ptTwo.first+ptOne.first)/2.0);
		}
	}

	if(numImbPts > 1)
	{
		for(size_t i = 1; i < numImbPts; ++i)
		{
			pair< double, double > ptOne(m_usbmDataImbibition[i-1]), ptTwo(m_usbmDataImbibition[i]);  ///. Sw - Pc
			areaImb += (ptTwo.second-ptOne.second)*(-(ptTwo.first+ptOne.first)/2.0);
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



