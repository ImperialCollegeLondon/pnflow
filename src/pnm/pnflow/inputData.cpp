#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <algorithm>
#include <map>
#include <set>

#include "inputData.h"


using namespace std;

/**
//
//  Available keywords are:
//
//X  TITLE               Title of run
//X  RAND_SEED           Seed for random number generator
//X  SAT_TARGET          Specifying saturation targets
//X  RELPERM_DEF         USe flowrate at single or residual saturation to calculate kr
//X  PRS_BDRS            What type of pressure boundaries should be used
//X  SAT_COMPRESS        Reduce saturation interval once kr drops below given threshold
//X  MAT_BAL             Calculates material balace of the fluids present
//X  APEX_PRS            Record tme minimum (drainage) pc observed in the last step
//X  POINT_SOURCE        Rather than inject across ine inlet face, inject from a single pore index
//X  SOLVER_TUNE         Tuning options for rel perm solver
//X  PORE_FILL_WGT       Weights for pore body filling mechnism
//X  PORE_FILL_ALG       Pore body filling mechnism algorithm
//X  MODIFY_RAD_DIST     Modify inscribed radii distribution
//X  MODIFY_G_DIST       Modify shape factor distribution
//X  MODIFY_PORO         Modify porosity
//X  MODIFY_CONN_NUM     Reduce connection number
//X  MODIFY_MOD_SIZE     Modify absolute model size
//X  SAT_COVERGENCE      Saturation covergence tuning
//X  TRAPPING            Trapping options
//X  FILLING_LIST        Create list for movie post processing
//X  OUTPUT              Create water sat map
//X  RES_FORMAT          Which output format to use
//X  GRAV_CONST          Vectorized gravity constant
//X  FLUID               Fluid properties
//X  CLAY_EDIT           Edit clay content
//X  A_CLOSE_SHAVE       Shave off boundaries to remove end-effects in net generation
//X  CALC_BOX            Where to set pressure boundaries for rel perm solver
//X  PRS_DIFF            Define the pressure differential across network for rel perm calculations
//X  NETWORK             The network to be used
//  NET_SERIES          Add networks in series
//  PERIODIC_PBC        Add periodic bc
//X  FRAC_CONT_ANG        Fractional wetting options
//X  INIT_CONT_ANG        Initial contact angle
//X  EQUIL_CON_ANG       Equilibration contact angle
//X  WRITE_NET           Write new network to file
//  SOLVER_DBG          Debugging options for solver
//
*/









const int       InputData::DUMMY_INDEX = -99;

InputData::InputData(const InputFile& inputFile)
:InputFile(inputFile,inputFile.name())
{
	averageThroatLength_ = 0.0;
	networkSeparation_ = 0.0;
	connectionsRemoved_ = 0;
	useAvrXOverThroatLen_ = false;
	addPeriodicBC_ = false;
	useAvrPbcThroatLen_ = false;
}





int InputData::randSeed() const
{	int seedNum;
	istringstream data;
	string keyword("RAND_SEED");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> seedNum;
		checkEndOfData(data,keyword);
	}
	else
		seedNum = (unsigned)time( NULL );
   
	return seedNum;
}


void InputData::prsBdrs(bool& usePrsBdr, bool& reportPrsBdr, int& numPlanes) const
{
	istringstream data;
	string keyword("PRS_BDRS");
	char usePrs, numPl;

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> usePrs >> numPl >> numPlanes;
		usePrsBdr = (usePrs == 'T' || usePrs == 't');
		reportPrsBdr = (numPl == 'T' || numPl == 't');

		if(!reportPrsBdr) numPlanes = 0;
		checkEndOfData(data, keyword);
	}
	else
	{
		usePrsBdr = false;
		reportPrsBdr = false;
		numPlanes = 0;
	}
}



void InputData::solverTune(double& eps, double& scaleFact, int& slvrOutput, bool& verbos, double& condCutOff) const
{
	istringstream data;
	char verb('F');
	string keyword("SOLVER_TUNE");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> eps >> scaleFact >> slvrOutput >> verb >> condCutOff;
		verbos = (verb == 'T' || verb == 't');
		checkEndOfData(data, keyword);
	}
	else
	{
		eps = 1.0E-15;
		scaleFact = 5;
		slvrOutput = 0;
		condCutOff = 0.0;
		verbos = false;
	}
}

/**
// Filling weights:
// Oren1/2 = 0.0, 0.5, 1.0, 5.0, 10.0, 50.0
// Blunt2 = 0.0, 15000, 15000, 15000, 15000, 15000
// Blunt1 = 0.0, 50E-6, 50E-6, 100E-6, 200E-6, 500E-6
*/
void InputData::poreFillWgt(vector< double >& weights) const
{
	istringstream data;
	string keyword("PORE_FILL_WGT");
	weights.resize(6);

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword ; cout.flush();
		data >> weights[0] >> weights[1] >> weights[2] >> weights[3] >> weights[4] >> weights[5];
		checkEndOfData(data, keyword);
	}
	else
	{
		weights[0] = 0.0;
		weights[1] = 15000.0;
		weights[2] = 15000.0;
		weights[3] = 15000.0;
		weights[4] = 15000.0;
		weights[5] = 15000.0;
		cout<< "Using default pore filling weights" ; cout.flush();

	}
   cout<< " " << endl;

}

void InputData::poreFillAlg(string& algorithm) const
{
	istringstream data;
	string keyword("PORE_FILL_ALG");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> algorithm;
		checkEndOfData(data, keyword);
	}
	else
	{
		algorithm = "blunt2";
	}
}


void InputData::getModifyRadDistOptions(int& throatModel, int& poreModel, string& throatOptions, string& poreOptions,
							  bool& maintainLtoR, bool& writeDistToFile, int& numPtsRDist)const
{
	istringstream data;
	char toFile, maintainAR;
	string keyword("MODIFY_RAD_DIST");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		getRadDist(data, throatModel, throatOptions);
		getRadDist(data, poreModel, poreOptions);
		data >> maintainAR >> toFile >> numPtsRDist;
		writeDistToFile = (toFile == 't' || toFile == 'T');
		maintainLtoR = (maintainAR == 't' || maintainAR == 'T');
		checkEndOfData(data, keyword);
	}
	else
	{
		throatModel = 0;
		poreModel = 0;
		writeDistToFile = false;
		maintainLtoR = false;
	}
}

void InputData::getModifyGDist(int& throatModel, int& poreModel, string& throatOptions, string& poreOptions,
							bool& writeDistToFile, int& numPts) const
{
	istringstream data;
	char toFile;
	string keyword("MODIFY_G_DIST");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		getRadDist(data, throatModel, throatOptions);
		getRadDist(data, poreModel, poreOptions);
		data >> toFile >> numPts;
		writeDistToFile = (toFile == 't' || toFile == 'T');
		checkEndOfData(data, keyword);
   }
	else
	{
		throatModel = 0;
		poreModel = 0;
		writeDistToFile = false;
	}
}

void InputData::getModifyPoro(double& netPoroTrgt, double& clayPoroTrgt) const
{
	istringstream data;
	string keyword("MODIFY_PORO");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> netPoroTrgt >> clayPoroTrgt;
		checkEndOfData(data, keyword);
	}
	//else
	//{
		//netPoroTrgt = -1.0;
		//clayPoroTrgt = -1.0;
	//}
}

void InputData::modifyConnNum(double& targetConnNum, string& model)const
{
	istringstream data;
	string keyword("MODIFY_CONN_NUM");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> targetConnNum >> model;
		checkEndOfData(data, keyword);
	}
	else
	{
		targetConnNum = -1.0;
	}
}

void InputData::getModifyModelSize(double& scaleFactor) const
{
	istringstream data;
	string keyword("MODIFY_MOD_SIZE");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> scaleFactor;
		checkEndOfData(data, keyword);
   }
	else
	{
		scaleFactor = -1.0;
	}
}


void InputData::getRadDist(istream& data, int& model, string& options) const
{
	ostringstream optionStr;
	data >> model;
	data.get (*(optionStr.rdbuf()));
	options = optionStr.str();
}


void InputData::fillingList(bool& isDrainage, bool& imbibition, bool& location) const
{
	istringstream data;
	string keyword("FILLING_LIST");
	char drain, imb, loc;

	if(getData(data, keyword))
	{
				cout<< "Reading " << keyword<<"\n ERROR disabled\n  ERROR \n  ERROR \n  ERROR \n  ERROR \n  ERROR \n " << endl;

		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> drain >> imb >> loc;
		isDrainage = (drain == 'T' || drain == 't');
		imbibition = (imb == 'T' || imb == 't');
		location = (loc == 'T' || loc == 't');
		checkEndOfData(data, keyword);
   }
	else
	{
	   isDrainage = false;
	   imbibition = false;
	   location = false;
	}
}



void InputData::resFormat(bool& matlabFormat, bool& excelFormat, bool& mcpFormat) const
{
	istringstream data;
	string keyword("RES_FORMAT");
	string resForm;

	if(getData(data, keyword))
	{
		//if (verbose)  cout<< "Reading " << keyword << endl;
		data >> resForm;
		matlabFormat = (resForm == "MATLAB" || resForm == "matlab" || resForm == "Matlab");
		excelFormat = (resForm == "EXCEL" || resForm == "excel" || resForm == "Excel");
		if (!excelFormat) excelFormat = (resForm == "excelAndMicroPorosity" || resForm == "EXCELANDMICROPOROSITY" || resForm == "ExcelAndMicroPorosity");
		if (excelFormat) 
			mcpFormat = (resForm == "excelAndMicroPorosity" || resForm == "EXCELANDMICROPOROSITY" || resForm == "ExcelAndMicroPorosity");
		else
			mcpFormat = (resForm == "upscaling" || resForm == "UPSCALING");

		checkEndOfData(data, keyword);
	}
}



void InputData::fluid(double& intfacTen, double& watVisc, double& oilVisc, double& watResist,
					  double& oilResist, double& watDens, double& oilDens) const
{
	istringstream data;
	string keyword("FLUID");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> intfacTen >> watVisc >> oilVisc >> watResist >> oilResist >> watDens >> oilDens;
		checkEndOfData(data, keyword);

		intfacTen *= 1.0E-3;                       // Get it into N/m
		watVisc *= 1.0E-3;                          // Get it into Pa.s
		oilVisc *= 1.0E-3;                          // Get it into Pa.s
	}
	//else
	//{
		//intfacTen = 30.0E-3;
		//watVisc = 1.0E-3;
		//oilVisc = 1.0E-3;
		//watResist = 1.0;
		//oilResist = 1000.0;
		//watDens = 1000.0;
		//oilDens = 1000.0;
	//}
}


void InputData::relPermCompression(bool& useComp, double& krThres,
		 double& deltaSw, bool& wettPhase, bool& nonWettPhase) const
{
	istringstream data;
	string keyword("SAT_COMPRESS");
	char wettP('F'), nonWettP('F');

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		useComp = true;
		data >> krThres >> deltaSw >> wettP >> nonWettP;
		wettPhase = (wettP == 'T' || wettP == 't');
		nonWettPhase = (nonWettP == 'T' || nonWettP == 't');
		checkEndOfData(data, keyword);
	}
	else
	{
		useComp = false;
		krThres = 0.0;
		deltaSw = 0.1;
	}
}


void InputData::calcBox(double& inletBdr, double& outletBdr) const
{
	istringstream data;
	string keyword("CALC_BOX");

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> inletBdr >> outletBdr;
		checkEndOfData(data, keyword);
   }
	else
	{
		inletBdr = 0.0;
		outletBdr = 1.0;
		cout<< "Using default  CALC_BOX: " << inletBdr << " " << outletBdr << endl;
	}
}

void InputData::prsDiff(double& inletPrs, double& outletPrs, bool& useGravInKr) const
{
	istringstream data;
	string keyword("PRS_DIFF");
	char grav;

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> inletPrs >> outletPrs >> grav;
		useGravInKr = (grav == 'T' || grav == 't');
		checkEndOfData(data, keyword);
	}
	else
	{
		inletPrs = 1.0;
		outletPrs = 0.0;
		useGravInKr = false;
	}
}


/**
// The nework data is supplied back as the pores and throats are created. The input
// streams are intially just opened and the headers read.
*/
void InputData::network(int& numPores, int& numThroats, double& xDim, double& yDim, double& zDim)
{
	istringstream data, dataNet, dataPbc;

	if(getData(dataPbc, "PERIODIC_BC"))
	{
		//if (verbose) cout<< "Reading " << "PERIODIC_BC" << endl;
		char usePbc, avrLen;
		dataPbc >> usePbc >> avrLen;
		useAvrPbcThroatLen_ = (avrLen == 'T' || avrLen == 't');
		addPeriodicBC_ = (usePbc == 'T' || usePbc == 't');
		checkEndOfData(dataPbc, "PERIODIC_BC");
	}
	else
	{
		addPeriodicBC_ = false;
		useAvrPbcThroatLen_ = false;
	}

	if(getData(dataNet, "NET_SERIES"))
	{
		//if (verbose) cout<< "Reading " << "NET_SERIES" << endl;
		char avrLen;
		dataNet >> numNetInSeries_ >> avrLen >> networkSeparation_;
		useAvrXOverThroatLen_ = (avrLen == 'T' || avrLen == 't');
		checkEndOfData(dataNet, "NET_SERIES");
	}
	else
	{
		numNetInSeries_ = 1;
		useAvrXOverThroatLen_ = false;
	}

	ststr     netNam;

	binaryFiles_ = false;
	if(!getVar(netNam, "networkFile"))
	 if(getData(data, "NETWORK",2))
	 {
		char binFile='F';
		std::string line;
		std::getline(data, line);
		std::istringstream iss(line);
		//if (verbose) cout<< "Reading " << "NETWORK" << endl;
		iss >> binFile >> netNam;
		binaryFiles_ = (binFile == 'T' || binFile == 't');
		checkEndOfData(data, "NETWORK");
	 }

	// rm "Net.xmf".   If we got here then new Net.xmf format is not supported
	if(netNam.size()>7 && netNam.substr(netNam.size()-7)=="Net.xmf") netNam=netNam.substr(0,netNam.size()-7);

	if(binaryFiles_)
	{
		string porePropFile(netNam + "_node.bin");
		poreProp_.open(porePropFile.c_str(), ios::binary);

		string throatPropFile(netNam + "_link.bin");
		throatProp_.open(throatPropFile.c_str(), ios::binary);

		poreProp_.read((char *)(&numPores), sizeof(int));
		poreProp_.read((char *)(&xDim), sizeof(double));
		poreProp_.read((char *)(&yDim), sizeof(double));
		poreProp_.read((char *)(&zDim), sizeof(double));
		throatProp_.read((char *)(&numThroats), sizeof(int));
	}
	else
	{
		string poreConnFile(netNam + "_node1.dat");          // Open file containing pore connection data
		poreConn_.open(poreConnFile.c_str());

		string porePropFile(netNam + "_node2.dat");          // Open file containing pore geometry data
		poreProp_.open(porePropFile.c_str());

		string throatConnFile(netNam + "_link1.dat");        // Open file containing throat connection data
		throatConn_.open(throatConnFile.c_str());

		string throatPropFile(netNam + "_link2.dat");        // Open file containing throat geometry data
		throatProp_.open(throatPropFile.c_str());


		poreConn_ >> numPores >> xDim >> yDim >> zDim;
		throatConn_ >> numThroats;
	}

	if (!poreProp_ || !throatProp_ || (!binaryFiles_ && (!poreConn_ || !throatConn_)))
	{
		std::cout<< "======================================== " << std::endl
			<<"Error: Unable to open network data files " << netNam << std::endl
			<<" failed on !"<< bool(poreProp_) <<"  || !"<< bool(throatProp_)  <<" || (!"<< bool(binaryFiles_) <<" && (!"<<bool(poreConn_)<<" || !"<<bool(throatConn_)<<"))" << std::endl
			<<"            "<< !poreProp_ <<"  ||  "<< !throatProp_  <<" || "<< (!binaryFiles_ && (!poreConn_ || !throatConn_)) << std::endl
			<< "======================================== " << endl;
		exit( -1 );
	}
	outD<<"network:" << netNam <<(binaryFiles_?" binary ":" ascii ")<<" Nt:"<<numThroats<<" Np:"<<numPores<< std::endl;
	origNumPores_ = numPores;
	origNumThroats_ = numThroats;
	origXDim_ = xDim;
	origYDim_ = yDim;
	origZDim_ = zDim;



	ensure(numNetInSeries_ >0, "Error: There need to be at least one net in series", -1 );


	loadPoreData();
	loadThroatData();

	numPores *= numNetInSeries_;
	numThroats = origNumThroats_*numNetInSeries_ - connectionsRemoved_*(numNetInSeries_-1);
	xDim = xDim*static_cast< double >(numNetInSeries_) + static_cast< double >(numNetInSeries_-1)*networkSeparation_;
}

/**
// As the pores are created the data is read from file and supplied back together with contact angles that are
// stored in vectors.
//
// The format of pore network files are:
// *_node1.dat:
// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
//
// *_node2.dat:
// index, volume, radius, shape factor, clay volume
*/
void InputData::loadPoreData()
{
	poreData_.resize(origNumPores_);
	for(int i = 0; i < origNumPores_; ++i)
	{
		if ((!binaryFiles_ && !poreConn_) || !poreProp_)
		{
			cerr << "=========================" << endl
				 << "Error while reading pores." << endl
				 << i << endl
				 << (!binaryFiles_ && !poreConn_) << " || "<< !poreProp_ << endl
				 << "==========================" << endl;
			exit( -1 );
		}

		PoreStruct *poreProp = new PoreStruct;
		int *connThroats, *connPores;

		if(binaryFiles_)
		{
			poreProp_.read((char *)(poreProp), sizeof(*poreProp));

			connThroats = new int[poreProp->connNum];
			connPores = new int[poreProp->connNum];
			poreProp_.read((char *)(connPores), poreProp->connNum*sizeof(int));
			poreProp_.read((char *)(connThroats), poreProp->connNum*sizeof(int));
		}
		else
		{
			int idx;
			bool isAtInletRes, isAtOutletRes;

			poreProp_ >> poreProp->index
				>> poreProp->volume
				>> poreProp->radius
				>> poreProp->shapeFact
				>> poreProp->clayVol;

			poreConn_ >> idx
				>> poreProp->x
				>> poreProp->y
				>> poreProp->z
				>> poreProp->connNum;

			ensure(idx == poreProp->index);

			connThroats = new int[poreProp->connNum];
			connPores = new int[poreProp->connNum];

			for(int k = 0; k < poreProp->connNum; ++k)
				poreConn_ >> connPores[k];

			poreConn_ >> isAtInletRes >> isAtOutletRes;

			for(int j = 0; j < poreProp->connNum; ++j)
				poreConn_ >> connThroats[j];
		}
		ThreeSome< PoreStruct*, int*, int* > elem(poreProp, connPores, connThroats);
		poreData_[i] = elem;

		for(int conn = 0; conn < poreProp->connNum; ++conn)
		{
			if(connPores[conn] == -1)
				inletPores_.push_back(elem);
			else if(connPores[conn] == 0)
				outletPores_.push_back(elem);
		}
	}

	if(addPeriodicBC_) findBoundaryPores();
	poreProp_.close();
	poreConn_.close();
}

void InputData::findBoundaryPores()
{
	int nDir = static_cast< int>(pow(origNumPores_, 1.0/3.0))+1;
	double xStep(origXDim_/nDir), yStep(origYDim_/nDir), zStep(origZDim_/nDir);
	for(int i = 0; i < nDir; ++i)
	{
		double xPos = xStep/2.0 + i*xStep;
		ensure(xPos < origXDim_);
		for(int j = 0; j < nDir; ++j)
		{
			double yPos = yStep/2.0 + j*yStep;
			double zPos = zStep/2.0 + j*zStep;
			ensure(xPos < origXDim_ && yPos < origYDim_ && zPos < origZDim_);
			vector< ThreeSome< double, double, double > > positonData;
			ThreeSome< double, double, double > xyMinus(xPos, yPos, 0.0);
			positonData.push_back(xyMinus);
			ThreeSome< double, double, double > xyPluss(xPos, yPos, origZDim_);
			positonData.push_back(xyPluss);
			ThreeSome< double, double, double > xzMinus(xPos, 0.0, zPos);
			positonData.push_back(xzMinus);
			ThreeSome< double, double, double > xzPluss(xPos, origYDim_, zPos);
			positonData.push_back(xzPluss);

			vector< int > poreIndecies;
			vector< double > p2BdrLength;
			findClosestPoreForPBC(positonData, poreIndecies, p2BdrLength);
			pair< int, int > xyConn(poreIndecies[0], poreIndecies[1]);
			pair< int, int > xzConn(poreIndecies[2], poreIndecies[3]);

			if(xyConn.first != xyConn.second && xyPbcConn_.count(xyConn) == 0)  // Detect 2D networks and duplicates
			{
				ThreeSome<int, int, double> pbcConn(xyConn.first, xyConn.second, p2BdrLength[0]+p2BdrLength[1]);
				pbcData_.push_back(pbcConn);
			}

			if(xzConn.first != xzConn.second && xzPbcConn_.count(xzConn) == 0)
			{
				ThreeSome<int, int, double> pbcConn(xzConn.first, xzConn.second, p2BdrLength[2]+p2BdrLength[3]);
				pbcData_.push_back(pbcConn);
			}

			xyPbcConn_.insert(xyConn);
			xzPbcConn_.insert(xzConn);
		}
	}
}

void InputData::findClosestPoreForPBC(const vector< ThreeSome< double, double, double > >& positonData,
									  vector< int >& poreIndecies, vector< double >& p2BdrLength) const
{
	poreIndecies.resize(positonData.size());
	p2BdrLength.resize(positonData.size(), 1.0E21);
	for(size_t i = 0; i < poreData_.size(); ++i)
	{
		for(size_t j = 0; j < positonData.size(); ++j)
		{
			double len = sqrt(pow(positonData[j].first()-poreData_[i].first()->x, 2.0)+
				pow(positonData[j].second()-poreData_[i].first()->y, 2.0)+
				pow(positonData[j].third()-poreData_[i].first()->z, 2.0));
			if(len < p2BdrLength[j])
			{
				p2BdrLength[j] = len;
				poreIndecies[j] = static_cast< int >(i)+1;
			}
		}
	}
}

void InputData::getOutletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const
{
	ThreeSome<int,int,double> entry = (*itr).second;
	thatIdx = entry.second() + origNumPores_*(numNetsInFront+1);
	bool throatToOutlet = outletThroat(entry.first());
	bool throatToInlet = inletThroat(entry.first());
	int hashedIdx(entry.first());
	if(numNetsInFront > 0 && throatToOutlet)
	{
		hashedIdx = reverseThroatHash_[hashedIdx-1] +
			numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
	}
	else if(throatToInlet)
	{
		ensure(reverseThroatHash_[hashedIdx-1] != DUMMY_INDEX);
		hashedIdx = reverseThroatHash_[hashedIdx-1] +
			(numNetsInFront+1)*origNumThroats_ - numNetsInFront*connectionsRemoved_;
	}
	throatIdx = hashedIdx;
}

void InputData::getInletData(MapItr itr, int numNetsInFront, int& throatIdx, int& thatIdx) const
{
	ThreeSome<int,int,double> entry = (*itr).second;
	thatIdx = entry.second() + origNumPores_*(numNetsInFront-1);
	bool throatToOutlet = outletThroat(entry.first());
	bool throatToInlet = inletThroat(entry.first());
	int hashedIdx(entry.first());
	if(numNetsInFront > 1 && throatToOutlet)
	{
		hashedIdx = reverseThroatHash_[hashedIdx-1] +
			(numNetsInFront-1)*origNumThroats_ - (numNetsInFront-2)*connectionsRemoved_;
	}
	else if(throatToInlet)
	{
		hashedIdx = reverseThroatHash_[hashedIdx-1] +
			numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
	}
	throatIdx = hashedIdx;
}

void InputData::poreLocation(int idx, double& xPos) const
{
	int stdIdx(idx), numNetsInFront(0);
	if(idx > origNumPores_)
	{
		cout<<"Error(?)";
		numNetsInFront = (idx-1)/origNumPores_;
		stdIdx = idx - numNetsInFront*origNumPores_;
	}

	PoreStruct *poreProp = poreData_[stdIdx-1].first();
	xPos = poreProp->x + numNetsInFront*(origXDim_+networkSeparation_);
}

/**
// There is a lot of memory assosiated with storing all network data.
// Best to clean up after ourself before proceeding
*/
void InputData::clearNetworkData()
{
	for(size_t i = 0; i < poreData_.size(); ++i)
	{
		delete[] poreData_[i].second();
		delete[] poreData_[i].third();
		delete poreData_[i].first();
	}

	for(size_t j = 0; j < throatData_.size(); ++j)
	{
		delete throatData_[j];
	}

	poreData_.clear();
	throatData_.clear();
	inletPores_.clear();
	outletPores_.clear();
	outletConnections_.clear();
	inletConnections_.clear();
	throatHash_.clear();
	reverseThroatHash_.clear();
}

/**
// As the throats are created the data is read from file and supplied back together with contact angles
// that are stored in vectors.
//
// The format of throat network files are:
// *_link1.dat:
// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
//
// *_link2.dat:
// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
*/
void InputData::loadThroatData()
{
	throatData_.resize(origNumThroats_);
	double lenSumPore(0.0), lenSumThroat(0.0);
	for(int i = 0; i < origNumThroats_; ++i)
	{
		if ((!binaryFiles_ && !throatConn_) || !throatProp_)
		{
			cerr << "=========================== " << endl
				<< "Error while reading throats." << endl
				<< "throat index: "<<i << endl
				<< bool(throatConn_)<<bool(throatProp_) << endl
				<< "============================ " << endl;
			exit( -1 );
		}

		ThroatStruct *throatProp = new ThroatStruct;

		if(binaryFiles_)
		{
			throatProp_.read((char *)(throatProp), sizeof(*throatProp));
		}
		else
		{
			int idx, tmp;

			throatConn_ >> throatProp->index
				>> throatProp->poreOne
				>> throatProp->poreTwo
				>> throatProp->radius
				>> throatProp->shapeFact
				>> throatProp->lenTot;

			throatProp_ >> idx
				>> tmp
				>> tmp
				>> throatProp->lenPoreOne
				>> throatProp->lenPoreTwo
				>> throatProp->lenThroat
				>> throatProp->volume
				>> throatProp->clayVol;

			ensure(idx == throatProp->index);
		}
		throatData_[i] = throatProp;
		lenSumThroat += throatProp->lenThroat;
		lenSumPore += (throatProp->lenPoreOne+throatProp->lenPoreTwo)/2.0;
	}
	throatConn_.close();
	throatProp_.close();
	averageThroatLength_ = lenSumThroat/origNumThroats_;
	averagePoreHalfLength_ = lenSumPore/origNumThroats_;

	if(addPeriodicBC_)
	{
		int index = origNumThroats_+1;
		for(size_t conn = 0; conn < pbcData_.size(); ++conn)
		{
			ThroatStruct *pbcThroat = new ThroatStruct;
			double randNum = double(rand()) / double(RAND_MAX);
			//double randNum = 0.5;   // delete me
			int randThroatIdx = (randNum*origNumThroats_);
			if(randThroatIdx >= origNumThroats_) randThroatIdx = origNumThroats_-1;
			pbcThroat->clayVol = 0.0;
			pbcThroat->index = index++;
			pbcThroat->poreOne = pbcData_[conn].first();
			pbcThroat->poreTwo = pbcData_[conn].second();
			pbcThroat->radius = throatData_[randThroatIdx]->radius;
			pbcThroat->shapeFact = throatData_[randThroatIdx]->shapeFact;
			pbcThroat->volume = 0.0;
			if(useAvrPbcThroatLen_ || pbcData_[conn].third() == 0.0)
			{
				pbcThroat->lenPoreOne = averagePoreHalfLength_;
				pbcThroat->lenPoreTwo = averagePoreHalfLength_;
				pbcThroat->lenThroat = averageThroatLength_;
				pbcThroat->lenTot = averageThroatLength_+2.0*averagePoreHalfLength_;
			}
			else
			{
				pbcThroat->lenTot = pbcData_[conn].third();
				pbcThroat->lenThroat = pbcData_[conn].third() *
					(averageThroatLength_/(averageThroatLength_+2.0*averagePoreHalfLength_));
				pbcThroat->lenPoreOne = (pbcThroat->lenTot-pbcThroat->lenThroat)/2.0;
				pbcThroat->lenPoreTwo = (pbcThroat->lenTot-pbcThroat->lenThroat)/2.0;
			}
			throatData_.push_back(pbcThroat);
			appendPoreData(pbcThroat->poreOne, pbcThroat->index, pbcThroat->poreTwo);
			appendPoreData(pbcThroat->poreTwo, pbcThroat->index, pbcThroat->poreOne);
		}
		origNumThroats_ += pbcData_.size();
	}

	set<int> addedInletThroats;

	if(numNetInSeries_ > 1)
	{
		for(int outT = 0; outT < origNumThroats_; ++outT)
		{
			if(throatData_[outT]->poreOne == 0 || throatData_[outT]->poreTwo == 0)
			{
				int outFacePore = throatData_[outT]->poreOne == 0 ? throatData_[outT]->poreTwo: throatData_[outT]->poreOne;
				double xPos = poreData_[outFacePore-1].first()->x - origXDim_ - networkSeparation_;
				double yPos = poreData_[outFacePore-1].first()->y;
				double zPos = poreData_[outFacePore-1].first()->z;
				double p2pLength(0.0);
				int inFacePore = findClosestPore(inletPores_, xPos, yPos, zPos, p2pLength);
				ThreeSome<int, int, double> outEntry(outT+1, inFacePore, p2pLength);
				ThreeSome<int, int, double> inEntry(outT+1, outFacePore, p2pLength);
				inletConnections_.insert(MultiConnValType(inFacePore, inEntry));
				outletConnections_.insert(MultiConnValType(outFacePore, outEntry));
			}
		}

		for(int inT = 0; inT < origNumThroats_; ++inT)
		{
			if(throatData_[inT]->poreOne == -1 || throatData_[inT]->poreTwo == -1)
			{
				int inFacePore = throatData_[inT]->poreOne == -1 ? throatData_[inT]->poreTwo: throatData_[inT]->poreOne;
				if(inletConnections_.count(inFacePore) == 0)
				{
					addedInletThroats.insert(inT+1);
					double xPos = poreData_[inFacePore-1].first()->x + origXDim_ + networkSeparation_;
					double yPos = poreData_[inFacePore-1].first()->y;
					double zPos = poreData_[inFacePore-1].first()->z;
					double p2pLength(0.0);
					int outFacePore = findClosestPore(outletPores_, xPos, yPos, zPos, p2pLength);
					ThreeSome<int, int, double> outEntry(inT+1, inFacePore, p2pLength);
					ThreeSome<int, int, double> inEntry(inT+1, outFacePore, p2pLength);
					inletConnections_.insert(MultiConnValType(inFacePore, inEntry));
					outletConnections_.insert(MultiConnValType(outFacePore, outEntry));
				}
				else
					++connectionsRemoved_;
			}
		}
	}

	int runningIndex(-1);
	throatHash_.resize(origNumThroats_-connectionsRemoved_);      // Since inlet throats do not exist in serial
	for(int j = 1; j <=  origNumThroats_; ++j)                       // nets, we end up having to take into account
	{                                                               // an offset for the indecies
		if((throatData_[j-1]->poreOne != -1 && throatData_[j-1]->poreTwo != -1) || addedInletThroats.count(j) > 0)
			throatHash_[++runningIndex] = j;
	}

	runningIndex = 1;
	reverseThroatHash_.resize(origNumThroats_);
	for(int k = 0; k < origNumThroats_; ++k)
	{
		int hashedIdx(DUMMY_INDEX);
		if((throatData_[k]->poreOne != -1 && throatData_[k]->poreTwo != -1) || addedInletThroats.count(k+1) > 0)
			hashedIdx = runningIndex++;

		reverseThroatHash_[k] = hashedIdx;
	}
}

int InputData::findClosestPore(const vector< ThreeSome< PoreStruct*, int*, int* > > pores, double xPos,
							   double yPos, double zPos, double& totalLen) const
{
	totalLen = 1.0E21;
	int index(DUMMY_INDEX);
	for(size_t i = 0; i < pores.size(); ++i)
	{
		double len = sqrt(pow(xPos-pores[i].first()->x, 2.0)+
			pow(yPos-pores[i].first()->y, 2.0)+
			pow(zPos-pores[i].first()->z, 2.0));
		if(len < totalLen)
		{
			totalLen = len;
			index = pores[i].first()->index;
		}
	}
	ensure(index != DUMMY_INDEX);
	return index;
}



void InputData::poreData(int idx, double& xPos, double& yPos, double& zPos, int& connNum, vector< int >& connThroats,
						 vector< int >& connPores, double& vol, double& volCl, double& rad, double& shapeFact)
{
	int stdIdx(idx), numNetsInFront(0);
	if(idx > origNumPores_)
	{
		numNetsInFront = (idx-1)/origNumPores_;
		stdIdx = idx - numNetsInFront*origNumPores_;
	}
	const PoreStruct *poreProp = poreData_[stdIdx-1].first();
	xPos = poreProp->x + numNetsInFront*(origXDim_+networkSeparation_);
	yPos = poreProp->y;
	zPos = poreProp->z;
	connNum = poreProp->connNum;
	vol = poreProp->volume;
	volCl = poreProp->clayVol;
	rad = poreProp->radius;
	shapeFact = poreProp->shapeFact;
	connThroats.resize(connNum);
	connPores.resize(connNum);
	bool outletPore(false), inletPore(false);
	for(int i = 0; i < connNum; ++i)
	{
		connPores[i] = poreData_[stdIdx-1].second()[i];
		connThroats[i] = poreData_[stdIdx-1].third()[i];

		if(connPores[i] != 0 && connPores[i] != -1 && numNetsInFront > 0)
		{
			connPores[i] += numNetsInFront*origNumPores_;
			int throatStdIdx = reverseThroatHash_[connThroats[i]-1];
			connThroats[i] = throatStdIdx + numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
		}
		else if(connPores[i] == 0)
			outletPore = true;
		else if(connPores[i] == -1)
			inletPore = true;
	}
	//vector< int > connPores; vector< int > connThroats; int connNum;
	if(numNetInSeries_ > 1 && (outletPore || inletPore))
	{
	//connNum = poreProp->connNum;
		pair<MapItr,MapItr> mapConns;
	//connThroats.resize(connNum);
		if(outletPore)
		{
	//connPores.resize(connNum);
			mapConns = outletConnections_.equal_range(stdIdx);
		}
		else
		{
	//bool outletPore(false), inletPore(false);
			mapConns = inletConnections_.equal_range(stdIdx);
		}
		MapItr itr = mapConns.first;
	//for(int i = 0; i < connNum; ++i)
		for(int j = 0; j < connNum; ++j)
		{
			if(numNetsInFront < numNetInSeries_-1 && connPores[j] == 0)
			{
				ensure(itr != mapConns.second);
				getOutletData(itr, numNetsInFront, connThroats[j], connPores[j]);
	//{
		//connPores[i] = poreData_[stdIdx-1].second()[i];
				++itr;
			}
			else if(numNetsInFront > 0 && connPores[j] == -1)
			{
				ensure(itr != mapConns.second);
				getInletData(itr, numNetsInFront, connThroats[j], connPores[j]);
		//connThroats[i] = poreData_[stdIdx-1].third()[i];
				++itr;
			}
			else if(numNetsInFront > 0 && connPores[j] == 0)
			{
				int throatStdIdx = reverseThroatHash_[connThroats[j]-1];
				connThroats[j] = throatStdIdx + numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
			}
		}

		while(itr != mapConns.second)
		{
			int poreIdx(DUMMY_INDEX), throatIdx(DUMMY_INDEX);
		//if(connPores[i] != 0 && connPores[i] != -1 && numNetsInFront > 0)
			if(numNetsInFront < numNetInSeries_-1 && outletPore)
				getOutletData(itr, numNetsInFront, throatIdx, poreIdx);
			else if(numNetsInFront > 0 && inletPore)
				getInletData(itr, numNetsInFront, throatIdx, poreIdx);
		//{
			if(poreIdx != DUMMY_INDEX)
			{
				connPores.push_back(poreIdx);
				connThroats.push_back(throatIdx);
				++connNum;
			}
			++itr;
			//connPores[i] += numNetsInFront*origNumPores_;
			//int throatStdIdx = reverseThroatHash_[connThroats[i]-1];
			//connThroats[i] = throatStdIdx + numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
		//}
		//else if(connPores[i] == 0)    outletPore = true;
		//else if(connPores[i] == -1)   inletPore = true;
	//}

	//if(numNetInSeries_ > 1 && (outletPore || inletPore))
	//{
		//pair<MapItr,MapItr> mapConns;

		//if(outletPore)
		//{
			//mapConns = outletConnections_.equal_range(stdIdx);
		//}
		//else
		//{
			//mapConns = inletConnections_.equal_range(stdIdx);
		//}
		//MapItr itr = mapConns.first;

		//for(int j = 0; j < connNum; ++j)
		//{	if(numNetsInFront < numNetInSeries_-1 && connPores[j] == 0)
			//{
				//ensure(itr != mapConns.second);
				//getOutletData(itr, numNetsInFront, connThroats[j], connPores[j]);

				//++itr;
			//}
			//else if(numNetsInFront > 0 && connPores[j] == -1)
			//{
				//ensure(itr != mapConns.second);
				//getInletData(itr, numNetsInFront, connThroats[j], connPores[j]);
				//++itr;
			//}
			//else if(numNetsInFront > 0 && connPores[j] == 0)
			//{
				//int throatStdIdx = reverseThroatHash_[connThroats[j]-1];
				//connThroats[j] = throatStdIdx + numNetsInFront*origNumThroats_ - (numNetsInFront-1)*connectionsRemoved_;
			//}
		//}

		//while(itr != mapConns.second)
		//{
			//int poreIdx(DUMMY_INDEX), throatIdx(DUMMY_INDEX);

			//if(numNetsInFront < numNetInSeries_-1 && outletPore)
				//getOutletData(itr, numNetsInFront, throatIdx, poreIdx);
			//else if(numNetsInFront > 0 && inletPore)
				//getInletData(itr, numNetsInFront, throatIdx, poreIdx);

			//if(poreIdx != DUMMY_INDEX)
			//{
				//connThroats.push_back(throatIdx);
				//++connNum;
			//}
			//++itr;
		//}
	//}
		}
	}
}


void InputData::throatData(int idx, int& pore1, int& pore2, double& vol, double& volCl, double& rad, double& shapeFact,
						   double& lenPore1, double& lenPore2, double& lenThroat, double& lenTot)
{
	int stdIdx(idx), numNetsInFront(0);
	if(idx > origNumThroats_)
	{
		numNetsInFront = 1+(idx-origNumThroats_-1)/(origNumThroats_-connectionsRemoved_);
		int tmpIdx = idx - numNetsInFront*origNumThroats_ + (numNetsInFront-1)*connectionsRemoved_;
		stdIdx = throatHash_[tmpIdx-1];
	}

	pore1 = throatData_[stdIdx-1]->poreOne;
	pore2 = throatData_[stdIdx-1]->poreTwo;
	vol = throatData_[stdIdx-1]->volume;
	volCl = throatData_[stdIdx-1]->clayVol;
	rad = throatData_[stdIdx-1]->radius;
	shapeFact = throatData_[stdIdx-1]->shapeFact;
	lenPore1 = throatData_[stdIdx-1]->lenPoreOne;
	lenPore2 = throatData_[stdIdx-1]->lenPoreTwo;
	lenThroat = throatData_[stdIdx-1]->lenThroat;
	lenTot = throatData_[stdIdx-1]->lenTot;

	if(idx > origNumThroats_)                      // Point to correct pores in the subsequent networks
	{
		if(pore1 != 0 && pore1 != -1)pore1 += numNetsInFront*origNumPores_;
		if(pore2 != 0 && pore2 != -1)pore2 += numNetsInFront*origNumPores_;
	}

	double totalLength(0.0);
	pair<MapItr,MapItr> mapPos;
	bool inlet(false), good2go(false);
	if(numNetsInFront != numNetInSeries_-1 && throatData_[stdIdx-1]->poreOne == 0)
	{
		good2go = true;
		mapPos = outletConnections_.equal_range(throatData_[stdIdx-1]->poreTwo);
	}
	else if(numNetsInFront != numNetInSeries_-1 && throatData_[stdIdx-1]->poreTwo == 0)
	{
		good2go = true;
		mapPos = outletConnections_.equal_range(throatData_[stdIdx-1]->poreOne);
	}
	else if(numNetsInFront > 0 && throatData_[stdIdx-1]->poreOne == -1)
	{
		good2go = true;
		inlet = true;
		mapPos = inletConnections_.equal_range(throatData_[stdIdx-1]->poreTwo);
	}
	else if(numNetsInFront > 0 && throatData_[stdIdx-1]->poreTwo == -1)
	{
		good2go = true;
		inlet = true;
		mapPos = inletConnections_.equal_range(throatData_[stdIdx-1]->poreOne);
	}

	if(good2go)
	{
		MapItr itr = mapPos.first;
		while(itr != mapPos.second)
		{
			if((*itr).second.first() == stdIdx)
			{
				int poreIdx = (*itr).second.second();

				if(inlet)
					poreIdx += (numNetsInFront-1)*origNumPores_;
				else
					poreIdx += (numNetsInFront+1)*origNumPores_;

				totalLength = (*itr).second.third();
				if(throatData_[stdIdx-1]->poreOne == 0 || throatData_[stdIdx-1]->poreOne == -1)
				{
					pore1 = poreIdx;
					lenPore1 = lenPore2;
				}
				else
				{
					pore2 = poreIdx;
					lenPore2 = lenPore1;
				}
			}
			++itr;
		}
		if(useAvrXOverThroatLen_)
		{
			lenThroat = averageThroatLength_;
			lenTot = 2.0*lenPore2 + lenThroat;
		}
		else if(totalLength > 2.0*lenPore2)
		{
			lenTot = totalLength;
			lenThroat = totalLength - 2.0*lenPore2;
		}
		else
		{
			lenPore1 = 1.0E-8;
			lenPore2 = 1.0E-8;
			lenThroat = 1.0E-8;
			lenTot = 3.0E-8;
		}
		double xSectArea(pow(rad, 2.0) / (4.0*shapeFact));  // Oren networks have 0 volume for in/outlet throats
		vol = lenThroat*xSectArea;                          // Assign volume based on a tube
	}
}




void InputData::solverDebug(bool& watMat, bool& oilMat, bool& resMat, bool& watVel, bool& oilVel, bool& resVel, bool& matlab, bool& initOnly)const
{

	istringstream data;
	string keyword("SOLVER_DBG");
	char writeWat, writeOil, writeRes, writeWatVel, writeOilVel, writeResVel, mat, forInit;

	if(getData(data, keyword))
	{
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> writeWat >> writeOil >> writeRes >> writeWatVel >> writeOilVel >> writeResVel >> mat >> forInit;
		watMat = (writeWat == 'T' || writeWat == 't');
		oilMat = (writeOil == 'T' || writeOil == 't');
		resMat = (writeRes == 'T' || writeRes == 't');
		watVel = (writeWatVel == 'T' || writeWatVel == 't');
		oilVel = (writeOilVel == 'T' || writeOilVel == 't');
		resVel = (writeResVel == 'T' || writeResVel == 't');
		matlab = (mat == 'T' || mat == 't');
		initOnly = (forInit == 'T' || forInit == 't');
		checkEndOfData(data, keyword);
		if(initOnly)
		{
			watMat = false;
			oilMat = false;
			resMat = false;
			watVel = false;
			oilVel = false;
			resVel = false;
		}
	}
	else
	{
		watMat = false;
		oilMat = false;
		resMat = false;
		watVel = false;
		oilVel = false;
		resVel = false;
		matlab = false;
		initOnly = false;
	}

}




