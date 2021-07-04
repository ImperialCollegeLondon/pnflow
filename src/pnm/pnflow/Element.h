#ifndef ELEMENT_CnM_H
#define ELEMENT_CnM_H

//! Pore and Throats and their connectivity 


#include <set>

#include "typses.h"
#include "fluid.h"

#define OutI  1  //! outlet pore index

 

class RockType;


class ElemModel;
class Pore;
class Throat;
class InOutBoundary;
class Polygon;
class CommonData;


enum TrappingCriteria {escapeToInlet = 0, escapeToOutlet, escapeToEither, escapeToBoth};
enum FluidBlob {filmBlob = 0, bulkBlob};



enum BounCond  : unsigned char
{
	NOTSET    = 0,
	INLET     = 1,
	INSIDE    = 2,
	OUTLET    = 4
};


class SlvrElm
{///.  connected path search for solver
public:
	inline void clearFlag() const { _passed = false; }
	inline void setAllFlags() const {	_passed = true;    _searched = true;   }

	void clearAllFlags() const { _passed = false;   _searched = false; };
	inline bool isPassed() const {return _passed;};
	inline bool isSearched() const {return _searched;};
	//bool isInSlvrBox() const {return _SIDE==INSIDE;}
	//bool isRightBC() const {return _SIDE==OUTLET;}
	//bool isLeftBC() const {return  _SIDE==INLET;}
	//void setBC(BounCond bc) {_SIDE=bc;}
	//const std::vector<Elem*>& neis() const;

	SlvrElm()
	{
		_SIDE = NOTSET;
		_SIDE = INSIDE;
		_passed = false;
		_searched = false;
		_condOrP=-10000.;
		//elem_=NULL;
		//nAzadPaths_=0;
	}


	//Elem*            elem_;
	BounCond            _SIDE;
	mutable bool        _passed;
	mutable bool        _searched;
	//short               nAzadPaths_;
	double              _condOrP;
	double              _condOrP2;
};



class ElemModel;
//////////////////////// BASE CLASS /////////////////////////////
class Apex
{
public:

	Apex() :  virgin(true),parentModel_(NULL), advancingPc_(1000.), receedingPc_(-1000.),
	entryPc_(0.), gravityCorrection_(0.), exists_(false),inited_(false) ,trappedCL_(-1, 0.) {}
	virtual ~Apex() {}

	void setConnections(ElemModel* parent, int subIndex){parentModel_=(parent);subIndex_=(subIndex);}

	inline double advancingPc() const { return advancingPc_; };
	inline double receedingPc() const { return receedingPc_; };
	const std::pair<int, double>& trappingCL() const {return trappedCL_;}
	int                           iTrap() const {return trappedCL_.first;}

	inline int subIndex() const {return subIndex_;};
	ElemModel*  parentModel() const {return parentModel_;};
	double gravCorrectedEntryPress() const {return entryPc_+gravityCorrection_;}
	double gravityCorrection() const {return gravityCorrection_;}
	double entryPc() const {return entryPc_;}

	bool isInWatFloodVec() const  {return isInWatFloodVec_;}
	void setInWatFloodVec(bool isIt) {isInWatFloodVec_ = isIt;}
	bool isInOilFloodVec() const  {return isInOilFloodVec_;}
	void setInOilFloodVec(bool isIt) {isInOilFloodVec_ = isIt;}

	bool exists() const {return exists_;}
	bool pinned() const {return inited_;}

	mutable bool virgin;
	double      					trapPcOld_;
	double      					creationPc;

	static int nErrors;


protected:


	static const double             PI;
	static const double             INF_NEG_NUM;
	static const double             LOWEST_LAYER_PC;
	static const double             INF_POS_NUM;
	static const double             EPSILON;
	static const double             SMALL_NUM;
	static const double             NEG_ALMOST_ZERO;
	static const double             POS_ALMOST_ZERO;
	static const double             MOLECULAR_LENGTH;
	static const int                MAX_ITR;

	ElemModel*                     parentModel_;
	int                            subIndex_;

	double                         advancingPc_;
	double                         receedingPc_;
	double                          entryPc_;
	double                          gravityCorrection_;

	bool                            isInWatFloodVec_;
	bool                            isInOilFloodVec_;

	bool                            exists_;
	bool                            inited_;
	std::pair<int, double>         trappedCL_;


};


class Elem : public Apex
{

	friend std::ostream& operator<< (std::ostream&, Elem&) ;
	friend class solverFlags ;

public:

	using TrotT = Throat;
	using PoreT = Pore;
	using BondT = InOutBoundary;
	using ShapT = Polygon;

	Elem(const CommonData&, int, dbl3, double, double, double, double, int, bool, int);
	virtual ~Elem();

	virtual void modifyLength(double scaleFactor) = 0;
	virtual bool crossesPlaneAt(double location) const = 0;
	virtual int indexOren() const = 0;
	virtual bool prevSolvrRes(const Fluid& fluid, double loc, double& res, double& flowRate) const = 0;
	virtual void prepare2() = 0;
	virtual void sortConnectingElems_DistToExit() = 0;
	virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const = 0;
	virtual void writeNetworkDataBinary(std::ostream& out) const = 0;
	virtual void updateLatticeIndex(int newIdx) = 0;
	virtual double snapOfLongitCurvature() const = 0 ;

	virtual int index() const { return index_; };

	inline void adjustVolume(double newNetVol, double newClayVol);

	void finalizeCopyConstruct(std::vector<Elem*>& throats) {cnctions_ = throats;}
	void identifyConnectedPoreElems();
	int removeFromNetwork();
	void severConnection(Elem* connection);
	static void set_useGravInKr(bool soAreWe) {USE_GRAV_IN_KR = soAreWe;}
	double rhogh(double rho, dbl3 dh) const;// {return comn_.rhogh(rho,dh);}

	void findMarkTrappedOilGanglia(double prs, std::vector<Elem*>& stor, double& elap, TrappingCriteria crit);
	void findMarkTrappedWaterGanglia(double prs, FluidBlob startPt, std::vector< std::pair<Elem*,FluidBlob> >& stor, double& elap, TrappingCriteria crit);




	//inline double rhogh(double density, dbl3 distance) const;
	double flowVolume() const {return flowVolume_;}
	double clayVolume() const {return clayVolume_;}
	double waterSaturation() const {return waterSaturation_;}
	double flowVolumeX() const {return flowVolume_+clayVolume_;} 
	double saturation() const  {return waterSaturation_;}

	Elem* neib(int conn) const {return cnctions_[conn];}
	const std::vector<Elem*>& connections() const {return cnctions_;}
	ElemModel* ChModel() {return model_;}
	const ElemModel* model() const {return model_;}

	int nCncts() const {return nCncts_;}



///. TO DELETE/REVISE:
	bool connectedToOutlet(const Fluid& fluid) const;

	bool isOnInletSlvrBdr() const {return isOnInletSlvrBdr_;}
	void setOnInletSlvrBdr(bool isIt) {isOnInletSlvrBdr_ = isIt;}
	bool isOnOutletSlvrBdr() const {return isOnOutletSlvrBdr_;}
	void setOnOutletSlvrBdr(bool isIt) {isOnOutletSlvrBdr_ = isIt;}
	bool isOnSlvrBdr() const {return isOnOutletSlvrBdr_ || isOnInletSlvrBdr_;}
	bool isInsideSolverBox() const {return isInsideSolverBox_;}
	bool isInCalcBox() const {return isInCalcBox_;}

	bool isExitRes() const {return isExitRes_;}
	bool isEntryRes() const {return isEntryRes_;}
	bool isEntryOrExitRes() const {return isExitRes_ || isEntryRes_;}
	bool connectedToEntryOrExit() const {return connectedToEntryOrExit_;}
	bool connectedToNetwork() const {return connectedToNetwork_;}
	void setConnectedToNetwork(bool isIt){connectedToNetwork_ = isIt;}
	bool isConnectedToExit() const {return isConnectedToExit_;}
	bool isConnectedToEntry() const {return isConnectedToEntry_;}
	void isConnectedToExit(bool isIt) {isConnectedToExit_ = isIt;}
	void isConnectedToEntry(bool isIt) {isConnectedToEntry_ = isIt;}
	bool convertToMicroPorosityForSven(bool entry1Exit0);
	bool iAmAPore() const {return iAmAPore_;}
	int rockIndex() const {return rockIndex_;}
	virtual const RockType* rockType() const { return nullptr; };



	void IncreaseNumOilCentreFeederNeis() {++numOilCentreFeederNeis_;}
	void ReduceNumOilCentreFeederNeis() {--numOilCentreFeederNeis_;}    
	void IncreaseNumWatCentreFeederNeis() {++numWatCentreFeederNeis_;}
	void ReduceNumWatCentreFeederNeis() {--numWatCentreFeederNeis_;}    
	int numOilCentreFeederNeis() const {return numOilCentreFeederNeis_;}
	int num_WatCentreFeederNeis() const {return numWatCentreFeederNeis_;}

	bool isTrappedOil() const {return trapIndexOil_.first > -1;}
	bool nonTrappedOil() const {return trapIndexOil_.first < 0;}
	void unTrapOil() {trapIndexOil_.first = -1;}
	void fillElemCentreWithOilRemoveLayers();
	int trapIndexOil() const {return trapIndexOil_.first;}

	inline bool isTrappedWat(FluidBlob blob) const { return (blob == filmBlob) ?	trapIndexWatFilm_.first > -1   :  trapIndexWatBulk_.first > -1;}
	inline bool notTrapdW(FluidBlob blob) const { return (blob == filmBlob) ?	trapIndexWatFilm_.first < 0   :  trapIndexWatBulk_.first < 0;}
	inline void unTrapWat(FluidBlob blob);
	void fillElemCentreWithWaterCreateLayers(bool snapOffOverRide = false);
	bool canBeAddedToEventVec(const Fluid& injectant) const;



	double updateSat_calcR(double cappPrs);



	const dbl3& node() const {return node_;}
	dbl3& node() {return node_;}


	///.  rare funcs
	const std::pair<int, double>& trappingOil() const {return trapIndexOil_;}    
	inline const std::pair<int, double>& trappingWat(FluidBlob blob) const;
	inline void setWatFilmTrappingFromBulk();
	inline int trapIndexWat(FluidBlob blob) const;
	bool addToLayerVec(const Fluid& injectant, Polygon* shyp, std::vector<int>& addCrn) const;
	const std::pair<int, double>& trappingWatBulk() const {return trapIndexWatBulk_;}
	const std::pair<int, double>& trappingWatFilm() const {return trapIndexWatFilm_;}


	void setGravityCorrection(const dbl3& node);
	void setGravityCorrection(const dbl3& nodeOne, const dbl3& nodeTwo);

	void calcCentreEntryPrsWatInj();
	void calcCentreEntryPrsOilInj();


	///. rubbish func
	void resetEventIndex() {eventI_ = -2;}
	char eventI() const {return eventI_;}

	static void conductanceCutOff(double cutoff){COND_CUT_OFF=cutoff;};
	static int             nErrs;

	double RRR() const;

	int ffaz() const;  ///. Viz only

	///valid for throats only
	void setPoreToPoreCond(fluidf ff, double cond) const {  WOISolver[ff-1]._condOrP=cond; }
	double poreToPoreCond(fluidf ff) const {return WOISolver[ff-1]._condOrP; }

	const SlvrElm& slvrCnct(fluidf ff) const {return WOISolver[ff-1];}
	SlvrElm& solverConectCh(fluidf ff) {return WOISolver[ff-1];}

protected:

	mutable SlvrElm    WOISolver[5];///    0: Pw,   1: Po,   2 Volt, not used yet, sync cnflow&gnflow

///.  trapping search
	bool foundEscapePathOil_trapOtherwise(double pc, std::vector<Elem*>& stor, TrappingCriteria crit);
		inline Elem* nextUntrappedOil(TrappingCriteria criteria);

		inline Elem* nextSuccessorWat(TrappingCriteria criteria, FluidBlob& blob);
	bool foundEscapePathWat_trapOtherwise(double pc, FluidBlob startPt, std::vector< std::pair<Elem*, FluidBlob> >& stor, TrappingCriteria criteria); 
		   
	inline void trapOil(double prs);
	inline void trapWat(double prs, FluidBlob blob);    


///.  connected path search for solver 
	bool markWaterElemForSolver(const Elem*elem, fluidf ff) const;
	inline const Elem* nextSuccSolvrWat(bool& outletFound, fluidf ff) const; // water and electricity
	inline const Elem* nextSuccSolvrOil(bool& outletFound) const;    


	void checkConnections() const; ///.  used in calcVolume_CheckIntegr

	typedef std::set<Elem*>::iterator ItrSet;

	static bool                 USE_GRAV_IN_KR;
	static double               COND_CUT_OFF;
	static const double         PI;

	const bool                  iAmAPore_;
	int                         rockIndex_;
	const CommonData&           comn_;
	double                      flowVolume_;
	double                      clayVolume_;
	int                         nCncts_;

	std::vector<Elem*>       cnctions_;
	ElemModel*                  model_;
	double                      waterSaturation_;

	bool                        isInsideSolverBox_;
	bool                        isInCalcBox_;
	bool                        isOnInletSlvrBdr_;
	bool                        isOnOutletSlvrBdr_;
	bool                        isEntryRes_,  isExitRes_;


	bool                        connectedToNetwork_;
	bool                        isConnectedToExit_;
	bool                        isConnectedToEntry_;
	bool                        connectedToEntryOrExit_;


	std::pair<int, double>      trapIndexOil_;
	std::pair<int, double>      trapIndexWatBulk_;
	std::pair<int, double>      trapIndexWatFilm_;
	int                         numOilCentreFeederNeis_;
	int                         numWatCentreFeederNeis_;


	int                        eventI_;

	int                        index_;
	dbl3                       node_;

};

///////////////////////////////// Pore Class //////////////////////////////////////
class Pore : public Elem
{
public:

	Pore(const CommonData&, int, dbl3, double, double, double, double, bool, bool, double, std::vector<Elem*>&, int);
	virtual ~Pore(){}

	virtual void modifyLength(double scaleFactor) {}
	virtual bool crossesPlaneAt(double location) const {return false;}
	virtual int indexOren() const { return index()-1; }
	virtual bool prevSolvrRes(const Fluid& fluid, double loc, double& res, double& flowRate) const;
	virtual void prepare2();
	virtual void sortConnectingElems_DistToExit();
	virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const;
	virtual void writeNetworkDataBinary(std::ostream& out) const;
	virtual void updateLatticeIndex(int newIdx) {}

	virtual double snapOfLongitCurvature() const;


	Elem* getConnectionPropSP(int conn, double& conductance, double& rhs_throat, const Fluid& fluid) const; // TODO: replace with Throat::calcR2SP

	void setSolverPrs(const Fluid& fluid, double res) const;

	Throat*       neiT(int i);
	const Throat* neiT(int i) const;
};




//////////////////////////////////////// Throat //////////////////////////////////////////
class Throat : public Elem
{
public:

	Throat(const CommonData&, int indx, dbl3, double, double, double, double, double, double, double, int);
	void addConnections(Elem* first, Elem* second, double inletBdr, double outletBdr, bool moveBoundary, bool setNod);

	virtual void modifyLength(double scaleFactor);
	virtual bool crossesPlaneAt(double location) const;
	virtual int index() const {return index_;}
	virtual int indexOren() const;// { return index_-comn_.numPores()-1; }
	virtual bool prevSolvrRes(const Fluid& fluid, double loc, double& res, double& flowRate) const;
	virtual void prepare2();
	virtual void sortConnectingElems_DistToExit();
	virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const;
	virtual void writeNetworkDataBinary(std::ostream& out) const;
	virtual void updateLatticeIndex(int newIdx) {index_ = newIdx;}

	void  neiSet(int i, Pore* por)  { neiPores_[i] = por;} // cnctions_[i] = por; 
	Pore*       neiPCh(int i)       { return neiPores_[i]; }
	const Pore*	neiP(int i) const   { return neiPores_[i]; }
	//int ihOf(const Elem* neip) const  {return (neiPores_[1]==neip)? 1 : ((neiPores_[0]==neip)? 0: 2);}

	void calcR2(const Fluid& fluid);
	const Pore* neighbouringPore(const Elem* callingPore) const;

	double snapOfLongitCurvature() const;


	double poreLength(int conn) const {return poreLength_[conn];}
	double throatLength() const {return length_;}

private:

	int                          index_;
	Pore*                        neiPores_[2]; 
	std::array<double,2>         poreLength_;
	double                       length_;
};



///////////////////////////////// Inlet & Outlet Pore /////////////////////////////////////

class InOutBoundary : public Pore
{
public:

	InOutBoundary(const CommonData&, int, const dbl3, std::vector<Elem*>&);
	~InOutBoundary();

	void prepare(const std::vector<Elem*>& connThroats, dbl2 calcBox);

	void fillElemCentreWithOilRemoveLayersIO(double pc);
	void fillElemCentreWithWaterCreateLayersIO(double pc);

	int   rockIndex() const {return -1;}


	virtual void prepare2();

	std::vector<Pore*>  miroredPores_;
};


//////////////////////////////////////////////////////////////////////////////////////////

inline Throat*       Pore::neiT(int i)       { return dynamic_cast<      Throat*>(cnctions_[i]); }
inline const Throat*	Pore::neiT(int i) const { return dynamic_cast<const Throat*>(cnctions_[i]); }
//inline const Pore*	Pore::neiP(int i) const { return dynamic_cast<const Pore*>( cnctions_[i]->neib(0)==this ? cnctions_[i]->neib(1) : cnctions_[i]->neib(0) ); }

/// adjust water, clay and pore volume ++
inline void Elem::adjustVolume(double newNetVol, double newClayVol)
{
	if(newNetVol>=0.) flowVolume_ = newNetVol;
	clayVolume_ = newClayVol;
}

inline const std::pair<int, double>& Elem::trappingWat(FluidBlob blob) const
{
	if(blob == filmBlob)  return trapIndexWatFilm_;
	else                  return trapIndexWatBulk_;
}






inline int Elem::trapIndexWat(FluidBlob blob) const
{
	return (blob==filmBlob) ? trapIndexWatFilm_.first : trapIndexWatBulk_.first;
}


inline void Elem::setWatFilmTrappingFromBulk()
{
	ensure(trapIndexWatBulk_.first == -1);
	trapIndexWatFilm_ = trapIndexWatBulk_;
}


inline void Elem::unTrapWat(FluidBlob blob)
{
	if(blob == bulkBlob)  trapIndexWatBulk_.first = -1;
	else                  trapIndexWatFilm_.first = -1;
}



#endif

