#ifndef ROCKELEM_H
#define ROCKELEM_H

//! Pore and Throats and their connectivity 


#include "assert.h"
#include <set>

//#include "sortedEvents.h"
//#include "elem_Model.h"
//
//#include "elem_porous.h"
//#include "elem_hetroPorous.h"
#include "apex.h"
class Fluid;
class ElemModel;
class Polygon;


#include "netsim_data.h"
 
class Fluid;
class commonData;

///////////////////////////////// Base Class //////////////////////////////////////





class Node
{
public:

    /**
	// Constructor for the node class. Inlet node has index == 0 whereas outlet
	// has the index has index == numPores + 1 If node is recognised to be inlet
	// or outlet reservoir the x location is moved outside the reservoir to prevent
	// them having the same x coord as the first line of pores (which usually are
	// at 0.0 or xSize. This is required for these to be recognised as in/outlet
	// if boxsize is set to be 1.0.
	*/
	Node(int index, int numPores, double xPos, double yPos, double zPos, double networkLength)
	: m_xPos(xPos), m_yPos(yPos), m_zPos(zPos),m_index(index),  m_networkLength(networkLength)	  //m_optimizedIndex(-1), //m_numPores(numPores),  //m_oldIndex(index),
	{
		m_isOutsideLattice = false;            
		m_isEntryRes = false;
		m_isExitRes = false;    

		if(m_index == 0)       					     m_isEntryRes = true;   // We are at inlet face
		else if(m_index == numPores + 1)             m_isExitRes = true;     // We are at outlet face
		else if(m_index < 0 || m_index > numPores+1) m_isOutsideLattice = true;// We are outside lattice
	};
	Node(){};
    bool operator==(const Node& rhs) const {return m_index == rhs.index();}
    bool operator>=(const Node& rhs) const {return m_index >= rhs.index();}
    bool operator<=(const Node& rhs) const {return m_index <= rhs.index();}
    bool operator>(const Node& rhs) const {return m_index > rhs.index();}
    bool operator<(const Node& rhs) const {return m_index < rhs.index();}

    inline bool optimizedIndex(int& optimizedIdx);
    int index() const {return m_index;}
    /**
	// In/outlet internally has indecies 0 and numPores+1. However when
	// writing the data to file this has to be changed to -1 and 0 to
	// be compatible with Oren's data format
	*/
	int indexOren() const	{
		if(m_isEntryRes)        return -1;
		else if(m_isExitRes)    return 0;
		else   			        return m_index;
	};
    double xPos() const {return m_xPos;}///. shifted to left by half of curtailed fraction
    double yPos() const {return m_yPos;}
    double zPos() const {return m_zPos;}
    inline void rePosition(double scaleFactor)
    {    m_xPos *= scaleFactor;    m_yPos *= scaleFactor;    m_zPos *= scaleFactor;    m_networkLength *= scaleFactor;	};

    bool isInsideBox(double boxStart, double boxEnd) const {
        return m_xPos >= m_networkLength*boxStart &&  m_xPos <=  m_networkLength*boxEnd;
    }

    bool isOutsideLattice() const {return m_isOutsideLattice;}
    bool isEntryRes() const {return m_isEntryRes;}
    bool isExitRes() const {return m_isExitRes;}
    bool isEntryOrExitRes() const {return m_isEntryRes || m_isExitRes;}

    double distToExit() const {return m_networkLength - m_xPos;}

    double                      m_xPos, m_yPos, m_zPos;     // The node coordinate ///. x_Pos shited to left by half of curtailed fraction
private:

    void initNode();

    int                         m_index;                    // Single consecutive index
    double                      m_networkLength;
    bool                        m_isEntryRes, m_isExitRes;  // Is the node inlet or outlet node
    bool                        m_isOutsideLattice;         // Is node outside lattice
    //int                         m_optimizedIndex;            // In order to reduce bandwidth

};


inline std::ostream& operator<< (std::ostream& out, const Node& node)
{
    out.flags(std::ios::showpoint);
    out.flags(std::ios::scientific);
    out << std::setprecision(4)          << std::setw(15) << node.m_xPos
        << std::setw(15) << node.m_yPos  << std::setw(15) << node.m_zPos;
    return out;
}




class Element : public Apex
{

    friend std::ostream& operator<< (std::ostream&, Element&) ;
    friend class solverFlags ;

public:

    Element(CommonData&, const Oil&, const Water&, double, double, double, double, int, bool, std::string);
    virtual ~Element();

    virtual double lenToRadRatio() const = 0;
    virtual void modifyLength(double scaleFactor) = 0;
    virtual const Node* node() const = 0;
    virtual Node* node() = 0;
    virtual bool crossesPlaneAt(double location) const = 0;
    virtual int latticeIndex() const = 0;
    virtual int orenIndex() const = 0;
    virtual int index2p() const = 0;
    virtual void addConnections(Element* first, Element* second, double inBdr, double outBdr,  bool moveBdr) = 0;
    virtual bool prevSolvrRes(const Fluid* fluid, int resistSolve, double loc, double& res, double& flowRate) const = 0;
    virtual void calcVolume_CheckIntegrity(double& vTot, double& vcTot, int& maxNonO, int& isoSum) const = 0;
    virtual void sortConnectingElems_DistToExit() = 0;
    virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const = 0;
    virtual void writeNetworkDataBinary(std::ostream& out) const = 0;
    virtual void updateLatticeIndex(int newIdx) = 0;
    
    virtual double snapOfLongitCurvature() const = 0 ;


    inline void adjustVolume(double newNetVol, double newClayVol, double& netVolSum, double& clayVolSum);

    void finalizeCopyConstruct(std::vector<Element*>& throats) {m_connections = throats;}
    void identifyConnectedPoreElems();
    int removeFromNetwork();
    unsigned int randomAssInt() const {return m_randomNum;}
    void severConnection(Element* connection);
    static void set_useGravInKr(bool soAreWe) {USE_GRAV_IN_KR = soAreWe;}

    void findMarkTrappedOilGanglia(double prs, std::vector<Element*>& stor, double& elap, TrappingCriteria crit);
    void findMarkTrappedWaterGanglia(double prs, FluidBlob startPt, std::vector< std::pair<Element*,FluidBlob> >& stor, double& elap, TrappingCriteria crit);
        



    inline double rhogh(double density, double xLen, double yLen, double zLen) const;
    double flowVolume() const {return m_flowVolume;}
    double clayVolume() const {return m_clayVolume;}
	double waterSaturation() const {return m_waterSaturation;}
	double volume() const {return m_flowVolume+m_clayVolume;}

    Element* connection(int conn) const {return m_connections[conn];}
    const std::vector<Element*>& connections() const {return m_connections;}
    ElemModel* ChModel() {return m_elemModel;}
    const ElemModel* model() const {return m_elemModel;}

    int connectionNum() const {return m_connectionNum;}



///. TO DELETE/REVISE:
    bool connectedToOutlet(const Fluid* fluid) const;
    void clearAllWSolverFlag() const;
    void clearAllOSolverFlag() const;
    inline bool canWBulkBePassedToSolver() const{return m_canWBulkBePassedToSolver;};
    inline bool canWFilmBePassedToSolver() const{return m_canWFilmBePassedToSolver;};
    inline bool canWBePassedToSolver() const {return m_canWFilmBePassedToSolver || m_canWBulkBePassedToSolver;};
	inline bool canBePassedToSolver(const Fluid* fluid) const {	return fluid->isOil() ? m_canOBulkBePassedToSolver :  m_canWBulkBePassedToSolver || m_canWFilmBePassedToSolver; }
    inline bool canOBePassedToSolver() const {return m_canOBulkBePassedToSolver;};
    bool isOnInletSlvrBdr() const {return m_isOnInletSlvrBdr;}
    void setOnInletSlvrBdr(bool isIt) {m_isOnInletSlvrBdr = isIt;}
    bool isOnOutletSlvrBdr() const {return m_isOnOutletSlvrBdr;}
    void setOnOutletSlvrBdr(bool isIt) {m_isOnOutletSlvrBdr = isIt;}
    bool isOnSlvrBdr() const {return m_isOnOutletSlvrBdr || m_isOnInletSlvrBdr;}
	bool isInsideSolverBox() const {return m_isInsideSolverBox;}
    bool isInsideSatBox() const {return m_isInsideSatBox;}

    bool isExitRes() const {return m_isExitRes;}
    bool isEntryRes() const {return m_isEntryRes;}
    bool isEntryOrExitRes() const {return m_isExitRes || m_isEntryRes;}
    bool connectedToEntryOrExit() const {return m_connectedToEntryOrExit;}
    bool connectedToNetwork() const {return m_connectedToNetwork;}
    void setConnectedToNetwork(bool isIt){m_connectedToNetwork = isIt;}
    bool isConnectedToExit() const {return m_isConnectedToExit;}
    bool isConnectedToEntry() const {return m_isConnectedToEntry;}
    void isConnectedToExit(bool isIt) {m_isConnectedToExit = isIt;}
    void isConnectedToEntry(bool isIt) {m_isConnectedToEntry = isIt;}
	bool convertToMicroPorosityForSven(bool entry1Exit0);
    bool iAmAPore() const {return m_iAmAPore;}
    int iRockType() const {return m_iRockType;}



    void IncreaseNumOilCentreFeederNeis() {++m_numOilCentreFeederNeis;}
    void ReduceNumOilCentreFeederNeis() {--m_numOilCentreFeederNeis;}    
    void IncreaseNumWatCentreFeederNeis() {++m_numWatCentreFeederNeis;}
    void ReduceNumWatCentreFeederNeis() {--m_numWatCentreFeederNeis;}    
    int numOilCentreFeederNeis() const {return m_numOilCentreFeederNeis;}
    int num_WatCentreFeederNeis() const {return m_numWatCentreFeederNeis;}

    bool isTrappedOil() const {return m_trappingIndexOil.first > -1;}
    bool nonTrappedOil() const {return m_trappingIndexOil.first < 0;}
    void unTrapOil() {m_trappingIndexOil.first = -1;}
    void fillElemCentreWithOilRemoveLayers();
    int trappingIndexOil() const {return m_trappingIndexOil.first;}

	inline bool isTrappedWat(FluidBlob blob) const { return (blob == filmBlob) ?	m_trappingIndexWatFilm.first > -1   :  m_trappingIndexWatBulk.first > -1;}
	inline bool nonTrappedWat(FluidBlob blob) const { return (blob == filmBlob) ?	m_trappingIndexWatFilm.first < 0   :  m_trappingIndexWatBulk.first < 0;}
    inline void unTrapWat(FluidBlob blob);
    void fillElemCentreWithWaterCreateLayers(bool snapOffOverRide = false);
    bool canBeAddedToEventVec(const Fluid* injectant) const;



    double updateSat_calcR(double cappPrs);





	///.  rare funcs
    const std::pair<int, double>& trappingOil() const {return m_trappingIndexOil;}    
    inline const std::pair<int, double>& trappingWat(FluidBlob blob) const;
    inline void setWatFilmTrappingFromBulk();
    inline int trappingIndexWat(FluidBlob blob) const;
    bool addToLayerVec(const Fluid* injectant, Polygon* polyShape, std::vector<int>& addCrn) const;
    const std::pair<int, double>& trappingWatBulk() const {return m_trappingIndexWatBulk;}
    const std::pair<int, double>& trappingWatFilm() const {return m_trappingIndexWatFilm;}


    inline void setGravityCorrection(const Node* node);
    inline void setGravityCorrection(const Node* nodeOne, const Node* nodeTwo);

	void calcCentreEntryPrsWatInj();
	void calcCentreEntryPrsOilInj();


    ///. rubbish func
    int fillingEventRecord() const {return m_fillingEventRecord;}
    void resetFillingEventRecord() {m_fillingEventRecord = -2;}
    
    static void conductanceCutOff(double cutoff){COND_CUT_OFF=cutoff;};
    static int             nErrs;

    double RRR() const;

	int ffaz() const;  ///. Viz only

	///valid for throats only
    mutable double   m_conductance[5]; ///. Viz only
	void setPoreToPoreCond(double cond, const Fluid *fluid, int resistiv) {m_poreToPoreCond = cond; m_conductance[int(fluid->isOil())+2*resistiv]=cond; }
	double saturation() const ;
protected:

    //virtual void printData(std::ostream& out) const = 0;


///.  trapping search
    bool foundEscapePathOil_trapOtherwise(double pc, std::vector<Element*>& stor, TrappingCriteria crit);
		inline Element* nextUntrappedOil(TrappingCriteria criteria);

		inline Element* nextSuccessorWat(TrappingCriteria criteria, FluidBlob& blob);
    bool foundEscapePathWat_trapOtherwise(double pc, FluidBlob startPt, std::vector< std::pair<Element*, FluidBlob> >& stor,
        TrappingCriteria criteria); 
           
    inline void trapOil(double prs);
    inline void trapWat(double prs, FluidBlob blob);    
    

///.  connected path search for solver 
    bool markWaterElemForSolver(std::pair<const Element*, FluidBlob> elem) const;
    inline const Element* nextSuccSolvrWat(bool& outletFound, FluidBlob& blob) const;
    inline const Element* nextSuccSolvrOil(bool& outletFound) const;    
    
    inline void clearWSolverFlag(FluidBlob blob = bulkBlob) const;
    inline void setWSolverFlag(FluidBlob blob = bulkBlob) const;    
    inline void clearOSolverFlag(FluidBlob blob = bulkBlob) const;
    inline void setOSolverFlag(FluidBlob blob = bulkBlob) const;    
    
    inline bool isWFilmTouchedInSearch() const {return m_isWFilmTouchedInSearch;};
    inline bool isWBulkTouchedInSearch() const {return m_isWBulkTouchedInSearch;};
    inline bool isWTouchedInSearch(FluidBlob blob) const {return blob == bulkBlob ? m_isWBulkTouchedInSearch :  m_isWFilmTouchedInSearch;};
    inline bool isOTouchedInSearch() const {return m_isOBulkTouchedInSearch;};



    void checkConnections() const; ///.  used in calcVolume_CheckIntegr

    typedef std::set<Element*>::iterator ItrSet;

    //static bool                             ERROR_STATE;
    static bool                             USE_GRAV_IN_KR;
    static double                           COND_CUT_OFF;
    static const double                     PI;

    const unsigned int                      m_randomNum;
    const bool                              m_iAmAPore;
    int                              		m_iRockType;
    const CommonData&                       m_comn;
    double                                  m_flowVolume;
    double                                  m_clayVolume;
    int                              		m_connectionNum;

    std::vector<Element*>                 m_connections;
    ElemModel*                              m_elemModel;
	double									m_waterSaturation;
    //double                                  m_averageAspectRatio;
    //double                                  m_minAspectRatio;
    //double                                  m_maxAspectRatio;
    
    bool                                    m_isInsideSolverBox;
    bool                                    m_isInsideSatBox;
    bool                                    m_isOnInletSlvrBdr;
    bool                                    m_isOnOutletSlvrBdr;
    bool                                    m_isExitRes;
    bool                                    m_isEntryRes;

    mutable bool                m_canWFilmBePassedToSolver;
	mutable bool                m_canWBulkBePassedToSolver;
    // mutable bool                m_canOFilmBePassedToSolver;
    mutable bool                m_canOBulkBePassedToSolver;
    mutable bool                m_isWFilmTouchedInSearch;
    mutable bool                m_isWBulkTouchedInSearch;
    // mutable bool                m_isOFilmTouchedInSearch;
    mutable bool                m_isOBulkTouchedInSearch;


    mutable double                          m_poreToPoreCond;///HASAN
       
    //bool                                    m_isInWatFloodVec;
    //bool                                    m_isInOilFloodVec;
    bool                                    m_connectedToNetwork;
    bool                                    m_isConnectedToExit;
    bool                                    m_isConnectedToEntry;
    bool                                    m_connectedToEntryOrExit;
    

    std::pair<int, double>                 m_trappingIndexOil;
    std::pair<int, double>                 m_trappingIndexWatBulk;
    std::pair<int, double>                 m_trappingIndexWatFilm;
    int                                    m_numOilCentreFeederNeis;
    int                                    m_numWatCentreFeederNeis;

///. rubish vars
    int                                     m_fillingEventRecord;
};

///////////////////////////////// Pore Class //////////////////////////////////////

class Pore : public Element
{
public:

    Pore(CommonData&, Node*, const Oil&, const Water&, double, double, double,
        double, bool, bool, double, std::vector<Element*>&, std::string);
    //Pore(CommonData&, const Fluid*, const Fluid*, const Pore&);
    virtual ~Pore(){delete m_node;}

    virtual double lenToRadRatio() const {return 0.0;}
    virtual void modifyLength(double scaleFactor) {}
    virtual const Node* node() const {return m_node;}
    virtual Node* node() {return m_node;}
    virtual bool crossesPlaneAt(double location) const {return false;}
    virtual int index() const {return m_node->index();}
    virtual int latticeIndex() const {return m_node->index();}
    virtual int orenIndex() const {return m_node->indexOren();}
    virtual int index2p() const {return m_node->indexOren()+1;}
    virtual void addConnections(Element* first, Element* second, double inBdr, double outBdr, bool moveBdr);
    virtual bool prevSolvrRes(const Fluid* fluid, int resistSolve, double loc, double& res, double& flowRate) const;
    virtual void calcVolume_CheckIntegrity(double& vTot, double& vcTot, int& maxNonO, int& isoSum) const;
    virtual void sortConnectingElems_DistToExit();
    virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const;
    virtual void writeNetworkDataBinary(std::ostream& out) const;
    virtual void updateLatticeIndex(int newIdx) {}
    //virtual void setRadiusFromAspectRatio(double aspectRadius);


    virtual double snapOfLongitCurvature() const;


	Element* getConnectionProp(int conn, double& conductance, double& deltaGrav, const Fluid *fluid, bool resistivitySolve) const;
	Element* getConnectionPropSP(int conn, double& conductance, double& rhs_throat, double& deltaGrav, const Fluid *fluid, bool resistivitySolve) const;

    void setSolverResults(const Fluid* fluid, int resistSolve, double res) const;

protected:

    //virtual void printData(std::ostream& out) const;

    Node                           *m_node;
    mutable double                  m_oilSolverPrs;
    mutable double                  m_watSolverPrs;
    mutable double                  m_oilSolverVolt;
    mutable double                  m_watSolverVolt;

};

///////////////////////////////// Inlet & Outlet Pore /////////////////////////////////////

class InOutBoundaryPore : public Pore
{
public:

    InOutBoundaryPore(CommonData&, Node*, const Oil&, const Water&, std::vector<Element*>&);
    //InOutBoundaryPore(CommonData& common, const Fluid* oil, const Fluid* water, const InOutBoundaryPore& endP) : Pore(common, oil, water, endP) {}
    
    void fillElemCentreWithOilRemoveLayersIO(double pc);
    void fillElemCentreWithWaterCreateLayersIO(double pc);

    
    virtual ~InOutBoundaryPore(){}

    virtual void calcVolume_CheckIntegrity(double& vTot, double& vcTot, int& maxNonO, int& isoSum) const;

};

//////////////////////////////////////// Throat //////////////////////////////////////////

class Throat : public Element
{
public:

    Throat(CommonData&, const Oil&, const Water&, double, double, double, double,
        double, double, double, int, std::string);
    //Throat(CommonData&, const Fluid*, const Fluid*, const Throat&);

    virtual double lenToRadRatio() const;
    virtual void modifyLength(double scaleFactor);
    virtual const Node* node() const;
    virtual Node* node();
    virtual bool crossesPlaneAt(double location) const;
    virtual int index() const {return m_latticeIndex;}
    virtual int index2p() const {return m_latticeIndex;}
    virtual int latticeIndex() const {return m_latticeIndex;}
    virtual int orenIndex() const {return m_latticeIndex-m_comn.numPores()-1;}
    virtual void addConnections(Element* first, Element* second, double inletBdr, double outletBdr, bool moveBoundary);
    virtual bool prevSolvrRes(const Fluid* fluid, int resistSolve, double loc, double& res, double& flowRate) const;
    virtual void calcVolume_CheckIntegrity(double& vTot, double& vcTot, int& maxNonO, int& isoSum) const;
    virtual void sortConnectingElems_DistToExit();
    virtual void writeNetworkData(std::ostream& outOne, std::ostream& outTwo) const;
    virtual void writeNetworkDataBinary(std::ostream& out) const;
    virtual void updateLatticeIndex(int newIdx) {m_latticeIndex = newIdx;}
    //virtual void setRadiusFromAspectRatio(double aspectRadius);
	

	double snapOfLongitCurvature() const;

	const Element* neighbouringPore(const Element* callingPore) const;
	double poreLength(const Element* callingPore) const;

	double poreLength(int conn) const {return m_poreLength[conn];}
	double length() const {return m_length;}

private:

    //virtual void printData(std::ostream& out) const;
    Node                           m_node;

    int                                 m_latticeIndex;
    double                             m_originalPoreLengths[3];
    std::vector< double >               m_poreLength;
    double                              m_length;
};

//////////////////////////////////////////////////////////////////////////////////////////

/**
 * adjust water, clay and pore volume ++
 */
inline void Element::adjustVolume(double newNetVol, double newClayVol,
                                   double& netVolSum, double& clayVolSum)
{
    m_flowVolume = newNetVol;
    m_clayVolume = newClayVol;

    if(m_isInsideSatBox)
    {
        netVolSum += m_flowVolume;
        clayVolSum += m_clayVolume;
    }
}

inline const std::pair<int, double>& Element::trappingWat(FluidBlob blob) const
{
    if(blob == filmBlob)
        return m_trappingIndexWatFilm;
    else
        return m_trappingIndexWatBulk;
}




inline double Element::rhogh(double density, double xLen, double yLen, double zLen) const
{
    return density*(m_comn.gravConstX()*xLen+m_comn.gravConstY()*yLen+m_comn.gravConstZ()*zLen);
}

inline void Element::setGravityCorrection(const Node* node)
{
    m_gravityCorrection = (m_comn.water().density()-m_comn.oil().density())*
        (m_comn.gravConstX()*node->xPos() +
        m_comn.gravConstY()*node->yPos() +
        m_comn.gravConstZ()*node->zPos());
        
        m_gravityCorrection = 0.0;///. ERROR

}

inline void Element::setGravityCorrection(const Node* nodeOne, const Node* nodeTwo)
{
    double xPos((nodeOne->xPos()+nodeTwo->xPos())/2.0);
    double yPos((nodeOne->yPos()+nodeTwo->yPos())/2.0);
    double zPos((nodeOne->zPos()+nodeTwo->zPos())/2.0);

    m_gravityCorrection = (m_comn.water().density()-m_comn.oil().density())*
        (m_comn.gravConstX()*xPos +
        m_comn.gravConstY()*yPos +
        m_comn.gravConstZ()*zPos);
        
        m_gravityCorrection = 0.0;///. ERROR
}








inline int Element::trappingIndexWat(FluidBlob blob) const
{
    if(blob == filmBlob)
        return m_trappingIndexWatFilm.first;
    else
        return m_trappingIndexWatBulk.first;
}


inline void Element::setWatFilmTrappingFromBulk()
{
    ensure(m_trappingIndexWatBulk.first == -1);
    m_trappingIndexWatFilm = m_trappingIndexWatBulk;
}


inline void Element::unTrapWat(FluidBlob blob)
{
    if(blob == bulkBlob)
		m_trappingIndexWatBulk.first = -1;
	else
		m_trappingIndexWatFilm.first = -1;
}

#endif
