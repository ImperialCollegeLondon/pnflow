#ifndef LAYERAPEX_H
#define LAYERAPEX_H


//! Oil layer connectivity and interface tracking 

#include <algorithm>
#undef max
#undef min
class Polygon;



/////////////////////// LAYER FILM /////////////////////////////////////////
class LayerApex : public Apex
{
public:

    LayerApex() : m_initedOLApexDist(-1.0) {};

    //LayerApex(CornerApex* innerCornerApex, Polygon* parent, int subIndex) 
		//: Apex(parent, subIndex),  m_parentShape(parent), /*m_lastStablePc(0.0),*/, m_innerCornerApex(innerCornerApex)


    void setLayerConnections(CornerApex* innerCornerApex, Polygon* parent, int subIndex)
    { m_innerCornerApex=(innerCornerApex);m_parentShape=(parent);setConnections(parent, subIndex);}

    virtual ~LayerApex() {}

    bool createOLayer(double pc, double recAng, double advAng, double maxSpontConAng, double halfAng, double intTen, bool isOilInj);

    bool initLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intTen, bool oilInj, bool silent=false);
    bool finitLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intTen, bool oilInj, bool overRideTrp);

    void getCAApexDistUntraped(double& apxDist, double& conA, const double& halfA, double pc, double ten, bool itr = false) const;
    void getCAApexDist(double& apxDist, double& conA, const double& halfA, double pc, double ten, bool itr = false) const;

    inline void markCentreLayerTrappings(const std::pair< int, double >&, double, double, double, double, double, bool);
    double layerCollPc() const {return m_entryPc;}
    inline void removeLayer();
    const std::pair<int, double>& trappedOLayer() const {return m_trappedCL;}
    inline void set_m_exists(bool isIt){ m_exists = isIt; }
    inline bool freeAtPrs(double pc) const;
    inline bool forcedSnapOff(double prs) const;
	inline int index() const {return m_subIndex;};

    void advConAng(double conAng) {m_advConAng = conAng;}
    inline bool stablePinnedInLastCycle(double minPcLastCycle) const {return (exists(/*st ab le*/) && /*pinned () &&*/ m_advancingPc == minPcLastCycle);}

    double pinnedApexDist() const {return m_initedOLApexDist;}
    double layerCollapsePc(double pc, double conAng, double halfAng, double intfacTen, bool injOil) const;
    double layerCollapsePc_fromEitherSide(double pc, double conAng, double halfAng, double intfacTen,bool debug=false) const;

    mutable int m_colType;
private:

    double layerCollapsePc_FromCentre(double outPc, double inPc, double conAng, double halfAng, double ten) const;
    double layerCollapsePc_FromCorner(double outPc, double inPc, double conAng, double halfAng, double ten) const;

    const Polygon*                	m_parentShape;

	double                          m_initedOLApexDist;

    CornerApex*                     m_innerCornerApex;
    double                          m_advConAng;


};



inline void LayerApex::markCentreLayerTrappings(const std::pair< int, double >& trap, double pc, double conAngRec, double conAngAdv,
                                      double halfAng, double intfacTen, bool injOil)
{
    if(!exists()) return;



    if(trap.first > -1 && m_trappedCL.first<0)    ///.  Becomming trapped
    {
        m_trappedCL.first = trap.first;
        m_trappedCL.second = trap.second;
        virgin = false;
    }
    else if(trap.first == -1)                         ///.  Becomming untrapped
    {
        m_trappedCL.first = -1;
    }
		//double conAng = injOil ? conAngRec: conAngAdv;
		//m_entryPc = layerCollapsePc(pc, conAng, halfAng, intfacTen, injOil);

    
}


/**
///////////////////////  LayerApex Inline Functions  ////////////////////////////////
*/
inline void LayerApex::removeLayer()
{
    m_exists = false;
    //LayerApex::set_m_stable(false);
    m_inited = false;
    m_trappedCL.first = -1;
    m_isInWatFloodVec = false;
	m_advancingPc=m_receedingPc+10000.0; ///.to affect unpinned calculations 
}



/**
// This really isn't completly correct. Since we don't know when
// the oil will become coalesced it becomes sligthly difficult.
*/
inline bool LayerApex::freeAtPrs(double pc) const { return (exists() && pc > m_entryPc) && m_trappedCL.first < 0; }/// only correct when waterInj





inline bool LayerApex::forcedSnapOff(double prs) const
{    return prs > m_receedingPc && prs > 0.0;	}




#endif
