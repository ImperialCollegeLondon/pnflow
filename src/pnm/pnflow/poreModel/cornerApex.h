#ifndef CORNERAPEX_H
#define CORNERAPEX_H

//! Water layer connectivity and interface tracking 

#include <algorithm>
#undef max
#undef min
class Polygon;

////////////////////////// CORNER FILM //////////////////////////////////////
class CornerApex : public Apex
{
public:

    CornerApex(): m_parentShape(NULL), m_initedApexDist(0.0), m_initOrMaxPcHist(-1.0e32), m_initOrMinApexDistHist(1.0e32){};
    void setCornerConnections(Apex* outerLayerApex, Polygon* parent, int subIndex)
    { m_outerLayerApex=(outerLayerApex);m_parentShape=(parent);setConnections(parent, subIndex);}

    virtual ~CornerApex() {}

    void initCornerApex(double pc, double recAng, double advAng, double halfAng, double intTen, bool oilInj);
    void finitCornerApex(double pc, double recAng, double advAng, double halfAng, double intTen, bool oilInj, bool overRideTrp);
    void createFilm(double pc, double recAng, double advAng, double halfAng, double intTen, bool isOilInj);
    inline void removeCorner();
    void getCApexDistConAng(double& apxDist, double& conA, double pc, double halfA, double ten, bool trapOveride = false, bool accurat = false, bool debug = false) const;
    inline void markTrappingCorner(const std::pair< int, double >& trpInside, double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool injOil);
	inline int index() const {return m_subIndex;};

    void updatePcsForDisconnectedOilLayer(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen);
    bool cornerExists() const {return m_exists;}

    const std::pair<int, double>& trappedCorner() const {return m_trappedCL;}
    const std::pair<int, double>& trappedCornerNOTTOBEUSED() const {return m_trappedCL;}
    bool pinnedInInitState() const {return (/*pinned () &&*/ m_initOrMinApexDistHist == m_initedApexDist);}
    inline void dump(double pc) const;
    
    double pinnedApexDist() const {return m_initedApexDist;}
    double initOrMinApexDistHist() const {return m_initOrMinApexDistHist;}


private:

    const Polygon*                	m_parentShape;

	double                          m_initedApexDist;
    double                          m_initOrMaxPcHist;
    double                          m_initOrMinApexDistHist;
    
    Apex*                     		m_outerLayerApex;

};




/**
///////////////////////  CornerApex Inline Functions  ///////////////////////////////
*/

inline void CornerApex::markTrappingCorner(const std::pair< int, double >& trpInside, double pc, double conAngRec, double conAngAdv,
                                      double halfAng, double intfacTen, bool injOil)
{
    if(!m_exists) return;

    if(trpInside.first > -1 && m_trappedCL.first<0)
    {
        m_trappedCL.first = trpInside.first;
        m_trappedCL.second = trpInside.second;
        virgin = false;
    }
    else if(trpInside.first == -1 && m_trappedCL.first>-1)
    {
        m_trappedCL.first = -1;    
        m_trapPcOld = m_trappedCL.second;
        m_trappedCL.second = 0.0;
                                              // If untrapping at a higher pressure
        //cout<<endl<<"B"<<endl;
    }
}


inline void CornerApex::removeCorner()
{
    m_exists = false;
    m_inited = false;
    m_trappedCL.first = -1;
	m_advancingPc=m_receedingPc+10000.0; ///.to affect unpinned calculations 

}



#endif
