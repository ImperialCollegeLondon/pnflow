#ifndef APEX_H
#define APEX_H


//! Layer and pore/throat connectivity and fluid occupancy tracking


#include <algorithm>
#undef max
#undef min
class ElemModel;
#include <string>
#include <iostream>
#include <sstream>

//////////////////////// BASE CLASS /////////////////////////////
class Apex
{
public:

    Apex() :  virgin(true),m_parentModel(NULL), m_advancingPc(1000.0), m_receedingPc(-1000.0),
    m_entryPc(0.0), m_gravityCorrection(0.0), m_exists(false),m_inited(false) ,m_trappedCL(-1, 0.0) {}
    //Apex(ElemModel* parent, int subIndex) :  virgin(true), m_parentModel(parent), m_subIndex(subIndex), // m_advancingPc(1000.0), m_receedingPc(0.0),
    //m_entryPc(0.0), m_gravityCorrection(0.0), m_exists(false) ,m_trappedCL(-1, 0.0) {}
    virtual ~Apex() {}

    void setConnections(ElemModel* parent, int subIndex){m_parentModel=(parent);m_subIndex=(subIndex);}

    //virtual void updatePcsForDisconnectedOilLayer(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen) = 0;
    inline double advancingPc() const { return m_advancingPc; };
    inline double receedingPc() const { return m_receedingPc; };
    //inline bool   pinned() const { return    m _pinned; };
    const std::pair<int, double>& trappingCL() const {return m_trappedCL;}

	inline int subIndex() const {return m_subIndex;};
	ElemModel*  parentModel() const {return m_parentModel;};
    double gravCorrectedEntryPress() const {return m_entryPc+m_gravityCorrection;}
    double gravityCorrection() const {return m_gravityCorrection;}
    double entryPc() const {return m_entryPc;}

    bool isInWatFloodVec() const  {return m_isInWatFloodVec;}
    void setInWatFloodVec(bool isIt) {m_isInWatFloodVec = isIt;}
    bool isInOilFloodVec() const  {return m_isInOilFloodVec;}
    void setInOilFloodVec(bool isIt) {m_isInOilFloodVec = isIt;}

    bool exists() const {return m_exists;}
    bool pinned() const {return m_inited;}

    //double appexDist() const {return m_VisAppexDist;}///. vtuwriter only
    //double conAng() const {return m_VisConAng;}///. vtuwriter only	
	mutable bool virgin;
    double      					m_trapPcOld;
    double      					creationPc;

	static int nErrors;
	static int debugMode;


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

    ElemModel*                	m_parentModel;
	int 							m_subIndex;

    double      					m_advancingPc;
    double       				    m_receedingPc;
    double                          m_entryPc;
    double                          m_gravityCorrection;

    bool                            m_isInWatFloodVec;
    bool                            m_isInOilFloodVec;

    bool                            m_exists;
    bool                            m_inited;
    std::pair<int, double>         m_trappedCL;
    
    
};






#endif
