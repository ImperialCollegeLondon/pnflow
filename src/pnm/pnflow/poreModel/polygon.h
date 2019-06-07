#ifndef POLYGON_H
#define POLYGON_H
/////////////////////////  BASE CLASS FOR POLYGON SHAPES  /////////////////////////////////
//! Pore and throat shapes 

class CornerApex;
class LayerApex;
//#include "apex.h"
#include "elem_Model.h"
#include "sortedEvents.h"
class PceImbCmp;
class PceDrainCmp;

class VoidElem : public ElemModel
{
 public:
    VoidElem(Element&, const CommonData&, double, double, int);
    virtual ~VoidElem(){};
    
    virtual void setContactAngle(double equilConAng, int wettClass, double modelTwoSepAng);
	virtual double shapeFactor() const {return m_shapeFactor;}
    virtual void setRadius(double newRad);

    double conAngEquil() const {return m_conAngEquil;}
    double conAngleAdv() const {return m_cntAngAdv;}
    double conAngleRec() const {return m_cntAngRec;}
    double minInitRecCntAng() const {return m_minInitRecCntAng;}
    virtual bool check(double pc) const {return true;}; 

protected:

    double                          m_shapeFactor;

     
    double                          m_minInitRecCntAng;
    double                          m_conAngEquil;
    double                          m_cntAngAdv;
    double                          m_cntAngRec;


};


class Polygon : public VoidElem
{
public:

    Polygon(Element&, const CommonData&, double, double, int, int);
    virtual ~Polygon();

    virtual double calcR(double pc);
    virtual void finitOilInjection(double pc);
    virtual void finitWaterInjection(double pc);
    virtual void initOilInjection(double pc);
    virtual void initWaterInjection(double pc);
    virtual void fillCentreWithWaterCreateLayers(bool snapOffOverRide = false);
    virtual void fillCentreWithOilRemoveLayers();
     double centreEntryPrsWatInj();
     double centreEntryPrsOilInj();
    virtual bool waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const;
    virtual bool hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const;
    //virtual void initialize();

    virtual bool hasOilLayer() const {return m_numLayers>0;}
    
    LayerApex* oilLayerCh() {return m_oilLayer;}
    const LayerApex* oilLayerConst() const {return m_oilLayer;}
    const CornerApex* waterInCorner() const {return m_waterInCorner;}
    int numCorners() const {return m_numCorners;}
    int numLayers() const {return m_numCorners;}
    double  cornerHalfAngles(int i) const { return  m_crnHafAngs[i];}

    void insertWatSnapEvent_IfSnapPcHgPc(SortedEvents< Apex*, PceImbCmp >& watEvents, double globalPc);
    void insertOilSnapEvent_IfSnapPcLgPc(SortedEvents< Apex*, PceDrainCmp >& oilEvents,       double globalPc);
        
    void calcOilLayerPc_syncTrappings(double pc);
    void calcOilLayerPc_markUntrappedFilms(double pc);
    
    //void setSnapOffPrs() {m_snapOffPrs = calcSnapOffPressureDrain();}

    bool Pc_growStableOilLayerDrain_UseLess(double Pc, int corner);///.  layer events
    double Pc_pin_disconnectOilLayer(int corner);  ///.  layer events

    virtual char displacementType() const {return m_displacementType;};


protected:

    const int                             m_numCorners;

    double calcSnapOffPressureDrain() const;
    double calcSnapOffPressureImb() const;
    virtual double snapOffPrsDrainSpontNR() const = 0;
    virtual double snapOffPrsImbNR() const = 0;

    double Pc_pistonType_Imbnww() const;
    double Pc_pistonType_ImbHingCLine() const;
    double Pc_pistonType_Drain(double conAng) const;
    double Pc_pistonType_DrainHing() const;

    void calcR_oilWithWaterInCorners(double cappPressure);
    void calcR_waterWithOilLayers(double cappPressure);

    std::vector< double >           m_crnHafAngs;
	CornerApex* 	  				m_waterInCorner;
	LayerApex* 						m_oilLayer;
    short                           m_numLayers;
    char 							m_displacementType;
    double                          m_maxConAngSpont;


};

//////////////////////////////  SQUARE PORE SHAPES  //////////////////////////////////////
class Square : public Polygon
{
public:

    Square(Element&, const CommonData&, double, int,istringstream&);
    //Square(Element&, CommonData*, const Fluid*, const Fluid*, const Square&);
    virtual ~Square(){}

    virtual void setShapeFactor(double shapeFact) {};

    virtual double SPConductance(double area, double visc) const;
   virtual int numCorners() const {return 4;};

private:

    //double snapOffPrsDrainFromCorner(bool& casePossible, int cor) const;
    virtual double snapOffPrsDrainSpontNR() const;
    virtual double snapOffPrsImbNR() const;

};

/////////////////////////////  TRIANGULAR PORE SHAPES  ////////////////////////////////////
class Triangle : public Polygon
{
public:

    Triangle(Element&, const CommonData&, double, double, int,istringstream&);
    //Triangle(Element&, CommonData*, const Fluid*, const Fluid*, const Triangle&);
    virtual ~Triangle(){}

    virtual void setShapeFactor(double shapeFact);

    virtual double SPConductance(double area, double visc) const;
   virtual int numCorners() const {return 3;};

private:

    virtual double snapOffPrsDrainSpontNR() const;
    virtual double snapOffPrsImbNR() const;
    void setHalfAngles();
    double snapOffPrsDrainFromCorner(bool& casePossible, int cor) const;

};

/////////////////////////////  CIRCULAR PORE SHAPES  ////////////////////////////////////

class Circle : public VoidElem
{
public:

    Circle(Element&, const CommonData&, double, int, istringstream&);
    virtual ~Circle(){}

    virtual double calcR(double pc);
    virtual void finitOilInjection(double pc);
    virtual void finitWaterInjection(double pc);
    virtual void initOilInjection(double pc);
    virtual void initWaterInjection(double pc);
    virtual void fillCentreWithWaterCreateLayers(bool snapOffOverRide = false);
    virtual void fillCentreWithOilRemoveLayers();
    virtual double centreEntryPrsWatInj();
    virtual double centreEntryPrsOilInj();
    virtual void setShapeFactor(double shapeFact) {};
    virtual bool waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const {return false;}
    virtual bool hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const {return false;}
    virtual char displacementType() const {return 'P';};

    virtual double SPConductance(double area, double visc) const;
    virtual bool check(double pc) const {return true;};
   virtual int numCorners() const {return 0;};


};



#endif

