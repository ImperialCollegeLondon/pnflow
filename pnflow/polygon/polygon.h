#ifndef POLYGON_H
#define POLYGON_H
/////////////////////////  BASE CLASS FOR POLYGON SHAPES  /////////////////////////////////
//! Pore and throat shapes
#include "../elem_Model.h"
#include "../sortedEvents.h"


class CornerApex;
class LayerApex;
class PceImbCmp;
class PceDrainCmp;

class VoidElem : public ElemModel
{
 public:
	VoidElem(Elem&, const CommonData&, double, double, int);
	virtual ~VoidElem(){};

	virtual void setContactAngle(double equilCA, int wettClass, double modelTwoSepAng);
	virtual double shapeFactor() const {return shapeFactor_;}
	virtual void setRadius(double newRad);

	double conAngleAdv() const {return cntAngAdv_;}
	double conAngleRec() const {return cntAngRec_;}
	double minInitRecCntAng() const {return minInitRecCntAng_;}
	virtual bool check(double pc) const {return true;};

protected:

	double                          shapeFactor_;


	double                          minInitRecCntAng_;
	double                          cntAngAdv_;
	double                          cntAngRec_;


};


class Polygon : public VoidElem
{
public:

	Polygon(Elem&, const CommonData&, double, double, int, int);
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

	virtual bool hasOilLayer() const {return numLayers_>0;}

	LayerApex* oilLayerCh() {return oilLayer_;}
	const LayerApex* oilLayerConst() const {return oilLayer_;}
	const CornerApex* waterInCorner() const {return waterInCorner_;}
	int numCorners() const {return numCorners_;}
	int numLayers() const {return numCorners_;}
	double  cornerHalfAngles(int i) const { return  crnHafAngs_[i];}

	void insertWatSnapEvent_IfSnapPcHgPc(Events< Apex*, PceImbCmp >& watEvents, double globalPc);
	void insertOilSnapEvent_IfSnapPcLgPc(Events< Apex*, PceDrainCmp >& oilEvents,       double globalPc);

	void calcOilLayerPc_syncTrappings(double pc);
	void calcOilLayerPc_markUntrappedFilms(double pc);

	//void setSnapOffPrs() {snapOffPrs_ = calcSnapOffPressureDrain();}

	bool Pc_growStableOilLayerDrain_UseLess(double Pc, int corner);///.  layer events
	double Pc_pin_disconnectOilLayer(int corner);  ///.  layer events

	virtual char displacementType() const {return displacementType_;};


protected:

	const int                             numCorners_;

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

	std::vector<double>     crnHafAngs_;
	CornerApex*               waterInCorner_;
	LayerApex*                oilLayer_;
	short                     numLayers_;
	char                      displacementType_;
	double                    maxConAngSpont_;


};

//////////////////////////////  SQUARE PORE SHAPES  //////////////////////////////////////
class Square : public Polygon
{
public:

	Square(Elem&, const CommonData&, double, int,int);

	void setShapeFactor(double shapeFact) {};

	double SPConductance(double area, double visc) const;
   int numCorners() const {return 4;};

private:

	double snapOffPrsDrainSpontNR() const;
	double snapOffPrsImbNR() const;

};

/////////////////////////////  TRIANGULAR PORE SHAPES  ////////////////////////////////////
class Triangle : public Polygon
{
public:

	Triangle(Elem&, const CommonData&, double, double, int,int);

	void setShapeFactor(double shapeFact);

	double SPConductance(double area, double visc) const;
   int numCorners() const {return 3;};

private:

	double snapOffPrsDrainSpontNR() const;
	double snapOffPrsImbNR() const;
	void setHalfAngles();
	double snapOffPrsDrainFromCorner(bool& casePossible, int cor) const;

};

/////////////////////////////  CIRCULAR PORE SHAPES  ////////////////////////////////////

class Circle : public VoidElem
{
  public:

	Circle(Elem&, const CommonData&, double, int, int);

	double calcR(double pc);
	void finitOilInjection(double pc);
	void finitWaterInjection(double pc);
	void initOilInjection(double pc);
	void initWaterInjection(double pc);
	void fillCentreWithWaterCreateLayers(bool snapOffOverRide = false);
	void fillCentreWithOilRemoveLayers();
	double centreEntryPrsWatInj();
	double centreEntryPrsOilInj();
	void setShapeFactor(double shapeFact) {};
	bool waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const {return false;}
	bool hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const {return false;}
	char displacementType() const {return 'P';};

	double SPConductance(double area, double visc) const;
	bool check(double pc) const {return true;};
	int numCorners() const {return 0;};


};



#endif
