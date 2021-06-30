#ifndef SHAPE_H
#define SHAPE_H




#include "Element.h"
#include "CommonData.h"



class CommonData;
////////////////////////////// BASE CLASS FOR PORE SHAPES ////////////////////////////////
class ElemModel
{

	friend std::ostream& operator<< (std::ostream&, ElemModel&);

public:

	ElemModel(Elem&, const CommonData&, double, int);
	virtual ~ElemModel(){}

	virtual double calcR(double pc) = 0;

	virtual void initOilInjection(double pc) = 0;
	virtual void initWaterInjection(double pc) = 0;
	virtual void finitOilInjection(double pc) = 0;
	virtual void finitWaterInjection(double pc) = 0;
	virtual void fillCentreWithWaterCreateLayers(bool snapOffOverRide = false) = 0;
	virtual void fillCentreWithOilRemoveLayers() = 0;
	virtual double centreEntryPrsWatInj() = 0;
	virtual double centreEntryPrsOilInj() = 0;
	virtual void setShapeFactor(double shapeFact) = 0;
	virtual bool waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const = 0;///. rare
	virtual bool hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const = 0;
	//virtual bool check(double cappPressure) const = 0 ;

	virtual void setContactAngle(double equilCA, int wettClass, double modelTwoSepAng) = 0;
	virtual void setRadius(double newRad) = 0;
	virtual double shapeFactor() const {return 0.0625;} ///. default shp
	virtual char displacementType() const {return 'p';};


	virtual int    rockIndex() const  {return 0; };
	int    index() const  {return elem_.index(); };
	const Elem* neib(int ind) { return elem_.neib(ind); };

	//virtual bool containCFluid(const Fluid* injectant) const  { return (bulkFluid_ == injectant); }
	virtual bool containCOil() const   {return (bulkFluid_->ff() == OIL);	}


	virtual bool canNOTReconfigure(const Fluid& injectant) const  
	{  ///. only dealing with centre reconfiguration
		if (bulkFluid_==&injectant) return true;

		if (injectant.isOil())   return elem_.trappingWatBulk().first>-1 && elem_.trappingWatFilm().first>-1;
		else                     return elem_.isTrappedOil();

	}

	virtual bool conductCOil() const  { return (bulkFluid_ == &comn_.oil() && elem_.nonTrappedOil());	}
	virtual bool conductCWater() const { return (bulkFluid_ == &comn_.water() && elem_.notTrapdW(bulkBlob));	}
	bool conductAnyWaterBlob(FluidBlob blob) const
	{
		if(blob == bulkBlob) return (bulkFluid_ == &comn_.water() && elem_.notTrapdW(bulkBlob));
		else                 return waterConnection_ && elem_.notTrapdW(filmBlob);
	}

	bool conductsAny(const Fluid& fluid) const
	{
		if(fluid.isOil())  return (oilConnection_  &&  elem_.nonTrappedOil());
		else               return (waterConnection_ && (elem_.notTrapdW(bulkBlob) || elem_.notTrapdW(filmBlob) ));
	}

	bool conductsAnyOil() const  { return oilConnection_  &&  elem_.nonTrappedOil(); }
	bool conductAnyWater() const { return waterConnection_ && (elem_.notTrapdW(bulkBlob) || elem_.notTrapdW(filmBlob)); }





	double RRR() const {return R_;}
	double area() const {return area_;}///.  statistics
	double porosity() const {return porosity_;}

	bool disConectedCentreWCornerW() const {return hasDisConectedCentreWCornerW_;}
	bool affectsAdjEntryPc(const Fluid& retreadingFluid) const {return bulkFluid_ == &retreadingFluid;}
	const Fluid* bulkFluid() const {return bulkFluid_;}

	const Elem* const eleman() const {return &elem_;}   
	Elem* ChParent() const {return &elem_;}   

	void setClusterIndex(int wtindex) {tetaClusterIndex_ = wtindex;}
	inline int clusterIndex() const {return tetaClusterIndex_;}


	double Pc_pistonTypeAdv() const {return Pc_pistonTypeAdv_;};
	double Pc_pistonTypeRec() const {return Pc_pistonTypeRec_;};
	void SetInOutletPc_pistonTypeAdv(double pc)  {Pc_pistonTypeAdv_=pc;};
	void SetInOutletPc_pistonTypeRec(double pc)  {Pc_pistonTypeRec_=pc;};

	//inline double gravCorrectedEntryPress() const {	return elem_.gravCorrectedEntryPress(); }

	double electricalConductance() const	{ return ElectricalConductance_+1e-260;	}

	inline double getConductance(const Fluid& fluid, bool neighbourToInOutlet) const;
	inline double getWaterConductance(FluidBlob blob, bool neighbourToInOutlet) const;



	virtual double SPConductance(double area, double visc) const = 0;
	double oilConductance() const {return conductanceOil_;};
	double wtrConductancePar() const {return conductanceWater_.first+conductanceWater_.second;};


	const CommonData& commonData() {return comn_;}


	double  K_E_SP, K_Q_SP;

	static const int                MAX_NEWT_ITR;
	static const double             EPSILON;
	static const double             PI;
	static const double             INF_NEG_NUM;


   virtual int numCorners() const =0;

///.tmp
	virtual void readKrPcSwData(std::istringstream&, double voxelSizeScaleFact) {std::cout<<"Error: not to be called this way"<<std::endl;};

	bool exists(fluidf ff) const { return bulkFluid_->ff()==ff; };

protected:

	Elem&                       elem_;
	const CommonData&              comn_;

	double                          R_;
	const Fluid*                    bulkFluid_;
	bool                            waterConnection_;
	bool                            oilConnection_;  

	bool                            hasDisConectedCentreWCornerW_;
	bool                            virginState_;
	int                             tetaClusterIndex_; ///. contact angle cluster, to remove


	double                          porosity_;
	double                          area_;
	double                          SatWater_;

	double                          Pc_pistonTypeRec_;
	double                          Pc_pistonTypeAdv_;


	double                          ElectricalConductance_;
	double                          conductanceOil_;
	std::pair<double, double>       conductanceWater_;
	std::pair<double, double>       conductanceWaterOldToCheck_;
	friend class InOutBoundary;

};


//////////////////////////// INLINE FUNCTION DEFINITIONS //////////////////////////////





inline double ElemModel::getWaterConductance(FluidBlob blob, bool neighbourToInOutlet = false) const
{
	if(neighbourToInOutlet)
		return SPConductance(area_, comn_.water().viscosity());/// TODO: such a good job,
	else 
	if(blob == filmBlob)
		return conductanceWater_.first*0.999999999+0.000000001*conductanceWater_.second;
	else
		return conductanceWater_.second*0.999999999+0.000000001*conductanceWater_.first;
}

inline double ElemModel::getConductance(const Fluid& fluid,  bool neighbourToInOutlet = false) const
{
	double flowConductance(0.);
	if(neighbourToInOutlet)
		flowConductance = SPConductance(area_, fluid.viscosity()); /// TODO: such a good job,
	else if(fluid.isOil())
		flowConductance = conductanceOil_;
	else 
		flowConductance = std::max(conductanceWater_.first,0.) + std::max(conductanceWater_.second,0.);

	dbgAsrt(flowConductance > 1e-50);
	//if(flowConductance <=  1e-50 && debugLevel>500)
	//{
		//std::cout<<"\nError: zero flow conductance "<<"\n" ;
		//if(fluid == &comn_.oil())        std::cout<<" for oil phase "<<"\n";
		//if(neighbourToInOutlet)        std::cout<<" which is neighbourToInOutlet"<<"\n";
		//std::cout<<"  SPConductance = "<<SPConductance(area_, fluid->viscosity())<<"  bulkOil"<<conductCOil()<<std::endl; 
	//}



	return flowConductance;
}



#endif

