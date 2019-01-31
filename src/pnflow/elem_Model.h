#ifndef SHAPE_H
#define SHAPE_H


class Fluid;
class Node;
class Element;
//class CoalesceWatFillCmp;
//class CoalesceOilFillCmp;

#include "netsim_data.h"
//#include "apex.h"
#include "threeSome.h"
#include "Element.h"
 

 

////////////////////////////// BASE CLASS FOR PORE SHAPES ////////////////////////////////
class ElemModel
{

    friend std::ostream& operator<< (std::ostream&, ElemModel&);

public:

    ElemModel(Element&, const CommonData&, double, int);
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

    virtual double conAngEquil() const = 0;///. duggy
    virtual void setContactAngle(double equilConAng, int wettClass, double modelTwoSepAng) = 0;
    virtual void setRadius(double newRad) = 0;
	virtual double shapeFactor() const {return 0.0625;} ///. default shp
    virtual char displacementType() const {return 'p';};
    



	//virtual bool containCFluid(const Fluid* injectant) const  { return (m_bulkFluid == injectant); }
	virtual bool containCOil() const  {return (m_bulkFluid == &m_comn.oil());	}


    virtual bool canNOTReconfigure(const Fluid* injectant) const  
    {  ///. only dealing with centre reconfiguration
		if (m_bulkFluid == injectant) return true;

		if (injectant == &m_comn.oil()) return m_elem.trappingWatBulk().first>-1 && m_elem.trappingWatFilm().first>-1;
		else 								  return m_elem.isTrappedOil();

	 }

    virtual bool conductCOil() const  { return (m_bulkFluid == &m_comn.oil() && m_elem.nonTrappedOil());	}
    virtual bool conductCWater() const { return (m_bulkFluid == &m_comn.water() && m_elem.nonTrappedWat(bulkBlob));	}
	bool conductAnyWaterBlob(FluidBlob blob) const
	{	 //double Warning = 1.0; 	
		if(blob == bulkBlob)	return (m_bulkFluid == &m_comn.water() && m_elem.nonTrappedWat(bulkBlob));
		else					return m_waterConnection && m_elem.nonTrappedWat(filmBlob);
	}
	
     bool conductsAny(const Fluid* fluid) const
	{
		if(fluid == &m_comn.oil())			return (m_oilConnection  &&  m_elem.nonTrappedOil());
		else									return (m_waterConnection && (m_elem.nonTrappedWat(bulkBlob) || m_elem.nonTrappedWat(filmBlob) ));
	}

    bool conductsAnyOil() const
     {return (m_oilConnection  &&  m_elem.nonTrappedOil());}
    bool conductAnyWater() const
     {return (m_waterConnection && (m_elem.nonTrappedWat(bulkBlob) || m_elem.nonTrappedWat(filmBlob) ));}





    double radius() const {return m_R;}
    double area() const {return m_area;}///.  statistics
	double porosity() const {return m_porosity;}

    bool disConectedCentreWCornerW() const {return m_hasDisConectedCentreWCornerW;}
    bool affectsNeiEntryPc(const Fluid& retreadingFluid) const {return m_bulkFluid == &retreadingFluid;}
    const Fluid* bulkFluid() const {return m_bulkFluid;}

    const Element* const eleman() const {return &m_elem;}   
    Element* ChParent() const {return &m_elem;}   

    void setClusterIndex(int wtindex) {m_tetaClusterIndex = wtindex;}
    inline int clusterIndex() const {return m_tetaClusterIndex;}


    double Pc_pistonTypeAdv() const {return m_Pc_pistonTypeAdv;};
    double Pc_pistonTypeRec() const {return m_Pc_pistonTypeRec;};
    void SetInOutletPc_pistonTypeAdv(double pc)  {m_Pc_pistonTypeAdv=pc;};
    void SetInOutletPc_pistonTypeRec(double pc)  {m_Pc_pistonTypeRec=pc;};

    double rhogh() const {return m_elem.gravityCorrection();}
    inline double gravCorrectedEntryPress() const {	return m_elem.gravCorrectedEntryPress(); }

	double electricalConductance(const Fluid* fluid, bool neighbourToExit = false) const	{ return m_ElectricalConductance+1.0e-260;	}

	inline double getConductance(const Fluid* fluid, bool neighbourToInOutlet) const;
	inline double getWaterConductance(FluidBlob blob, bool neighbourToInOutlet) const;
    
    
    
    virtual double SPConductance(double area, double visc) const = 0;
    double oilConductance() const {return m_conductanceOil;};


    const CommonData& commonData() {return m_comn;}


	double  K_E_SP, K_Q_SP;

    static const int                MAX_NEWT_ITR;
    static const double             EPSILON;
    static const double             PI;
    static const double             INF_NEG_NUM;
    
    
   virtual int numCorners() const =0;

///.tmp
	virtual void readKrPcSwData(istringstream&, double voxelSizeScaleFact) {cout<<"Error: not to be called this way"<<endl;};

    
protected:

    Element&                       m_elem;
    const CommonData&              m_comn;

    double                          m_R;
    const Fluid*                    m_bulkFluid;

    bool                            m_hasDisConectedCentreWCornerW;
    bool                            m_waterConnection;
    bool                            m_oilConnection;  
    bool                            m_virginState;
    int                             m_tetaClusterIndex; ///. contact angle cluster, to remove


    double                          m_porosity;
    double                          m_area;
    double                          m_SatWater;
    
    double                          m_Pc_pistonTypeRec;
    double                          m_Pc_pistonTypeAdv;

    
	double                          m_ElectricalConductance;
	double                          m_conductanceOil;
    std::pair<double, double>       m_conductanceWater;
    std::pair<double, double>       m_conductanceWaterOldToCheck;

};


//////////////////////////// INLINE FUNCTION DEFINITIONS //////////////////////////////





inline double ElemModel::getWaterConductance(FluidBlob blob, bool neighbourToInOutlet = false) const
{
    if(neighbourToInOutlet)
        return SPConductance(m_area, m_comn.water().viscosity());/// TODO: such a good job, tabarakallah
    else 
    if(blob == filmBlob)
        return m_conductanceWater.first*0.999999999+0.000000001*m_conductanceWater.second;
    else
        return m_conductanceWater.second*0.999999999+0.000000001*m_conductanceWater.first;
}

inline double ElemModel::getConductance(const Fluid* fluid,  bool neighbourToInOutlet = false) const
{/*bool filmFlow, bool bulkFlow,*/
    double flowConductance(0.0);
    if(neighbourToInOutlet)
        flowConductance = SPConductance(m_area, fluid->viscosity()); /// TODO: such a good job, tabarakallah
    else if(fluid == &m_comn.oil())
        flowConductance = m_conductanceOil;
    else 
        flowConductance = max(m_conductanceWater.first,0.0) + max(m_conductanceWater.second,0.0);

	softAssert(m_comn.debugMode<1 || flowConductance > 1.0e-50);
	if(flowConductance <=  1.0e-50 && m_comn.debugMode>500)
	{
		cout<<"\nError: zero flow conductance "<<"\n" ;
		if(fluid == &m_comn.oil())        cout<<" for oil phase "<<"\n";
		//if(fluid == &m_comn.water())        cout<<" for water in bulk"<<bulkFlow<<" film"<<filmFlow<<" conductance = "<<flowConductance<<"\n";
		//printInfo(*m_elem);
		if(neighbourToInOutlet)        cout<<" which is neighbourToInOutlet"<<"\n";
		cout<<"  SPConductance = "<<SPConductance(m_area, fluid->viscosity())<<"  bulkOil"<<conductCOil()<<endl; 
	}



	//if (timestep>50 && !fluid->isOil() && flowConductance < watCondOld*0.9999 && 
	//m_elem.trappingWatBulk().first<0 &&  m_elem.trappingWatFilm().first<0)
	//{
		//cout
		//<<" "<<m_elem.canWBulkBePassedToSolver()<<m_elem.canWFilmBePassedToSolver()  
		//<<containCOil()<<conductCOil()<<"   "<<cCondOld<<conductCWater()<<conductAnyWater()<<"  "<<flowConductance<<" "<<watCondOld<<"\\\n"
		//<<"  "<<m_conductanceWater.first<<" "<<m_conductanceWater.second<<"\\\n";
		//
		//exit(-1); 
	//}
	//timestep++;
	//if ( !fluid->isOil())
	//{
		//watCondOld=flowConductance;
		//cCondOld=conductAnyWater();
	//}




	
	return flowConductance;
}



#endif

