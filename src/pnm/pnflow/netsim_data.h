#ifndef COMMONDATA_H
#define COMMONDATA_H

class Element;
class InputData;
#undef max 
#undef min
#include "inputFile.h"
#include "fluid.h"
#include <algorithm>
#include "inputData.h"

enum TrappingCriteria {escapeToInlet = 0, escapeToOutlet, escapeToEither, escapeToBoth};
enum FluidBlob {filmBlob = 0, bulkBlob};






class CommonData
{
public:

    CommonData(const InputData& input)		   :
			debugMode(0),
			m_maxEverCappPress(0.0),
			m_minEverCappPress(0.0),		   
		   m_trappedRegionsOil(m_trappedRegionsOil_),
		   m_trappedRegionsWat(m_trappedRegionsWat_)

	{
		cout<<"creating flowData"<<endl;

		m_circWatCondMultFact = input.getOr(0.0,"SURFACE_FILM_COND_FACT");
		input.gravityConst(m_gravConstX, m_gravConstY, m_gravConstZ);
		input.poreFillWgt(m_poreFillWeights);
		input.poreFillAlg(m_poreFillAlg);
		input.getVar(debugMode,"debugLevel");
		m_floodingCycle = 0;
		//m_imbibitionCycle = 1;
		m_maxPcLastDrainCycle = 0.0;
		m_minPcLastImbCycle = 0.0;
		m_numSquares = 0;
		m_numTriangles = 0;
		m_numCircles = 0;
		m_injectant = 0;
		m_guessPc = 0;
		m_drainagePhase = false;
		
		double interfacTen, watDens(1000.0), oilDens(1000.0) , watVisc, oilVisc, watResist, oilResist;
		input.fluid(interfacTen, watVisc, oilVisc, watResist, oilResist, watDens, oilDens);
		double clayResistivity(watResist), surfaceResistivity(1.0e12);
		input.getVar(clayResistivity, "clayResistivity");
		input.getVar(surfaceResistivity, "surfaceResistivity");

		m_solveSP=false; input.getVar(m_solveSP, "solveSP");
		m_SwPowerStreamCond = 1.0;
		double epsilonW(80.1*8.854e-12), zetaW(0.061), zetaW_OW(0.061), epsilonO(80.1*8.854e-12*0.0001), zetaOW(-0.061), zetaO(0.0);
		
		if (m_solveSP)
		{
			if (!input.getVar(epsilonW, "epsilonW")) { cout<< "Error missing keyword epsilonW";    }
			if (!input.getVar(zetaW, "zetaW")) { cout<< "Error missing keyword zetaW";   }
			if (!input.getVar(zetaW_OW, "zetaW_OW")) { cout<< "Error missing keyword zetaW_OW";   }
			if (!input.getVar(epsilonO, "epsilonO")) { cout<< "Error missing keyword epsilonO";   }
			if (!input.getVar(zetaOW, "zetaOW")) { cout<< "Error missing keyword zetaOW";   }
			if (!input.getVar(m_SwPowerStreamCond, "SwPowerStreamCond")) { cout<< "Error missing keyword SwPowerStreamCond";   }

			cout<< endl << endl;
			cout<< " epsilonW " << epsilonW << endl;
			cout<< " zetaW    " << zetaW << endl;
			cout<< " epsilonO " << epsilonO << endl;
			cout<< " zetaW_OW    " << zetaW_OW << endl;
			cout<< " zetaOW    " << zetaOW << endl;
			cout<< " watResist    " << watResist << endl;
			cout<< " oilResist    " << oilResist << endl;
			cout<< " clayResistivity    " << clayResistivity << endl;
			cout<< " surfaceResistivity    " << surfaceResistivity << endl;
			cout<< " SwPowerStreamCond    " << m_SwPowerStreamCond << endl;
			cout<< endl << endl;
		}

		m_water.setFluidProps(watVisc, interfacTen, watResist, watDens, epsilonW, zetaW);
		m_clay.setFluidProps(watVisc*1.0e12, 0, clayResistivity, watDens, epsilonW, zetaW);
		m_rockWaterSurface.setFluidProps(watVisc*1.0e12, 0, surfaceResistivity, watDens, epsilonW, zetaW);
		m_rockWaterSurface_OW.setFluidProps(watVisc*1.0e12, 0, surfaceResistivity, watDens, epsilonW, zetaW_OW);
		m_oilWaterInterface.setFluidProps(watVisc*1.0e12, 0, surfaceResistivity, watDens, epsilonW, zetaOW); ///. only zetaOW is used
		m_oil.setFluidProps(oilVisc, interfacTen, oilResist, oilDens, epsilonO, zetaO); ///. zetaO is not used

	}


    ~CommonData(){}

    void finalizeCopyConstruct(const CommonData& data, const std::vector<Element*>& elems);

    double poreFillWeights(int i) const {return m_poreFillWeights[i];}
    const std::string& poreFillAlg() const {return m_poreFillAlg;}



    int floodingCycle() const {return m_floodingCycle;}
    //int imbibitionCycle() const {return m_imbibitionCycle;}
    void incrementFloodCycle() {++m_floodingCycle;}
    //void incrementImbCycle() {++m_imbibitionCycle;}



    void setDrainageCycle(bool isIt) {m_drainagePhase = isIt;}
    bool isDrainageCycle() const {return m_drainagePhase;}

    double maxEverPc() const {return m_maxEverCappPress;}
    double maxPcLastDrainCycle() const {return m_maxPcLastDrainCycle;}
    inline void setMaxPcLastDrainCycle(double pc);
    double minEverCappPress() const {return m_minEverCappPress;}
    double minPcLastImbCycle() const {return m_minPcLastImbCycle;}
    inline void setMinPcLastCycle(double pc);
    void GuessCappPress(double pc) {m_guessPc = pc;}
    double GuessCappPress() const {return m_guessPc;}
    double gravConstX() const {return m_gravConstX;}
    double gravConstY() const {return m_gravConstY;}
    double gravConstZ() const {return m_gravConstZ;}

    double KrwatcornAtSw0() const {return m_circWatCondMultFact;}

    void injectant(const Fluid* injFluid) {m_injectant = injFluid;}
    const Fluid* injectant() const {return m_injectant;}
    const Oil& oil() const {return m_oil;}
	const Water& water() const { return m_water; }
	const Clay& clay() const { return m_clay; }
	const Surface& surface() const { return m_rockWaterSurface; }
	const Surface& surfaceOW() const { return m_rockWaterSurface_OW; }
	const Surface& oilWaterInterface() const { return m_oilWaterInterface; }
	double SwPowerStreamCond() const { return m_SwPowerStreamCond; }

    size_t newOilTrappingIndex() const {return m_trappedRegionsOil.size();}
    size_t newWatTrappingIndex() const {return m_trappedRegionsWat.size();}

    inline void removeTrappedOilElem(int idx, Element* elem);
    inline void removeTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem);
           void addTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem) {    m_trappedRegionsWat_[idx].push_back(elem); }

    void addTrappedRegionOil(const std::vector<Element*>& elems) {m_trappedRegionsOil_.push_back(elems);}
    void addTrappedRegionWat(const std::vector< std::pair<Element*, FluidBlob> >& elems) {m_trappedRegionsWat_.push_back(elems);}
    void removeTrappedRegionOil(int region) {m_trappedRegionsOil_[region].clear();}
    void removeTrappedRegionWat(int region) {m_trappedRegionsWat_[region].clear();}

    const std::vector<Element*>& trappedRegionsOil(int region) const {return m_trappedRegionsOil[region];}
    const std::vector< std::pair<Element*,FluidBlob> >& trappedRegionsWat(int region) const {return m_trappedRegionsWat[region];}

	/// for writting statistics
    const std::vector< std::vector<Element*> >	&				        trappedOilRegions() const {return m_trappedRegionsOil;}
    const std::vector< std::vector< std::pair<Element*,FluidBlob> > > &    trappedWatRegions() const {return m_trappedRegionsWat;}



    int numPores() const {return m_numPores;}
    int numThroats() const {return m_numThroats;}
    void setNumElem(int numP, int numT) {m_numPores = numP; m_numThroats = numT;}

    inline size_t numTrappedOilRegions() const;
    inline size_t numTrappedOilElems() const;
    inline size_t numTrappedWatRegions() const;
    inline size_t numTrappedWatElems() const;    
    
    void countCircle() {++m_numCircles;}
    void countTriangle() {++m_numTriangles;}
    void countSquare() {++m_numSquares;}
    void countPorous() {++m_numPorous;}
    int numCircles() const {return m_numCircles;}
    int numTriangles() const {return m_numTriangles;}
    int numSquares() const {return m_numSquares;}
    int numPorous() const {return m_numPorous;}    
    
    std::ofstream                        dbgOut;
    int                        debugMode;
    bool                    		m_solveSP; ///. TODO make const
    
private:

    void readDBGangles(const std::string& fileName);

	
    Oil                    		m_oil; ///. TODO make const
	Water                    	m_water;
	Clay                    	m_clay;
	Surface                    	m_rockWaterSurface;
	Surface                    	m_rockWaterSurface_OW;
	Surface                    	m_oilWaterInterface;

	double                      m_SwPowerStreamCond;
	const Fluid*                                        m_injectant;

    double										        m_maxEverCappPress;
    double                                              m_maxPcLastDrainCycle;
	double										        m_minEverCappPress;
    double                                              m_minPcLastImbCycle;
    double                                              m_guessPc;
    double										        m_circWatCondMultFact;
    double                                        m_gravConstX;
    double                                        m_gravConstY;
    double                                        m_gravConstZ;
    int									                m_numPores;
    int									                m_numThroats;
    bool                                                m_drainagePhase;

    int                                                 m_floodingCycle;
    //int                                                 m_imbibitionCycle;
    std::vector< double >							        m_poreFillWeights;
    std::string										        m_poreFillAlg;
    std::vector< std::vector<Element*> >					        m_trappedRegionsOil_;
    std::vector< std::vector< std::pair<Element*,FluidBlob> > >	    m_trappedRegionsWat_;
    const std::vector< std::vector<Element*> >	&				        m_trappedRegionsOil;
    const std::vector< std::vector< std::pair<Element*,FluidBlob> > >	&    m_trappedRegionsWat;

    ///.  not used data
    int											        m_numSquares;
    int											        m_numTriangles;
    int											        m_numCircles;
    int											        m_numPorous;
 };

inline size_t CommonData::numTrappedOilElems() const
{
    size_t numTrapped(0);
    for(size_t i = 0; i < m_trappedRegionsOil.size(); ++i)
        numTrapped += m_trappedRegionsOil[i].size();

    return numTrapped;
}

inline size_t CommonData::numTrappedWatElems() const
{
    size_t numTrapped(0);
    for(size_t i = 0; i < m_trappedRegionsWat.size(); ++i)
        numTrapped += m_trappedRegionsWat[i].size();

    return numTrapped;
}

inline size_t CommonData::numTrappedOilRegions() const
{
    size_t numTrpRegions(0);
    for(size_t i = 0; i < m_trappedRegionsOil.size(); ++i)
    {
        if(!m_trappedRegionsOil[i].empty())
            ++numTrpRegions;
    }
    return numTrpRegions;
}

inline size_t CommonData::numTrappedWatRegions() const
{
    size_t numTrpRegions(0);
    for(size_t i = 0; i < m_trappedRegionsWat.size(); ++i)
    {
        if(!m_trappedRegionsWat[i].empty())
            ++numTrpRegions;
    }
    return numTrpRegions;
}

inline void CommonData::setMaxPcLastDrainCycle(double pc)
{
    m_maxEverCappPress = std::max(pc, m_maxEverCappPress);
    m_maxPcLastDrainCycle = pc;
}

inline void CommonData::setMinPcLastCycle(double pc)
{
    m_minEverCappPress = std::min ( pc, m_minEverCappPress);
    m_minPcLastImbCycle = pc;
}

///. deletes elem from trappedRegionsOil //called from checkUntrapOilIfUnstableConfigsImb
inline void CommonData::removeTrappedOilElem(int idx, Element* elem)
{
    m_trappedRegionsOil_[idx].erase(
        std::remove(m_trappedRegionsOil_[idx].begin(),m_trappedRegionsOil_[idx].end(),elem),
        m_trappedRegionsOil_[idx].end());
}

///. deletes elem blob from trappedRegionsWat //called from checkUntrapWaterIfUnstableConfigsDrain
inline void CommonData::removeTrappeWatElem(int idx, std::pair<Element*, FluidBlob> elem)
{
    m_trappedRegionsWat_[idx].erase(
        std::remove(m_trappedRegionsWat_[idx].begin(),m_trappedRegionsWat_[idx].end(),elem),
        m_trappedRegionsWat_[idx].end());
}




#endif

