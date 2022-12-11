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

	LayerApex() : initedOLApexDist_(-1.) {};

	//LayerApex(CornerApex* innerCornerApex, Polygon* parent, int subIndex)
		//: Apex(parent, subIndex),  parentShape_(parent), /*lastStablePc_(0.),*/, innerCornerApex_(innerCornerApex)


	void setLayerConnections(CornerApex* innerCornerApex, Polygon* parent, int subIndex)
	{ innerCornerApex_=(innerCornerApex);parentShape_=(parent);setConnections(parent, subIndex);}

	virtual ~LayerApex() {}

	bool createOLayer(double pc, double recAng, double advAng, double maxSpontConAng, double halfAng, double intTen, bool isOilInj);

	bool initLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intTen, bool oilInj, bool silent=false);
	bool finitLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intTen, bool oilInj, bool overRideTrp);

	void getCAApexDistUntraped(double& apxDist, double& conA, const double& halfA, double pc, double ten, bool itr = false) const;
	void getCAApexDist(double& apxDist, double& conA, const double& halfA, double pc, double ten, bool itr = false) const;

	inline void markCentreLayerTrappings(const std::pair<int,double>&, double, double, double, double, double, bool);
	double layerCollPc() const {return entryPc_;}
	inline void removeLayer();
	const std::pair<int, double>& trappedOLayer() const {return trappedCL_;}
	inline void set_m_exists(bool isIt){ exists_ = isIt; }
	inline bool freeAtPrs(double pc) const;
	inline bool forcedSnapOff(double prs) const;
	//inline int index() const {return subIndex_;};

	void advConAng(double conAng) {advConAng_ = conAng;}
	inline bool stablePinnedInLastCycle(double minPcLastCycle) const {return (exists(/*st ab le*/) && /*pinned () &&*/ advancingPc_ == minPcLastCycle);}

	double pinnedApexDist() const {return initedOLApexDist_;}
	double layerCollapsePc(double pc, double conAng, double halfAng, double intfacTen, bool injOil) const;
	double layerCollapsePc_fromEitherSide(double pc, double conAng, double halfAng, double intfacTen,bool debug=false) const;

	mutable int colType_;
private:

	double layerCollapsePc_FromCentre(double outPc, double inPc, double conAng, double halfAng, double ten) const;
	double layerCollapsePc_FromCorner(double outPc, double inPc, double conAng, double halfAng, double ten) const;

	const Polygon*                	parentShape_;

	double                          initedOLApexDist_;

	CornerApex*                     innerCornerApex_;
	double                          advConAng_;


};



inline void LayerApex::markCentreLayerTrappings(const std::pair<int,double>& trap, double pc, double conAngRec, double conAngAdv,
									  double halfAng, double intfacTen, bool injOil)
{
	if(!exists()) return;



	if(trap.first > -1 && trappedCL_.first<0)    ///.  Becomming trapped
	{
		trappedCL_.first = trap.first;
		trappedCL_.second = trap.second;
		virgin = false;
	}
	else if(trap.first == -1)                         ///.  Becomming untrapped
	{
		trappedCL_.first = -1;
	}
		//double conAng = injOil ? conAngRec: conAngAdv;
		//entryPc_ = layerCollapsePc(pc, conAng, halfAng, intfacTen, injOil);


}


/**
///////////////////////  LayerApex Inline Functions  ////////////////////////////////
*/
inline void LayerApex::removeLayer()
{
	exists_ = false;
	//LayerApex::set_m_stable(false);
	inited_ = false;
	trappedCL_.first = -1;
	isInWatFloodVec_ = false;
	advancingPc_=receedingPc_+10000.; ///.to affect unpinned calculations
}



/**
// This really isn't completly correct. Since we don't know when
// the oil will become coalesced it becomes sligthly difficult.
*/
inline bool LayerApex::freeAtPrs(double pc) const { return (exists() && pc > entryPc_) && trappedCL_.first < 0; }/// only correct when waterInj





inline bool LayerApex::forcedSnapOff(double prs) const
{    return prs > receedingPc_ && prs > 0.;	}




#endif
