#ifndef CORNERAPEX_H
#define CORNERAPEX_H

//! Water layer connectivity and interface tracking 

#include <algorithm>




class Polygon;

////////////////////////// CORNER FILM //////////////////////////////////////
class CornerApex : public Apex
{
public:

	CornerApex(): parentShape_(NULL), initedApexDist_(0.), initOrMaxPcHist_(-1e32), initOrMinApexDistHist_(1e32){};
	void setCornerConnections(Apex* outerLayerApex, Polygon* parent, int subIndex)
	{ outerLayerApex_=(outerLayerApex);parentShape_=(parent);setConnections(parent, subIndex);}

	virtual ~CornerApex() {}

	void initCornerApex(double pc, double recAng, double advAng, double halfAng, double intTen, bool oilInj);
	void finitCornerApex(double pc, double recAng, double advAng, double halfAng, double intTen, bool oilInj, bool overRideTrp);
	void createFilm(double pc, double recAng, double advAng, double halfAng, double intTen, bool isOilInj);
	inline void removeCorner();
	void getCApexDistConAng(double& apxDist, double& conA, double pc, double halfA, double ten, bool trapOveride = false, bool accurat = false, bool debug = false) const;
	inline void markTrappingCorner(const std::pair< int, double >& trpInside, double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool injOil);

	void updatePcsForDisconnectedOilLayer(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen);
	bool cornerExists() const {return exists_;}

	const std::pair<int, double>& trappedCorner() const {return trappedCL_;}
	const std::pair<int, double>& trappedCornerNOTTOBEUSED() const {return trappedCL_;}
	bool pinnedInInitState() const {return (/*pinned () &&*/ initOrMinApexDistHist_ == initedApexDist_);}
	inline void dump(double pc) const;

	double pinnedApexDist() const {return initedApexDist_;}
	double initOrMinApexDistHist() const {return initOrMinApexDistHist_;}


private:

	const Polygon*                	parentShape_;

	double                          initedApexDist_;
	double                          initOrMaxPcHist_;
	double                          initOrMinApexDistHist_;

	Apex*                     		outerLayerApex_;

};




/**
///////////////////////  CornerApex Inline Functions  ///////////////////////////////
*/

inline void CornerApex::markTrappingCorner(const std::pair< int, double >& trpInside, double pc, double conAngRec, double conAngAdv,
									  double halfAng, double intfacTen, bool injOil)
{
	if(!exists_) return;

	if(trpInside.first > -1 && trappedCL_.first<0)
	{
		trappedCL_.first = trpInside.first;
		trappedCL_.second = trpInside.second;
		virgin = false;
	}
	else if(trpInside.first == -1 && trappedCL_.first>-1)
	{
		trappedCL_.first = -1;    
		trapPcOld_ = trappedCL_.second;
		trappedCL_.second = 0.;
											  // If untrapping at a higher pressure
		//cout<<endl<<"B"<<endl;
	}
}


inline void CornerApex::removeCorner()
{
	exists_ = false;
	inited_ = false;
	trappedCL_.first = -1;
	advancingPc_=receedingPc_+10000.; ///.to affect unpinned calculations 

}



#endif
