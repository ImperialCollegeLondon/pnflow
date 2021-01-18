#ifndef APEX_H
#define APEX_H


//! Layer and pore/throat connectivity and fluid occupancy tracking


#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include "globals.h"

class ElemModel;
//////////////////////// BASE CLASS /////////////////////////////
class Apex
{
public:

	Apex() :  virgin(true),parentModel_(NULL), advancingPc_(1000.0), receedingPc_(-1000.0),
	entryPc_(0.0), gravityCorrection_(0.0), exists_(false),inited_(false) ,trappedCL_(-1, 0.0) {}
	virtual ~Apex() {}

	void setConnections(ElemModel* parent, int subIndex){parentModel_=(parent);subIndex_=(subIndex);}

	inline double advancingPc() const { return advancingPc_; };
	inline double receedingPc() const { return receedingPc_; };
	const std::pair<int, double>& trappingCL() const {return trappedCL_;}

	inline int subIndex() const {return subIndex_;};
	ElemModel*  parentModel() const {return parentModel_;};
	double gravCorrectedEntryPress() const {return entryPc_+gravityCorrection_;}
	double gravityCorrection() const {return gravityCorrection_;}
	double entryPc() const {return entryPc_;}

	bool isInWatFloodVec() const  {return isInWatFloodVec_;}
	void setInWatFloodVec(bool isIt) {isInWatFloodVec_ = isIt;}
	bool isInOilFloodVec() const  {return isInOilFloodVec_;}
	void setInOilFloodVec(bool isIt) {isInOilFloodVec_ = isIt;}

	bool exists() const {return exists_;}
	bool pinned() const {return inited_;}

	mutable bool virgin;
	double      					trapPcOld_;
	double      					creationPc;

	static int nErrors;


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

	ElemModel*                	parentModel_;
	int 							subIndex_;

	double      					advancingPc_;
	double       					receedingPc_;
	double                          entryPc_;
	double                          gravityCorrection_;

	bool                            isInWatFloodVec_;
	bool                            isInOilFloodVec_;

	bool                            exists_;
	bool                            inited_;
	std::pair<int, double>         trappedCL_;


};






#endif
