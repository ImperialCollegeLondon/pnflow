
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <set>
#include <map>

#include "fluid.h"
#include "Element.h"
#include "elem_Model.h"

using namespace std;

const int       ElemModel::MAX_NEWT_ITR = 1000;
const double    ElemModel::INF_NEG_NUM = -1.0E21;
const double    ElemModel::EPSILON = 1.0E-6;
const double    ElemModel::PI = acos(-1.0);

/// Base Class Constructor
ElemModel::ElemModel(Element& parent, const CommonData& common, double radius, int connNum):
	 elem_(parent), comn_(common), R_(radius) ,
	 bulkFluid_(&common.water()), 	waterConnection_ (true),  oilConnection_(false),
	 virginState_(true), porosity_(1.0)
{

	hasDisConectedCentreWCornerW_ = false;
	tetaClusterIndex_ = 0;
	ElectricalConductance_ = 1000.0; ///for IinletOutlet
	Pc__pistonTypeAdv = -20000.0*comn_.sigmaOW()*(-2.0) / R_;     ///. TOBE initialised properly later, delete to check for errors
	Pc__pistonTypeRec = -20000.0*comn_.sigmaOW()*(2.0) / R_;     ///. TOBE initialised properly later, delete to check for errors
}



