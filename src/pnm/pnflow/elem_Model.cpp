
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
const double    ElemModel::PI = acos(-1.);

/// Base Class Constructor
ElemModel::ElemModel(Elem& parent, const CommonData& common, double radius, int connNum):
	 elem_(parent), comn_(common), R_(radius) ,
	 bulkFluid_(&common.water()), 	waterConnection_ (true),  oilConnection_(false),
	 virginState_(true), porosity_(1.)  {

	hasDisConectedCentreWCornerW_ = false;
	tetaClusterIndex_ = 0;
	ElectricalConductance_ = 1000.; ///for IinletOutlet
	Pc_pistonTypeAdv_ = -20000.*comn_.sigmaOW()*(-2.) / R_;     ///. TOBE initialised properly later, delete to check for errors
	Pc_pistonTypeRec_ = -20000.*comn_.sigmaOW()*(2.) / R_;     ///. TOBE initialised properly later, delete to check for errors
}



