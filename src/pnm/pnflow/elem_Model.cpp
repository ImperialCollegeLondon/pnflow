#ifdef WIN32
#pragma warning(disable:4786)
#endif

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
using namespace std;

#include "fluid.h"
#include "Element.h"
#include "elem_Model.h"

const int       ElemModel::MAX_NEWT_ITR = 1000;
const double    ElemModel::INF_NEG_NUM = -1.0E21;
const double    ElemModel::EPSILON = 1.0E-6;
const double    ElemModel::PI = acos(-1.0);


/**
// Base Class Constructors
*/
ElemModel::ElemModel(Element& parent, const CommonData& common, double radius, int connNum)://, m_oil(oil),m_water(water),
	         m_elem(parent), m_comn(common), m_R(radius), m_bulkFluid(&common.water()),
             m_virginState(true),
             m_porosity(1.0) 
             //m_entryPc(0.0),
             //m_gravityCorrection(0.0)			
{
    //m_KrwatcornAtSw0 = m_comn.KrwatcornAtSw0();
    m_waterConnection = true;
    m_oilConnection = false;

    m_hasDisConectedCentreWCornerW = false;
    m_tetaClusterIndex = 0;
	m_ElectricalConductance = 1000.0; ///for IinletOutlet
	//m_Rc_pistonTypeRec = 1.0*m_R; /// set to highest possible value
    //m_Rc_pistonTypeAdv = -1.0*m_R; /// set to lowest possible value
	m_Pc_pistonTypeAdv = -20000.0*m_comn.oil().interfacialTen()*(-2.0) / m_R;     ///. TOBE initialised properly later, delete to check for errors
	m_Pc_pistonTypeRec = -20000.0*m_comn.oil().interfacialTen()*(2.0) / m_R;     ///. TOBE initialised properly later, delete to check for errors
}



//#include "elem_hetroPorous.cpp"
//#include "elem_porous.cpp"
//#include "polygon.cpp"
//#include "polygonDrain.cpp"
//#include "polygonImb.cpp"

