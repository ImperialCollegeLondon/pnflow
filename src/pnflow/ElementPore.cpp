#ifdef WIN32
#pragma warning(disable:4786)
#endif

//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include <cmath>
//#include <cstdlib>
//#include <algorithm>
//#include <functional>
//#include <typeinfo>
#include <vector>
//#include <ctime>
//#include <stack>
//#include <set>
#include <cassert>
//#include <utility>
//#include <map>
using namespace std;

//#include "f2c.h"
//#include "sortedEvents.h"
//#include "threeSome.h"
#include "fluid.h"
//#include "elem_Model.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"
#include "inputData.h"

#include "Element.h"


Pore::Pore(CommonData& common, Node *node, const Oil& oil, const Water& water, double radius, double volume,
           double volumeClay, double shapeFactor, bool insideSlvBox, bool insideSatBox,
           double initSolverPrs, vector<Element*>& connThroats, std::string rockType
          ) 
   : Element(common, oil, water,
   radius, volume, volumeClay, shapeFactor, static_cast< int >(connThroats.size()), true,rockType), m_node(node)
{
    m_isInsideSatBox = insideSatBox;
    m_isInsideSolverBox = insideSlvBox;
    m_connections = connThroats;
    m_oilSolverPrs = initSolverPrs;
    m_watSolverPrs = initSolverPrs;
    m_watSolverVolt = initSolverPrs;
    m_oilSolverVolt = initSolverPrs;

    /*m_elemModel->*/setGravityCorrection(m_node);

    //double minRad(1.0e+30), maxRad(1.0e-30), radSum(0.0);
    //for(size_t i = 0; i < m_connections.size(); ++i)
    //{
        //double rad = m_connections[i]->model()->radius();
        //minRad = min(minRad, rad);
        //maxRad = max(minRad, rad);
        //radSum += rad;
    //}

    //m_averageAspectRatio = m_elemModel->radius()*m_connections.size()/radSum;
    //m_maxAspectRatio = m_elemModel->radius()/maxRad;
    //m_minAspectRatio = m_elemModel->radius()/minRad;
}


/**
// Endpores are assumed to be triangular and have zero volume and contact angles.
// Radius and interfacial tension is set to some arbitrary value to prevent any
// divison by zero that may occur. They are further assumed to be triangular so that
// they wont suddenly not contain water.
*/
InOutBoundaryPore::InOutBoundaryPore(CommonData& common, Node *node, const Oil& oil, const Water& water,
                 vector<Element*>& connThroats)
	 : Pore(common, node, oil, water,
	 1.0E-9, 0.0, 0.0, sqrt(3.0)/36.0, false, false, 0.0, connThroats,"0")///.  warning hardcode pore type, last entry
{
    //m_isInOilFloodVec = true;
    //m_isInWatFloodVec = true;
    Polygon* polyShape = dynamic_cast< Polygon* >(m_elemModel);
    if(polyShape)
    {
      for(int i = 0; i < polyShape->numCorners(); ++i)
      {
        polyShape->oilLayerCh()[i].setInWatFloodVec(true);
        polyShape->oilLayerCh()[i].setInOilFloodVec(true);
      }
        //polySh ape->calc R(0.0);
	}

    m_waterSaturation = 0.5;

    if(node->isExitRes())
        m_isExitRes = true;
    else
        m_isEntryRes = true;

    m_isConnectedToExit = m_isExitRes;
    m_isConnectedToEntry = m_isEntryRes;
    m_connectedToNetwork = false;
    //polyShape->m_m_area = 1.0
	//m_isInOilFloodVec=true;
	//m_isInWatFloodVec=true;
}


void InOutBoundaryPore::fillElemCentreWithOilRemoveLayersIO(double pc)
{
	if ( model()->bulkFluid() != &(m_comn.oil()) )
		fillElemCentreWithOilRemoveLayers();
	m_elemModel->SetInOutletPc_pistonTypeRec(pc);
	m_elemModel->SetInOutletPc_pistonTypeAdv(pc);
	//m_isInOilFloodVec=true;
	//m_isInWatFloodVec=true;
    m_connectedToNetwork = false;

}
void InOutBoundaryPore::fillElemCentreWithWaterCreateLayersIO(double pc)
{
	if ( model()->bulkFluid() != &(m_comn.water()) )
		fillElemCentreWithWaterCreateLayers();
	m_elemModel->SetInOutletPc_pistonTypeRec(pc);
	m_elemModel->SetInOutletPc_pistonTypeAdv(pc);
	//m_isInOilFloodVec=true;
	//m_isInWatFloodVec=true;
    m_connectedToNetwork = false;

}


/**
// In order to reach the outlet quickly during the trapping routine, we sort the
// connecting elements according to the distance to outlet
*/
void Pore::sortConnectingElems_DistToExit()
{
    sort(m_connections.begin(), m_connections.end(), DistToExitCompareThroats());
}











// NOT IMPORTATNTS:


/*
//void Pore::printData(ostream& out) const
//{
//// Pore data is written to output stream, containing following data:
//// x pos, y pos, z pos, radius, area, shape factor, volume, clay volume, connection number
    //m_node->printInfo(out);
    //out << *m_elemModel;
//
    //out.flags(ios::showpoint);
    //out.flags(ios::scientific);
//
    //double radAspectRatioSum(0.0);
    //for(int i = 0; i < m_connectionNum; ++i)
        //radAspectRatioSum += m_connections[i]->model()->radius();
    //double aspectRat(m_elemModel->radius()*m_connectionNum/radAspectRatioSum);
//
    //out << setprecision(4)
        //<< setw(15) << m_flowVolume*1.0E18
        //<< setw(15) << m_clayVolume*1.0E18
        //<< setw(4)  << m_isInsideSolverBox
        //<< setw(4)  << m_connectionNum
        //<< setw(4)  << m_isOnInletSlvrBdr
        //<< setw(4)  << m_isOnOutletSlvrBdr
        //<< setw(15) << aspectRat;
//}*/


/**
// These functions checks the integrity on the network, ie that pores and throats
// agree that they point to each other and that there are no NULL pointers lurking
// about. The total volume and number of connections contained within the network
// is also counted here.
*/
void Pore::calcVolume_CheckIntegrity(double& totNetVolume, double& totClayVolume, int& maxNonZeros, int& isolatedSum) const
{
    if(!m_connectedToNetwork) ++isolatedSum;
    if(m_isInsideSatBox)
    {
        totNetVolume += m_flowVolume;
        totClayVolume += m_clayVolume;
    }
    maxNonZeros += m_connectionNum + 1;
    if(index()<10)
    {
		 //outD<<"p:"<<(connection(0))->index()<<"-"<<connection(1)->index()<<":"<<index()<<"\n";
		outD<<"R_p"<<index2p()<<":"<<RRR()<<" "<<"  X_p"<<index2p()<<":"<<m_node->xPos()<<" "<<m_node->yPos()<<" "<<m_node->zPos()<<" \n";
	 }
    checkConnections();
}

void InOutBoundaryPore::calcVolume_CheckIntegrity(double& totNetVolume, double& totClayVolume, int& maxNonZeros, int& isolatedSum) const
{
    checkConnections();
}


/**
// The pore data is written to file in following format:
//
// *_node1.dat (outOne):
// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
//
// *_node2.dat (outTwo):
// index, volume, radius, shape factor, clay volume
*/
void Pore::writeNetworkData(ostream& outOne, ostream& outTwo) const
{
    size_t i;
    bool connToIn(false), connToOut(false);
    outOne.flags(ios::showpoint);
    outOne.flags(ios::scientific);
    outTwo.flags(ios::showpoint);
    outTwo.flags(ios::scientific);

    outOne << setw(7)  << orenIndex()
        << *m_node
        << setw(5) << m_connectionNum;

    for(i = 0; i < m_connections.size(); ++i)
    {
        const Throat* throat = dynamic_cast< Throat* >(m_connections[i]);
        assert(throat);
        const Element* nextPore = throat->neighbouringPore(this);
        if(nextPore->isEntryRes()) connToIn = true;
        if(nextPore->isExitRes()) connToOut = true;

        outOne << setw(7) << nextPore->orenIndex();                                 // Connecting nodes
    }

    outOne << setw(7) << connToIn << setw(7) << connToOut;                          // In and outlet?

    for(i = 0; i < m_connections.size(); ++i)
        outOne << setw(7) << m_connections[i]->orenIndex();                     // Connecting throats

    outOne << endl;
	
    outTwo << setw(7) << orenIndex()
        << setw(15) << m_flowVolume* (m_iRockType>0 ? 1.0/m_elemModel->porosity() : 1.0)
        << setw(15) << m_elemModel->radius()
        << setw(15) << m_elemModel->shapeFactor()
        << setw(15) << m_clayVolume
        << endl;
}

void Pore::writeNetworkDataBinary(ostream& out) const
{
    PoreStruct Prop;
    Prop.index = orenIndex();
    Prop.x = m_node->xPos();
    Prop.y = m_node->yPos();
    Prop.z = m_node->zPos();
    Prop.connNum = m_connectionNum;

    Prop.radius = m_elemModel->radius();
    Prop.shapeFact = m_elemModel->shapeFactor();
    Prop.volume = m_flowVolume* (m_iRockType>0 ? 1.0/m_elemModel->porosity() : 1.0);
    Prop.clayVol = m_clayVolume;
    out.write((char *)(&Prop), sizeof(Prop));
    for(size_t i = 0; i < m_connections.size(); ++i)
    {
        const Throat* throat = dynamic_cast< Throat* >(m_connections[i]);
        assert(throat);
        int idx = throat->neighbouringPore(this)->orenIndex();
        out.write((char *)(&idx), sizeof(int));
    }
    for(size_t j = 0; j < m_connections.size(); ++j)
    {
        int idx = m_connections[j]->orenIndex();
        out.write((char *)(&idx), sizeof(int));
    }

} 
 

/**

// All connecting throats for a pore should already have been set in the
// pore constructor.
*/ 
void Pore::addConnections(Element* first, Element* second, double inBdr, double outBdr, bool moveBdr)
{
    cerr << "Illegal operation" << endl; exit(-1);
}


double Pore::snapOfLongitCurvature() const
{ ///. to be implemented later
	//double length = 0.0;
	//double length = 0.0;
    //for(i = 0; i < m_connections.size(); ++i)
    //{
        //const Throat* throat = dynamic_cast< Throat* >(m_connections[i]);
        //assert(throat);
        //const Element* nextPore = throat->neighbouringPore(this);
        //if(nextPore->isEntryRes()) connToIn = true;
        //if(nextPore->isExitRes()) connToOut = true;
        //outOne << setw(7) << nextPore->orenIndex();                                 // Connecting nodes
    //}
    
		return 0.0;
}
