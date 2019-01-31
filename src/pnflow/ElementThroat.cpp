#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include <cassert>

using namespace std;


#include "fluid.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"
#include "inputData.h"

#include "Element.h"



Throat::Throat(CommonData& common, const Oil& oil, const Water& water, double radius, double volume, double volumeClay,
               double shapeFactor, double length, double lengthPore1, double lengthPore2, int index,
               std::string rockType) : Element(common, oil,water, radius, volume, volumeClay, shapeFactor, 2, false, rockType), 
              m_latticeIndex(index)
{
    m_length = length;
    m_poreLength.push_back(lengthPore1);
    m_poreLength.push_back(lengthPore2);
}










void Throat::calcVolume_CheckIntegrity(double& totNetVolume, double& totClayVolume, int& maxNonZeros, int& isolatedSum) const
{
    if(!m_connectedToNetwork) ++isolatedSum;
    if(m_isInsideSatBox)
    {
        totNetVolume += m_flowVolume;
        totClayVolume += m_clayVolume;
    }
    if(index2p()<10) outD<<"g:"<<connection(0)->index2p()<<"-"<<connection(1)->index2p()<<":"<<index2p()<<" R:"<<RRR()<<" nCor_"<<index2p()<<":"<<model()->numCorners()<<"\n";
    checkConnections();

}



double Throat::snapOfLongitCurvature() const
{
	//double lengthEff = (m_poreLength[1]+m_poreLength[0]+m_length);
	double delRSqr = 0.5*(m_connections[1]->model()->radius()+m_connections[0]->model()->radius())-model()->radius();
	//cout<<lengthEff<<"  "<<m_connections[1]->model()->radius()<<"  "<<m_connections[0]->model()->radius()<<"  "<<model()->radius()<<"  "<<endl;
	if (delRSqr<0.0) return 0.0; ///. Errror
	delRSqr *= delRSqr;
	
	return 0.0; 
	//return 4.0*std::sqrt(delRSqr/(delRSqr+lengthEff*lengthEff*0.025330296))/lengthEff; //lengthEff/lengthEff;
	///. return 2*sigma*sin(45)*sqrt(2.0*delRSqr/(delRSqr+(L_poreToPore/2PI)^2)/(lengthEff/2.0) //lengthEff/lengthEff;

}


double Throat::poreLength(const Element* callingPore) const
{
    if(m_connections.size() != 2)
    {
        cerr << endl
            << "============================================" << endl
            << "For optimized network to be written to file " << endl
            << "the option to drain singlets must be enabled" << endl
            << "============================================" << endl;        exit(-1);
    }
    if(m_connections[1] == callingPore)
        return m_poreLength[1];
    else
        return m_poreLength[0];
}


const Element* Throat::neighbouringPore(const Element* callingPore) const
{
    if(m_connections.size() != 2)
    {
        cerr << endl
            << "============================================" << endl
            << "For optimized network to be written to file " << endl
            << "the option to drain singlets must be enabled" << endl
            << "============================================" << endl;        exit(-1);
    }
    if(m_connections[0] == callingPore)
        return m_connections[1];
    else
        return m_connections[0];
}





/**
// Does a throat cross any given interior plane
*/
bool Throat::crossesPlaneAt(double location) const
{
    double ptOne(m_connections[0]->node()->xPos()), ptTwo(m_connections[1]->node()->xPos());

    if(ptTwo < ptOne)
    {
        double tmp = ptTwo;         // Ensure that that the two points we pick up acyually are within
        ptTwo = ptOne;              // the inteior of the model
        ptOne = tmp;
    }

    return ptOne < location && ptTwo >= location;
}


/**
// When sorting the connecting pores in the throats we have to make sure that we also sort
// the associated pore lengths ;-).
*/
void Throat::sortConnectingElems_DistToExit()
{
    Element *oldFirst = m_connections[0];
    sort(m_connections.begin(), m_connections.end(), DistToExitComparePores());

    assert(m_connections.size() == 2 && m_poreLength.size() == 2);
    if(oldFirst != m_connections[0])    // Pores got switched
    {
        double tmp(m_poreLength[0]);
        m_poreLength[0] = m_poreLength[1];
        m_poreLength[1] = tmp;

        if(m_originalPoreLengths)
        {
            double tmpLen(m_originalPoreLengths[0]);
            m_originalPoreLengths[0] = m_originalPoreLengths[2];
            m_originalPoreLengths[2] = tmpLen;
        }
    }
}


/**
// The throat data is written to file in following format:
//
// *_link1.dat (outOne):
// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
//
// *_link2.dat (outTwo):
// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
*/
void Throat::writeNetworkData(ostream& outOne, ostream& outTwo) const
{
    outOne.flags(ios::showpoint);
    outOne.flags(ios::scientific);
    outTwo.flags(ios::showpoint);
    outTwo.flags(ios::scientific);
    double lenPoreOne(m_poreLength[0]), lenPoreTwo(m_poreLength[1]), lenThroat(m_length);
    double lenTotal(m_poreLength[0]+m_poreLength[1]+m_length);

    if(m_originalPoreLengths)
    {
        lenPoreOne = m_originalPoreLengths[0];  // The pore lengths were modified when moving
        lenThroat = m_originalPoreLengths[1];   // pressure boundaries
        lenPoreTwo = m_originalPoreLengths[2];
        lenTotal = lenPoreOne+lenThroat+lenPoreTwo;
    }

    if(m_connections.size() != 2)
    {
        cerr << endl
            << "============================================" << endl
            << "For optimized network to be written to file " << endl
            << "the option to drain singlets must be enabled" << endl
            << "============================================" << endl;        exit(-1);
    }

    outOne << setw(7)   << orenIndex()
        << setw(7)      << m_connections[0]->orenIndex()
        << setw(7)      << m_connections[1]->orenIndex()
        << setw(15)     << m_elemModel->radius()
        << setw(15)     << m_elemModel->shapeFactor()
        << setw(15)     << lenTotal
        << endl;

    outTwo << setw(7)   << orenIndex()
        << setw(7)      << m_connections[0]->orenIndex()
        << setw(7)      << m_connections[1]->orenIndex()
        << setw(15)     << lenPoreOne
        << setw(15)     << lenPoreTwo
        << setw(15)     << lenThroat
        << setw(15)     << m_flowVolume* (m_iRockType>0 ? 1.0/m_elemModel->porosity() : 1.0)
        << setw(15)     << m_clayVolume
        << endl;
}

void Throat::writeNetworkDataBinary(ostream& out) const
{
    double lenPoreOne(m_poreLength[0]), lenPoreTwo(m_poreLength[1]), lenThroat(m_length);
    double lenTotal(m_poreLength[0]+m_poreLength[1]+m_length);
    if(m_originalPoreLengths)
    {
        lenPoreOne = m_originalPoreLengths[0];  // The pore lengths were modified when moving
        lenThroat = m_originalPoreLengths[1];   // pressure boundaries
        lenPoreTwo = m_originalPoreLengths[2];
        lenTotal = lenPoreOne+lenThroat+lenPoreTwo;
    }

    ThroatStruct Prop;
    Prop.index = orenIndex();
    Prop.poreOne = m_connections[0]->orenIndex();
    Prop.poreTwo = m_connections[1]->orenIndex();

    Prop.radius = m_elemModel->radius();
    Prop.shapeFact = m_elemModel->shapeFactor();
    Prop.lenPoreOne = lenPoreOne;
    Prop.lenPoreTwo = lenPoreTwo;
    Prop.lenThroat = lenThroat;
    Prop.lenTot = lenTotal;
    Prop.volume = m_flowVolume* (m_iRockType>0 ? 1.0/m_elemModel->porosity() : 1.0);
    Prop.clayVol = m_clayVolume;
    out.write((char *)(&Prop), sizeof(Prop));

} 
 
 
 


/**

// The connecting pores are added to the throat. From the pores it is also
// possible to determine if the throat is inside the calculation box. If part of
// the connection (pore-throat-pore) is outside the calculation box that lenght
// is set to zero. This has the effect when solving the pressure field, no
// pressure loss occurs outside the box. However if the pore inside the box itself
// is on the boundary, pressure losses do occur within that pore. This constraint
// is needed in the case we're using the whole network for rel perm calculations.
// The first pores are usually located at x position 0.0. These are still within
// the box. There has to be some length to the outlet to make the problem
// solveable, hence we allow pressure drop to occur within that pore.
*/ 
void Throat::addConnections(Element* first, Element* second, double inletBdr, double outletBdr, bool moveBoundary)
{
	m_node.m_xPos=0.5*(second->node()->m_xPos + first->node()->m_xPos);
	m_node.m_yPos=0.5*(second->node()->m_yPos + first->node()->m_yPos);
	m_node.m_zPos=0.5*(second->node()->m_zPos + first->node()->m_zPos);

	
	
    if(first->isEntryRes() || first->isExitRes())
        /*m_elemModel->*/setGravityCorrection(second->node());
    else if(second->isEntryRes() || second->isExitRes())
        /*m_elemModel->*/setGravityCorrection(first->node());
    else
        /*m_elemModel->*/setGravityCorrection(first->node(), second->node());

    if(first->isEntryRes() || second->isEntryRes())
    {
        m_isConnectedToEntry = true;
        first->isConnectedToEntry(first->isEntryRes());
        second->isConnectedToEntry(second->isEntryRes());
    }
    else if(first->isExitRes() || second->isExitRes())
    {
        m_isConnectedToExit = true;
        first->isConnectedToExit(first->isExitRes());
        second->isConnectedToExit(second->isExitRes());
    }

    m_isInsideSolverBox = (first->isInsideSolverBox() || second->isInsideSolverBox());
    m_isInsideSatBox = (first->isInsideSatBox() || second->isInsideSatBox());
    m_connectedToEntryOrExit = (first->isEntryOrExitRes() || second->isEntryOrExitRes());

    double oldPOneLen(m_poreLength[0]), oldPTwoLen(m_poreLength[1]), oldThrLen(m_length);

    if(m_isInsideSolverBox && !second->isInsideSolverBox())
    {
        second->setOnInletSlvrBdr(second->node()->xPos() < inletBdr);
        second->setOnOutletSlvrBdr(second->node()->xPos() > outletBdr);

        if(moveBoundary)
        {
            double scaleFact = (second->node()->xPos() - first->node()->xPos()) / (m_length+m_poreLength[0]+m_poreLength[1]);
            if(second->isEntryOrExitRes())      // We don't know position of exit and entry res.
                scaleFact = (second->node()->xPos()-first->node()->xPos()) / fabs(second->node()->xPos()-first->node()->xPos());;

            double bdr = (second->node()->xPos() < inletBdr) ? inletBdr: outletBdr;
            double throatStart = first->node()->xPos() + m_poreLength[0]*scaleFact;
            double throatEnd = throatStart + m_length*scaleFact;
            m_originalPoreLengths[0] = m_poreLength[0];
            m_originalPoreLengths[1] = m_length;
            m_originalPoreLengths[2] = m_poreLength[1];


            if(second->isEntryOrExitRes())                                          // Keep throat lengths if whole model is being used
                m_poreLength[1] = 0.0;                                              // for calculations
            else if(throatEnd > inletBdr && throatEnd < outletBdr)                  // Both pore1 and throat are within the box
                m_poreLength[1] *= (bdr - throatEnd)/(m_poreLength[1]*scaleFact);
            else if(throatStart > inletBdr && throatStart < outletBdr)              // Onle pore 1 is fully within box
            {
                m_poreLength[1] = 0.0;
                m_length *= (bdr - throatStart)/(m_length*scaleFact);
            }
            else                                                                    // Pore 1 is only partially within box
            {
                m_poreLength[1] = 0.0;
                m_length = 0.0;
            }
        }
    }
    else if(m_isInsideSolverBox && !first->isInsideSolverBox())             // Pore 1 is outside box
    {
        first->setOnInletSlvrBdr(first->node()->xPos() < inletBdr);
        first->setOnOutletSlvrBdr(first->node()->xPos() > outletBdr);

        if(moveBoundary)
        {
            double scaleFact = (first->node()->xPos() - second->node()->xPos()) / (m_length+m_poreLength[0]+m_poreLength[1]);
            if(first->isEntryOrExitRes())       // We don't know position of exit and entry res.
                scaleFact = (first->node()->xPos()-second->node()->xPos()) / fabs(first->node()->xPos()-second->node()->xPos());

            double bdr = (first->node()->xPos() < inletBdr) ? inletBdr: outletBdr;
            double throatStart = second->node()->xPos() + m_poreLength[1]*scaleFact;
            double throatEnd = throatStart + m_length*scaleFact;
            m_originalPoreLengths[0] = m_poreLength[0];
            m_originalPoreLengths[1] = m_length;
            m_originalPoreLengths[2] = m_poreLength[1];

            if(first->isEntryOrExitRes())
                m_poreLength[0] = 0.0;
            else if(throatEnd > inletBdr && throatEnd < outletBdr)               // Both pore 2 and throat are within the box
                m_poreLength[0] *= (bdr - throatEnd)/(m_poreLength[0]*scaleFact);
            else if(throatStart > inletBdr && throatStart < outletBdr)      // Only pore 2 is fully within box
            {
                m_poreLength[0] = 0.0;
                m_length *= (bdr - throatStart)/(m_length*scaleFact);
            }
            else                                                            // Pore 2 is only partially within box
            {
                m_poreLength[0] = 0.0;
                m_length = 0.0;
            }
        }
    }

    if(m_poreLength[0] > 1.1*oldPOneLen || m_poreLength[1] > 1.1*oldPTwoLen || m_length > 1.1*oldThrLen)
    {
        cout<< endl
            << "==============================================="    << endl
            << "Warning: The new lengths for elements connected"    << endl
            << "to the pressure boundary are larger than the"       << endl
            << "original ones. The lengths should be smaller"       << endl
            << "since we do not want pressure drops occuring"       << endl
            << "outside the box across which we're calculating"     << endl
            << "relative permeability."                             << endl
            << "============================================== "     << endl
            << endl;
    }

    m_connections.push_back(first);         // Add pore connections
    m_connections.push_back(second);

    //double minRad(100), maxRad(0.0), radSum(0.0);
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








void Throat::modifyLength(double scaleFactor)
{
    m_poreLength[0] *= scaleFactor;
    m_poreLength[1] *= scaleFactor;
    m_length *= scaleFactor;

    if(m_originalPoreLengths)
    {
        m_originalPoreLengths[0] *= scaleFactor;
        m_originalPoreLengths[1] *= scaleFactor;
        m_originalPoreLengths[2] *= scaleFactor;
    }
}

double Throat::lenToRadRatio() const
{
    assert(m_poreLength.size() == 2);
    return (m_poreLength[0]+m_poreLength[1]+m_length)/m_elemModel->radius();
}


const Node* Throat::node() const
{
    //cerr << "No node associated with throats" << endl;    //exit(-1);
    //return m_connections[0]->node();
    return &m_node;
}

Node* Throat::node()
{
    //cerr << "No node associated with throats" << endl;   // exit(-1);
    return &m_node;
}


