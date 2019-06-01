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
#include <set>
#include <utility>
#include <cassert>
using namespace std;


#include "Element.h"
#include "elem_Model.h"
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"


const double    Apex::PI = acos(-1.0);
const double    Apex::EPSILON = 1.0E-10;
const double    Apex::INF_NEG_NUM = -1.0E21;
const double    Apex::LOWEST_LAYER_PC = -1.0E20;//TOCHANGE lowest possible oil in corner PC
const double    Apex::INF_POS_NUM = 1.0E21;
const double    Apex::NEG_ALMOST_ZERO = -1.0E-12;
const double    Apex::POS_ALMOST_ZERO = 1.0E-12;
const double    Apex::MOLECULAR_LENGTH = 1.0E-10;
const double    Apex::SMALL_NUM = 1.0E-7;
const int       Apex::MAX_ITR = 5000;

int Apex::nErrors = 0;
int Apex::debugMode = 0;



/**
 * sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
 */
void CornerApex::createFilm(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool isOilInj)
{
    ensure(!m_exists);
    ensure(!m_inited);

	if(m_exists) return;
		
	double conAng = isOilInj ? conAngRec: conAngAdv;
    m_exists = conAng < PI/2.0 - halfAng;

	
    if(!m_exists) return;



	m_initedApexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
	ensure(m_initedApexDist > 0.0);
	if (m_initedApexDist > m_parentShape->radius()/tan(halfAng)*2) 	{ cout<<" xkg ";	}
	creationPc=pc;


    m_advancingPc = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(m_initedApexDist*sin(halfAng));
    m_receedingPc = intfacTen*cos( min(PI, conAngRec)+halfAng )/(m_initedApexDist*sin(halfAng));
	m_inited = true;


    if(pc > m_initOrMaxPcHist)
    {
        m_initOrMinApexDistHist = m_initedApexDist;
        m_initOrMaxPcHist = pc;
    }

}
 

/**
 * sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
 */
void CornerApex::finitCornerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool overwriteTrapping)
{
    if(!m_exists) return;

	int nTraps=m_parentShape->eleman()->isTrappedWat(filmBlob) + (m_parentShape->eleman()->isTrappedOil() || m_outerLayerApex->trappingCL().first>-1);
	if( m_inited/*(overwriteTrapping && nTraps < 2)*/ || nTraps < 1)
	{
		double apxDistOld=m_initedApexDist;
		double apxDist;
		double conAng(oilInj ? conAngRec : conAngAdv); 
		getCApexDistConAng(apxDist, conAng, pc, halfAng, intfacTen, true); ///. warnign don't send m_initedApexDist directly
		if (apxDist <=0.0 || apxDist > m_parentShape->radius()/tan(halfAng)*2) 	{ cout<<" xkf ";	}

		m_initedApexDist = apxDist;
		
		m_advancingPc = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(m_initedApexDist*sin(halfAng));
		m_receedingPc = intfacTen*cos( min(PI, conAngRec)+halfAng )/(m_initedApexDist*sin(halfAng));

		if(pc > m_initOrMaxPcHist)
		{
			m_initOrMinApexDistHist = m_initedApexDist;
			m_initOrMaxPcHist = pc;
		}
		m_inited = false;
	}
	


}
/**
 * sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
 */
void CornerApex::initCornerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj)
{
    if(!m_exists || m_trappedCL.first > -1 || m_parentShape->eleman()->isTrappedOil() || m_outerLayerApex->trappingCL().first>-1) return;
	m_inited = true;

    
		double  recPc = intfacTen*cos( min(PI, conAngRec+halfAng) )/(m_initedApexDist*sin(halfAng));
		if (m_receedingPc < recPc) m_receedingPc = recPc;  ///. to handle strange wettability change

		m_advancingPc = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(m_initedApexDist*sin(halfAng));


}





void CornerApex::getCApexDistConAng(double & apexDist, double & conAng, double pc, double halfAng, double intfacTen, bool overidetrapping, bool accurat, bool debug) const
{
    double delta = accurat ? 0.0: SMALL_NUM;
		bool Warning_ensure_deactivated=True;
		ensure(True || m_inited || m_trappedCL.first>-1 || m_parentShape->eleman()->isTrappedOil() || m_outerLayerApex->trappingCL().first>-1);
		ensure(m_exists);

	if(!m_exists)	{ apexDist = MOLECULAR_LENGTH;  return; }

    if(!overidetrapping)
    {
		
		apexDist=m_initedApexDist;
		if( m_trappedCL.first>-1)
		{
			pc=m_trappedCL.second;
			double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			if (debug) cout<< "d8sAA";
			return;
		}
		else if (m_parentShape->eleman()->trappingOil().first>-1)
		{
			pc=m_parentShape->eleman()->trappingOil().second;
			double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			if (debug) cout<< "d8sBB";
			return;
		}
		else if (m_outerLayerApex->trappingCL().first>-1)
		{
			pc=m_outerLayerApex->trappingCL().second;
			double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			if (debug) cout<< "d8sCC";

			return;
		}
	}
    
    
    if(m_advancingPc-delta <= pc && pc <= m_receedingPc+delta)
    {
        double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
        part = std::min(part, 0.999999);
        part = std::max(part, -0.999999);
        double hingAng(acos(part)-halfAng);
        ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
        conAng = (std::min(std::max(hingAng, 0.0), PI));
		apexDist=m_initedApexDist;
    }
    else if(pc < m_advancingPc)
    {

        conAng = m_parentShape->conAngleAdv();
			if (debug) 	cout<<"   efs1 ";

        apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
		if (apexDist < m_initedApexDist)
		{

			double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			apexDist = m_initedApexDist;
			if (debug) 	cout<<"   efs2 ";

		}
        else if (apexDist<0) cout<<"ErrApL<0";

	}
	else if(pc > m_initOrMaxPcHist)
    {   ///. recover pinning receeding contact angle

		conAng=m_parentShape->minInitRecCntAng();
        apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
			if (debug) 	cout<<"   efs3 ";

	}
    else  if(pc > m_receedingPc)
    {
        conAng = m_parentShape->conAngleRec();

        apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
		if (apexDist > m_initedApexDist)
		{

			
			double part(pc*m_initedApexDist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			
			if (debug) 	cout<<"   pcrecPcusds "<<pc<<" "<<m_receedingPc<<" "<<apexDist<<" <= "<<m_initedApexDist<<" "<<conAng<<"         ";

			apexDist = m_initedApexDist;			
			if (debug) 	cout<< apexDist / m_parentShape->radius()*tan(halfAng) <<"    ";
			if (debug) 	cout<<"   efs4 ";

		}
		else if(apexDist<m_initOrMinApexDistHist)
		{
			double part(pc*m_initOrMinApexDistHist*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng+SMALL_NUM >= m_parentShape->minInitRecCntAng() && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.0), PI));
			if (debug) cout<<"  pcrecPcuadsdsds "<<pc<<" "<<m_receedingPc<<"   ";		
			apexDist=m_initOrMinApexDistHist;
			if (debug) 	cout<<"   efs5 ";

		}
	}
	else
	{
        conAng = conAng;
        apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
			if (debug) 	cout<<"   efs6 ";

		ensure(conAng >= 0.0 && conAng <=  PI);
	}
	
	if (debug) if (apexDist > m_parentShape->radius()/tan(halfAng)*2) 
	{
		cout<<" syp ";
	}
	 if (debug) cout<<" pc"<<pc<<",conAng"<<conAng<<" "<<endl;

}


