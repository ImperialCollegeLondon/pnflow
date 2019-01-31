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
#include "layerApex.h"








void LayerApex::getCAApexDistUntraped(double& apxDist, double& conAng, const double& halfAng, double pc, double intfacTen, bool itr) const
{
    double delta = itr ? 0.0: SMALL_NUM;
    if(pc > m_advancingPc-delta && pc < m_receedingPc+delta)
    {
        conAng = (acos(pc*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        conAng = std::max(std::min(PI, conAng), 0.0);
        apxDist = m_initedOLApexDist;
    }
    else
    {
		apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
		if(itr && apxDist <=  0.0) apxDist=1.0E-15;
		softAssert(apxDist > 0.0);
		softAssert(conAng >= 0.0 && conAng <=  PI);
	}
}


/*double LayerApex::hingingConAngUntraped(double pc, double conAng, double halfAng, double intfacTen, bool accurat) const
{
    double delta = accurat ? 0.0: SMALL_NUM;
    if(pc > m_advancingPc-delta && pc < m_receedingPc+delta)
    {
        double hingAng(acos(pc*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        hingAng = std::max(std::min(PI, hingAng), 0.0);
        return hingAng;
    }

    softAssert(conAng >= 0.0 && conAng <=  PI, "k");
    return conAng;
}

double LayerApex::getApexDistanceUntraped(double pc, double conAng, double halfAng, double intfacTen, bool itrRoutine) const
{
	if(pc > m_advancingPc-SMALL_NUM && pc < m_receedingPc+SMALL_NUM)
        return m_initedOLApexDist;
	
    double apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
    if(itrRoutine && apxDist <=  0.0) return 1.0E-15;
    softAssert(apxDist > 0.0, "iu");
    return apxDist;
}
*/



void LayerApex::getCAApexDist(double& apxDist, double& conAng, const double& halfAng, double pc, double intfacTen, bool itr) const
{

	softAssert(m_inited || m_trappedCL.first>-1 || m_parentShape->eleman()->trappingWatBulk().first>-1);	

    const double delta = itr ? 0.0: SMALL_NUM;
    if(m_trappedCL.first>-1)
    {
        conAng = (acos(m_trappedCL.second*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        conAng = std::max(std::min(PI, conAng), 0.0);
        apxDist = m_initedOLApexDist;
        //softAssert(conAng >= 0.0 && conAng <=  PI, "j");
        return;
    }
    else if (m_parentShape->eleman()->trappingWatBulk().first>-1)
    {
        conAng = (acos(m_parentShape->eleman()->trappingWatBulk().second*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        conAng = std::max(std::min(PI, conAng), 0.0);
        //softAssert(conAng >= 0.0 && conAng <=  PI, "j");
        apxDist = m_initedOLApexDist;
        return;
	}
    else if (pc > m_advancingPc-delta && pc < m_receedingPc+delta)
    {
        conAng  = (acos(pc*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        conAng = std::max(std::min(PI, conAng), 0.0);
        //softAssert(conAng >= 0.0 && conAng <=  PI, "j");
        apxDist = m_initedOLApexDist;
        return;
    }
    else
	{
		softAssert(conAng >= 0.0 && conAng <=  PI);
		
		apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
    

		if(itr && apxDist <=  0.0) apxDist = 1.0E-15;
		softAssert(apxDist > 0.0);
		if(apxDist < 0.0)
		{
			cout<< apxDist<<"  "<< pc<<"  "<< m_receedingPc<<"  "<< m_advancingPc<<"  "<<conAng<<"  "<<halfAng<<"  "<<intfacTen<<"  "<< itr<<endl;
		}
	}






}


/*double LayerApex::hingingConAng(double pc, double conAng, double halfAng, double intfacTen, bool accurat) const
{
    const double delta = accurat ? 0.0: SMALL_NUM;
    if(m_trappedCL.first>-1)
    {
        double hingAng(acos(m_trappedCL.second*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        hingAng = std::max(std::min(PI, hingAng), 0.0);
        //softAssert(hingAng >= 0.0 && hingAng <=  PI, "j");
        return hingAng;
    }
    else if (m_parentShape->eleman()->trappingWatBulk().first>-1)
    {
        double hingAng(acos(m_parentShape->eleman()->trappingWatBulk().second*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        hingAng = std::max(std::min(PI, hingAng), 0.0);
        //softAssert(hingAng >= 0.0 && hingAng <=  PI, "j");
        return hingAng;
	}
    else if (pc > m_advancingPc-delta && pc < m_receedingPc+delta)
    {
        double hingAng(acos(pc*m_initedOLApexDist*sin(halfAng)/intfacTen)+halfAng);
        hingAng = std::max(std::min(PI, hingAng), 0.0);
        //softAssert(hingAng >= 0.0 && hingAng <=  PI, "j");
        return hingAng;
    }
    else
	{
		softAssert(conAng >= 0.0 && conAng <=  PI, "k");
		return conAng;
	}
}

double LayerApex::getApexDistance(double pc, double conAng, double halfAng, double intfacTen, bool itrRoutine) const
{
	softAssert(m_inited || m_trappedCL.first>-1 || m_parentShape->eleman()->trappingWatBulk().first>-1, " qljt ");	

	if( m_trappedCL.first>-1 || m_parentShape->eleman()->trappingWatBulk().first>-1 || (pc > m_advancingPc-SMALL_NUM && pc < m_receedingPc+SMALL_NUM) )
        return m_initedOLApexDist;
	
    double apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
    

    if(itrRoutine && apxDist <=  0.0) return 1.0E-15;
    softAssert(apxDist > 0.0, "i");
    if(apxDist < 0.0)
    {
		cout<< apxDist<<"  "<< pc<<"  "<< m_receedingPc<<"  "<< m_advancingPc<<"  "<<conAng<<"  "<<halfAng<<"  "<<intfacTen<<"  "<< itrRoutine<<endl;
	}
    return apxDist;
}
*/








//bool LayerApex::create Film(double pc, double conAng, double maxSpontConAng, double halfAng, double intfacTen)
bool LayerApex::createOLayer(double pc, double conAngRec, double conAngAdv, double maxSpontConAng, double halfAng, double intfacTen, bool oilInj)
{
    softAssert(!m_exists);
	double conAng = oilInj ? conAngRec: conAngAdv;
    m_exists = conAng > PI/2.0 + halfAng  && conAng > maxSpontConAng;
    
    m_gravityCorrection = m_parentShape->rhogh();
    

    if(!exists()) return false;


	m_initedOLApexDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
		
		
	if(m_initedOLApexDist < m_innerCornerApex->pinnedApexDist())
	{
		 m_initedOLApexDist= -1.0; ///. used in re-creation
		 m_exists = false;
		return false;                               
	}

	softAssert(m_initedOLApexDist>0);
  

    m_advancingPc = intfacTen*cos(conAngAdv-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set
    m_receedingPc = intfacTen*cos(conAngRec-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set

    m_inited = true; ///. warning this is used by lay erCollapsePc_ fromEitherSide


        m_entryPc = layerCollapsePc(pc, conAng, halfAng, intfacTen,oilInj);

	m_exists = pc > m_entryPc; 
	if (!m_exists) { m_inited=false; m_initedOLApexDist= -1.0; }
	softAssert(m_initedOLApexDist>0);





    return m_exists;
}

bool LayerApex::finitLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool overwriteTrapping)
{   ///. called from initOilInjection/initWaterInjection
    if(!exists()) return true;
    softAssert(m_inited ||  m_trappedCL.first > -1 || m_parentShape->eleman()->trappingWatBulk().first>-1); //,"opin"
    softAssert(!m_parentShape->containCOil());//,"finx"

	double conAng = oilInj ? conAngRec : conAngAdv;

	int nTraps= m_parentShape->eleman()->isTrappedOil() + (m_parentShape->eleman()->trappingWatBulk().first>-1); 
	if( m_inited || /*(overwriteTrapping && nTraps  < 2) ||*/ nTraps < 1 )
	{ 

		///. m_initedOLApexDist = 
		getCAApexDist(m_initedOLApexDist, conAng, halfAng, pc, intfacTen);//const
		
		softAssert(m_initedOLApexDist > m_innerCornerApex->pinnedApexDist());
		m_inited = false;


	}

    m_advancingPc = intfacTen*cos(conAngAdv-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set
    m_receedingPc = intfacTen*cos(conAngRec-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set  

	if(m_initedOLApexDist < m_innerCornerApex->pinnedApexDist())
	{
		 cout<<" fdraw "; cout.flush();
		 cout<<"  "<<m_initedOLApexDist<<"   " << m_innerCornerApex->pinnedApexDist(); cout.flush();
		 //removeLayer();
		return false; 
	}




  return true;
}



bool LayerApex::initLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool silent)
{   ///. called from initOilInjection/initWaterInjection


    if( (!oilInj && !exists()) || m_trappedCL.first > -1 || m_parentShape->eleman()->trappingWatBulk().first>-1 )
    {
		return true;
	}
	else if (oilInj && !exists())
	{ ///. calc layer creation Pc
		if (conAngRec <= PI/2.0 + halfAng  || m_initedOLApexDist <= 0.0 ) return true; ///. will not form
    
		m_inited = true;

		m_initedOLApexDist = m_parentShape->radius()/tan(halfAng);//(intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);

		m_advancingPc = intfacTen*cos(conAngAdv-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set
		m_receedingPc = intfacTen*cos(conAngRec-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set

		m_entryPc = layerCollapsePc(pc, conAngRec, halfAng, intfacTen, oilInj); ///.  set
		m_inited = false;

		return true;
	}
	else
	{


		int WarningReInitLayers;


		bool stable = true;


		if(m_initedOLApexDist < m_innerCornerApex->pinnedApexDist())
		{
			 cout<<"\nWarning118   ";cout.flush();
			//removeLayer();

			stable = false;                               
		}


		m_advancingPc = intfacTen*cos(conAngAdv-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set
		m_receedingPc = intfacTen*cos(conAngRec-halfAng)/(m_initedOLApexDist*sin(halfAng)); ///.  set
		m_inited = true;

		double currConAng = oilInj ? conAngRec : conAngAdv;
		m_entryPc = layerCollapsePc(pc, currConAng, halfAng, intfacTen, oilInj); ///.  set

		if(pc < m_entryPc)
		{
			 cout<<"\nWarning168   ";cout.flush();
			 //removeLayer();
		      stable = false;                               
		}

	    return stable;
	}
}




double LayerApex::layerCollapsePc_FromCentre(double outPc, double inPc, double conAng, double halfAng,
                                           double intfacTen) const
{
    softAssert(m_trappedCL.first<0|| m_parentShape->eleman()->trappingWatBulk().first<0);

    double inerApexDist, innerConAng(m_parentShape->minInitRecCntAng()); m_innerCornerApex->getCApexDistConAng(inerApexDist, innerConAng, inPc, halfAng, intfacTen);
    softAssert(innerConAng != conAng);
    double ki = (intfacTen/inPc)*(cos(innerConAng)/sin(halfAng)-1.0);

    if(m_inited)
    {
        double slipPrs = m_advancingPc;
        double ko = (intfacTen/slipPrs)*(cos(m_advConAng)/sin(halfAng)+1.0);
        if(ko - ki < NEG_ALMOST_ZERO)
        {
            double oldRad(intfacTen/(slipPrs+SMALL_NUM));    // Don't get caught by rounding errors
            for(int itr = 0; itr < MAX_ITR; ++itr)
            {
                double outConAng = acos(m_initedOLApexDist*sin(halfAng)/oldRad)+halfAng;
                double outConAngDr = m_initedOLApexDist*sin(halfAng)/(sin(outConAng-halfAng)*oldRad*oldRad);
                double func = oldRad*(cos(outConAng)/sin(halfAng)+1.0)-ki;
                double funcDr = (cos(outConAng)-oldRad*outConAngDr*sin(outConAng))/sin(halfAng) + 1.0;
                double newRad = oldRad - func/funcDr;
                if(fabs((newRad-oldRad)/newRad) < EPSILON)
                {	m_colType =30;

                    double pc(intfacTen/newRad);
                    softAssert(fabs(ki-newRad*(cos(outConAng)/sin(halfAng)+1.0)) < EPSILON);
                    return pc;
                }
                oldRad = newRad;
            }

			std::cerr << std::endl
				<< "============================================"   << std::endl
				<< "Error: failed to converge to a solution for"    << std::endl
				<< "layer collapsing pressure.   "                     << std::endl
				<< "============================================"   << std::endl;    exit(-1);
			return 0.0;

        }
    }
    else if(m_exists) 
    {
		cout<<"b";cout.flush();
		//m_parentShape->ChParent()[1000000].setInOilFloodVec(true);

	}

	m_colType =31;


    //double inerApexDist = m_innerCornerApex->getCApexDistance(inPc, innerConAng, halfAng, intfacTen);
    double collPc1 = intfacTen*cos(m_advConAng-halfAng)/(m_innerCornerApex->pinnedApexDist()*sin(halfAng)); 

    double collPc2 = intfacTen*(cos(m_advConAng)+sin(halfAng))/(ki*sin(halfAng));
    double collPc = max(collPc1, collPc2);
    softAssert(collPc < 0.0);
    if (collPc > 0.0) cout<<"CA"<<m_advConAng<<" "<<conAng<<" ";
    return collPc;
}

double LayerApex::layerCollapsePc_FromCorner(double outPc, double inPc, double conAng, double halfAng, double intfacTen) const
{
    softAssert(m_trappedCL.first<0);



    double outerConAng(conAng), apxDistDummy;// = hingingConAng(outPc, conAng, halfAng, intfacTen);
    getCAApexDist(apxDistDummy, outerConAng, halfAng, outPc, intfacTen);
    double ko = (intfacTen/outPc)*(cos(outerConAng)/sin(halfAng)+1.0);
     //m_inited ? acos(m_initedOLApexDist*sin(halfAng)*outPc/intfacTen)+halfAng: conAng;
    double innerSlipPrs = m_innerCornerApex->advancingPc();
    double maxConAng = m_advConAng;
    double ki = (intfacTen/innerSlipPrs)*(cos(maxConAng)/sin(halfAng)-1.0);

    if(ko - ki < NEG_ALMOST_ZERO  && inPc < m_innerCornerApex->advancingPc())
    {
        double oldRad(intfacTen/(innerSlipPrs+SMALL_NUM));
        double bi = m_innerCornerApex->pinnedApexDist();
        for(int itr = 0; itr < MAX_ITR; ++itr)
        {
            softAssert((m_innerCornerApex->pinnedApexDist()*sin(halfAng)/oldRad) >= -1.0 &&
                (m_innerCornerApex->pinnedApexDist()*sin(halfAng)/oldRad) <=  1.0);

            double inConAng = acos(bi*sin(halfAng)/oldRad)-halfAng;
            double inConAngDr = bi*sin(halfAng)/(sin(inConAng+halfAng)*oldRad*oldRad);
            double func = oldRad*(cos(inConAng)/sin(halfAng)-1.0)-ko;
            double funcDr = (cos(inConAng)-oldRad*inConAngDr*sin(inConAng))/sin(halfAng) - 1.0;

            double newRad = oldRad - func/funcDr;
            if(fabs((newRad-oldRad)/newRad) < EPSILON)
            {				m_colType =10;

                double pc(intfacTen/newRad);
                softAssert(fabs(ko-(intfacTen/pc)*(cos(inConAng)/sin(halfAng)-1.0)) < EPSILON);//, "outTrap"
                return pc;
            }
            oldRad = newRad;
        }
        //return Error_layer CollapsePc();

			std::cerr << std::endl
				<< "============================================"   << std::endl
				<< "Error: could not converge to a solution for"    << std::endl
				<< "layer collapsing pressure from corner."                     << std::endl
				<< outPc<<"  "<< inPc << "  " << m_advancingPc                << std::endl
				<< bi << "  "                << std::endl
				<< conAng << "  " << m_advConAng                << std::endl
				<< "  " << halfAng                << std::endl
				<< m_innerCornerApex->advancingPc() << "  " << m_innerCornerApex->receedingPc()                 << std::endl
				//<< m_innerCornerApex->trappingCL().first << "  " << m_innerCornerApex->getCApexDistance(inPc, conAng, halfAng, intfacTen)      << std::endl 
				<< "============================================"   << std::endl;    exit(-1);
				
        int warning;
        
		std::cerr << std::endl
			<< "============================================"   << std::endl
			<< "Error: Failed to converge to a solution     "    << std::endl
			<< "for layer collapsing pressure.      "                     << std::endl
			<< "============================================"   << std::endl;    exit(-1);
		return 0.0;
    }

	m_colType =11;
    return innerSlipPrs;
}





inline double dimLessCornerArea(double halfAng, double contactAng) 
{
    if(fabs(contactAng + halfAng - ElemModel::PI/2.0) < 0.01)
    {
        return sin(halfAng)*cos(halfAng);
    }
    else
    {
        return pow(sin(halfAng)/cos(halfAng+contactAng), 2.0)
        * (cos(contactAng)*cos(halfAng+contactAng)/sin(halfAng) + halfAng + contactAng - ElemModel::PI/2.0);
    }
}

double LayerApex::layerCollapsePc_fromEitherSide(double pc, double conAng, double halfAng, double intfacTen,bool debug) const
{
    softAssert(m_trappedCL.first<0);

//if (debug)
//{
    if((m_exists && !m_inited) || (!m_exists && m_inited))
    {
    softAssert(m_inited); 
    //softAssert(m_exists,"kfh");
		//m_parentShape->ChParent()[1000000].setInOilFloodVec(true);
	}
//}
    //if(m_inited)
    //{
	double slipPrs = m_inited ? m_advancingPc : pc;
	double innerSlipPrs = m_innerCornerApex->cornerExists() ? m_innerCornerApex->advancingPc() : pc;
	softAssert(slipPrs < 0.0);
	//double innerHingConAng = m_innerCornerApex->cornerExists() ? m_innerCornerApex->getCA pexDistConAng(innerSlipPrs, conAng, halfAng, intfacTen) : m_parentShape->minInitRecCntAng();
    double bi, innerHingConAng(m_parentShape->minInitRecCntAng()); 
    if (m_innerCornerApex->cornerExists()) m_innerCornerApex->getCApexDistConAng(bi, innerHingConAng, innerSlipPrs, halfAng, intfacTen,true);


	double /*ki*/innerApexDistThroughCentre = (intfacTen/innerSlipPrs)*(cos(innerHingConAng)/sin(halfAng)-1.0);
	double ko = (intfacTen/slipPrs)*(cos(m_advConAng)/sin(halfAng)+1.0);
softAssert(ko>0);
	if(ko - innerApexDistThroughCentre < NEG_ALMOST_ZERO   && pc < innerSlipPrs +0.1)//
	{
		double oldRad(intfacTen/(slipPrs+SMALL_NUM));
		
		//m_innerCornerApex->getCApexDistConAng(bi, innerHingConAng, innerSlipPrs, halfAng, intfacTen,true);
		
		double bo = m_initedOLApexDist;
		for(int itr = 0; itr < MAX_ITR; ++itr)
		{
			double inConAng = acos(bi*sin(halfAng)/oldRad)-halfAng;
			double outConAng = acos(bo*sin(halfAng)/oldRad)+halfAng;
			double inConAngDr = bi*sin(halfAng)/(sin(inConAng+halfAng)*oldRad*oldRad);
			double outConAngDr = bo*sin(halfAng)/(sin(outConAng-halfAng)*oldRad*oldRad);

			double func = (cos(inConAng)-cos(outConAng))/sin(halfAng)-2.0;
			double funcDr = (outConAngDr*sin(outConAng)-inConAngDr*sin(inConAng))/sin(halfAng);
			double newRad = oldRad - func/funcDr;
			if(fabs((newRad-oldRad)/newRad) < EPSILON)
			{
				softAssert(fabs(newRad*(cos(innerHingConAng)/sin(halfAng)-1.0) -
					newRad*(cos(m_advConAng)/sin(halfAng)+1.0)) < EPSILON);
				softAssert(intfacTen/newRad < 0.0);


			m_colType =20;
				
				return intfacTen/newRad;
			}
			oldRad = newRad;
		}
		//return Error _layer CollapsePc();


		std::cerr << std::endl
			<< "============================================"   << std::endl
			<< "Error: could not converge to a solution for"    << std::endl
			<< "layer collapsing pressure.         "     << std::endl
			<< pc << "  " << m_advancingPc                << std::endl
			<< bi << "  " << bo                << std::endl
			<< conAng << "  " << m_advConAng                << std::endl
			<< innerHingConAng << "  " << halfAng                << std::endl
			<< m_innerCornerApex->advancingPc() << "  " << m_innerCornerApex->receedingPc()                 << std::endl
			<< m_innerCornerApex->trappingCL().first << "  " << bi             << std::endl
			<< "============================================"   << std::endl;    exit(-1);
			
							
			
			
			
		return 0.0;
	}
	//else if(pc >= innerSlipPrs)
	{	
	double collPc1;


		innerApexDistThroughCentre=/*m_innerCornerApex->getCApexDistanceUntraped(pc, innerHingConAng, halfAng, intfacTen)**/
		bi*(cos(halfAng)+sin(halfAng)*sin(0.5*(innerHingConAng+halfAng-PI/2.0))/cos(0.5*(innerHingConAng+halfAng-PI/2.0)));
		//double outertouchApexDist=innerApexDistThroughCentre/(cos(halfAng)+sin(halfAng)*cos(PI-conAng+halfAng)/sin(PI-conAng+halfAng))WRONG;




		double outertouchApexRad=innerApexDistThroughCentre/(cos(m_advConAng)/sin(halfAng)+1.0);
		collPc1= intfacTen/outertouchApexRad;//(outertouchApexDist*sin(halfAng)); 
		





		double pinnedOLApexDist = outertouchApexRad*cos(m_advConAng-halfAng)/(sin(halfAng));
		
		softAssert (dimLessCornerArea(halfAng,PI-m_advConAng)*pinnedOLApexDist*pinnedOLApexDist > dimLessCornerArea(halfAng,innerHingConAng)*bi*bi-1.0e-18);
		
		if (debug || !(dimLessCornerArea(halfAng,PI-m_advConAng)*pinnedOLApexDist*pinnedOLApexDist > dimLessCornerArea(halfAng,innerHingConAng)*bi*bi-1.0e-18))
		{
			cout<< "  layerCollapsePc_fromEitherSide:"<<"  inCorExists:"<< m_innerCornerApex->cornerExists()<<"  hAng:"<<halfAng<<"  CAAdv"<<m_advConAng<<"  ";
			//exit(-1);
		}
		int warningtocheck;

		collPc1=max(collPc1,m_innerCornerApex->advancingPc());
		

			m_colType = 21;


		double collPc2 = (m_innerCornerApex->cornerExists()) ? intfacTen*cos(m_advConAng-halfAng)/(bi*sin(halfAng)) : collPc1-0.1;





		if (debug)
		{
			
			double outertouchApexDist=outertouchApexRad*cos(conAng-halfAng)/sin(halfAng);
	double outerApexDistThroughCentre = outertouchApexRad*(cos(m_advConAng)/sin(halfAng)+1.0);

			cout<<m_colType<<"\n  innerHingConAng "<<innerHingConAng<<" "//<<m_innerCornerApex->getCApex DistConAng(innerSlipPrs, conAng, halfAng, intfacTen)<<" "
				<<"  m_advConAng "<<m_advConAng<<"  conAng"<<conAng
				<<"  halfAng "<<halfAng<<"  outertouchApexRad"<<outertouchApexRad<<"   pc " <<pc<<"  inSlipPc"<<innerSlipPrs<<endl
				<<"    innerApexDistThroughCentre "<<innerApexDistThroughCentre<<"    koColapse "<<outerApexDistThroughCentre
				<<"    inApexDist "//<<m_innerCornerApex->getCApexDistanceUntraped(innerSlipPrs, innerHingConAng, halfAng, intfacTen)<<"  outertouchApexDist"<<outertouchApexDist
				<<" pc:"<<pc<<" pci:"<<(m_innerCornerApex->cornerExists() ? m_innerCornerApex->advancingPc() : slipPrs)<<" pco:"<<slipPrs<<" Pc1:"<<collPc1<<" Pc2:"<<collPc2<<" max:"<<max(collPc1,collPc2)<<" \n";
		//exit(-1);
		}
			
		softAssert(collPc2 < 0.0);
		if (collPc1<collPc2)
		{
			m_colType = 22;
			//cout<<"colPcvk"<<collPc2<<"  "<<conAng<<"  "<<m_advConAng<<"  "<<halfAng<<endl;
			return collPc2;
		}		
		
    return collPc1;
	}
	
    
}





double LayerApex::layerCollapsePc(double pc, double conAng, double halfAng, double intfacTen, bool injOil) const
{
    double collPrs;
	m_colType=0;
    if(!m_innerCornerApex->cornerExists())
    {
        //softAssert(exists(), "quacki di quack");
        int WarningFiEnablesoftAssert;
        collPrs = LOWEST_LAYER_PC;                      // No water in corner  => never to collapse
    }
	else if(injOil && exists())                      // Oil Inj, stable layer  => no update needed
		collPrs = INF_NEG_NUM;
	else if(!injOil && !exists())               // Wat Inj, collapsed layer  => no update needed
		collPrs = INF_POS_NUM;
	else if(m_trappedCL.first>-1)                          // Wat Inj, both trapped  => never collapse
		collPrs = INF_NEG_NUM;
	else
    {
		bool freeWBulk=m_parentShape->eleman()->trappingWatBulk().first<0;
		bool freeWFilm=m_innerCornerApex->cornerExists() && m_innerCornerApex->trappedCorner().first<0;
		
		
		
		int Warning_unsynced_Trappings;
		 //softAssert(m_parentShape->eleman()->trappingWatFilm().first<0 && m_innerCornerApex->trappedCorner().first>-1, " Fq1");
		 //softAssert(m_parentShape->eleman()->trappingWatFilm().first>-1 && m_innerCornerApex->trappedCorner().first<0, " Fq2");
//
	//if(!(m_parentShape->eleman()->trappingWatFilm().first<0 && m_innerCornerApex->trappedCorner().first>-1))
	//{
		//conAng=m_parentShape[10000].conAngleRec();
	//}




		if(freeWBulk && freeWFilm) // No trapping
			collPrs = layerCollapsePc_fromEitherSide(pc, conAng, halfAng, intfacTen);
		else if(freeWBulk)  ///. Only Corner trapped
			collPrs = layerCollapsePc_FromCentre(pc, m_innerCornerApex->trappedCornerNOTTOBEUSED().second, conAng, halfAng, intfacTen);
		else if(freeWFilm)  ///. Only Centre trapped
			collPrs = layerCollapsePc_FromCorner(m_trappedCL.second, pc, conAng, halfAng, intfacTen);
		else if(injOil)                 // Oil Inj, both trapped  => never grow
			collPrs = INF_POS_NUM;
		else                            // Wat Inj, both trapped  => never collapse
		   collPrs = INF_NEG_NUM;
	}
       
    //if(!m_innerCornerApex->cornerExists())
    //{
        //softAssert(exists(), "quacki di quack");
        //collPrs = LOWEST_LAYER_PC;                      // No water in corner  => never to collapse
    //}
    //else if(injOil && exists())                      // Oil Inj, stable layer  => no update needed
        //collPrs = INF_NEG_NUM;
    //else if(!injOil && !exists())               // Wat Inj, collapsed layer  => no update needed
        //collPrs = INF_POS_NUM;
    //else if(m_innerCornerApex->trappedCornerNOTTOBEUSED().first<0 && m_trappedCL.first<0) // No trapping
        //collPrs = layerCollapsePc_fromEitherSide(pc, conAng, halfAng, intfacTen);
    //else if(m_innerCornerApex->trappedCornerNOTTOBEUSED().first>-1 && m_trappedCL.first<0)  ///. Only Corner trapped
        //collPrs = layerCollapsePc_FromCentre(pc, m_innerCornerApex->trappedCornerNOTTOBEUSED().second, conAng, halfAng, intfacTen);
    //// else if(m_innerCornerApex->trappedCornerNOTTOBEUSED().first<0 && m_trappedCL.first>-1)  ///. Only Centre trapped
        //// collPrs = layerCollapsePc_FromCorner(m_trappedCL.second, pc, conAng, halfAng, intfacTen);
    //else if(injOil)                 // Oil Inj, both trapped  => never grow
        //collPrs = INF_POS_NUM;
    //else                            // Wat Inj, both trapped  => never collapse
       //collPrs = INF_NEG_NUM;

    return collPrs;
}




