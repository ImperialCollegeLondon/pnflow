
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "sortedEvents.h"
#include "compareFuncs.h"



/**
* Before invoking the solver and calculating saturations we must update the state of the
* element, ie calculate saturation and conductance. Since the configuration in imbibition can
* be more complex than in drainage, the calculations are separate for each event
*/
double Polygon::calcR(double pc)
{

		//m_sumOilLayerAreas = 0.0;
	if(m_bulkFluid == &m_comn.water() && m_numLayers == 0)
	{
        m_SatWater = 1.0;
        m_conductanceOil = 0.0;
        m_conductanceWater.first = 0.0;
        m_conductanceWater.second = SPConductance(m_area, m_comn.water().viscosity());
		

	}
	else if(m_bulkFluid == &m_comn.oil() && !m_waterInCorner[0].cornerExists())
	{
        m_SatWater = 0.0;
        m_conductanceOil = SPConductance(m_area, m_comn.oil().viscosity());
        m_conductanceWater.first = m_comn.KrwatcornAtSw0() * m_conductanceOil; ///. m_KrwatcornAtSw0 is read from TRAPPING keyword
        m_conductanceWater.second = 0.0;


	}
	else if(m_bulkFluid == &m_comn.oil() && m_waterInCorner[0].cornerExists())
    {
		calcR_oilWithWaterInCorners(pc-rhogh()-m_elem.snapOfLongitCurvature()*m_comn.oil().interfacialTen());


    }
	else
    {
		calcR_waterWithOilLayers(pc-rhogh()-m_elem.snapOfLongitCurvature()*m_comn.oil().interfacialTen());
    }

    softAssert(m_conductanceOil >= 0.0 && (m_conductanceWater.first >= 0.0 || m_conductanceWater.second >= 0.0));
    if ( !(m_conductanceOil >= 0.0 && (m_conductanceWater.first >= 0.0 || m_conductanceWater.second >= 0.0) ) )
		cout<<"Ks = "<< m_conductanceOil << "  " << m_conductanceWater.first << "  " << m_conductanceWater.second <<endl;
    ///. TODO: sensitivity analysis, add non-linearity by decoupling Aw_dl from Sw
    

	m_ElectricalConductance = m_SatWater*m_area/m_comn.water().resistivity() + sqrt(m_SatWater*m_area/m_shapeFactor) / m_comn.surface().resistivity();
	m_ElectricalConductance += m_area*(1.0-m_SatWater)/m_comn.oil().resistivity();





    if (m_comn.debugMode > 500 && m_comn.injectant() == &m_comn.water() && ( m_conductanceWater.first +m_conductanceWater.second  <  m_conductanceWaterOldToCheck.first +m_conductanceWaterOldToCheck.second  || m_conductanceWater.second < m_conductanceWaterOldToCheck.second) )
    {
		cout<< m_SatWater <<"  "<<   m_elem.trappingWatBulk().first <<"  "<< m_conductanceWater.second <<"  "<< m_conductanceWaterOldToCheck.second<<"     ";
		cout<< m_conductanceWater.first <<"  "<<   m_conductanceWaterOldToCheck.first <<"  "<< m_conductanceWater.second <<"  "<< m_conductanceWaterOldToCheck.second
		<<"  kw\\"<<endl; //exit(-1);
	}
    m_conductanceWaterOldToCheck.first = m_conductanceWater.first;
    m_conductanceWaterOldToCheck.second = m_conductanceWater.second;


    return m_SatWater; ///. saturation calculation
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
    //{
		//double cosTerm=cos(halfAng+contactAng);
        //return pow(sin(halfAng)/cos(halfAng+contactAng), 2.0)
        //* (cos(halfAng+contactAng)*sin(halfAng+contactAng)+cos(halfAng+contactAng)*cos(halfAng+contactAng)/tan(halfAng) + halfAng + contactAng - ElemModel::PI/2.0);
    //}
    //{
		//double cosTerm=cos(halfAng+contactAng);
        //return pow(sin(halfAng)/cosTerm, 2.0)
        //* (cosTerm*sin(halfAng+contactAng)+cosTerm*cosTerm/tan(halfAng) + halfAng + contactAng - ElemModel::PI/2.0);
    //}
}


/**
// The corner conductance version by Paal-Eric
*/
double cornerConductance(double dimLessCornerA, double meniscusApexDist, double halfAng,
                                  double contactAngle, double viscosity) 
{
    double cornerGstar((cos(halfAng)*sin(halfAng)) / (4.0*pow(1+sin(halfAng), 2.0)));
    double cornerG(cornerGstar);

    if(fabs(contactAngle + halfAng - ElemModel::PI/2.0) > 0.01)
    {
        cornerG = dimLessCornerA /
        (4.0 * pow(1.0 - (sin(halfAng)/cos(halfAng+contactAngle)) * (halfAng+contactAngle-ElemModel::PI/2.0), 2.0));
    }


    double cFactor = 0.364 + 0.28 * cornerGstar / cornerG;
    double dimLessCornerConduct = cFactor * dimLessCornerA * dimLessCornerA * cornerG;

    softAssert(dimLessCornerA != 0.0 && halfAng != 0.0 && cornerG != 0.0);

    return dimLessCornerConduct * pow(meniscusApexDist, 4.0) / viscosity;
}

/**
* My very own layer conductance routine. Wohoo :-)
*/
double layerConductance(double insideBC, double outsideBC, double layerArea, double cornerArea, double outerApexDist,
                                 double innerCornerApexDist, double innerConAng, double outerConAng, double halfAng, double viscosity)
{
    double coeff[3];
    coeff[0] = -2.4010e-002;
    coeff[1] = 2.8402e-001;
    coeff[2] = -2.9531;

    if(halfAng < 10.0*ElemModel::PI/180.0)
    {
        coeff[0] = -1.0610E-002;
        coeff[1] = 5.1608E-001;
        coeff[2] = -2.0645;
    }
    else if(halfAng < 20.0*ElemModel::PI/180.0)
    {
        coeff[0] = -2.6810E-002;
        coeff[1] = 1.8672E-001;
        coeff[2] = -3.5977;
    }
    else if(halfAng < 30.0*ElemModel::PI/180.0)
    {
        coeff[0] = -4.4021E-002;
        coeff[1] = -6.3195E-002;
        coeff[2] = -4.3749;
    }
    else if(halfAng < 40.0*ElemModel::PI/180.0)
    {
        coeff[0] = -3.1523E-002;
        coeff[1] = 1.6948E-001;
        coeff[2] = -3.3600;
    }
    else if(halfAng < 50.0*ElemModel::PI/180.0)
    {
        coeff[0] = -3.1367E-002;
        coeff[1] = 1.9327E-001;
        coeff[2] = -3.2673;
    }
    else if(halfAng < 60.0*ElemModel::PI/180.0)
    {
        coeff[0] = -2.3201E-002;
        coeff[1] = 3.1178E-001;
        coeff[2] = -2.9251;
    }
    else
    {
        coeff[0] = -3.5760E-002;
        coeff[1] = -6.5370E-004;
        coeff[2] = -4.7019;
    }

    double dimlessAo(layerArea / (outerApexDist*outerApexDist));
    double dimlessBi(innerCornerApexDist/outerApexDist);

    double outerOWlen = - (2.0*sin(halfAng)*(halfAng-outerConAng+ElemModel::PI/2.0))
        / cos(ElemModel::PI-outerConAng+halfAng);

    double innerOWLen = - (2.0*dimlessBi*sin(halfAng)*(halfAng+innerConAng-ElemModel::PI/2.0))
        / cos(innerConAng+halfAng);

    double perimeterLen = outerOWlen + innerOWLen + 2.0*(1.0-dimlessBi);
    double shapeFactorOil = dimlessAo / (perimeterLen*perimeterLen);
    double xGroup = log(pow(dimlessAo, 3.0)*shapeFactorOil);

    double dimlessOilCond = exp(coeff[0]*xGroup*xGroup + coeff[1]*xGroup + coeff[2]);
    double dimCond = pow(outerApexDist, 4.0) * dimlessOilCond / viscosity;
    return dimCond;
}



/**
* Water will remain in some of the corners of the polygon, Conductance and area is calculated
* for each corner.
*/
void Polygon::calcR_oilWithWaterInCorners(double cappPressure)
{
    double conAng(0.0), cornerAreas(0.0), cornerCond(0.0);
    //bool overRideTrapping(false);
	double pcSuggested = cappPressure;
	if(m_comn.injectant() == &m_comn.water())
    {
        conAng = m_cntAngAdv;
    }
	else
	{
		conAng = m_cntAngRec;
	}
	double cappPcWL = cappPressure;
    if( m_elem.trappingOil().first>-1 ) cappPcWL = m_elem.trappingOil().second;


    for(int i = 0; i < m_numCorners; ++i)
    {
		double apexDist(0.0), conAngCurr(conAng); 
        if(m_waterInCorner[i].cornerExists())
        {
			if( m_waterInCorner[i].trappingCL().first>-1 ) cappPcWL = m_waterInCorner[i].trappingCL().second;

			m_waterInCorner[i].getCApexDistConAng(apexDist, conAngCurr, cappPressure, m_crnHafAngs[i], m_comn.oil().interfacialTen());
			if (apexDist >radius()/tan(m_crnHafAngs[0])*2) 	{	cout<<" sxp ";	}
            
            if(m_comn.debugMode>1 &&(apexDist < 0.0))         cout<<" apexD"<<apexDist<<" ";

            double dlCornerArea = dimLessCornerArea(m_crnHafAngs[i], conAngCurr);
			double conductance = cornerConductance(dlCornerArea, apexDist, m_crnHafAngs[i],
				conAngCurr, m_comn.water().viscosity());

            softAssert(conductance > 0.0);
            cornerCond += conductance;
			cornerAreas += apexDist * apexDist * dlCornerArea;


        }

    }
    m_SatWater = min(cornerAreas/m_area,1.0);
    softAssert( m_SatWater < 1.0 && cornerCond > 0.0);
    
	if ( m_comn.debugMode>10 && !(cornerAreas/m_area < 1.0 && cornerCond > 0.0))
	{	cout<< "============ !(cornerAreas/m_area < 1.0 && cornerCond > 0.0  =================";
		exit(-1);
	}
	
    double oilSat(max(1.0-m_SatWater,0.0));
    m_conductanceOil = SPConductance(m_area, m_comn.oil().viscosity()) * oilSat;
    m_conductanceWater.first = cornerCond;
    m_conductanceWater.second = 0.0;
    


}


/**
* If oil is present at the same time that water occupies the center there must be
* oil layer sandwisched between water in corners and center.
*/
void Polygon::calcR_waterWithOilLayers(double cappPressure)
{
    double conAng(0.0);
    double cornerAreaWat(0.0), cornerCondWat(0.0);
    double cornerAreaOil(0.0), cornerCondOil(0.0);
    double layerAreaOil(0.0), layerCondOil(0.0);
    //bool overRideTrapping(false);

    if(m_comn.injectant() == &m_comn.water())		
        conAng = m_cntAngAdv;
	else	
    {
        conAng = m_cntAngRec;
    }
	//if(cappPressure > m_snapOffPrs && m_oilLayer[0].trappedOLayer().first<0)
		//cappPressure = m_snapOffPrs-EPSILON;

	double cappPcOL = cappPressure;
    if( m_elem.trappingWatBulk().first>-1 ) cappPcOL = m_elem.trappingWatBulk().second;

    double baseLength(m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1])));
    double oilApexBaseLengths(0.0), watApexBaseLengths(0.0);

    for(int i = 0; i < m_numCorners; ++i)
    {

      if(m_oilLayer[i].LayerApex::exists(/*st ab le*/))
      {
		double cappPcOL = cappPressure;
		if( m_oilLayer[i].trappingCL().first>-1 ) cappPcOL = m_oilLayer[i].trappingCL().second;
		
        if(m_waterInCorner[i].cornerExists())
        {
			double cappPcWL = cappPcOL;
			if( m_waterInCorner[i].trappingCL().first>-1 ) cappPcWL = m_waterInCorner[i].trappingCL().second;

            double layerPC(cappPressure); ///. trapping pc handled in calcR _getIFDistAng
            //if(!m_oilLayer[i].LayerApex::freeAtPrs(layerPC))
                //layerPC = m_oilLayer[i].LayerApex::layerCollPc()+EPSILON;

            //double outerConAng(conAng), outerApexDist(0.0);
            //m_oilLayer[i].LayerApex::calcR_getLIFDistAng(outerConAng, outerApexDist, layerPC, m_crnHafAngs[i], m_comn.oil().interfacialTen(), m_cntAngAdv, m_cntAngRec);
			double outerConAng(conAng);// = m_oilLayer[i].hingingConAng(cappPressure,conAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
			double outerApexDist;// = m_oilLayer[i].getApexDistance(cappPressure,outerConAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
			m_oilLayer[i].getCAApexDist(outerApexDist, outerConAng,m_crnHafAngs[i], cappPressure, m_comn.oil().interfacialTen());

			softAssert(outerApexDist>0);

			softAssert(outerApexDist < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
			softAssert(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));

            double dlTotCornerArea = dimLessCornerArea(m_crnHafAngs[i], PI-outerConAng);
            double areaCornerTot = outerApexDist * outerApexDist * dlTotCornerArea;

			oilApexBaseLengths += outerApexDist;    // A permanaent check for not overstepping oil volumes
            
            if(i < 2) 
            { 
				softAssert( oilApexBaseLengths <=  baseLength*1.01);
			}

            //double innerConAng(conAng), innerCornerApexDist(0.0);
            //m_waterInCorner[i].calcR_getIFDistAng(innerConAng, innerCornerApexDist, layerPC, m_crnHafAngs[i],  m_comn.oil().interfacialTen(), m_cntAngAdv, m_cntAngRec);
            
			//double innerConAng = m_waterInCorner[i].hingingCConAng(cappPressure,conAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
			//double innerCornerApexDist = m_waterInCorner[i].getCApexDistance(cappPressure,innerConAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
				double innerCornerApexDist, innerConAng(conAng); 
				m_waterInCorner[i].getCApexDistConAng(innerCornerApexDist, innerConAng, cappPressure, m_crnHafAngs[i], m_comn.oil().interfacialTen());

			watApexBaseLengths+=innerCornerApexDist;
			
            double dlWatCornerArea = dimLessCornerArea(m_crnHafAngs[i], innerConAng);
            double areaWater = innerCornerApexDist * innerCornerApexDist * dlWatCornerArea;
            softAssert(outerApexDist > innerCornerApexDist);
            if (outerApexDist <=  innerCornerApexDist)
            {
				cout<<int(m_elem.isTrappedWat(bulkBlob))<<"  ";
				cout<<" eXiT "<<endl;			//	exit(-1);
			}



            double areaOil = max(areaCornerTot - areaWater,1.0e-50);

			double waterCond = cornerConductance(dlWatCornerArea, innerCornerApexDist, m_crnHafAngs[i],
                innerConAng, m_comn.water().viscosity());
            double oilCond = layerConductance(1.0, 1.0, areaOil, areaCornerTot, outerApexDist, innerCornerApexDist,
                innerConAng, outerConAng, m_crnHafAngs[i], m_comn.oil().viscosity());

            if (m_comn.debugMode>1 && ((areaWater*waterCond*(areaCornerTot - areaWater)*oilCond <=  0.0) || (oilApexBaseLengths >  baseLength*1.01 && i < 2)))
            {
				cout<<" if (m_comn.debugMode>1 && ((areaWater*waterCond*(areaCornerTot - areaWater)*oilCond <=  0.0) || (oilApexBaseLengths >  baseLength*1.01 && i < 2)))   "<< m_oilLayer[i].m_colType<<" ";
				exit(-1);
			}
            
			softAssert(waterCond > 0.0);

            cornerAreaWat += areaWater;
            cornerCondWat += waterCond;
            layerAreaOil += areaOil;
            layerCondOil += oilCond;
        }
        else 
        {
            //double currentConAng(conAng), apexDist(0.0);
            //m_oilLayer[i].LayerApex::calcR_getLIFDistAng(currentConAng, apexDist, cappPressure, m_crnHafAngs[i],  m_comn.oil().interfacialTen(), m_cntAngAdv, m_cntAngRec);
			double currentConAng(conAng);// = m_oilLayer[i].hingingConAng(cappPressure,conAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
			double apexDist;// = m_oilLayer[i].getApexDistance(cappPressure,currentConAng,m_crnHafAngs[i],  m_comn.oil().interfacialTen());
			m_oilLayer[i].getCAApexDist(apexDist,currentConAng,m_crnHafAngs[i], cappPressure, m_comn.oil().interfacialTen());
			
			softAssert(apexDist>0);


            double dlCornerArea = dimLessCornerArea(m_crnHafAngs[i], PI-currentConAng);
            double conductance = cornerConductance(dlCornerArea, apexDist, m_crnHafAngs[i],
                PI-currentConAng, m_comn.oil().viscosity());

            cornerCondOil += conductance;
            cornerAreaOil += apexDist * apexDist * dlCornerArea;
		if ( m_comn.debugMode>10 && !(layerAreaOil >= 0.0))
		{
			cout<<cornerAreaOil<<endl;
		}

        }
	  }
	  else
	  {
		  ///. everything water, 
	  }
    }

    softAssert(layerAreaOil > 0.0 || cornerAreaOil > 0.0);


    m_SatWater = 1.0 - (cornerAreaOil + layerAreaOil)/m_area;
    m_conductanceOil = layerCondOil + cornerCondOil;

    double watSat((m_area-cornerAreaOil-layerAreaOil-cornerAreaWat) / m_area);

    double centerWatCond = SPConductance(m_area, m_comn.water().viscosity()) * watSat;
    softAssert(m_comn.debugMode==0 || (watSat > 0.0 && watSat < 1.0));

    if(m_hasDisConectedCentreWCornerW)
    {
        softAssert(cornerCondWat > 0.0);
        
        m_conductanceWater.first = cornerCondWat;
        m_conductanceWater.second = centerWatCond;
    }
    else
    {
        m_conductanceWater.first = 0.0;
        m_conductanceWater.second = centerWatCond+cornerCondWat;
    }
    
}




/**
* Circular elements can only hold a single fluid.
*/
double Circle::calcR(double pc)
{

    if(m_bulkFluid == &m_comn.water())
    {
        m_SatWater = 1.0;
        m_conductanceOil = 0.0;
        m_conductanceWater.first = 0.0;
        m_conductanceWater.second = SPConductance(m_area, m_comn.water().viscosity());
		m_ElectricalConductance = m_area/m_comn.water().resistivity() + sqrt(m_area / m_shapeFactor) / m_comn.surface().resistivity();
    }
    else
    {
		m_SatWater = 0.0;
		m_conductanceOil = SPConductance(m_area, m_comn.oil().viscosity());
		m_conductanceWater.first = m_comn.KrwatcornAtSw0() * m_conductanceOil;
		m_conductanceWater.second = 0.0;
		m_ElectricalConductance = m_area/m_comn.oil().resistivity();
	}

	m_ElectricalConductance = (m_SatWater*m_area)/m_comn.water().resistivity() + sqrt(m_SatWater*m_area/m_shapeFactor) / m_comn.surface().resistivity();
	m_ElectricalConductance += m_area*(1.0-m_SatWater)/m_comn.oil().resistivity();

	return m_SatWater;
}


/**
 * Computed the conductance of a fluid occupying the center of the element
 *
 * There are some choices when it comes to compute the conducance of the
 * centre of an element.
 *
 * g_p(A) = 3r^2A / 20mju   (1)     or      g_p(A) = 3A^2G / 5mju   (2)
 *
 * g = g_p(A_t) * S_p       (A)     or      g = g_p(A_p)            (B)
 *
 * Most combinations seem to give similar results (with the exception of
 * 2B). Lets run with 2A....
 */
double Square::SPConductance(double area, double visc) const
{
    // return (m_R * m_R * area) / (7.1136 * visc);   // 1
    return (0.5623 * area * area * m_shapeFactor) / visc;       // 2
}

/**
 * Computed the conductance of a fluid occupying the center of the element
*/
double Triangle::SPConductance(double area, double visc) const
{
    // return (3.0 * m_R * m_R * area) / (20.0 * visc);   // 1
    return (3.0 * area * area * m_shapeFactor) / (5.0 * visc);      // 2
}

/**
// Computed the conductance of a fluid occupying the center of the element
*/
double Circle::SPConductance(double area, double visc) const
{
    return (m_R * m_R * area) / (8.0 * visc);
}







