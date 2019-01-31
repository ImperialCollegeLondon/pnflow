#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

void Circle::finitWaterInjection(double cappPrs)
{
    //softAssert(m_elem.connectionNum() > 0, "23");
    //m_elem.resetFillingEventRecord();
    //m_elem.setInWatFloodVec(false);
    //m_elem.setInOilFloodVec(false);
 
    //if(m_bulkFluid == &m_comn.water()) return;
 //
    //m _maxConAngSpont = PI/2.0;
    //m _Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;
  
    //c alcCentreEntryPrsWatInj();
} 

void Circle::initWaterInjection(double cappPrs)
{
    softAssert(m_elem.connectionNum() > 0);
    m_elem.resetFillingEventRecord();
    m_elem.setInWatFloodVec(false);
    //m_elem.setInOilFloodVec(false);
 
    if(m_bulkFluid == &m_comn.water()) return;
 
    //m _maxConAngSpont = PI/2.0;
    m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;
  	//timestep = -1;

    //calcCentreEntryPrsWatInj();
} 
 
///. move to pore/throat classes (?)
double Circle::centreEntryPrsWatInj()
{ 
    int num_WatCentreFeederNeis(m_elem.num_WatCentreFeederNeis());
    double entryPres= m_Pc_pistonTypeAdv;

    if(m_elem.iAmAPore() && m_cntAngAdv < PI/2.0 && num_WatCentreFeederNeis > 0)
    {
        double radSum(0.0);
        int iEvent(0);
        string poreBodyFillAlg(m_comn.poreFillAlg());
        for(int i = 0; i < m_elem.connectionNum(); ++i)
        {
            if(m_elem.connection(i)->model()->affectsNeiEntryPc(m_comn.oil()))
            {
                if(poreBodyFillAlg[0] == 'o' || poreBodyFillAlg[0] == 'O')
                {
                    radSum += m_comn.poreFillWeights(min(iEvent, 5))*
                        m_elem.connection(i)->model()->radius()*
                        double(rand())/double(RAND_MAX);
                }
                else
                {
                    radSum += m_comn.poreFillWeights(min(iEvent, 5))*
                        double(rand())/double(RAND_MAX);
                }
                ++iEvent;
            }
        }
        //softAssert(iEvent == m_elem.connectionNum()-num_WatCentreFeederNeis, "Failed on circle imb I Events");
        if(poreBodyFillAlg == "blunt2")
            entryPres = m_comn.oil().interfacialTen()*
                                  (2.0*cos(m_cntAngAdv)/m_R - radSum);
        else
            entryPres = 2.0*m_comn.oil().interfacialTen()*
                                 cos(m_cntAngAdv)/(m_R+radSum);
    }
    
	double maxNeiPistonEntryPrs=-1.0e26;
	for(int i = 0; i < m_elem.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( m_elem.connection(i)->model()->conductCWater())
			maxNeiPistonEntryPrs = max( maxNeiPistonEntryPrs,
				m_elem.connection(i)->model()->Pc_pistonTypeAdv());
	}
	if (maxNeiPistonEntryPrs<-1.0e25 && !conductAnyWater()) 
	{
		cout<<"'";cout.flush();
		//cout<<m_elem.connection(0)->model()[100000000].Pc_pistonTypeAdv();
	}
	//else
	entryPres = min(maxNeiPistonEntryPrs*0.999+0.001*entryPres,entryPres); ///. serial: the hardest of the pore and throat
	return entryPres;
	
}




/**
// Calculate the various entry pressures for imbibition, This function is used for all
// displacement cycles
*/
void Polygon::finitWaterInjection(double cappPrs)
{
	double tension =  m_comn.oil().interfacialTen();
    for(int i = 0; i < m_numCorners; ++i)
    {
        m_waterInCorner[i].CornerApex::finitCornerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, true, false);

        if(m_bulkFluid == &m_comn.water())
        {
            if (!m_oilLayer[i].LayerApex::finitLayerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, false, false))
			{
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" FBPD "<<endl;
									((double*)&m_R)[10000000]=0.0;

			}
			else if (m_oilLayer[i].exists())
			{
				softAssert(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
			}

        }
    }

}

/**
// Calculate the various entry pressures for imbibition, This function is used for all
// displacement cycles
*/
void Polygon::initWaterInjection(double cappPrs)
{
    softAssert(m_elem.connectionNum() > 0);
    m_elem.resetFillingEventRecord();
    m_elem.setInWatFloodVec(false);
    for(int j = 0; j < m_numCorners; ++j)
    {
        m_oilLayer[j].setInWatFloodVec(false);
    }

	double tension =  m_comn.oil().interfacialTen();
    for(int i = 0; i < m_numCorners; ++i)
    {
        m_waterInCorner[i].CornerApex::initCornerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, false);

        if(m_bulkFluid == &m_comn.water())
        {
            if (!m_oilLayer[i].LayerApex::initLayerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, false))
			{
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" ji ";
			};

        }
    }

    if(m_bulkFluid == &m_comn.water()) 
    {
        m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;     
        if (cappPrs < m_Pc_pistonTypeRec)
        {
			m_Pc_pistonTypeAdv = cappPrs * cos(m_cntAngAdv) / cos(m_cntAngRec); ///. Warning risk of division by zero 
		}
		return;
	}


    m_virginState = false;

    {
		//double recConAng(m_waterInCorner[0].cornerExists() ? m_waterInCorner[0].getCApex DistConAng(cappPrs, m_cntAngAdv, m_crnHafAngs[0], m_comn.oil().interfacialTen()) : m_minInitRecCntAng);
		double apexDist, recConAng(m_minInitRecCntAng); 
		m_waterInCorner[0].getCApexDistConAng(apexDist, recConAng, cappPrs, m_crnHafAngs[0], m_comn.oil().interfacialTen(),true);
		int WarningDontKnowWhatsGoingOn;

		double maxLocalPcLastCycle(m_comn.maxPcLastDrainCycle()-m_elem.gravityCorrection());
		double maxLocalPc(m_comn.maxEverPc()-m_elem.gravityCorrection());
		double max_Pc = m_waterInCorner[0].pinnedInInitState() ? maxLocalPc: maxLocalPcLastCycle;
		double angSum(0.0), normThresPress((m_R*max_Pc )/m_comn.oil().interfacialTen());//, curvRad(0.0);
		for(int i = 0; i < m_numCorners; ++i)
		{
			if(m_waterInCorner[i].cornerExists()) 
			{
				angSum += cos(recConAng + m_crnHafAngs[i]);
			}
		}

		double rhsMaxAdvConAng = (-4.0*m_shapeFactor*angSum)                      
			/ (normThresPress-cos(recConAng)+12.0*m_shapeFactor*sin(recConAng));
		rhsMaxAdvConAng = max(rhsMaxAdvConAng, -1.0);    // Prevent falling out of range [1, -1].  This is only applicable when r is very small
		rhsMaxAdvConAng = min(rhsMaxAdvConAng, 1.0);  // and these elems are probably not drained.
		m_maxConAngSpont = acos(rhsMaxAdvConAng); 
	}
 
    if(m_cntAngAdv < m_maxConAngSpont)
        m_Pc_pistonTypeAdv = Pc_pistonType_ImbHingCLine();
    else if(m_cntAngAdv <=  PI/2.0 + m_crnHafAngs[0])            // If oil layers are not possible in all
        m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;     // corners, we do not use MSP procedure
    else                                                                //  => do the check with max(beta)
        m_Pc_pistonTypeAdv = Pc_pistonType_Imbnww();

}


/**
// When the number of oil filled neighbours change the imbibition entry pressure
// will need to be updated. E.g. piston type displacement might now be possible
// or we move from an I3 to I2 pore body filling event. If water injection is
// forced, coopertive pore body filling will not occur.
*/
double Polygon::centreEntryPrsWatInj()
{

    int num_WatCentreFeederNeis(m_elem.num_WatCentreFeederNeis());
    double pistonEntryPrs;

    if(m_elem.iAmAPore() && m_cntAngAdv < PI/2.0 && num_WatCentreFeederNeis > 0 )
    {
        double radSum(0.0);
        int iEvent(0);
        string poreBodyFillAlg(m_comn.poreFillAlg());
        for(int i = 0; i < m_elem.connectionNum(); ++i)
        {
            if( m_elem.connection(i)->model()->affectsNeiEntryPc( m_comn.oil() ) && !m_elem.iRockType() )
            {
                if(poreBodyFillAlg[0] == 'o' || poreBodyFillAlg[0] == 'O')
                {
                    radSum += m_comn.poreFillWeights(min(iEvent, 5))*
                        m_elem.connection(i)->model()->radius()*
                        double(rand())/double(RAND_MAX);
                }
                else
                {
                    radSum += m_comn.poreFillWeights(min(iEvent, 5))*
                        double(rand())/double(RAND_MAX);
                }
                ++iEvent;
            }
        }

        if(poreBodyFillAlg == "blunt2")
            pistonEntryPrs = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)/m_R - radSum );
        else if(poreBodyFillAlg == "blunt1" || poreBodyFillAlg == "oren1")
            pistonEntryPrs = m_comn.oil().interfacialTen()*2.0*cos(m_cntAngAdv)/(m_R+radSum);
        else
            pistonEntryPrs = m_comn.oil().interfacialTen()*(1.0+2.0*sqrt(PI*m_shapeFactor))*cos(m_cntAngAdv)/(m_R+radSum);
    }
    else
		pistonEntryPrs = m_Pc_pistonTypeAdv;

	//double mySlfpistonEntryPrs = pistonEntryPrs;


	double maxNeiPistonEntryPrs=-1.0e26;
	for(int i = 0; i < m_elem.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( m_elem.connection(i)->model()->conductCWater())
			maxNeiPistonEntryPrs = max( maxNeiPistonEntryPrs,
				m_elem.connection(i)->model()->Pc_pistonTypeAdv());
	}
	if (maxNeiPistonEntryPrs>-1.0e25) ///. Not necessary, but anyway
		pistonEntryPrs = min(maxNeiPistonEntryPrs*0.999+0.001*pistonEntryPrs,pistonEntryPrs); ///. serial: the hardest of the pore and throat


	double snapOffPrs = m_waterInCorner[0].trappedCorner().first<0 ? calcSnapOffPressureImb() : -1.0e27;
    if(num_WatCentreFeederNeis > 0 && pistonEntryPrs > snapOffPrs)
    {
		m_displacementType = 'P';
		 if (m_elem.iAmAPore()) outD<<" ep"<<m_elem.index2p()<<":"<<pistonEntryPrs<<" ";
		 else
		 {
			outD<<" e"<<m_displacementType<<m_elem.index2p()<<":"<<pistonEntryPrs<<" ";
			outD<<" e"<<'P'<<m_elem.index2p()<<":"<<snapOffPrs<<" ";
		 }
       return pistonEntryPrs;
	}
    else
    {
        m_displacementType = 'S';
		 if (m_elem.iAmAPore()) outD<<" ep"<<m_elem.index2p()<<":"<<pistonEntryPrs<<" ";
		 else 
			{
				outD<<" e"<<m_displacementType<<m_elem.index2p()<<":"<<snapOffPrs<<" ";
				outD<<" e"<<'P'<<m_elem.index2p()<<":"<<pistonEntryPrs<<" ";
			}
        return snapOffPrs;
	}
  
}



/**
// The threshold pressures are based on the Mayer-Stowe-Princen method for calculating entry pressure.
// The contact angle used in the calculations is PI-contactAngle when the contactangle is > PI/2.
// Initally we assume that oil layers are stable in all coners. For the then given entry pressure we
// make sure that the collapsing pressure in the most oblique has not been reached. If it has we
// compare that pressure to the entry pressure for two layers. If the entry pressure for oil in two
// coners
*/
double Polygon::Pc_pistonType_Imbnww() const
{
    double contactAngle = PI - m_cntAngAdv;
    vector< double > potentialCurveRad;

    double sOne(0.0), sTwo(0.0), sThree(0.0);
    for(int i = 0; i < m_numCorners; ++i)
    {
        if(m_cntAngAdv > PI/2.0 + m_crnHafAngs[i])  // Only cases where oil layers/oil in corner might exist
        {
            sOne += cos(contactAngle)*cos(contactAngle+m_crnHafAngs[i])/sin(m_crnHafAngs[i])
                - (PI/2.0-contactAngle-m_crnHafAngs[i]);
            sTwo += cos(contactAngle+m_crnHafAngs[i])/sin(m_crnHafAngs[i]);
            sThree += 2.0 * (PI/2.0 - contactAngle - m_crnHafAngs[i]);

            double dFact = sOne - 2.0*sTwo*cos(contactAngle) + sThree;
            double rootFact = 1.0 + 4.0*m_shapeFactor*dFact/(cos(contactAngle)*cos(contactAngle));

            double radOne = m_R*cos(contactAngle)*(1.0-sqrt(rootFact))/(4.0*m_shapeFactor*dFact);
            double radTwo = m_R*cos(contactAngle)*(1.0+sqrt(rootFact))/(4.0*m_shapeFactor*dFact);
            potentialCurveRad.push_back(max(radOne, radTwo));
        }
    }


    for(int j = potentialCurveRad.size()-1; j >= 0; --j)
    {
        double tension(m_comn.oil().interfacialTen());
        double pc(tension / potentialCurveRad[j]);
        double layerPc(INF_NEG_NUM);
        if(m_waterInCorner[j].cornerExists()) 
        {///. Warning layers should be initialised for WatInj befor this
            layerPc = m_oilLayer[j].entryPc();//layerCollapsePc(pc, m_cntAngAdv, m_crnHafAngs[j], tension,false);
		}

        if(pc > layerPc) return   m_comn.oil().interfacialTen() / potentialCurveRad[j];
    }

    return  m_comn.oil().interfacialTen() * (2.0*cos(m_cntAngAdv)) / m_R;
}



/**
// Piston type displacement (also referred to as event 1 displacement) for advancing contact angles less
// than the critical is obtained by solving a set of non-linear equations using the Newton-Raphson tech.
// The hinging contact angle will vary between the receding and advancing, fixed at a position bi. This
// procedure is identical to that used by Patzek and Oren.
*/
double Polygon::Pc_pistonType_ImbHingCLine() const
{
    double tension(m_comn.oil().interfacialTen()), err(1000.0);
	double newPc(1.1*tension*(2.0*cos(m_cntAngAdv)) / m_R);
	
	double oldPc;

    int itr;
	double newPc2 = newPc;
    for(itr = 0; itr < MAX_NEWT_ITR+1; ++itr)
    {
        double sumOne(0.0), sumTwo(0.0), sumThree(0.0), sumFour(0.0);
		oldPc = newPc;
        for(int i = 0; i < m_numCorners; ++i)
        { 
            if(m_waterInCorner[i].cornerExists())
            { 
				
				double meniscusApexDist, hingConAng(m_cntAngAdv); 
				m_waterInCorner[i].getCApexDistConAng(meniscusApexDist, hingConAng, oldPc, m_crnHafAngs[i], m_comn.oil().interfacialTen(),true,true);

				double partus(meniscusApexDist * sin(m_crnHafAngs[i]) * oldPc/tension);
				softAssert(partus >= -1.0 && partus <=  1.0);
				if(!(partus >= -1.0 && partus <=  1.0))
				{
					cout<<partus<<" = "<<meniscusApexDist<<" * sin(" <<m_crnHafAngs[i]<<") / "<<tension/oldPc<<endl;
					cout<<"  oldPc "<<oldPc<<endl;
					cout<<"  m_maxConAngSpont "<<m_maxConAngSpont<<endl;
					cout<<"  tension/oldPc "<<tension/oldPc<<endl;
				}

                sumOne += meniscusApexDist*cos(hingConAng);
                sumTwo += PI/2.0-hingConAng-m_crnHafAngs[i];
                sumThree += asin(partus);
                sumFour += meniscusApexDist;
            }
        }


            double a= 2.0*sumThree - sumTwo;
            double b=cos(m_cntAngAdv)*(m_R/(2.0*m_shapeFactor)) -2.0*sumFour + sumOne;
            double c=-m_R*m_R/(4.0*m_shapeFactor);
            if (b*b-4.0*a*c>0) 
            {
				newPc2=tension*(2.0*a)/
				    ( (-b+sqrt(b*b-4.0*a*c)) );
			} 
            else
            {
				newPc2=tension*(2.0*a)/
				    ( (-b) );
		
			}
            

		newPc = newPc2;

        err = 2.0*fabs((newPc - oldPc)/(abs(oldPc)+abs(newPc)+1.0e-3));
        if(err < EPSILON)  break;

    }



		
    if(err < 0.0001)         
    {
		return newPc;
	}
	else if(err < 0.1)
    {
        cout<< "Problem in Polygon::mspCurveRadHingImb error:" <<err <<" %of " <<newPc <<"  teta " <<m_cntAngAdv << endl;
        return newPc;
    }


    {cerr << endl
        << "=================================================" << endl
        << "Error: Failed to obtain valid value for threshold" << endl
        << "radius of curvature in piston type displacement  " << endl
        << "during water injection."                           << endl
        << "Trapped water: " << m_waterInCorner[0].trappedCorner().first << endl
        << "Err  " << err << endl
        << "Con ang " << m_cntAngAdv*180.0/PI << endl
        << "Con ang " << m_cntAngAdv << endl
        << "Con ang " << m_cntAngRec << endl
        << "Radius " << m_R << endl
        << "oldPc " << oldPc << endl
        << "newPc " << newPc << endl
        << "G " << m_shapeFactor << endl
        << "Iteration " << itr << endl
        << "adPc " << m_waterInCorner[0].advancingPc() << endl
        << "recPc " << m_waterInCorner[0].receedingPc() << endl
        << "recPc " << m_waterInCorner[0].pinnedApexDist() << endl
        << "=================================================" << endl;   
					//((double*)&m_R)[10000000]=0.0;
		//exit(-1);
	}
    return 0.0;
}


///. TO PERMANENTLY DELETE
bool Polygon::waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const
{  ///. checks if we need to untrap oil, insertReCalcImbibeEntryPrs ...         never does anything
	if (!m_waterInCorner[0].cornerExists()) return false;
	if (m_waterInCorner[0].trappedCorner().first>-1 && m_elem.trappingWatFilm().first<0) 
	{
		cout<<" * Error: unsynced Trapping * "<<endl;
	}
	if (m_waterInCorner[0].trappedCorner().first>-1 || eleman()->isTrappedOil() || m_oilLayer->trappingCL().first>-1) return false;
    if (m_comn.injectant() == &m_comn.oil()) return false;
    if (conductCOil())
    {


		for (int i=0;i<m_numCorners;++i)
		{
			if (m_waterInCorner[i].cornerExists() && m_waterInCorner[i].trappedCorner().first<0 && !m_waterInCorner[i].pinned()) 
			{
				m_waterInCorner[i].initCornerApex(cappPrs,m_cntAngRec,m_cntAngAdv,m_crnHafAngs[i],m_comn.oil().interfacialTen(),false);
				cout<<"j";
			}
		}
		
		double snapPc = calcSnapOffPressureImb();

		return  cappPrs < snapPc;
    
	}
    
	int WARninG;
	return false;
}




///. Warning also called during Oil Injection, when decreasing oil pressure due to coalescence 
double Polygon::Pc_pin_disconnectOilLayer(int cor)  ///rare
{
    softAssert(m_oilConnection && m_numLayers>0);    
	LayerApex& oilLayer = m_oilLayer[cor];
    double layerPc(oilLayer.layerCollPc());
	double tension =  m_comn.oil().interfacialTen();
    if (m_waterInCorner[cor].cornerExists())
    {
		
		m_waterInCorner[cor].finitCornerApex(layerPc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[cor],tension, false, false);
		m_waterInCorner[cor].initCornerApex(layerPc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[cor],tension, false);

		m_hasDisConectedCentreWCornerW = false;
		bool Warn;
	}

    oilLayer.removeLayer();
    --m_numLayers;
	softAssert(m_numLayers >= 0 && m_numLayers < m_numCorners);

    if(!m_numLayers)
    {
        softAssert(!m_hasDisConectedCentreWCornerW);
        m_oilConnection = false;
        if (!containCOil() && m_elem.isTrappedOil()) {cout<<" axsj ";m_elem.unTrapOil();}
    }
    return layerPc;
}






