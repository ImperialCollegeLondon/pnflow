#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

void Circle::finitOilInjection(double cappPrs)
{
	//softAssert(m_elem.connectionNum() > 0, "20");
	//softAssert(m_elem.connectionNum() > 0, "20");
	m_elem.resetFillingEventRecord();
	//m_elem.setInWatFloodVec(false);
	m_elem.setInOilFloodVec(false);
 
	//if(m_bulkFluid == &m_comn.oil()) return;
 //
	//m _maxConAngSpont = PI/2.0;
	//m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R ;
  //
	//calcCentreEntryPrsOilInj();
	//timestep = -1;

} 
 
void Circle::initOilInjection(double cappPrs)
{
	softAssert(m_elem.connectionNum() > 0);
	//m_elem.resetFillingEventRecord();
	//m_elem.setInWatFloodVec(false);
	 softAssert(!m_elem.isInOilFloodVec());

	if(m_bulkFluid == &m_comn.oil()) return;
 
	//m _maxConAngSpont = PI/2.0;
	m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R ;
  
	//calcCentreEntryPrsOilInj();
} 

///. move to pore/throat classes (?)
double Circle::centreEntryPrsOilInj()
{ 
	int num_OilCentreFeederNeis(m_elem.numOilCentreFeederNeis());
	double conAng(m_cntAngRec);
	double entryPres = m_Pc_pistonTypeRec;
 
	if(m_elem.iAmAPore() && conAng > PI/2.0 && num_OilCentreFeederNeis != 0)
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(m_comn.poreFillAlg());
		for(int i = 0; i < m_elem.connectionNum(); ++i)
		{
			if(m_elem.connection(i)->model()->affectsNeiEntryPc(m_comn.water()))
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
		softAssert(iEvent == m_elem.connectionNum()-num_OilCentreFeederNeis);
		if(poreBodyFillAlg == "blunt2")
			entryPres = m_comn.oil().interfacialTen()*(2.0*cos(conAng)/m_R - radSum);
		else
			entryPres = 2.0*m_comn.oil().interfacialTen()*cos(conAng)/(m_R+radSum);
	}
	
	double minNeiPistonEntryPrs(1.0e26);
	for(int i = 0; i < m_elem.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( m_elem.connection(i)->model()->conductCOil())
			minNeiPistonEntryPrs = min( minNeiPistonEntryPrs,
				m_elem.connection(i)->model()->Pc_pistonTypeRec());
	}
	//if (minNeiPistonEntryPrs<-1.0e35) cout<<"Rc";
	//else
	if (minNeiPistonEntryPrs<1.0e25) ///. Not necessary, but anyway
		entryPres = max(minNeiPistonEntryPrs*0.999+0.001*entryPres,entryPres); ///. serial: the hardest of the pore and throat
	return entryPres;
}




/**
// Calculate the various entry pressures for drainage, This function is used for all
// displacement cycles
*/
void Polygon::finitOilInjection(double cappPrs)
{
	double tension =  m_comn.oil().interfacialTen();
	m_elem.setInOilFloodVec(false);

	for(int i = 0; i < m_numCorners; ++i)
	{
		m_waterInCorner[i].CornerApex::finitCornerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, false, false);

		if(m_bulkFluid == &m_comn.water() && m_oilLayer[i].exists() &&  m_oilLayer[i].trappedOLayer().first < 0 )
		{
			if ( !m_oilLayer[i].LayerApex::finitLayerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, true, false) )
			{
				Pc_pin_disconnectOilLayer(i);
				cout<<" KFXO  ";
			}
			else
			{
				softAssert(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
			}

		}
	}

}


/**
// Calculate the various entry pressures for drainage, This function is used for all
// displacement cycles
*/
void Polygon::initOilInjection(double cappPrs)
{
	softAssert(m_elem.connectionNum() > 0);
	m_elem.resetFillingEventRecord();
	//m_elem.setInWatFloodVec(false);
	 softAssert(!m_elem.isInOilFloodVec());

	for(int j = 0; j < m_numCorners; ++j)
	{
		m_oilLayer[j].setInWatFloodVec(false);
		m_oilLayer[j].setInOilFloodVec(false);
	}

	double tension =  m_comn.oil().interfacialTen();
	for(int i = 0; i < m_numCorners; ++i)
	{
		m_waterInCorner[i].CornerApex::initCornerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, true);

		if(m_bulkFluid == &m_comn.water())
		{
			if ( !m_oilLayer[i].LayerApex::initLayerApex(cappPrs, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], tension, true) )
			{
				Pc_pin_disconnectOilLayer(i);
				cout<<" jd ";
			}
		}
	}

	double conAng(m_cntAngRec);

	if(m_bulkFluid == &m_comn.oil()) 
	{
		m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R;	 
		if (cappPrs > m_Pc_pistonTypeAdv)
		{
			m_Pc_pistonTypeRec = cappPrs * cos(m_cntAngRec) / cos(m_cntAngRec); ///. Warning risk of division by zero 
		}
		return;
	}


 
	double minLocalPcLastCycle(m_comn.minPcLastImbCycle()-m_elem.gravityCorrection());
	double minLocalPc(m_comn.minEverCappPress()-m_elem.gravityCorrection());
	double min_Pc = m_oilLayer[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc;
	double angSum(0.0), normThresPress((m_R*min_Pc )/m_comn.oil().interfacialTen());
	for(int i = 0; i < m_numCorners; ++i)
	{
		if(m_oilLayer[i].exists(/*st ab le*/)) angSum += cos(m_cntAngAdv - m_crnHafAngs[i]);
	}

	double rhsMaxRecConAng = (4.0*m_shapeFactor*angSum)								 
		/ (normThresPress-cos(m_cntAngAdv)-12.0*m_shapeFactor*sin(m_cntAngAdv));
	rhsMaxRecConAng = min(max(rhsMaxRecConAng, -1.0), 1.0);  //Prevent falling out of range [1, -1]. This is only applicable when r is very small  and these elems are probably not drained.
	m_maxConAngSpont = acos(rhsMaxRecConAng);
 
	if(conAng > m_maxConAngSpont && m_oilLayer[0].exists())
		m_Pc_pistonTypeRec = Pc_pistonType_DrainHing();
	else if(conAng >= PI/2.0 - m_crnHafAngs[0])
		m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(conAng)) / m_R;
	else
		m_Pc_pistonTypeRec = Pc_pistonType_Drain(conAng);

}


double Polygon::centreEntryPrsOilInj()
{
	int num_OilCentreFeederNeis(m_elem.numOilCentreFeederNeis());
	double conAng(m_cntAngRec);
	double pistonEntryPrs;

	if(m_elem.iAmAPore() && conAng > PI/2.0 && num_OilCentreFeederNeis != 0)
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(m_comn.poreFillAlg());
		for(int i = 0; i < m_elem.connectionNum(); ++i)
		{
			if( m_elem.connection(i)->model()->affectsNeiEntryPc( m_comn.water() ) && !m_elem.iRockType() )
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

		softAssert(iEvent == m_elem.connectionNum()-num_OilCentreFeederNeis);
		if(poreBodyFillAlg == "blunt2")
			pistonEntryPrs = m_comn.oil().interfacialTen()*(2.0*cos(conAng)/m_R - radSum);
		else if(poreBodyFillAlg == "blunt1" || poreBodyFillAlg == "oren1")
			pistonEntryPrs = m_comn.oil().interfacialTen()*2.0*cos(conAng)/(m_R+radSum);
		else
			pistonEntryPrs = m_comn.oil().interfacialTen()*(1.0+2.0*sqrt(PI*m_shapeFactor))*cos(conAng)/(m_R+radSum);
	}
	else
		pistonEntryPrs = m_Pc_pistonTypeRec;

	double mySlfpistonEntryPrs = pistonEntryPrs;

	double minNeiPistonEntryPrs(1.0e26);
	for(int i = 0; i < m_elem.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( m_elem.connection(i)->model()->conductCOil())
			minNeiPistonEntryPrs = min( minNeiPistonEntryPrs,
				m_elem.connection(i)->model()->Pc_pistonTypeRec());
	}
	if (minNeiPistonEntryPrs<1.0e25) ///. Not necessary, but anyway
	pistonEntryPrs = max(minNeiPistonEntryPrs*0.999+0.001*pistonEntryPrs,pistonEntryPrs); ///. serial: the hardest of the pore and throat

	double entryPres = pistonEntryPrs;
	m_displacementType = 'P';

	if(m_numLayers && m_oilLayer[0].trappedOLayer().first < 0)
	{
	  double snapOffPrs = calcSnapOffPressureDrain();
	  
	  if(	m_oilLayer[0].freeAtPrs(snapOffPrs)
		 &&( num_OilCentreFeederNeis == 0 || pistonEntryPrs > snapOffPrs )
		)
	  {
		 entryPres = snapOffPrs;
		 m_displacementType = 'S';
		 if (m_elem.iAmAPore()) outD<<" ep"<<m_elem.index2p()<<":"<<mySlfpistonEntryPrs<<" ";
		 else outD<<" e"<<m_displacementType<<m_elem.index2p()<<":"<<snapOffPrs<<" ";
		}
		else
		{
			 if (m_elem.iAmAPore()) outD<<" ep"<<m_elem.index2p()<<":"<<mySlfpistonEntryPrs<<" ";
			 else outD<<" e"<<m_displacementType<<m_elem.index2p()<<":"<<mySlfpistonEntryPrs<<" ";
		}
	}
	else
	{
		 if (m_elem.iAmAPore()) outD<<" ep"<<m_elem.index2p()<<":"<<mySlfpistonEntryPrs<<" ";
		 else outD<<" e"<<m_displacementType<<m_elem.index2p()<<":"<<mySlfpistonEntryPrs<<" ";
	}
	return entryPres;
  
}



/**
// The threshold pressures are based on the Mayer-Stowe-Princen method for calculating entry pressure.
// The contact angle used in the calculations is PI-contactAngle when the contactangle is > PI/2.
// Initally we assume that oil layers are stable in all coners. For the then given entry pressure we
// make sure that the collapsing pressure in the most oblique has not been reached. If it has we
// compare that pressure to the entry pressure for two layers. If the entry pressure for oil in two
// coners
*/
double Polygon::Pc_pistonType_Drain(double conAng) const
{
	double ACornsDL(0.0);

	for(int i = 0; i < m_numCorners; ++i)
	{
	   if(conAng < PI/2.0 - m_crnHafAngs[i])
		  ACornsDL += cos(conAng) * cos(conAng+m_crnHafAngs[i])/sin(m_crnHafAngs[i]) - (PI/2.0 - conAng - m_crnHafAngs[i]);
	}

	double funcF = (1.0 + sqrt(1.0 - 4.0*m_shapeFactor*ACornsDL/(cos(conAng)*cos(conAng)))) / (1.0 + 2.0*sqrt(PI * m_shapeFactor));

	return m_comn.oil().interfacialTen()*((1.0 + 2.0*sqrt(PI*m_shapeFactor)) * cos(conAng) * funcF) / m_R;
}


double Polygon::Pc_pistonType_DrainHing() const
{
	double tension(m_comn.oil().interfacialTen());
	double minLocalPcLastCycle(m_comn.minPcLastImbCycle()-m_elem.gravityCorrection());
	double minLocalPc(m_comn.minEverCappPress()-m_elem.gravityCorrection());

	double err, oldRad( tension / ( m_oilLayer[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc ) );
	int itr(0);
 
	if(m_oilLayer[0].trappedOLayer().first>-1)
		oldRad = tension/m_oilLayer[0].trappedOLayer().second;
 
	int numCorn(m_numCorners);
	while(numCorn >= 0) 
	{
  
		for(itr = 0; itr < MAX_NEWT_ITR; ++itr)
		{
			double sumOne(0.0), sumTwo(0.0), sumThree(0.0), sumFour(0.0);
 
			for(int j = 0; j < numCorn; ++j)
			{
				if(m_oilLayer[j].freeAtPrs(tension/oldRad) && !m_oilLayer[j].forcedSnapOff(tension/oldRad))
				{ 
					double hingConAng(m_cntAngRec);// = m_oilLayer[j].hingingConAngUntraped(tension/oldRad, m_cntAngRec, m_crnHafAngs[j], tension, true);
					double meniscusApexDist;// = m_oilLayer[j].getApexDistanceUntraped(tension/oldRad, hingConAng, m_crnHafAngs[j], tension, true);
					m_oilLayer[j].getCAApexDistUntraped(meniscusApexDist,hingConAng, m_crnHafAngs[j], tension/oldRad, tension, true);
							softAssert(meniscusApexDist>0);

					double partus(-meniscusApexDist * sin(m_crnHafAngs[j])/oldRad);
					softAssert(partus >= -1.0 && partus <=  1.0);
  
					sumOne += meniscusApexDist*cos(hingConAng);
					sumTwo += hingConAng-m_crnHafAngs[j]-PI/2.0;
					sumThree += asin(partus);
					sumFour += meniscusApexDist;
				}
			} 

			double newRad = (m_R*m_R/(4.0*m_shapeFactor) - oldRad*sumOne + oldRad*oldRad*sumTwo)
				/ (2.0*oldRad*sumThree + cos(m_cntAngRec)*(m_R/(2.0*m_shapeFactor) - 2.0*sumFour));

			err = fabs((newRad - oldRad)/oldRad);
			if(err < EPSILON) return  m_comn.oil().interfacialTen() / newRad;

			oldRad = newRad;
		}
		--numCorn;
	}

	cerr << "\n Error: failed to obtain valid value for threshold radius of curvature" 
		 << "\n   in piston type displacement during drainage."   << "\n   Err  " << err << "		Con ang " << m_cntAngRec*180.0/PI 
		 << "\n   Radius " << m_R << "   G " << m_shapeFactor  << "\n   Iteration " << itr << endl ;	exit(-1);


	return 0.0;
}




bool Polygon::hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const
{ ///. checks if we need to untrap water, reCalcDrainEntryPrs ...
	bool oilFlood(m_comn.injectant() == &m_comn.oil());
  
	return oilFlood && m_bulkFluid == &m_comn.water() && m_numLayers>0 && !m_oilLayer[0].trappedOLayer().first /*&& oilTrp.second > m_snapOffPrs*/ && cappPrs > calcSnapOffPressureDrain();
}




bool Polygon::Pc_growStableOilLayerDrain_UseLess(double Pc, int corner)
{///. oil film growth
	LayerApex& oilLayer = m_oilLayer[corner];
	softAssert(!containCOil());
	softAssert(!oilLayer.exists());
	

	

	
	if (m_oilLayer[corner].createOLayer(Pc,m_cntAngRec,m_cntAngAdv, m_maxConAngSpont, m_crnHafAngs[corner],  m_comn.oil().interfacialTen(),true))
	{
		softAssert(m_oilLayer[corner].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
		
		
		
		m_oilConnection = true;
		++m_numLayers;
		softAssert(m_numLayers > 0 && m_numLayers <=  m_numCorners);
		
		m_hasDisConectedCentreWCornerW = m_numLayers == m_numCorners;
	}
	return m_oilLayer[corner].exists();
}






