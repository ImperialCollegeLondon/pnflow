#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"



using namespace std;

void Circle::finitOilInjection(double cappPrs)
{
	//ensure(elem_.connectionNum() > 0, "20");
	//ensure(elem_.connectionNum() > 0, "20");
	elem_.resetEventIndex();
	//elem_.setInWatFloodVec(false);
	elem_.setInOilFloodVec(false);
 
	//if(bulkFluid_ == &comn_.oil()) return;
 //
	//m _maxConAngSpont = PI/2.0;
	//Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_ ;
  //
	//calcCentreEntryPrsOilInj();
	//timestep = -1;

} 
 
void Circle::initOilInjection(double cappPrs)
{
	ensure(elem_.connectionNum() > 0);
	//elem_.resetEventIndex();
	//elem_.setInWatFloodVec(false);
	 ensure(!elem_.isInOilFloodVec());

	if(bulkFluid_ == &comn_.oil()) return;
 
	//m _maxConAngSpont = PI/2.0;
	Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_ ;
  
	//calcCentreEntryPrsOilInj();
} 

///. move to pore/throat classes (?)
double Circle::centreEntryPrsOilInj()
{ 
	int num_OilCentreFeederNeis(elem_.numOilCentreFeederNeis());
	double conAng(cntAngRec_);
	double entryPres = Pc__pistonTypeRec;
 
	if(elem_.iAmAPore() && conAng > PI/2.0 && num_OilCentreFeederNeis != 0)
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(comn_.poreFillAlg());
		for(int i = 0; i < elem_.connectionNum(); ++i)
		{
			if(elem_.connection(i)->model()->affectsNeiEntryPc(comn_.water()))
			{
				if(poreBodyFillAlg[0] == 'o' || poreBodyFillAlg[0] == 'O')
				{
					radSum += comn_.poreFillWeights(min(iEvent, 5))*
						elem_.connection(i)->model()->RRR()*
						double(rand())/double(RAND_MAX);
				}
				else
				{
					radSum += comn_.poreFillWeights(min(iEvent, 5))*
						double(rand())/double(RAND_MAX);
				}
				++iEvent;
			}
		}
		ensure(iEvent == elem_.connectionNum()-num_OilCentreFeederNeis);
		if(poreBodyFillAlg == "blunt2")
			entryPres = comn_.sigmaOW()*(2.0*cos(conAng)/R_ - radSum);
		else
			entryPres = 2.0*comn_.sigmaOW()*cos(conAng)/(R_+radSum);
	}

	double minNeiPistonEntryPrs(1.0e26);
	for(int i = 0; i < elem_.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( elem_.connection(i)->model()->conductCOil())
			minNeiPistonEntryPrs = min( minNeiPistonEntryPrs,
				elem_.connection(i)->model()->Pc_pistonTypeRec());
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
	double tension =  comn_.sigmaOW();
	elem_.setInOilFloodVec(false);

	for(int i = 0; i < numCorners_; ++i)
	{
		waterInCorner_[i].CornerApex::finitCornerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, false, false);

		if(bulkFluid_ == &comn_.water() && oilLayer_[i].exists() &&  oilLayer_[i].trappedOLayer().first < 0 )
		{
			if ( !oilLayer_[i].LayerApex::finitLayerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, true, false) )
			{
				Pc_pin_disconnectOilLayer(i);
				cout<<" KFXO  ";
			}
			else
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
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
	ensure(elem_.connectionNum() > 0);
	elem_.resetEventIndex();
	//elem_.setInWatFloodVec(false);
	 ensure(!elem_.isInOilFloodVec());

	for(int j = 0; j < numCorners_; ++j)
	{
		oilLayer_[j].setInWatFloodVec(false);
		oilLayer_[j].setInOilFloodVec(false);
	}

	double tension =  comn_.sigmaOW();
	for(int i = 0; i < numCorners_; ++i)
	{
		waterInCorner_[i].CornerApex::initCornerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, true);

		if(bulkFluid_ == &comn_.water())
		{
			if ( !oilLayer_[i].LayerApex::initLayerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, true) )
			{
				Pc_pin_disconnectOilLayer(i);
				cout<<" jd ";
			}
		}
	}

	double conAng(cntAngRec_);

	if(bulkFluid_ == &comn_.oil()) 
	{
		Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_;
		if (cappPrs > Pc__pistonTypeAdv)
		{
			Pc__pistonTypeRec = cappPrs * cos(cntAngRec_) / cos(cntAngRec_); ///. Warning risk of division by zero 
		}
		return;
	}


 
	double minLocalPcLastCycle(comn_.minPcLastImbCycle()-elem_.gravityCorrection());
	double minLocalPc(comn_.minEverCappPress()-elem_.gravityCorrection());
	double min_Pc = oilLayer_[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc;
	double angSum(0.0), normThresPress((R_*min_Pc )/comn_.sigmaOW());
	for(int i = 0; i < numCorners_; ++i)
	{
		if(oilLayer_[i].exists(/*st ab le*/)) angSum += cos(cntAngAdv_ - crnHafAngs_[i]);
	}

	double rhsMaxRecConAng = (4.0*shapeFactor_*angSum)
		/ (normThresPress-cos(cntAngAdv_)-12.0*shapeFactor_*sin(cntAngAdv_));
	rhsMaxRecConAng = min(max(rhsMaxRecConAng, -1.0), 1.0);  //Prevent falling out of range [1, -1]. This is only applicable when r is very small  and these elems are probably not drained.
	maxConAngSpont_ = acos(rhsMaxRecConAng);
 
	if(conAng > maxConAngSpont_ && oilLayer_[0].exists())
		Pc__pistonTypeRec = Pc_pistonType_DrainHing();
	else if(conAng >= PI/2.0 - crnHafAngs_[0])
		Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(conAng)) / R_;
	else
		Pc__pistonTypeRec = Pc_pistonType_Drain(conAng);

}


double Polygon::centreEntryPrsOilInj()
{
	int num_OilCentreFeederNeis(elem_.numOilCentreFeederNeis());
	double conAng(cntAngRec_);
	double pistonEntryPrs;

	if(elem_.iAmAPore() && conAng > PI/2.0 && num_OilCentreFeederNeis != 0)
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(comn_.poreFillAlg());
		for(int i = 0; i < elem_.connectionNum(); ++i)
		{
			if( elem_.connection(i)->model()->affectsNeiEntryPc( comn_.water() ) && !elem_.rockIndex() )
			{
				if(poreBodyFillAlg[0] == 'o' || poreBodyFillAlg[0] == 'O')
				{
					radSum += comn_.poreFillWeights(min(iEvent, 5))*
						elem_.connection(i)->model()->RRR()*
						double(rand())/double(RAND_MAX);
				}
				else
				{
					radSum += comn_.poreFillWeights(min(iEvent, 5))*
						double(rand())/double(RAND_MAX);
				}
				++iEvent;
			}
		}

		ensure(iEvent == elem_.connectionNum()-num_OilCentreFeederNeis);
		if(poreBodyFillAlg == "blunt2")
			pistonEntryPrs = comn_.sigmaOW()*(2.0*cos(conAng)/R_ - radSum);
		else if(poreBodyFillAlg == "blunt1" || poreBodyFillAlg == "oren1")
			pistonEntryPrs = comn_.sigmaOW()*2.0*cos(conAng)/(R_+radSum);
		else
			pistonEntryPrs = comn_.sigmaOW()*(1.0+2.0*sqrt(PI*shapeFactor_))*cos(conAng)/(R_+radSum);
	}
	else
		pistonEntryPrs = Pc__pistonTypeRec;

	double mySlfpistonEntryPrs = pistonEntryPrs;

	double minNeiPistonEntryPrs(1.0e26);
	for(int i = 0; i < elem_.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( elem_.connection(i)->model()->conductCOil())
			minNeiPistonEntryPrs = min( minNeiPistonEntryPrs,
				elem_.connection(i)->model()->Pc_pistonTypeRec());
	}
	if (minNeiPistonEntryPrs<1.0e25) ///. Not necessary, but anyway
	pistonEntryPrs = max(minNeiPistonEntryPrs*0.999+0.001*pistonEntryPrs,pistonEntryPrs); ///. serial: the hardest of the pore and throat

	double entryPres = pistonEntryPrs;
	displacementType_ = 'P';

	if(numLayers_ && oilLayer_[0].trappedOLayer().first < 0)
	{
	  double snapOffPrs = calcSnapOffPressureDrain();
	  
	  if(    oilLayer_[0].freeAtPrs(snapOffPrs)
		 &&( num_OilCentreFeederNeis == 0 || pistonEntryPrs > snapOffPrs )
		)
	  {
		 entryPres = snapOffPrs;
		 displacementType_ = 'S';
		 if (elem_.iAmAPore()) outD<<" ep"<<elem_.index()<<":"<<mySlfpistonEntryPrs<<" ";
		 else outD<<" e"<<displacementType_<<elem_.index()<<":"<<snapOffPrs<<" ";
		}
		else
		{
			 if (elem_.iAmAPore()) outD<<" ep"<<elem_.index()<<":"<<mySlfpistonEntryPrs<<" ";
			 else outD<<" e"<<displacementType_<<elem_.index()<<":"<<mySlfpistonEntryPrs<<" ";
		}
	}
	else
	{
		 if (elem_.iAmAPore()) outD<<" ep"<<elem_.index()<<":"<<mySlfpistonEntryPrs<<" ";
		 else outD<<" e"<<displacementType_<<elem_.index()<<":"<<mySlfpistonEntryPrs<<" ";
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

	for(int i = 0; i < numCorners_; ++i)
	{
	   if(conAng < PI/2.0 - crnHafAngs_[i])
		  ACornsDL += cos(conAng) * cos(conAng+crnHafAngs_[i])/sin(crnHafAngs_[i]) - (PI/2.0 - conAng - crnHafAngs_[i]);
	}

	double funcF = (1.0 + sqrt(1.0 - 4.0*shapeFactor_*ACornsDL/(cos(conAng)*cos(conAng)))) / (1.0 + 2.0*sqrt(PI * shapeFactor_));

	return comn_.sigmaOW()*((1.0 + 2.0*sqrt(PI*shapeFactor_)) * cos(conAng) * funcF) / R_;
}


double Polygon::Pc_pistonType_DrainHing() const
{
	double tension(comn_.sigmaOW());
	double minLocalPcLastCycle(comn_.minPcLastImbCycle()-elem_.gravityCorrection());
	double minLocalPc(comn_.minEverCappPress()-elem_.gravityCorrection());

	double err, oldRad( tension / ( oilLayer_[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc ) );
	int itr(0);
 
	if(oilLayer_[0].trappedOLayer().first>-1)
		oldRad = tension/oilLayer_[0].trappedOLayer().second;
 
	int numCorn(numCorners_);
	while(numCorn >= 0) 
	{
  
		for(itr = 0; itr < MAX_NEWT_ITR; ++itr)
		{
			double sumOne(0.0), sumTwo(0.0), sumThree(0.0), sumFour(0.0);
 
			for(int j = 0; j < numCorn; ++j)
			{
				if(oilLayer_[j].freeAtPrs(tension/oldRad) && !oilLayer_[j].forcedSnapOff(tension/oldRad))
				{ 
					double hingConAng(cntAngRec_);// = oilLayer_[j].hingingConAngUntraped(tension/oldRad, cntAngRec_, crnHafAngs_[j], tension, true);
					double meniscusApexDist;// = oilLayer_[j].getApexDistanceUntraped(tension/oldRad, hingConAng, crnHafAngs_[j], tension, true);
					oilLayer_[j].getCAApexDistUntraped(meniscusApexDist,hingConAng, crnHafAngs_[j], tension/oldRad, tension, true);
							ensure(meniscusApexDist>0);

					double partus(-meniscusApexDist * sin(crnHafAngs_[j])/oldRad);
					ensure(partus >= -1.0 && partus <=  1.0);
  
					sumOne += meniscusApexDist*cos(hingConAng);
					sumTwo += hingConAng-crnHafAngs_[j]-PI/2.0;
					sumThree += asin(partus);
					sumFour += meniscusApexDist;
				}
			} 

			double newRad = (R_*R_/(4.0*shapeFactor_) - oldRad*sumOne + oldRad*oldRad*sumTwo)
				/ (2.0*oldRad*sumThree + cos(cntAngRec_)*(R_/(2.0*shapeFactor_) - 2.0*sumFour));

			err = fabs((newRad - oldRad)/oldRad);
			if(err < EPSILON) return  comn_.sigmaOW() / newRad;

			oldRad = newRad;
		}
		--numCorn;
	}

	cerr << "\n Error: failed to obtain valid value for threshold radius of curvature" 
		 << "\n   in piston type displacement during drainage."   << "\n   Err  " << err << "        Con ang " << cntAngRec_*180.0/PI 
		 << "\n   Radius " << R_ << "   G " << shapeFactor_  << "\n   Iteration " << itr << endl ;    exit(-1);


	return 0.0;
}




bool Polygon::hasOilLayer_TrappedOutside_PcHsnapPc(double cappPrs) const
{ ///. checks if we need to untrap water, reCalcDrainEntryPrs ...
	bool oilFlood(comn_.injectant() == &comn_.oil());
  
	return oilFlood && bulkFluid_ == &comn_.water() && numLayers_>0 && !oilLayer_[0].trappedOLayer().first /*&& oilTrp.second > snapOffPrs_*/ && cappPrs > calcSnapOffPressureDrain();
}




bool Polygon::Pc_growStableOilLayerDrain_UseLess(double Pc, int corner)
{///. oil film growth
	LayerApex& oilLayer = oilLayer_[corner];
	ensure(!containCOil());
	ensure(!oilLayer.exists());





	if (oilLayer_[corner].createOLayer(Pc,cntAngRec_,cntAngAdv_, maxConAngSpont_, crnHafAngs_[corner],  comn_.sigmaOW(),true))
	{
		ensure(oilLayer_[corner].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));



		oilConnection_ = true;  ++numLayers_;
		
		ensure(numLayers_ > 0 && numLayers_ <=  numCorners_);

		hasDisConectedCentreWCornerW_ = numLayers_ == numCorners_;
	}
	return oilLayer_[corner].exists();
}






