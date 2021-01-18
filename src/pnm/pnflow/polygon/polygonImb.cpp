#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"
using namespace std;


void Circle::finitWaterInjection(double cappPrs)
{
	//ensure(elem_.connectionNum() > 0, "23");
	//elem_.resetEventIndex();
	//elem_.setInWatFloodVec(false);
	//elem_.setInOilFloodVec(false);
 
	//if(bulkFluid_ == &comn_.water()) return;
 //
	//m _maxConAngSpont = PI/2.0;
	//m _Pc_pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;
  
	//c alcCentreEntryPrsWatInj();
} 

void Circle::initWaterInjection(double cappPrs)
{
	ensure(elem_.connectionNum() > 0);
	elem_.resetEventIndex();
	elem_.setInWatFloodVec(false);
	//elem_.setInOilFloodVec(false);
 
	if(bulkFluid_ == &comn_.water()) return;
 
	//m _maxConAngSpont = PI/2.0;
	Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;
  	//timestep = -1;

	//calcCentreEntryPrsWatInj();
} 
 
///. move to pore/throat classes (?)
double Circle::centreEntryPrsWatInj()
{ 
	int num_WatCentreFeederNeis(elem_.num_WatCentreFeederNeis());
	double entryPres= Pc__pistonTypeAdv;

	if(elem_.iAmAPore() && cntAngAdv_ < PI/2.0 && num_WatCentreFeederNeis > 0)
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(comn_.poreFillAlg());
		for(int i = 0; i < elem_.connectionNum(); ++i)
		{
			if(elem_.connection(i)->model()->affectsNeiEntryPc(comn_.oil()))
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
		//ensure(iEvent == elem_.connectionNum()-num_WatCentreFeederNeis, "Failed on circle imb I Events");
		if(poreBodyFillAlg == "blunt2")
			entryPres = comn_.sigmaOW()*
								  (2.0*cos(cntAngAdv_)/R_ - radSum);
		else
			entryPres = 2.0*comn_.sigmaOW()*
								 cos(cntAngAdv_)/(R_+radSum);
	}

	double maxNeiPistonEntryPrs=-1.0e26;
	for(int i = 0; i < elem_.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( elem_.connection(i)->model()->conductCWater())
			maxNeiPistonEntryPrs = max( maxNeiPistonEntryPrs,
				elem_.connection(i)->model()->Pc_pistonTypeAdv());
	}
	if (maxNeiPistonEntryPrs<-1.0e25 && !conductAnyWater()) 
	{
		cout<<"'";cout.flush();
		//cout<<elem_.connection(0)->model()[100000000].Pc_pistonTypeAdv();
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
	double tension =  comn_.sigmaOW();
	for(int i = 0; i < numCorners_; ++i)
	{
		waterInCorner_[i].CornerApex::finitCornerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, true, false);

		if(bulkFluid_ == &comn_.water())
		{
			if (!oilLayer_[i].LayerApex::finitLayerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, false, false))
			{
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" FBPD "<<endl;
									((double*)&R_)[10000000]=0.0;

			}
			else if (oilLayer_[i].exists())
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
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
	ensure(elem_.connectionNum() > 0);
	elem_.resetEventIndex();
	elem_.setInWatFloodVec(false);
	for(int j = 0; j < numCorners_; ++j)
	{
		oilLayer_[j].setInWatFloodVec(false);
	}

	double tension =  comn_.sigmaOW();
	for(int i = 0; i < numCorners_; ++i)
	{
		waterInCorner_[i].CornerApex::initCornerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, false);

		if(bulkFluid_ == &comn_.water())
		{
			if (!oilLayer_[i].LayerApex::initLayerApex(cappPrs, cntAngRec_, cntAngAdv_, crnHafAngs_[i], tension, false))
			{
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" ji ";
			};

		}
	}

	if(bulkFluid_ == &comn_.water()) 
	{
		Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;     
		if (cappPrs < Pc__pistonTypeRec)
		{
			Pc__pistonTypeAdv = cappPrs * cos(cntAngAdv_) / cos(cntAngRec_); ///. Warning risk of division by zero 
		}
		return;
	}


	virginState_ = false;

	{
		//double recConAng(waterInCorner_[0].cornerExists() ? waterInCorner_[0].getCApex DistConAng(cappPrs, cntAngAdv_, crnHafAngs_[0], comn_.sigmaOW()) : minInitRecCntAng_);
		double apexDist, recConAng(minInitRecCntAng_); 
		waterInCorner_[0].getCApexDistConAng(apexDist, recConAng, cappPrs, crnHafAngs_[0], comn_.sigmaOW(),true);
		//*************  WarningDontKnowWhatsGoingOn; **********

		double maxLocalPcLastCycle(comn_.maxPcLastDrainCycle()-elem_.gravityCorrection());
		double maxLocalPc(comn_.maxEverPc()-elem_.gravityCorrection());
		double max_Pc = waterInCorner_[0].pinnedInInitState() ? maxLocalPc: maxLocalPcLastCycle;
		double angSum(0.0), normThresPress((R_*max_Pc )/comn_.sigmaOW());//, curvRad(0.0);
		for(int i = 0; i < numCorners_; ++i)
		{
			if(waterInCorner_[i].cornerExists()) 
			{
				angSum += cos(recConAng + crnHafAngs_[i]);
			}
		}

		double rhsMaxAdvConAng = (-4.0*shapeFactor_*angSum)                      
			/ (normThresPress-cos(recConAng)+12.0*shapeFactor_*sin(recConAng));
		rhsMaxAdvConAng = max(rhsMaxAdvConAng, -1.0);    // Prevent falling out of range [1, -1].  This is only applicable when r is very small
		rhsMaxAdvConAng = min(rhsMaxAdvConAng, 1.0);  // and these elems are probably not drained.
		maxConAngSpont_ = acos(rhsMaxAdvConAng); 
	}
 
	if(cntAngAdv_ < maxConAngSpont_)
		Pc__pistonTypeAdv = Pc_pistonType_ImbHingCLine();
	else if(cntAngAdv_ <=  PI/2.0 + crnHafAngs_[0])            // If oil layers are not possible in all
		Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;     // corners, we do not use MSP procedure
	else                                                                //  => do the check with max(beta)
		Pc__pistonTypeAdv = Pc_pistonType_Imbnww();

}


/**
// When the number of oil filled neighbours change the imbibition entry pressure
// will need to be updated. E.g. piston type displacement might now be possible
// or we move from an I3 to I2 pore body filling event. If water injection is
// forced, coopertive pore body filling will not occur.
*/
double Polygon::centreEntryPrsWatInj()
{

	int num_WatCentreFeederNeis(elem_.num_WatCentreFeederNeis());
	double pistonEntryPrs;

	if(elem_.iAmAPore() && cntAngAdv_ < PI/2.0 && num_WatCentreFeederNeis > 0 )
	{
		double radSum(0.0);
		int iEvent(0);
		string poreBodyFillAlg(comn_.poreFillAlg());
		for(int i = 0; i < elem_.connectionNum(); ++i)
		{
			if( elem_.connection(i)->model()->affectsNeiEntryPc( comn_.oil() ) && !elem_.rockIndex() )
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

		if(poreBodyFillAlg == "blunt2")
			pistonEntryPrs = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)/R_ - radSum );
		else if(poreBodyFillAlg == "blunt1" || poreBodyFillAlg == "oren1")
			pistonEntryPrs = comn_.sigmaOW()*2.0*cos(cntAngAdv_)/(R_+radSum);
		else
			pistonEntryPrs = comn_.sigmaOW()*(1.0+2.0*sqrt(PI*shapeFactor_))*cos(cntAngAdv_)/(R_+radSum);
	}
	else
		pistonEntryPrs = Pc__pistonTypeAdv;

	//double mySlfpistonEntryPrs = pistonEntryPrs;


	double maxNeiPistonEntryPrs=-1.0e26;
	for(int i = 0; i < elem_.connectionNum(); ++i)///. parallel invasion; find the easiest throat
	{
		if( elem_.connection(i)->model()->conductCWater())
			maxNeiPistonEntryPrs = max( maxNeiPistonEntryPrs,
				elem_.connection(i)->model()->Pc_pistonTypeAdv());
	}
	if (maxNeiPistonEntryPrs>-1.0e25) ///. Not necessary, but anyway
		pistonEntryPrs = min(maxNeiPistonEntryPrs*0.999+0.001*pistonEntryPrs,pistonEntryPrs); ///. serial: the hardest of the pore and throat


	double snapOffPrs = waterInCorner_[0].trappedCorner().first<0 ? calcSnapOffPressureImb() : -1.0e27;
	if(num_WatCentreFeederNeis > 0 && pistonEntryPrs > snapOffPrs)
	{
		displacementType_ = 'P';
		 if (elem_.iAmAPore()) outD<<" ep"<<elem_.index()<<":"<<pistonEntryPrs<<" ";
		 else
		 {
			outD<<" e"<<displacementType_<<elem_.index()<<":"<<pistonEntryPrs<<" ";
			outD<<" e"<<'P'<<elem_.index()<<":"<<snapOffPrs<<" ";
		 }
	   return pistonEntryPrs;
	}
	else
	{
		displacementType_ = 'S';
		 if (elem_.iAmAPore()) outD<<" ep"<<elem_.index()<<":"<<pistonEntryPrs<<" ";
		 else 
			{
				outD<<" e"<<displacementType_<<elem_.index()<<":"<<snapOffPrs<<" ";
				outD<<" e"<<'P'<<elem_.index()<<":"<<pistonEntryPrs<<" ";
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
	double contactAngle = PI - cntAngAdv_;
	vector< double > potentialCurveRad;

	double sOne(0.0), sTwo(0.0), sThree(0.0);
	for(int i = 0; i < numCorners_; ++i)
	{
		if(cntAngAdv_ > PI/2.0 + crnHafAngs_[i])  // Only cases where oil layers/oil in corner might exist
		{
			sOne += cos(contactAngle)*cos(contactAngle+crnHafAngs_[i])/sin(crnHafAngs_[i])
				- (PI/2.0-contactAngle-crnHafAngs_[i]);
			sTwo += cos(contactAngle+crnHafAngs_[i])/sin(crnHafAngs_[i]);
			sThree += 2.0 * (PI/2.0 - contactAngle - crnHafAngs_[i]);

			double dFact = sOne - 2.0*sTwo*cos(contactAngle) + sThree;
			double rootFact = 1.0 + 4.0*shapeFactor_*dFact/(cos(contactAngle)*cos(contactAngle));

			double radOne = R_*cos(contactAngle)*(1.0-sqrt(rootFact))/(4.0*shapeFactor_*dFact);
			double radTwo = R_*cos(contactAngle)*(1.0+sqrt(rootFact))/(4.0*shapeFactor_*dFact);
			potentialCurveRad.push_back(max(radOne, radTwo));
		}
	}


	for(int j = potentialCurveRad.size()-1; j >= 0; --j)
	{
		double tension(comn_.sigmaOW());
		double pc(tension / potentialCurveRad[j]);
		double layerPc(INF_NEG_NUM);
		if(waterInCorner_[j].cornerExists()) 
		{///. Warning layers should be initialised for WatInj befor this
			layerPc = oilLayer_[j].entryPc();//layerCollapsePc(pc, cntAngAdv_, crnHafAngs_[j], tension,false);
		}

		if(pc > layerPc) return   comn_.sigmaOW() / potentialCurveRad[j];
	}

	return  comn_.sigmaOW() * (2.0*cos(cntAngAdv_)) / R_;
}



/**
// Piston type displacement (also referred to as event 1 displacement) for advancing contact angles less
// than the critical is obtained by solving a set of non-linear equations using the Newton-Raphson tech.
// The hinging contact angle will vary between the receding and advancing, fixed at a position bi. This
// procedure is identical to that used by Patzek and Oren.
*/
double Polygon::Pc_pistonType_ImbHingCLine() const
{
	double tension(comn_.sigmaOW()), err(1000.0);
	double newPc(1.1*tension*(2.0*cos(cntAngAdv_)) / R_);

	double oldPc;

	int itr;
	double newPc2 = newPc;
	for(itr = 0; itr < MAX_NEWT_ITR+1; ++itr)
	{
		double sumOne(0.0), sumTwo(0.0), sumThree(0.0), sumFour(0.0);
		oldPc = newPc;
		for(int i = 0; i < numCorners_; ++i)
		{ 
			if(waterInCorner_[i].cornerExists())
			{ 

				double meniscusApexDist, hingConAng(cntAngAdv_); 
				waterInCorner_[i].getCApexDistConAng(meniscusApexDist, hingConAng, oldPc, crnHafAngs_[i], comn_.sigmaOW(),true,true);

				double partus(meniscusApexDist * sin(crnHafAngs_[i]) * oldPc/tension);
				ensure(partus >= -1.0 && partus <=  1.0);
				if(!(partus >= -1.0 && partus <=  1.0))
				{
					cout<<partus<<" = "<<meniscusApexDist<<" * sin(" <<crnHafAngs_[i]<<") / "<<tension/oldPc<<endl;
					cout<<"  oldPc "<<oldPc<<endl;
					cout<<"  maxConAngSpont_ "<<maxConAngSpont_<<endl;
					cout<<"  tension/oldPc "<<tension/oldPc<<endl;
				}

				sumOne += meniscusApexDist*cos(hingConAng);
				sumTwo += PI/2.0-hingConAng-crnHafAngs_[i];
				sumThree += asin(partus);
				sumFour += meniscusApexDist;
			}
		}


			double a= 2.0*sumThree - sumTwo;
			double b=cos(cntAngAdv_)*(R_/(2.0*shapeFactor_)) -2.0*sumFour + sumOne;
			double c=-R_*R_/(4.0*shapeFactor_);
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
		cout<< "Problem in Polygon::mspCurveRadHingImb error:" <<err <<" %of " <<newPc <<"  teta " <<cntAngAdv_ << endl;
		return newPc;
	}


	{cerr << endl
		<< "=================================================" << endl
		<< "Error: Failed to obtain valid value for threshold" << endl
		<< "radius of curvature in piston type displacement  " << endl
		<< "during water injection."                           << endl
		<< "Trapped water: " << waterInCorner_[0].trappedCorner().first << endl
		<< "Err  " << err << endl
		<< "Con ang " << cntAngAdv_*180.0/PI << endl
		<< "Con ang " << cntAngAdv_ << endl
		<< "Con ang " << cntAngRec_ << endl
		<< "Radius " << R_ << endl
		<< "oldPc " << oldPc << endl
		<< "newPc " << newPc << endl
		<< "G " << shapeFactor_ << endl
		<< "Iteration " << itr << endl
		<< "adPc " << waterInCorner_[0].advancingPc() << endl
		<< "recPc " << waterInCorner_[0].receedingPc() << endl
		<< "recPc " << waterInCorner_[0].pinnedApexDist() << endl
		<< "=================================================" << endl;   
					//((double*)&R_)[10000000]=0.0;
		//exit(-1);
	}
	return 0.0;
}


///. TO PERMANENTLY DELETE
bool Polygon::waterLayer_UntrappedCorner_PcLsnapPc(double cappPrs) const
{  ///. checks if we need to untrap oil, insertReCalcImbibeEntryPrs ...         never does anything
	if (!waterInCorner_[0].cornerExists()) return false;
	if (waterInCorner_[0].trappedCorner().first>-1 && elem_.trappingWatFilm().first<0) 
	{
		cout<<" * Error: unsynced Trapping * "<<endl;
	}
	if (waterInCorner_[0].trappedCorner().first>-1 || eleman()->isTrappedOil() || oilLayer_->trappingCL().first>-1) return false;
	if (comn_.injectant() == &comn_.oil()) return false;
	if (conductCOil())
	{


		for (int i=0;i<numCorners_;++i)
		{
			if (waterInCorner_[i].cornerExists() && waterInCorner_[i].trappedCorner().first<0 && !waterInCorner_[i].pinned()) 
			{
				waterInCorner_[i].initCornerApex(cappPrs,cntAngRec_,cntAngAdv_,crnHafAngs_[i],comn_.sigmaOW(),false);
				cout<<"j";
			}
		}

		double snapPc = calcSnapOffPressureImb();

		return  cappPrs < snapPc;

	}


	return false;/// WARninG;
}




///. Warning also called during Oil Injection, when decreasing oil pressure due to coalescence 
double Polygon::Pc_pin_disconnectOilLayer(int cor)  ///rare
{
	ensure(oilConnection_ && numLayers_>0);    
	LayerApex& oilLayer = oilLayer_[cor];
	double layerPc(oilLayer.layerCollPc());
	double tension =  comn_.sigmaOW();
	if (waterInCorner_[cor].cornerExists())
	{

		waterInCorner_[cor].finitCornerApex(layerPc, cntAngRec_, cntAngAdv_, crnHafAngs_[cor],tension, false, false);
		waterInCorner_[cor].initCornerApex(layerPc, cntAngRec_, cntAngAdv_, crnHafAngs_[cor],tension, false);

		hasDisConectedCentreWCornerW_ = false;
		/// Warn;
	}

	oilLayer.removeLayer();
	--numLayers_;
	ensure(numLayers_ >= 0 && numLayers_ < numCorners_);

	if(!numLayers_)
	{
		ensure(!hasDisConectedCentreWCornerW_);
		oilConnection_ = false;// if(!numLayers_)
		if (!containCOil() && elem_.isTrappedOil()) {cout<<" axsj ";elem_.unTrapOil();}
	}
	return layerPc;
}






