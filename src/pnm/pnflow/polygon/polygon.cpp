
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "sortedEvents.h"
#include "compareFuncs.h"


using namespace std;

/**
// Polygon Class Constructors
*/
VoidElem::VoidElem(Element& parent, const CommonData& common, 
	double radius, double shapeFactor, int connNum) 
	: ElemModel(parent, common, radius, connNum), shapeFactor_(shapeFactor)
				 
{
	area_ = pow(R_, 2.0) / (4.0 * shapeFactor_);     //From Oren; This is ~correct for all shape
	SatWater_ = 1.0;
	minInitRecCntAng_ = 180.0;
	cntAngAdv_ = -1.0;
   cntAngRec_ = -1.0;


	ensure(area_ > 0.0);

}


void Polygon::insertWatSnapEvent_IfSnapPcHgPc(SortedEvents<Apex*,PceImbCmp>& watEvents,  double globalPc)
{
	if(bulkFluid_ == &comn_.oil())  // Water snap off when the bulk is filled by
	{                                                   // oil. Make sure there is no trapping in
		if (waterInCorner_[0].cornerExists() && (waterInCorner_[0].trappedCorner().first<0) && !(elem_.isTrappedOil()))
			if(elem_.entryPc() > globalPc) 
				watEvents.quickInsert(&elem_);
	}
	else if( oilLayer_[0].exists() && (oilLayer_[0].trappedOLayer().first<0) &&
		(waterInCorner_[0].trappedCorner().first<0 || elem_.trappingWatBulk().first<0))   // Collapsing oil layers. Again make sure that
	{                                                                                           // trapping allows us to collapse it
		for(int i = 0; i < numCorners_; ++i)
		 if(oilLayer_[i].LayerApex::exists(/*st ab le*/))
			if(oilLayer_[i].layerCollPc() > globalPc)  ///. Why such check
			{
				pair<Polygon*, int> newEvent(this, i);//, collPc+elem_.gravityCorrection()
				watEvents.quickInsert(&oilLayer_[i]);
			}
	}
}


void Polygon::insertOilSnapEvent_IfSnapPcLgPc(SortedEvents<Apex*,PceDrainCmp>& oilEvents, double globalPc)
{
	if(bulkFluid_ == &comn_.water() && oilLayer_[0].exists(/*st ab le*/) && oilLayer_[0].trappedOLayer().first<0)
	{   
		ensure(numLayers_==0);

		double snapOffPrs = calcSnapOffPressureDrain();

		 cout<<"  coel Pcs:"<< elem_.entryPc() << "  " <<snapOffPrs<<endl;
		if(snapOffPrs < globalPc)
		{
			//pair<Polygon*, int> newEvent(this, -1);///, /*_*/snapOffPrs+elem_.gravityCorrection()
			oilEvents.quickInsert(&elem_);
		}
	}
}




/**
 * set radius and update area and waterInCentre conductivity accordingly
 */
void VoidElem::setRadius(double newRad)
{
	R_ = newRad;
	area_ = pow(newRad, 2.0) / (4.0 * shapeFactor_);
	ensure(area_ > 0.0);
}




///  Contact angles are assigned based on the principle of equilibrium contacr angles
/// DAng  is (modelTwo/5) Seperation angle,
void VoidElem::setContactAngle(double refCAng, int wettClass, double DAng)
{


	if(wettClass == 1)
	{
		cntAngRec_ = refCAng;
		cntAngAdv_ = refCAng;
	}
	else if(wettClass == 2)
	{
		double growthExp((PI+DAng)/PI);

		cntAngRec_ = max(0.0, growthExp*refCAng-DAng);
		cntAngAdv_ = min(PI, growthExp*refCAng);
	}
	else  if(wettClass == 3)
	{
		if     (refCAng < 0.38349) cntAngRec_ = 0.0;
		else if(refCAng < 1.5289) cntAngRec_ = (0.5*exp(0.05*refCAng*180.0/PI)-1.5)*PI/180.0;
		else if(refCAng < 2.7646) cntAngRec_ = 2.0*(refCAng-1.19680);
		else cntAngRec_ = PI;

		if     (refCAng < 0.38349) cntAngAdv_ = 0.0;
		else if(refCAng < 1.61268) cntAngAdv_ = 2.0*(refCAng-0.38349);
		else if(refCAng < 2.75805) cntAngAdv_ = (181.5 - 4051.0*exp(-0.05*refCAng*180.0/PI))*PI/180.0;
		else                       cntAngAdv_ = PI;
	}
	else  if(wettClass == 4)
	{
		cntAngAdv_ = refCAng;
		cntAngRec_ = pow(PI-1.3834263-pow(PI-cntAngAdv_+0.004,0.45), 1.0/0.45)-0.004;
	}
	else
	{
		double PlusCoef = PI-(0.1171859 *DAng*DAng*DAng -0.6614868 *DAng*DAng + 1.632065 *DAng) ;
		double powrCoef = 1.0-(0.01502745 *DAng*DAng*DAng -0.1015349 *DAng*DAng + 0.4734059 *DAng) ;

		cntAngAdv_ = refCAng;
		cntAngRec_ = pow(PlusCoef-pow(PI-cntAngAdv_+0.004,powrCoef), 1.0/powrCoef)-0.004;
	}

	Polygon *polyShape = dynamic_cast< Polygon* >(this);
	if(polyShape)
	{
		for(int i = 0; i < polyShape->numCorners(); ++i)
			polyShape->oilLayerCh()[i].advConAng(cntAngAdv_);
	}

	minInitRecCntAng_=min(cntAngRec_,minInitRecCntAng_);

	//Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;     // TOBE initialised properly later
	//Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_;     // TOBE initialised properly later

}




/**
// Polygon Class Constructors
*/
Polygon::Polygon(Element& parent, const CommonData& common, double radius, 
				 double shapeFactor, int nCorners, int connNum)
		: VoidElem(parent, common, radius, shapeFactor, connNum),
	numCorners_(nCorners),displacementType_('X')
{
	numLayers_ = 0;
	maxConAngSpont_ = PI/2;

	///. Warning should be constructed as array to allow comparison 
	waterInCorner_ = new CornerApex[nCorners];
	oilLayer_ = new LayerApex[nCorners];
	//waterInCorner_.reserve(nCorners);
	//oilLayer_.reserve(nCorners);    
	for(int i = 0; i < nCorners; ++i)
	{
		waterInCorner_[i].setCornerConnections(oilLayer_+i,this,i);
		oilLayer_[i].setLayerConnections(waterInCorner_+i,this,i);
	}
}






/**
// Triangle Class Constructors
*/
Triangle::Triangle(Element& parent, const CommonData& common,  double radius,
				   double shapeFactor, int connNum, int) 
	  : Polygon(parent, common, radius, shapeFactor, 3, connNum)
{
	crnHafAngs_.resize(numCorners_);
	setHalfAngles();
}

void Triangle::setShapeFactor(double shapeFact)
{
	ensure(shapeFact <=  sqrt(3.0)/36.0);
	shapeFactor_ = max(shapeFact,1.0e-6);
	area_ = pow(R_, 2.0) / (4.0 * shapeFactor_);
	ensure(area_ > 0.0);
	//areaWater_ = area_;
	setHalfAngles();
	//conductanceWater_.second = SPConductance(areaWater_, comn_.water().viscosity());
}

/**
// This function evaluates the three half angles that make up the triangular pore. The routine
// follows the outline that was described by Patzek.
*/
void Triangle::setHalfAngles()
{
	double beta_2_min = atan((2.0/sqrt(3.0))*cos(acos(-12.0*sqrt(3.0)*shapeFactor_)/3.0+4.0*PI/3.0));
	double beta_2_max = atan((2.0/sqrt(3.0))*cos(acos(-12.0*sqrt(3.0)*shapeFactor_)/3.0));
	double randNum = double(rand()) / double(RAND_MAX);
	randNum = 0.5*(randNum+0.5);///.25-.75

	crnHafAngs_[1] = beta_2_min + (beta_2_max - beta_2_min)*randNum;
	crnHafAngs_[0] = 0.5*(asin((tan(crnHafAngs_[1])+4.0*shapeFactor_)
		* sin(crnHafAngs_[1]) / (tan(crnHafAngs_[1])-4.0*shapeFactor_))-crnHafAngs_[1]);
	crnHafAngs_[2] = PI/2.0 - crnHafAngs_[1] - crnHafAngs_[0];
}


/**
// Square Class Constructors
*/
Square::Square(Element& parent, const CommonData& common, double radius, int connNum, int) 
   : Polygon(parent, common, radius, 0.0625, 4, connNum)
{
	crnHafAngs_.resize(numCorners_, PI/4.0);
}


/**
// Circle Class Constructors
*/
Circle::Circle(Element& parent, const CommonData& common, double radius, int connNum, int)
	  : VoidElem(parent, common, radius, 1.0/(4.0*PI), connNum)
{
}




/**
// Polygon Class Destructor
*/
Polygon::~Polygon()
{
	//for(int i = 0; i < numCorners_; ++i)
	//{
		delete[] waterInCorner_;
		delete[] oilLayer_;
	//}
}







///. WARNING: does not check from which throat water injection is taking place
void Polygon::fillCentreWithWaterCreateLayers(bool snapOffOverRide)
{
	bulkFluid_ = &comn_.water(); 	waterConnection_ = true;	oilConnection_ = false;
	//ensure(waterInCorner_[0].trappedCorner().first<0," flxpq ");
	//double snapOffPrs = ;
	double pc = comn_.GuessCappPress()-elem_.gravityCorrection();

	for(int i = 0; i < numCorners_; ++i)
	{
		double snapOffPrs = waterInCorner_[0].trappedCorner().first<0 ? calcSnapOffPressureImb() : -1.0e27;
		if( !snapOffOverRide && elem_.entryPc() >  snapOffPrs)
		{
			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i], comn_.sigmaOW(),false);
			//LayerApex* oillayeri = ;
			if ( oilLayer_[i].createOLayer(pc, cntAngRec_, cntAngAdv_, maxConAngSpont_, crnHafAngs_[i],  comn_.sigmaOW(),false))
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
				++numLayers_;
				//this->addFillingEventToHistory(5)
			}
			else
			{
				waterInCorner_[i].removeCorner();///Warning;  inconsistent film trapping and connectivity handling
			}
		}

	}

	if(numLayers_)   oilConnection_ = true; 
	hasDisConectedCentreWCornerW_ = numLayers_ == numCorners_;

	Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;

}


void Polygon::fillCentreWithOilRemoveLayers()
{
	 oilConnection_ = true; bulkFluid_ = &comn_.oil();
	 ensure(cntAngRec_>=0.0,"contact angle not yet set",2);
   

	for(int i = 0; i < numCorners_; ++i)
	{
	  
		oilLayer_[i].LayerApex::removeLayer();//_ exists


		if (!waterInCorner_[i].cornerExists())
			waterInCorner_[i].createFilm(elem_.entryPc(), cntAngRec_, cntAngAdv_, crnHafAngs_[i], comn_.sigmaOW(), true);
	}
	numLayers_ = 0;
	hasDisConectedCentreWCornerW_ = false;

	waterConnection_ = waterInCorner_[0].cornerExists();

	Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_;

}







void Circle::fillCentreWithWaterCreateLayers(bool snapOffOverRide)
{
	bulkFluid_ = &comn_.water();  oilConnection_ = false;  waterConnection_ = true;
	Pc__pistonTypeAdv = comn_.sigmaOW()*(2.0*cos(cntAngAdv_)) / R_;

}


void Circle::fillCentreWithOilRemoveLayers()
{

	oilConnection_ = true; bulkFluid_ = &comn_.oil();	waterConnection_ = false;///K

	Pc__pistonTypeRec = comn_.sigmaOW()*(2.0*cos(cntAngRec_)) / R_;

}


void Polygon::calcOilLayerPc_syncTrappings(double pc)
{///. pc is used for error checking only

	const pair<int, double>& oilTrp = elem_.trappingOil();
	const pair<int, double>& watTrpBulk = elem_.trappingWatBulk();
	const pair<int, double>& watTrpFilm = elem_.trappingWatFilm();

	//ensure (elem_->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1)

	///Warning double ComToFix = 1.0;

	bool oilInj =  comn_.injectant() == &comn_.oil();
	double tension =  comn_.sigmaOW();
	if(oilTrp.first > -1)
	{
		if (containCOil()&&watTrpBulk.first>0) cout<<"\n * Wrong trapping flags * \n";



		for(int i = 0; i < numCorners_; ++i)
		{				double oldLPc=oilLayer_[i].layerCollPc();



			waterInCorner_[i].finitCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj, true);
			if (!oilLayer_[i].finitLayerApex(oilTrp.second, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj, true))
			{
				if (debugLevel>100) 
				{	vector<int> addCrns;
					cout<<endl<<" q3jd "; 
					cout
					<<elem_.isInWatFloodVec()
					<<elem_.canBeAddedToEventVec(comn_.water())
					<<elem_.addToLayerVec(comn_.water(), this, addCrns)
					<<(elem_.trappingWatFilm().first>-1)
					<<(elem_.trappingWatBulk().first>-1)
					<<(oilLayer_[i].trappedOLayer().first>-1) 
					<<(oilLayer_[i].isInWatFloodVec())
					<<canNOTReconfigure(comn_.water())
					<<(elem_.trappingOil().first>-1)

					<<"\n   "<<oilTrp.second<<" <?> "<<oilLayer_[i].layerCollPc()<<" <?> "<<oldLPc<<" <?> "<<pc<< "    "
					<<"   "<<waterInCorner_[i].advancingPc()<<" <?> "<<oilLayer_[i].advancingPc()<<" <?> "<<0<<" <?> "<<pc<< "    "
					<<endl;
					//crnHafAngs_[100000]=0;
					//exit(-1);
				}

				Pc_pin_disconnectOilLayer(i);


			}
			else  if (oilLayer_[i].exists())
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
			}

			oilLayer_[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			comn_.sigmaOW(), oilInj);

			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj);
		}

		elem_.virgin = false;
	}



	if(watTrpFilm.first > -1)
	{

		for(int i = 0; i < numCorners_; ++i)
		{
			waterInCorner_[i].finitCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj, true);
			if (!oilLayer_[i].finitLayerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj,false))
			{ 
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" q4jd ";

			}
  			waterInCorner_[i].CornerApex::markTrappingCorner(watTrpFilm, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			   tension, oilInj); ///. 

			oilLayer_[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			   tension, oilInj); ///. 


			if (!oilLayer_[i].initLayerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj))
					if (debugLevel>0) cout<<" q4jd2";


		}
	}

	if(watTrpBulk.first > -1)
	{

		for(int i = 0; i < numCorners_; ++i)
		{
			waterInCorner_[i].finitCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj,false);
			if (!oilLayer_[i].finitLayerApex(watTrpBulk.second, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj, true))
			{
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" q8jd ";

			}
			else if (oilLayer_[i].exists())
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
				if(!(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1])))))
				{
					//elem_.unTrapWat(bulkBlob);
					cout<<endl
					<<elem_.canBeAddedToEventVec(comn_.oil())
					<<elem_.isInOilFloodVec()
					<<(elem_.trappingWatFilm().first>-1)
					<<canNOTReconfigure(comn_.oil())
					<<(oilLayer_[i].trappedOLayer().first>-1) 
					<< (elem_.trappingOil().first>-1)
					<< (elem_.trappingWatBulk().first>-1)<<endl;

					cout<<snapOffPrsDrainSpontNR()<<" <?> "
					<<"   "<<watTrpBulk.second<<" <?> "<<elem_.entryPc()<< "    ";
					cout<<oilLayer_[i].pinnedApexDist() <<"  <  "<< R_<<" * (1.0/tan( "<<crnHafAngs_[0]<<" ) + 1.0/tan( "<<crnHafAngs_[1]<<endl;
					cout<<oilLayer_[i].pinnedApexDist() <<"  <  "<< (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1])))<<endl;
					//((double*)&R_)[100000]=0.0;

					//exit(-1);
				}
			}

			//oilLayer_[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			   //tension, oilInj); ///. 
			   
			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj);
			if (!oilLayer_[i].initLayerApex(watTrpBulk.second, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj))
					if (debugLevel>0) cout<<" q4jd3";

		}
		elem_.virgin = false;
	}



}




void Polygon::calcOilLayerPc_markUntrappedFilms(double pc)
{///. pc is used for error checking only

	const pair<int, double>& oilTrp = elem_.trappingOil();
	const pair<int, double>& watTrpBulk = elem_.trappingWatBulk();
	const pair<int, double>& watTrpFilm = elem_.trappingWatFilm();

	//ensure (elem_->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1)

	///Warningdouble ComToFix = 1.0;

	bool oilInj =  comn_.injectant() == &comn_.oil();
	double tension =  comn_.sigmaOW();

	//if(oilTrp.first > -1)
	{

		for(int i = 0; i < numCorners_; ++i)
		{				double oldLPc=oilLayer_[i].layerCollPc();

  			waterInCorner_[i].CornerApex::markTrappingCorner(watTrpFilm, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			   tension, oilInj); ///. 
			oilLayer_[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			comn_.sigmaOW(), oilInj);

			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj);
			if (!oilLayer_[i].initLayerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj, true))
			{
				if (debugLevel>100) 
				{	vector<int> addCrns;
					cout<<" q3jdu "; 
					cout
					<<elem_.isInWatFloodVec()
					<<elem_.canBeAddedToEventVec(comn_.water())
					<<elem_.addToLayerVec(comn_.water(), this, addCrns)
					<<(elem_.trappingWatFilm().first>-1)
					<<(elem_.trappingWatBulk().first>-1)
					<<(oilLayer_[i].trappedOLayer().first>-1) 
					<<(oilLayer_[i].isInWatFloodVec())
					<<canNOTReconfigure(comn_.water())
					<<(elem_.trappingOil().first>-1)

					<<"\n   "<<oilTrp.second<<" <?> "<<oilLayer_[i].layerCollPc()<<" <?> "<<oldLPc<<" <?> "<<pc<< "    "
					<<"   "<<waterInCorner_[i].advancingPc()<<" <?> "<<oilLayer_[i].advancingPc()<<" <?> "<<0<<" <?> "<<pc<< "    "
					<<endl;
					//crnHafAngs_[100000]=0;
					//exit(-1);
				}

				Pc_pin_disconnectOilLayer(i);
			}
			else  if (oilLayer_[i].exists())
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
			}

			//waterInCorner_[i].finitCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj,false);
			//if (!oilLayer_[i].finitLayerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj,false))
			//{
					//if (debugLevel>0) cout<<" q4jd3";
				//Pc_pin_disconnectOilLayer(i);
//
			//}
		}
	}

	return;

	//if(watTrpFilm.first > -1)
	{

		for(int i = 0; i < numCorners_; ++i)
		{
 			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj);
			if (!oilLayer_[i].initLayerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj))
			{ 
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" q4jdu ";

			}

		}

	}

	//if(watTrpBulk.first > -1)
	{

		for(int i = 0; i < numCorners_; ++i)
		{
			//oilLayer_[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],
			   //tension, oilInj); ///. 

			waterInCorner_[i].initCornerApex(pc, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj);
			if (!oilLayer_[i].initLayerApex(watTrpBulk.second, cntAngRec_, cntAngAdv_, crnHafAngs_[i],tension, oilInj))
			{
				Pc_pin_disconnectOilLayer(i);
				if (debugLevel>0) cout<<" q8jdu ";
			}
			else if (oilLayer_[i].exists())
			{
				ensure(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1]))));
				if(!(oilLayer_[i].pinnedApexDist() < (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1])))))
				{
					//elem_.unTrapWat(bulkBlob);
					cout<<endl
					<<elem_.canBeAddedToEventVec(comn_.oil())
					<<elem_.isInOilFloodVec()
					<<(elem_.trappingWatFilm().first>-1)
					<<canNOTReconfigure(comn_.oil())
					<<(oilLayer_[i].trappedOLayer().first>-1) 
					<< (elem_.trappingOil().first>-1)
					<< (elem_.trappingWatBulk().first>-1)<<endl;

					cout<<snapOffPrsDrainSpontNR()<<" <?> "
					<<"   "<<watTrpBulk.second<<" <?> "<<elem_.entryPc()<< "    ";
					cout<<oilLayer_[i].pinnedApexDist() <<"  <  "<< R_<<" * (1.0/tan( "<<crnHafAngs_[0]<<" ) + 1.0/tan( "<<crnHafAngs_[1]<<endl;
					cout<<oilLayer_[i].pinnedApexDist() <<"  <  "<< (R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[1])))<<endl;
					//((double*)&R_)[100000]=0.0;

					//exit(-1);
				}
			}


		}
	}



}




/**
// Spontaneous snap off (Contact ang < PI/2 - beta_min):
// Snap off ocuurs when the base of the menisci in the two sharpest corners meet.
//
// Forced snap off:
// Snap off occurs when the contact angle in the sharpest corner reach advancing
// contact angle.
*/
double Polygon::calcSnapOffPressureImb() const
{
	double snapOffPrs;
	//int * sasasas;
	if (waterInCorner_[0].trappedCorner().first>-1 || eleman()->isTrappedOil() || oilLayer_->trappingCL().first>-1) 
	{
		if (debugLevel>100) cout<<" nkdps ";
		return -1.0e28;
	} 

	if(cntAngAdv_ < PI/2.0 - crnHafAngs_[0])              // Spontaneous
		snapOffPrs = snapOffPrsImbNR();
	else                                                   // Forced
		snapOffPrs = waterInCorner_[0].cornerExists() ? waterInCorner_[0].advancingPc() : -1.0e64;
   
	return snapOffPrs-elem_.snapOfLongitCurvature()*comn_.sigmaOW();

}


double Polygon::calcSnapOffPressureDrain() const
{
	double snapOffPrs;
	if(cntAngRec_ > PI/2.0 + crnHafAngs_[0])      // Spontaneous
	{
		snapOffPrs = snapOffPrsDrainSpontNR();
		ensure(snapOffPrs!=0.0);
	}
	else
	{                                                       // Forced
		snapOffPrs = oilLayer_[0].receedingPc();
		ensure(snapOffPrs>0.0);
		if (snapOffPrs <=  0.0)
		{
			cout<<oilLayer_[0].receedingPc()<<endl;
		}
	}


	return snapOffPrs-elem_.snapOfLongitCurvature()*comn_.sigmaOW();
}



/**
// Spontaneous snap off can occur in two ways. If advancing contact angle is small
// all the menisci will move and snap off occurs when the menicii in the two
// sharpest corners meet. For larger contact angles one or two of the meniscii may
// remain pinned. Snap off occurs when the sharpest coner meets any of the two
// others. The event actually occuring is the one at the highest pressure.
*/
double Triangle::snapOffPrsImbNR() const
{

//not converging

	double tension(comn_.sigmaOW());
	//- Tom Bultreys bug-report & fix R-> 1/R
	double snapOffPrsOne = tension * cos(cntAngAdv_+crnHafAngs_[0] )/((R_/tan(crnHafAngs_[0]) + R_/tan(crnHafAngs_[2]) - waterInCorner_[2].pinnedApexDist())*sin(crnHafAngs_[0]));

	double oldPc, errorVal(1.0);
	if(waterInCorner_[0].trappedCorner().first>-1)
		oldPc = waterInCorner_[0].trappedCorner().second;
	else if(comn_.isDrainage())
		oldPc = comn_.GuessCappPress()-elem_.gravityCorrection();
	else
		oldPc = comn_.maxPcLastDrainCycle()-elem_.gravityCorrection();

	const double L0pL2 = (R_/tan(crnHafAngs_[0]) + R_/tan(crnHafAngs_[1]))   ;
	for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
	{

		double apexDist2(0.0), teta2(cntAngAdv_); 
		waterInCorner_[1].getCApexDistConAng(apexDist2, teta2, oldPc, crnHafAngs_[1], tension, true, true);

		double rL2 = -apexDist2*sin(crnHafAngs_[1])/(tension*sin(teta2+crnHafAngs_[1])); 

		double func = oldPc - tension*(cos(cntAngAdv_)/tan(crnHafAngs_[0]) - sin(cntAngAdv_) + cos(teta2)/tan(crnHafAngs_[1]) - sin(teta2)) / L0pL2;

		double funcDpc = 1 + tension*(rL2*sin(teta2)/tan(crnHafAngs_[1]) +  rL2*cos(teta2)) /L0pL2;

		double newPc = (abs(funcDpc)>1.0e-32) ? oldPc - func / funcDpc  : oldPc - func;
		if (i>MAX_NEWT_ITR/2 ) newPc = 0.5*(newPc+oldPc);

		errorVal = fabs(newPc-oldPc)/(abs(oldPc)+1.0e-8);
		if(errorVal < EPSILON)
		{
			if(teta2 <= cntAngAdv_+0.000001)
				return max(newPc, snapOffPrsOne)+0.0001;     // Return the pressure occuring at the highest Pc
			else
				return snapOffPrsOne+0.0001;
		}

		oldPc = newPc;
	}

/*    double tension(comn_.sigmaOW());
	double snapOffPrsOne = (tension/R_) * (cos(cntAngAdv_) -
			(2.0*sin(cntAngAdv_))/(1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[1])));    // Small Contact Ang

	double oldPc, errorVal(1.0);
	if(waterInCorner_[0].trappedCorner().first>-1)
		oldPc = waterInCorner_[0].trappedCorner().second;
	else if(comn_.isDrainage())
		oldPc = comn_.GuessCappPress()-elem_.gravityCorrection();
	else
		oldPc = comn_.maxPcLastDrainCycle()-elem_.gravityCorrection();

	for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
	{

		double pinnedApexDistance, hingAng(cntAngAdv_); 
		waterInCorner_[2].getCApexDistConAng(pinnedApexDistance, hingAng, oldPc, crnHafAngs_[2], comn_.water().interfacialTen(), true, true);


		double hingAngDpc = -pinnedApexDistance*sin(crnHafAngs_[2])/(tension*sin(hingAng+crnHafAngs_[2]));

		double func = oldPc - (comn_.sigmaOW()/R_)*((cos(cntAngAdv_)/tan(crnHafAngs_[0]) -
			sin(cntAngAdv_) + cos(hingAng)/tan(crnHafAngs_[2]) - sin(hingAng)) /
			(1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[2])));

		double funcDpc = 1 + (comn_.sigmaOW()/R_)*((hingAngDpc*sin(hingAng)/tan(crnHafAngs_[2]) +
			hingAngDpc*cos(hingAng)) / (1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[2])));

		double newPc = oldPc - func / funcDpc;

		errorVal = fabs((newPc-oldPc)/oldPc);
		if(errorVal < EPSILON)
		{
			if(hingAng < cntAngAdv_)
				return max(newPc, snapOffPrsOne);     // Return the pressure occuring at the highest Pc
			else
				return snapOffPrsOne;
		}

		oldPc = newPc;
	}*/


	cerr <<"\n============  snap Off Prs Imb NR ================ " << endl
		<< "Error: Failed to obtain valid value for threshold" << endl
		<< "capillary pressure during snap off.              " << endl
		<< "cntAngAdv_: " << cntAngAdv_                               << endl
		<< "trapping ind" << elem_.trappingOil().first << endl
		<< "error " << errorVal<<oldPc << endl
		<< "=================================================" << endl;   

		 //exit(-1);

	return oldPc;
}

double Square::snapOffPrsImbNR() const///. to check
{
		return (comn_.sigmaOW()/R_) * (cos(cntAngAdv_) - sin(cntAngAdv_));
}

/**
// The problem with these NR approaches are that they're a bit sensitive
// to the first guess.
*/
double Triangle::snapOffPrsDrainFromCorner(bool& casePossible, int cor) const
{
	if(!oilLayer_[cor].exists()) return 0.0;

	if (cntAngRec_ <= PI/2.0 + crnHafAngs_[cor])
	{ 
		casePossible=true;
		return oilLayer_[cor].receedingPc();
	}

	double minLocalPcLastCycle(comn_.minPcLastImbCycle()-elem_.gravityCorrection());
	double minLocalPc(comn_.minEverCappPress()-elem_.gravityCorrection());

	double min_Pc = oilLayer_[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc;
	double tension(comn_.sigmaOW()), errorVal(1.0);
	vector< double > initGuesses(3);

	initGuesses[2] = min_Pc ;                                                 // First attempt
	initGuesses[1] = comn_.GuessCappPress()-elem_.gravityCorrection();    // Second attempt
	if(oilLayer_[cor].trappedOLayer().first>-1)                          // Getting desparete now
		initGuesses[0] = oilLayer_[cor].trappedOLayer().second;
	else
		initGuesses[0] = elem_.trappingOil().second;

	while(!initGuesses.empty())
	{
		double oldPc = initGuesses.back();
		initGuesses.pop_back();

		for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
		{
			double tetaCor(cntAngRec_);// = oilLayer_[cor].hingingConAng(oldPc,cntAngRec_,crnHafAngs_[cor],tension,true);// acos(part)+crnHafAngs_[cor];
			double pinnApexThree,halfAng(crnHafAngs_[cor]);// = oilLayer_[cor].getApexDistance(oldPc,tetaCor,crnHafAngs_[cor],tension);// acos(part)+crnHafAngs_[cor];
			oilLayer_[cor].getCAApexDist(pinnApexThree,tetaCor,halfAng,oldPc,tension);// acos(part)+crnHafAngs_[cor];

			ensure(pinnApexThree>0);

			double tetaCorDpc = -pinnApexThree*sin(crnHafAngs_[cor])/(tension*sin(tetaCor-crnHafAngs_[cor]));

			double func = oldPc - (tension/R_)*(
				(cos(cntAngRec_)/tan(crnHafAngs_[0]) + sin(cntAngRec_) +  cos(tetaCor)/tan(crnHafAngs_[cor]) + sin(tetaCor)) /
				(1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[cor]))
				);

			double funcDpc = 1 - (tension/R_)*((tetaCorDpc*sin(tetaCor)/tan(crnHafAngs_[cor]) -
				tetaCorDpc*cos(tetaCor)) / (1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[cor])));

			double newPc = oldPc - func / funcDpc;
			errorVal = fabs((newPc - oldPc)/oldPc);

			if(errorVal < EPSILON)
			{
				ensure(tetaCor >= cntAngRec_ || debugLevel < 1000);
				casePossible = tetaCor >= cntAngRec_;//  && oilLayer_[cor].freeAtPrs(newPc);

				return newPc;
			}

			oldPc = newPc;
		}
	}

	if ( debugLevel > 0 )
	cerr <<"==============oil Apex Meet Corner =============" << endl
		<< "Error: failed to obtain valid value for threshold" << endl
		<< "       capillary pressure during snap off.       " << endl
		<< "Error: " << errorVal                               << endl
		<< "cntAngRec_: " << cntAngRec_                               << endl
		<< "corner: " << cor                               << endl
		<< "minLocalPcLastCycle: " << minLocalPcLastCycle                               << endl
		<< "minLocalPc: " << minLocalPc                               << endl
		<< "min_Pc: " << min_Pc                               << endl
		<< "halfAng: " << crnHafAngs_[cor]                               << endl
		<< "=================================================" << endl;   ///. exit(-1);

	return 0.0;
}

/**
 * ///. revision needed
// 4 Different scenarios exists for oil snap off:
//
//  Case 1: There is only oil layer in sharpest corner. Snap off occurs when
//          interface meets opposing wall.
//  Case 2: Interfaces in the two sharpest corners slip. Snap off occurs when
//          they meet at the base. ///. Error
//  Case 3: Interfaces exist in all corners. Only the sharpest corner slip.
//          Snap off occurs when it meets the pinned interface in the most
//          oblique corner.
//  Case 4: Interfaces exist only in the two sharpest corners. Only the
//          sharpest corner slip. Snap off occurs when it meets the pinned
//          interface.
*/
double Triangle::snapOffPrsDrainSpontNR() const
{
	double tension(comn_.sigmaOW());
	double caseOneSnapOff = tension*cos(cntAngRec_-crnHafAngs_[0])
	/(sin(crnHafAngs_[0])* R_*(1.0/tan(crnHafAngs_[0])+1.0/tan(crnHafAngs_[2])));

	double caseTwoSnapOff = (tension/R_) * (cos(cntAngRec_) +
			(2.0*sin(cntAngRec_))/(1.0/tan(crnHafAngs_[0]) + 1.0/tan(crnHafAngs_[1])));

	bool case3Possible(false), case4Possible(false);

	double caseThreeSnapOff = snapOffPrsDrainFromCorner(case3Possible, 2);
	double caseFourSnapOff = snapOffPrsDrainFromCorner(case4Possible, 1);

	//if(case3Possible)
		//return min(caseThreeSnapOff, caseTwoSnapOff);
	//else if(case4Possible)
		//return min(caseFourSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
	//else if(oilLayer_[1].freeAtPrs(caseTwoSnapOff))
		//return min(caseTwoSnapOff, caseOneSnapOff);
	//else
		//return caseOneSnapOff;
	if(case3Possible && case4Possible)
		return min(caseFourSnapOff, min(caseThreeSnapOff, min(caseTwoSnapOff, caseOneSnapOff)));
	else if(case3Possible)
		return min(caseThreeSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
	else if(case4Possible)
		return min(caseFourSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
	else if(oilLayer_[1].exists() && cntAngRec_ > PI/2.0 + crnHafAngs_[0] )  //|| trappedCL_.first < 0;
		return min(caseTwoSnapOff, caseOneSnapOff);
	else
		return caseOneSnapOff;
}

double Square::snapOffPrsDrainSpontNR() const
{
		return (comn_.sigmaOW()/R_) * (cos(cntAngRec_) + sin(cntAngRec_));
}



