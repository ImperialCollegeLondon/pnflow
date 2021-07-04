
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "sortedEvents.h"
#include "compareFuncs.h"


using namespace std;


/**
* Before invoking the solver and calculating saturations we must update the state of the
* element, ie calculate saturation and conductance. Since the configuration in imbibition can
* be more complex than in drainage, the calculations are separate for each event
*/
double Polygon::calcR(double pc)
{

		//sumOilLayerAreas_ = 0.;
	if(bulkFluid_ == &comn_.water() && numLayers_ == 0)
	{
		SatWater_ = 1.;
		conductanceOil_ = 0.;
		conductanceWater_.first = 0.;
		conductanceWater_.second = SPConductance(area_, comn_.water().viscosity());
		

	}
	else if(bulkFluid_ == &comn_.oil() && !waterInCorner_[0].cornerExists())
	{
		SatWater_ = 0.;
		conductanceOil_ = SPConductance(area_, comn_.oil().viscosity());
		conductanceWater_.first = comn_.KrwatcornAtSw0() * conductanceOil_; ///. KrwatcornAtSw0_ is read from TRAPPING keyword
		conductanceWater_.second = 0.;


	}
	else if(bulkFluid_ == &comn_.oil() && waterInCorner_[0].cornerExists())
	{
		calcR_oilWithWaterInCorners(pc-eleman()->gravityCorrection()-elem_.snapOfLongitCurvature()*comn_.sigmaOW());


	}
	else
	{
		calcR_waterWithOilLayers(pc-eleman()->gravityCorrection()-elem_.snapOfLongitCurvature()*comn_.sigmaOW());
	}

	ensure(conductanceOil_ >= 0. && (conductanceWater_.first >= 0. || conductanceWater_.second >= 0.));
	if ( !(conductanceOil_ >= 0. && (conductanceWater_.first >= 0. || conductanceWater_.second >= 0.) ) )
		cout<<"Ks = "<< conductanceOil_ << "  " << conductanceWater_.first << "  " << conductanceWater_.second <<endl;
	///. TODO: sensitivity analysis, add non-linearity by decoupling Aw_dl from Sw
	

	ElectricalConductance_ = SatWater_*area_/comn_.water().resistivity() + sqrt(SatWater_*area_/shapeFactor_) / comn_.surface().resistivity();
	ElectricalConductance_ += area_*(1.-SatWater_)/comn_.oil().resistivity();
	ensure(ElectricalConductance_<1e64);





	if (debugLevel > 500 && comn_.injectant() == &comn_.water() && ( conductanceWater_.first +conductanceWater_.second  <  conductanceWaterOldToCheck_.first +conductanceWaterOldToCheck_.second  || conductanceWater_.second < conductanceWaterOldToCheck_.second) )
	{
		cout<< SatWater_ <<"  "<<   elem_.trappingWatBulk().first <<"  "<< conductanceWater_.second <<"  "<< conductanceWaterOldToCheck_.second<<"     ";
		cout<< conductanceWater_.first <<"  "<<   conductanceWaterOldToCheck_.first <<"  "<< conductanceWater_.second <<"  "<< conductanceWaterOldToCheck_.second
		<<"  kw\\"<<endl; //exit(-1);
	}
	conductanceWaterOldToCheck_.first = conductanceWater_.first;
	conductanceWaterOldToCheck_.second = conductanceWater_.second;


	return SatWater_; ///. saturation calculation
}







inline double dimLessCornerArea(double halfAng, double contactAng) 
{
	if(fabs(contactAng + halfAng - ElemModel::PI/2.) < 0.01)
	{
		return sin(halfAng)*cos(halfAng);
	}
	else
	{
		return pow(sin(halfAng)/cos(halfAng+contactAng), 2.)
		* (cos(contactAng)*cos(halfAng+contactAng)/sin(halfAng) + halfAng + contactAng - ElemModel::PI/2.);
	}
	//{
		//double cosTerm=cos(halfAng+contactAng);
		//return pow(sin(halfAng)/cos(halfAng+contactAng), 2.)
		//* (cos(halfAng+contactAng)*sin(halfAng+contactAng)+cos(halfAng+contactAng)*cos(halfAng+contactAng)/tan(halfAng) + halfAng + contactAng - ElemModel::PI/2.);
	//}
	//{
		//double cosTerm=cos(halfAng+contactAng);
		//return pow(sin(halfAng)/cosTerm, 2.)
		//* (cosTerm*sin(halfAng+contactAng)+cosTerm*cosTerm/tan(halfAng) + halfAng + contactAng - ElemModel::PI/2.);
	//}
}


/**
// The corner conductance version by Paal-Eric
*/
double cornerConductance(double dimLessCornerA, double meniscusApexDist, double halfAng,
								  double contactAngle, double viscosity) 
{
	double cornerGstar((cos(halfAng)*sin(halfAng)) / (4.*pow(1+sin(halfAng), 2.)));
	double cornerG(cornerGstar);

	if(fabs(contactAngle + halfAng - ElemModel::PI/2.) > 0.01)
	{
		cornerG = dimLessCornerA /
		(4. * pow(1. - (sin(halfAng)/cos(halfAng+contactAngle)) * (halfAng+contactAngle-ElemModel::PI/2.), 2.));
	}


	double cFactor = 0.364 + 0.28 * cornerGstar / cornerG;
	double dimLessCornerConduct = cFactor * dimLessCornerA * dimLessCornerA * cornerG;

	ensure(dimLessCornerA != 0. && halfAng != 0. && cornerG != 0.);

	return dimLessCornerConduct * pow(meniscusApexDist, 4.) / viscosity;
}

// My very own layer conductance routine. Wohoo :-)
double layerConductance(double insideBC, double outsideBC, double layerArea, double cornerArea, double outerApexDist,
						 double innerCornerApexDist, double innerConAng, double outerConAng, double halfAng, double viscosity)
{
	double coeff[3];
	coeff[0] = -2.4010e-002;
	coeff[1] = 2.8402e-001;
	coeff[2] = -2.9531;

	if(halfAng < 10.*ElemModel::PI/180.)
	{
		coeff[0] = -1.0610E-002;
		coeff[1] = 5.1608E-001;
		coeff[2] = -2.0645;
	}
	else if(halfAng < 20.*ElemModel::PI/180.)
	{
		coeff[0] = -2.6810E-002;
		coeff[1] = 1.8672E-001;
		coeff[2] = -3.5977;
	}
	else if(halfAng < 30.*ElemModel::PI/180.)
	{
		coeff[0] = -4.4021E-002;
		coeff[1] = -6.3195E-002;
		coeff[2] = -4.3749;
	}
	else if(halfAng < 40.*ElemModel::PI/180.)
	{
		coeff[0] = -3.1523E-002;
		coeff[1] = 1.6948E-001;
		coeff[2] = -3.3600;
	}
	else if(halfAng < 50.*ElemModel::PI/180.)
	{
		coeff[0] = -3.1367E-002;
		coeff[1] = 1.9327E-001;
		coeff[2] = -3.2673;
	}
	else if(halfAng < 60.*ElemModel::PI/180.)
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

	double outerOWlen = - (2.*sin(halfAng)*(halfAng-outerConAng+ElemModel::PI/2.))
		/ cos(ElemModel::PI-outerConAng+halfAng);

	double innerOWLen = - (2.*dimlessBi*sin(halfAng)*(halfAng+innerConAng-ElemModel::PI/2.))
		/ cos(innerConAng+halfAng);

	double perimeterLen = outerOWlen + innerOWLen + 2.*(1.-dimlessBi);
	double shapeFactorOil = dimlessAo / (perimeterLen*perimeterLen);
	double xGroup = log(pow(dimlessAo, 3.)*shapeFactorOil);

	double dimlessOilCond = exp(coeff[0]*xGroup*xGroup + coeff[1]*xGroup + coeff[2]);
	double dimCond = pow(outerApexDist, 4.) * dimlessOilCond / viscosity;
	return dimCond;
}



/**
* Water will remain in some of the corners of the polygon, Conductance and area is calculated
* for each corner.
*/
void Polygon::calcR_oilWithWaterInCorners(double cappPressure)
{
	double conAng(0.), cornerAreas(0.), cornerCond(0.);
	//bool overRideTrapping(false);
	//double pcSuggested = cappPressure;
	if(comn_.injectant() == &comn_.water())
	{
		conAng = cntAngAdv_;
	}
	else
	{
		conAng = cntAngRec_;
	}


	for(int i = 0; i < numCorners_; ++i)
	{
		double apexDist(0.), conAngCurr(conAng); 
		if(waterInCorner_[i].cornerExists())
		{

			waterInCorner_[i].getCApexDistConAng(apexDist, conAngCurr, cappPressure, crnHafAngs_[i], comn_.sigmaOW());
			if (apexDist >RRR()/tan(crnHafAngs_[0])*2) 	{	cout<<" sxp ";	}
			
			if(debugLevel>1 &&(apexDist < 0.))         cout<<" apexD"<<apexDist<<" ";

			double dlCornerArea = dimLessCornerArea(crnHafAngs_[i], conAngCurr);
			double conductance = cornerConductance(dlCornerArea, apexDist, crnHafAngs_[i], conAngCurr, comn_.water().viscosity());

			ensure(conductance > 0.);
			ensure(dlCornerArea < 1.);
			cornerCond += conductance;
			cornerAreas += apexDist * apexDist * dlCornerArea;

		}


	}
	SatWater_ = min(cornerAreas/area_,1.);
	ensure( SatWater_ < 1. && cornerCond > 0.);
	
	if ( debugLevel>10 && !(cornerAreas/area_ < 1. && cornerCond > 0.))
	{	cout<< "============ !(cornerAreas/area_ < 1. && cornerCond > 0.  =================";
		exit(-1);
	}
	
	double oilSat(max(1.-SatWater_,0.));
	conductanceOil_ = SPConductance(area_, comn_.oil().viscosity()) * oilSat;
	conductanceWater_.first = cornerCond;
	conductanceWater_.second = 0.;
	


}


/**
* If oil is present at the same time that water occupies the center there must be
* oil layer sandwisched between water in corners and center.
*/
void Polygon::calcR_waterWithOilLayers(double cappPressure)
{
	double conAng(0.);
	double cornerAreaWat(0.), cornerCondWat(0.);
	double cornerAreaOil(0.), cornerCondOil(0.);
	double layerAreaOil(0.), layerCondOil(0.);
	//bool overRideTrapping(false);

	if(comn_.injectant() == &comn_.water())		
		conAng = cntAngAdv_;
	else	
	{
		conAng = cntAngRec_;
	}
	//if(cappPressure > snapOffPrs_ && oilLayer_[0].trappedOLayer().first<0)
		//cappPressure = snapOffPrs_-EPSILON;

	//double cappPcOL = cappPressure;
	//if( elem_.trappingWatBulk().first>-1 ) cappPcOL = elem_.trappingWatBulk().second;

	double baseLength(R_*(1./tan(crnHafAngs_[0])+1./tan(crnHafAngs_[1])));
	double oilApexBaseLengths(0.), watApexBaseLengths(0.);

	for(int i = 0; i < numCorners_; ++i)
	{

	  if(oilLayer_[i].LayerApex::exists(/*st ab le*/))
	  {
		
		if(waterInCorner_[i].cornerExists())
		{

			//double layerPC(cappPressure); ///. trapping pc handled in calcR _getIFDistAng
			//if(!oilLayer_[i].LayerApex::freeAtPrs(layerPC))
				//layerPC = oilLayer_[i].LayerApex::layerCollPc()+EPSILON;

			//double outerConAng(conAng), outerApexDist(0.);
			//oilLayer_[i].LayerApex::calcR_getLIFDistAng(outerConAng, outerApexDist, layerPC, crnHafAngs_[i], comn_.sigmaOW(), cntAngAdv_, cntAngRec_);
			double outerConAng(conAng);// = oilLayer_[i].hingingConAng(cappPressure,conAng,crnHafAngs_[i],  comn_.sigmaOW());
			double outerApexDist;// = oilLayer_[i].getApexDistance(cappPressure,outerConAng,crnHafAngs_[i],  comn_.sigmaOW());
			oilLayer_[i].getCAApexDist(outerApexDist, outerConAng,crnHafAngs_[i], cappPressure, comn_.sigmaOW());

			ensure(outerApexDist>0);

			ensure(outerApexDist < (R_*(1./tan(crnHafAngs_[0])+1./tan(crnHafAngs_[1]))));
			ensure(oilLayer_[i].pinnedApexDist() < (R_*(1./tan(crnHafAngs_[0])+1./tan(crnHafAngs_[1]))));

			double dlTotCornerArea = dimLessCornerArea(crnHafAngs_[i], PI-outerConAng);
			double areaCornerTot = outerApexDist * outerApexDist * dlTotCornerArea;

			oilApexBaseLengths += outerApexDist;    // A permanaent check for not overstepping oil volumes
			
			if(i < 2) 
			{ 
				ensure( oilApexBaseLengths <=  baseLength*1.01);
			}

			//double innerConAng(conAng), innerCornerApexDist(0.);
			//waterInCorner_[i].calcR_getIFDistAng(innerConAng, innerCornerApexDist, layerPC, crnHafAngs_[i],  comn_.sigmaOW(), cntAngAdv_, cntAngRec_);
			
			//double innerConAng = waterInCorner_[i].hingingCConAng(cappPressure,conAng,crnHafAngs_[i],  comn_.sigmaOW());
			//double innerCornerApexDist = waterInCorner_[i].getCApexDistance(cappPressure,innerConAng,crnHafAngs_[i],  comn_.sigmaOW());
				double innerCornerApexDist, innerConAng(conAng); 
				waterInCorner_[i].getCApexDistConAng(innerCornerApexDist, innerConAng, cappPressure, crnHafAngs_[i], comn_.sigmaOW());

			watApexBaseLengths+=innerCornerApexDist;
			
			double dlWatCornerArea = dimLessCornerArea(crnHafAngs_[i], innerConAng);
			double areaWater = innerCornerApexDist * innerCornerApexDist * dlWatCornerArea;
			ensure(outerApexDist > innerCornerApexDist);
			if (outerApexDist <=  innerCornerApexDist)
			{
				cout<<int(elem_.isTrappedWat(bulkBlob))<<"  ";
				cout<<" eXiT "<<endl;			//	exit(-1);
			}



			double areaOil = max(areaCornerTot - areaWater,1e-50);

			double waterCond = cornerConductance(dlWatCornerArea, innerCornerApexDist, crnHafAngs_[i],
				innerConAng, comn_.water().viscosity());
			double oilCond = layerConductance(1., 1., areaOil, areaCornerTot, outerApexDist, innerCornerApexDist,
				innerConAng, outerConAng, crnHafAngs_[i], comn_.oil().viscosity());

			if (debugLevel>1 && ((areaWater*waterCond*(areaCornerTot - areaWater)*oilCond <=  0.) || (oilApexBaseLengths >  baseLength*1.01 && i < 2)))
			{
				cout<<" if (debugLevel>1 && ((areaWater*waterCond*(areaCornerTot - areaWater)*oilCond <=  0.) || (oilApexBaseLengths >  baseLength*1.01 && i < 2)))   "<< oilLayer_[i].colType_<<" ";
				exit(-1);
			}
			
			ensure(waterCond > 0.);

			cornerAreaWat += areaWater;
			cornerCondWat += waterCond;
			layerAreaOil += areaOil;
			layerCondOil += oilCond;
		}
		else 
		{
			//double currentConAng(conAng), apexDist(0.);
			//oilLayer_[i].LayerApex::calcR_getLIFDistAng(currentConAng, apexDist, cappPressure, crnHafAngs_[i],  comn_.sigmaOW(), cntAngAdv_, cntAngRec_);
			double currentConAng(conAng);// = oilLayer_[i].hingingConAng(cappPressure,conAng,crnHafAngs_[i],  comn_.sigmaOW());
			double apexDist;// = oilLayer_[i].getApexDistance(cappPressure,currentConAng,crnHafAngs_[i],  comn_.sigmaOW());
			oilLayer_[i].getCAApexDist(apexDist,currentConAng,crnHafAngs_[i], cappPressure, comn_.sigmaOW());
			
			ensure(apexDist>0);


			double dlCornerArea = dimLessCornerArea(crnHafAngs_[i], PI-currentConAng);
			double conductance = cornerConductance(dlCornerArea, apexDist, crnHafAngs_[i],
				PI-currentConAng, comn_.oil().viscosity());

			cornerCondOil += conductance;
			cornerAreaOil += apexDist * apexDist * dlCornerArea;
		if ( debugLevel>10 && !(layerAreaOil >= 0.))
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

	ensure(layerAreaOil > 0. || cornerAreaOil > 0.);


	SatWater_ = 1. - (cornerAreaOil + layerAreaOil)/area_;
	conductanceOil_ = layerCondOil + cornerCondOil;

	double watSat((area_-cornerAreaOil-layerAreaOil-cornerAreaWat) / area_);

	double centerWatCond = SPConductance(area_, comn_.water().viscosity()) * watSat;
	dbgAsrt(watSat > 0. && watSat < 1.);

	if(hasDisConectedCentreWCornerW_)
	{
		ensure(cornerCondWat > 0.);
		
		conductanceWater_.first = cornerCondWat;
		conductanceWater_.second = centerWatCond;
	}
	else
	{
		conductanceWater_.first = 0.;
		conductanceWater_.second = centerWatCond+cornerCondWat;
	}
	
}




/**
* Circular elements can only hold a single fluid.
*/
double Circle::calcR(double pc)
{

	if(bulkFluid_ == &comn_.water())
	{
		SatWater_ = 1.;
		conductanceOil_ = 0.;
		conductanceWater_.first = 0.;
		conductanceWater_.second = SPConductance(area_, comn_.water().viscosity());
		ElectricalConductance_ = area_/comn_.water().resistivity() + sqrt(area_ / shapeFactor_) / comn_.surface().resistivity();
	}
	else
	{
		SatWater_ = 0.;
		conductanceOil_ = SPConductance(area_, comn_.oil().viscosity());
		conductanceWater_.first = comn_.KrwatcornAtSw0() * conductanceOil_;
		conductanceWater_.second = 0.;
		ElectricalConductance_ = area_/comn_.oil().resistivity();
	}

	ElectricalConductance_ = (SatWater_*area_)/comn_.water().resistivity() + sqrt(SatWater_*area_/shapeFactor_) / comn_.surface().resistivity();
	ElectricalConductance_ += area_*(1.-SatWater_)/comn_.oil().resistivity();

	return SatWater_;
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
	// return (R_ * R_ * area) / (7.1136 * visc);   // 1
	return (0.5623 * area * area * shapeFactor_) / visc;       // 2
}

/**
 * Computed the conductance of a fluid occupying the center of the element
*/
double Triangle::SPConductance(double area, double visc) const
{
	// return (3. * R_ * R_ * area) / (20. * visc);   // 1
	return (3. * area * area * shapeFactor_) / (5. * visc);      // 2
}

/**
// Computed the conductance of a fluid occupying the center of the element
*/
double Circle::SPConductance(double area, double visc) const
{
	return (R_ * R_ * area) / (8. * visc);
}







