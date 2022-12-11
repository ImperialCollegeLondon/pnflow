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

#include "Element.h"
#include "elem_Model.h"
#include "polygon.h"

#include "cornerApex.h"
#include "layerApex.h"



using namespace std;





void LayerApex::getCAApexDistUntraped(double& apxDist, double& conAng, const double& halfAng, double pc, double intfacTen, bool itr) const
{
	double delta = itr ? 0.: SMALL_NUM;
	if(pc > advancingPc_-delta && pc < receedingPc_+delta)  {
		conAng = (acos(pc*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		conAng = std::max(std::min(PI, conAng), 0.);
		apxDist = initedOLApexDist_;
	}
	else
	{
		apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
		if(itr && apxDist <=  0.) apxDist=1.0E-15;
		ensure(apxDist > 0.);
		ensure(conAng >= 0. && conAng <=  PI);
	}
}


/*double LayerApex::hingingConAngUntraped(double pc, double conAng, double halfAng, double intfacTen, bool accurat) const
{
	double delta = accurat ? 0.: SMALL_NUM;
	if(pc > advancingPc_-delta && pc < receedingPc_+delta)  {
		double hingAng(acos(pc*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		hingAng = std::max(std::min(PI, hingAng), 0.);
		return hingAng;
	}

	ensure(conAng >= 0. && conAng <=  PI, "k");
	return conAng;
}

double LayerApex::getApexDistanceUntraped(double pc, double conAng, double halfAng, double intfacTen, bool itrRoutine) const
{
	if(pc > advancingPc_-SMALL_NUM && pc < receedingPc_+SMALL_NUM)
		return initedOLApexDist_;

	double apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);
	if(itrRoutine && apxDist <=  0.) return 1.0E-15;
	ensure(apxDist > 0., "iu");
	return apxDist;
}
*/



void LayerApex::getCAApexDist(double& apxDist, double& conAng, const double& halfAng, double pc, double intfacTen, bool itr) const
{

	ensure(inited_ || trappedCL_.first>-1 || parentShape_->eleman()->trappingWatBulk().first>-1);

	const double delta = itr ? 0.: SMALL_NUM;
	if(trappedCL_.first>-1)  {
		conAng = (acos(trappedCL_.second*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		conAng = std::max(std::min(PI, conAng), 0.);
		apxDist = initedOLApexDist_;
		//ensure(conAng >= 0. && conAng <=  PI, "j");
		return;
	}
	else if (parentShape_->eleman()->trappingWatBulk().first>-1)  {
		conAng = (acos(parentShape_->eleman()->trappingWatBulk().second*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		conAng = std::max(std::min(PI, conAng), 0.);
		//ensure(conAng >= 0. && conAng <=  PI, "j");
		apxDist = initedOLApexDist_;
		return;
	}
	else if (pc > advancingPc_-delta && pc < receedingPc_+delta)  {
		conAng  = (acos(pc*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		conAng = std::max(std::min(PI, conAng), 0.);
		//ensure(conAng >= 0. && conAng <=  PI, "j");
		apxDist = initedOLApexDist_;
		return;
	}
	else
	{
		ensure(conAng >= 0. && conAng <=  PI);

		apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);


		if(itr && apxDist <=  0.) apxDist = 1.0E-15;
		ensure(apxDist > 0.);
		if(apxDist < 0.)  {
			cout<< apxDist<<"  "<< pc<<"  "<< receedingPc_<<"  "<< advancingPc_<<"  "<<conAng<<"  "<<halfAng<<"  "<<intfacTen<<"  "<< itr<<endl;
		}
	}






}


/*double LayerApex::hingingConAng(double pc, double conAng, double halfAng, double intfacTen, bool accurat) const
{
	const double delta = accurat ? 0.: SMALL_NUM;
	if(trappedCL_.first>-1)  {
		double hingAng(acos(trappedCL_.second*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		hingAng = std::max(std::min(PI, hingAng), 0.);
		//ensure(hingAng >= 0. && hingAng <=  PI, "j");
		return hingAng;
	}
	else if (parentShape_->eleman()->trappingWatBulk().first>-1)  {
		double hingAng(acos(parentShape_->eleman()->trappingWatBulk().second*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		hingAng = std::max(std::min(PI, hingAng), 0.);
		//ensure(hingAng >= 0. && hingAng <=  PI, "j");
		return hingAng;
	}
	else if (pc > advancingPc_-delta && pc < receedingPc_+delta)  {
		double hingAng(acos(pc*initedOLApexDist_*sin(halfAng)/intfacTen)+halfAng);
		hingAng = std::max(std::min(PI, hingAng), 0.);
		//ensure(hingAng >= 0. && hingAng <=  PI, "j");
		return hingAng;
	}
	else
	{
		ensure(conAng >= 0. && conAng <=  PI, "k");
		return conAng;
	}
}

double LayerApex::getApexDistance(double pc, double conAng, double halfAng, double intfacTen, bool itrRoutine) const
{
	ensure(inited_ || trappedCL_.first>-1 || parentShape_->eleman()->trappingWatBulk().first>-1, " qljt ");

	if( trappedCL_.first>-1 || parentShape_->eleman()->trappingWatBulk().first>-1 || (pc > advancingPc_-SMALL_NUM && pc < receedingPc_+SMALL_NUM) )
		return initedOLApexDist_;

	double apxDist = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);


	if(itrRoutine && apxDist <=  0.) return 1.0E-15;
	ensure(apxDist > 0., "i");
	if(apxDist < 0.)  {
		cout<< apxDist<<"  "<< pc<<"  "<< receedingPc_<<"  "<< advancingPc_<<"  "<<conAng<<"  "<<halfAng<<"  "<<intfacTen<<"  "<< itrRoutine<<endl;
	}
	return apxDist;
}
*/








//bool LayerApex::create Film(double pc, double conAng, double maxSpontConAng, double halfAng, double intfacTen)
bool LayerApex::createOLayer(double pc, double conAngRec, double conAngAdv, double maxSpontConAng, double halfAng, double intfacTen, bool oilInj)  {
	ensure(!exists_);
	double conAng = oilInj ? conAngRec: conAngAdv;
	exists_ = conAng > PI/2. + halfAng  && conAng > maxSpontConAng;

	gravityCorrection_ = parentShape_->eleman()->gravityCorrection();


	if(!exists()) return false;


	initedOLApexDist_ = (intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);


	if(initedOLApexDist_ < innerCornerApex_->pinnedApexDist())  {
		 initedOLApexDist_= -1.; ///. used in re-creation
		 exists_ = false;
		return false;
	}

	ensure(initedOLApexDist_>0);


	advancingPc_ = intfacTen*cos(conAngAdv-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set
	receedingPc_ = intfacTen*cos(conAngRec-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set

	inited_ = true; ///. warning this is used by lay erCollapsePc_ fromEitherSide


		entryPc_ = layerCollapsePc(pc, conAng, halfAng, intfacTen,oilInj);

	exists_ = pc > entryPc_;
	if (!exists_) { inited_=false; initedOLApexDist_= -1.; }
	ensure(initedOLApexDist_>0);





	return exists_;
}

bool LayerApex::finitLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool overwriteTrapping)  {   ///. called from initOilInjection/initWaterInjection
	if(!exists()) return true;
	ensure(inited_ ||  trappedCL_.first > -1 || parentShape_->eleman()->trappingWatBulk().first>-1); //,"opin"
	ensure(!parentShape_->containCOil());//,"finx"

	double conAng = oilInj ? conAngRec : conAngAdv;

	int nTraps= parentShape_->eleman()->isTrappedOil() + (parentShape_->eleman()->trappingWatBulk().first>-1);
	if( inited_ || /*(overwriteTrapping && nTraps  < 2) ||*/ nTraps < 1 )  {

		///. initedOLApexDist_ =
		getCAApexDist(initedOLApexDist_, conAng, halfAng, pc, intfacTen);//const

		ensure(initedOLApexDist_ > innerCornerApex_->pinnedApexDist());
		inited_ = false;


	}

	advancingPc_ = intfacTen*cos(conAngAdv-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set
	receedingPc_ = intfacTen*cos(conAngRec-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set

	if(initedOLApexDist_ < innerCornerApex_->pinnedApexDist())  {
		 cout<<" fdraw "; cout.flush();
		 cout<<"  "<<initedOLApexDist_<<"   " << innerCornerApex_->pinnedApexDist(); cout.flush();
		 //removeLayer();
		return false;
	}




  return true;
}



bool LayerApex::initLayerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool silent)  {   ///. called from initOilInjection/initWaterInjection


	if( (!oilInj && !exists()) || trappedCL_.first > -1 || parentShape_->eleman()->trappingWatBulk().first>-1 )  {
		return true;
	}
	else if (oilInj && !exists())  { ///. calc layer creation Pc
		if (conAngRec <= PI/2. + halfAng  || initedOLApexDist_ <= 0. ) return true; ///. will not form

		inited_ = true;

		initedOLApexDist_ = parentShape_->RRR()/tan(halfAng);//(intfacTen/pc)*cos(conAng-halfAng)/sin(halfAng);

		advancingPc_ = intfacTen*cos(conAngAdv-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set
		receedingPc_ = intfacTen*cos(conAngRec-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set

		entryPc_ = layerCollapsePc(pc, conAngRec, halfAng, intfacTen, oilInj); ///.  set
		inited_ = false;

		return true;
	}
	else
	{


		//int WarningReInitLayers;


		bool stable = true;


		if(initedOLApexDist_ < innerCornerApex_->pinnedApexDist())  {
			 cout<<"\nWarning118   ";cout.flush();
			//removeLayer();

			stable = false;
		}


		advancingPc_ = intfacTen*cos(conAngAdv-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set
		receedingPc_ = intfacTen*cos(conAngRec-halfAng)/(initedOLApexDist_*sin(halfAng)); ///.  set
		inited_ = true;

		double currConAng = oilInj ? conAngRec : conAngAdv;
		entryPc_ = layerCollapsePc(pc, currConAng, halfAng, intfacTen, oilInj); ///.  set

		if(pc < entryPc_)  {
			 cout<<"\nWarning168   ";cout.flush();
			 //removeLayer();
			  stable = false;
		}

		return stable;
	}
}




double LayerApex::layerCollapsePc_FromCentre(double outPc, double inPc, double conAng, double halfAng, double intfacTen) const
{
	ensure(trappedCL_.first<0|| parentShape_->eleman()->trappingWatBulk().first<0);

	double inerApexDist, innerConAng(parentShape_->minInitRecCntAng()); innerCornerApex_->getCApexDistConAng(inerApexDist, innerConAng, inPc, halfAng, intfacTen);
	ensure(innerConAng != conAng);
	double ki = (intfacTen/inPc)*(cos(innerConAng)/sin(halfAng)-1.);

	if(inited_)  {
		double slipPrs = advancingPc_;
		double ko = (intfacTen/slipPrs)*(cos(advConAng_)/sin(halfAng)+1.);
		if(ko - ki < NEG_ALMOST_ZERO)  {
			double oldRad(intfacTen/(slipPrs+SMALL_NUM));    // Don't get caught by rounding errors
			for(int itr = 0; itr < MAX_ITR; ++itr)  {
				double outConAng = acos(initedOLApexDist_*sin(halfAng)/oldRad)+halfAng;
				double outConAngDr = initedOLApexDist_*sin(halfAng)/(sin(outConAng-halfAng)*oldRad*oldRad);
				double func = oldRad*(cos(outConAng)/sin(halfAng)+1.)-ki;
				double funcDr = (cos(outConAng)-oldRad*outConAngDr*sin(outConAng))/sin(halfAng) + 1.;
				double newRad = oldRad - func/funcDr;
				if(fabs((newRad-oldRad)/newRad) < EPSILON)  {
					colType_ =30;

					double pc(intfacTen/newRad);
					ensure(fabs(ki-newRad*(cos(outConAng)/sin(halfAng)+1.)) < EPSILON);
					return pc;
				}
				oldRad = newRad;
			}

			std::cerr << std::endl
				<< "============================================"   << std::endl
				<< "Error: failed to converge to a solution for"    << std::endl
				<< "layer collapsing pressure.   "                     << std::endl
				<< "============================================"   << std::endl;    exit(-1);
			return 0.;

		}
	}
	else if(exists_) {cout<<"b";cout.flush(); }//parentShape_->ChParent()[1000000].setInOilFloodVec(true);

	colType_ =31;


	//double inerApexDist = innerCornerApex_->getCApexDistance(inPc, innerConAng, halfAng, intfacTen);
	double collPc1 = intfacTen*cos(advConAng_-halfAng)/(innerCornerApex_->pinnedApexDist()*sin(halfAng));

	double collPc2 = intfacTen*(cos(advConAng_)+sin(halfAng))/(ki*sin(halfAng));
	double collPc = max(collPc1, collPc2);
	ensure(collPc < 0.);
	if (collPc > 0.) cout<<"CA"<<advConAng_<<" "<<conAng<<" ";
	return collPc;
}

double LayerApex::layerCollapsePc_FromCorner(double outPc, double inPc, double conAng, double halfAng, double intfacTen) const
{
	ensure(trappedCL_.first<0);



	double outerConAng(conAng), apxDistDummy;// = hingingConAng(outPc, conAng, halfAng, intfacTen);
	getCAApexDist(apxDistDummy, outerConAng, halfAng, outPc, intfacTen);
	double ko = (intfacTen/outPc)*(cos(outerConAng)/sin(halfAng)+1.);
	 //inited_ ? acos(initedOLApexDist_*sin(halfAng)*outPc/intfacTen)+halfAng: conAng;
	double innerSlipPrs = innerCornerApex_->advancingPc();
	double ki = (intfacTen/innerSlipPrs)*(cos(advConAng_)/sin(halfAng)-1.);

	if(ko - ki < NEG_ALMOST_ZERO  && inPc < innerCornerApex_->advancingPc())  {
		double oldRad(intfacTen/(innerSlipPrs+SMALL_NUM));
		double bi = innerCornerApex_->pinnedApexDist();
		for(int itr = 0; itr < MAX_ITR; ++itr)  {
			ensure((innerCornerApex_->pinnedApexDist()*sin(halfAng)/oldRad) >= -1. &&
				(innerCornerApex_->pinnedApexDist()*sin(halfAng)/oldRad) <=  1.);

			double inConAng = acos(bi*sin(halfAng)/oldRad)-halfAng;
			double inConAngDr = bi*sin(halfAng)/(sin(inConAng+halfAng)*oldRad*oldRad);
			double func = oldRad*(cos(inConAng)/sin(halfAng)-1.)-ko;
			double funcDr = (cos(inConAng)-oldRad*inConAngDr*sin(inConAng))/sin(halfAng) - 1.;

			double newRad = oldRad - func/funcDr;
			if(fabs((newRad-oldRad)/newRad) < EPSILON)  {				colType_ =10;

				double pc(intfacTen/newRad);
				ensure(fabs(ko-(intfacTen/pc)*(cos(inConAng)/sin(halfAng)-1.)) < EPSILON);//, "outTrap"
				return pc;
			}
			oldRad = newRad;
		}
		//return Error_layer CollapsePc();

			std::cerr << std::endl
				<< "============================================"   << std::endl
				<< "Error: could not converge to a solution for"    << std::endl
				<< "layer collapsing pressure from corner."                     << std::endl
				<< outPc<<"  "<< inPc << "  " << advancingPc_                << std::endl
				<< bi << "  "                << std::endl
				<< conAng << "  " << advConAng_                << std::endl
				<< "  " << halfAng                << std::endl
				<< innerCornerApex_->advancingPc() << "  " << innerCornerApex_->receedingPc()                 << std::endl
				//<< innerCornerApex_->trappingCL().first << "  " << innerCornerApex_->getCApexDistance(inPc, conAng, halfAng, intfacTen)      << std::endl
				<< "============================================"   << std::endl;    exit(-1);

		//int warning;

		std::cerr << std::endl
			<< "============================================"   << std::endl
			<< "Error: Failed to converge to a solution     "    << std::endl
			<< "for layer collapsing pressure.      "                     << std::endl
			<< "============================================"   << std::endl;    exit(-1);
		return 0.;
	}

	colType_ =11;
	return innerSlipPrs;
}





inline double dimLessCornerArea(double halfAng, double contactAng)
{
	if(fabs(contactAng + halfAng - ElemModel::PI/2.) < 0.01)  {
		return sin(halfAng)*cos(halfAng);
	}
	else
	{
		return pow(sin(halfAng)/cos(halfAng+contactAng), 2.)
		* (cos(contactAng)*cos(halfAng+contactAng)/sin(halfAng) + halfAng + contactAng - ElemModel::PI/2.);
	}
}

double LayerApex::layerCollapsePc_fromEitherSide(double pc, double conAng, double halfAng, double intfacTen,bool debug) const  {
	ensure(trappedCL_.first<0);

//if (debug)
//{
	if((exists_ && !inited_) || (!exists_ && inited_))  {
	ensure(inited_);
	//ensure(exists_,"kfh");
		//parentShape_->ChParent()[1000000].setInOilFloodVec(true);
	}
//}
	//if(inited_)
	//{
	double slipPrs = inited_ ? advancingPc_ : pc;
	double innerSlipPrs = innerCornerApex_->cornerExists() ? innerCornerApex_->advancingPc() : pc;
	ensure(slipPrs < 0.);
	//double innerHingConAng = innerCornerApex_->cornerExists() ? innerCornerApex_->getCA pexDistConAng(innerSlipPrs, conAng, halfAng, intfacTen) : parentShape_->minInitRecCntAng();
	double bi, innerHingConAng(parentShape_->minInitRecCntAng());
	if (innerCornerApex_->cornerExists()) innerCornerApex_->getCApexDistConAng(bi, innerHingConAng, innerSlipPrs, halfAng, intfacTen,true);


	double /*ki*/innerApexDistThroughCentre = (intfacTen/innerSlipPrs)*(cos(innerHingConAng)/sin(halfAng)-1.);
	double ko = (intfacTen/slipPrs)*(cos(advConAng_)/sin(halfAng)+1.);
		ensure(ko>0);
	if(ko - innerApexDistThroughCentre < NEG_ALMOST_ZERO   && pc < innerSlipPrs +0.1)  {
		double oldRad(intfacTen/(slipPrs+SMALL_NUM));

		//innerCornerApex_->getCApexDistConAng(bi, innerHingConAng, innerSlipPrs, halfAng, intfacTen,true);

		double bo = initedOLApexDist_;
		for(int itr = 0; itr < MAX_ITR; ++itr)  {
			double inConAng = acos(bi*sin(halfAng)/oldRad)-halfAng;
			double outConAng = acos(bo*sin(halfAng)/oldRad)+halfAng;
			double inConAngDr = bi*sin(halfAng)/(sin(inConAng+halfAng)*oldRad*oldRad);
			double outConAngDr = bo*sin(halfAng)/(sin(outConAng-halfAng)*oldRad*oldRad);

			double func = (cos(inConAng)-cos(outConAng))/sin(halfAng)-2.;
			double funcDr = (outConAngDr*sin(outConAng)-inConAngDr*sin(inConAng))/sin(halfAng);
			double newRad = oldRad - func/funcDr;
			if(fabs((newRad-oldRad)/newRad) < EPSILON)  {
				ensure(fabs(newRad*(cos(innerHingConAng)/sin(halfAng)-1.) -
					newRad*(cos(advConAng_)/sin(halfAng)+1.)) < EPSILON);
				ensure(intfacTen/newRad < 0.);


			colType_ =20;

				return intfacTen/newRad;
			}
			oldRad = newRad;
		}
		//return Error _layer CollapsePc();


		std::cerr << std::endl
			<< "============================================"   << std::endl
			<< "Error: could not converge to a solution for"    << std::endl
			<< "layer collapsing pressure.         "     << std::endl
			<< pc << "  " << advancingPc_                << std::endl
			<< bi << "  " << bo                << std::endl
			<< conAng << "  " << advConAng_                << std::endl
			<< innerHingConAng << "  " << halfAng                << std::endl
			<< innerCornerApex_->advancingPc() << "  " << innerCornerApex_->receedingPc()                 << std::endl
			<< innerCornerApex_->trappingCL().first << "  " << bi             << std::endl
			<< "============================================"   << std::endl;    exit(-1);





		return 0.;
	}
	//else if(pc >= innerSlipPrs)
	{
		double collPc1;


		innerApexDistThroughCentre=/*innerCornerApex_->getCApexDistanceUntraped(pc, innerHingConAng, halfAng, intfacTen)**/
		bi*(cos(halfAng)+sin(halfAng)*sin(0.5*(innerHingConAng+halfAng-PI/2.))/cos(0.5*(innerHingConAng+halfAng-PI/2.)));
		//double outertouchApexDist=innerApexDistThroughCentre/(cos(halfAng)+sin(halfAng)*cos(PI-conAng+halfAng)/sin(PI-conAng+halfAng))WRONG;




		double outertouchApexRad=innerApexDistThroughCentre/(cos(advConAng_)/sin(halfAng)+1.);
		collPc1= intfacTen/outertouchApexRad;//(outertouchApexDist*sin(halfAng));






		double pinnedOLApexDist = outertouchApexRad*cos(advConAng_-halfAng)/(sin(halfAng));

		ensure (dimLessCornerArea(halfAng,PI-advConAng_)*pinnedOLApexDist*pinnedOLApexDist > dimLessCornerArea(halfAng,innerHingConAng)*bi*bi-1e-18);

		if (debug || !(dimLessCornerArea(halfAng,PI-advConAng_)*pinnedOLApexDist*pinnedOLApexDist > dimLessCornerArea(halfAng,innerHingConAng)*bi*bi-1e-18))  {
			cout<< "  layerCollapsePc_fromEitherSide:"<<"  inCorExists:"<< innerCornerApex_->cornerExists()<<"  hAng:"<<halfAng<<"  CAAdv"<<advConAng_<<"  ";
			//exit(-1);
		}
		///warningtocheck;

		collPc1=max(collPc1,innerCornerApex_->advancingPc());


			colType_ = 21;


		double collPc2 = (innerCornerApex_->cornerExists()) ? intfacTen*cos(advConAng_-halfAng)/(bi*sin(halfAng)) : collPc1-0.1;





		if (debug)  {

			//double outertouchApexDist=outertouchApexRad*cos(conAng-halfAng)/sin(halfAng);
	double outerApexDistThroughCentre = outertouchApexRad*(cos(advConAng_)/sin(halfAng)+1.);

			cout<<colType_<<"\n  innerHingConAng "<<innerHingConAng<<" "//<<innerCornerApex_->getCApex DistConAng(innerSlipPrs, conAng, halfAng, intfacTen)<<" "
				<<"  advConAng_ "<<advConAng_<<"  conAng"<<conAng
				<<"  halfAng "<<halfAng<<"  outertouchApexRad"<<outertouchApexRad<<"   pc " <<pc<<"  inSlipPc"<<innerSlipPrs<<endl
				<<"    innerApexDistThroughCentre "<<innerApexDistThroughCentre<<"    koColapse "<<outerApexDistThroughCentre
				<<"    inApexDist "//<<innerCornerApex_->getCApexDistanceUntraped(innerSlipPrs, innerHingConAng, halfAng, intfacTen)<<"  outertouchApexDist"<<outertouchApexDist
				<<" pc:"<<pc<<" pci:"<<(innerCornerApex_->cornerExists() ? innerCornerApex_->advancingPc() : slipPrs)<<" pco:"<<slipPrs<<" Pc1:"<<collPc1<<" Pc2:"<<collPc2<<" max:"<<max(collPc1,collPc2)<<" \n";
		//exit(-1);
		}

		ensure(collPc2 < 0.);
		if (collPc1<collPc2)  {
			colType_ = 22;
			//cout<<"colPcvk"<<collPc2<<"  "<<conAng<<"  "<<advConAng_<<"  "<<halfAng<<endl;
			return collPc2;
		}

	return collPc1;
	}


}





double LayerApex::layerCollapsePc(double pc, double conAng, double halfAng, double intfacTen, bool injOil) const  {
	double collPrs;
	colType_=0;
	if(!innerCornerApex_->cornerExists())  {
		//ensure(exists(), "quacki di quack");
		///WarningFiEnableensure;
		collPrs = LOWEST_LAYER_PC;                      // No water in corner  => never to collapse
	}
	else if(injOil && exists())                      // Oil Inj, stable layer  => no update needed
		collPrs = INF_NEG_NUM;
	else if(!injOil && !exists())               // Wat Inj, collapsed layer  => no update needed
		collPrs = INF_POS_NUM;
	else if(trappedCL_.first>-1)                          // Wat Inj, both trapped  => never collapse
		collPrs = INF_NEG_NUM;
	else
	{
		bool freeWBulk=parentShape_->eleman()->trappingWatBulk().first<0;
		bool freeWFilm=innerCornerApex_->cornerExists() && innerCornerApex_->trappedCorner().first<0;



		///Warning_unsynced_Trappings;
		 //ensure(parentShape_->eleman()->trappingWatFilm().first<0 && innerCornerApex_->trappedCorner().first>-1, " Fq1");
		 //ensure(parentShape_->eleman()->trappingWatFilm().first>-1 && innerCornerApex_->trappedCorner().first<0, " Fq2");
//
	//if(!(parentShape_->eleman()->trappingWatFilm().first<0 && innerCornerApex_->trappedCorner().first>-1))
	//{
		//conAng=parentShape_[10000].conAngleRec();
	//}




		if(freeWBulk && freeWFilm) // No trapping
			collPrs = layerCollapsePc_fromEitherSide(pc, conAng, halfAng, intfacTen);
		else if(freeWBulk)  ///. Only Corner trapped
			collPrs = layerCollapsePc_FromCentre(pc, innerCornerApex_->trappedCornerNOTTOBEUSED().second, conAng, halfAng, intfacTen);
		else if(freeWFilm)  ///. Only Centre trapped
			collPrs = layerCollapsePc_FromCorner(trappedCL_.second, pc, conAng, halfAng, intfacTen);
		else if(injOil)                 // Oil Inj, both trapped  => never grow
			collPrs = INF_POS_NUM;
		else                            // Wat Inj, both trapped  => never collapse
		   collPrs = INF_NEG_NUM;
	}

	//if(!innerCornerApex_->cornerExists())
	//{
		//ensure(exists(), "quacki di quack");
		//collPrs = LOWEST_LAYER_PC;                      // No water in corner  => never to collapse
	//}
	//else if(injOil && exists())                      // Oil Inj, stable layer  => no update needed
		//collPrs = INF_NEG_NUM;
	//else if(!injOil && !exists())               // Wat Inj, collapsed layer  => no update needed
		//collPrs = INF_POS_NUM;
	//else if(innerCornerApex_->trappedCornerNOTTOBEUSED().first<0 && trappedCL_.first<0) // No trapping
		//collPrs = layerCollapsePc_fromEitherSide(pc, conAng, halfAng, intfacTen);
	//else if(innerCornerApex_->trappedCornerNOTTOBEUSED().first>-1 && trappedCL_.first<0)  ///. Only Corner trapped
		//collPrs = layerCollapsePc_FromCentre(pc, innerCornerApex_->trappedCornerNOTTOBEUSED().second, conAng, halfAng, intfacTen);
	//// else if(innerCornerApex_->trappedCornerNOTTOBEUSED().first<0 && trappedCL_.first>-1)  ///. Only Centre trapped
		//// collPrs = layerCollapsePc_FromCorner(trappedCL_.second, pc, conAng, halfAng, intfacTen);
	//else if(injOil)                 // Oil Inj, both trapped  => never grow
		//collPrs = INF_POS_NUM;
	//else                            // Wat Inj, both trapped  => never collapse
	   //collPrs = INF_NEG_NUM;

	return collPrs;
}
