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
using namespace std;


const double    Apex::PI = acos(-1.);
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



//! sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
void CornerApex::createFilm(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool isOilInj)  {
	ensure(!exists_);
	ensure(!inited_);

	if(exists_) return;

	double conAng = isOilInj ? conAngRec: conAngAdv;
	exists_ = conAng < PI/2. - halfAng;


	if(!exists_) return;



	initedApexDist_ = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng); // TODO: valgrind error uninitialized halfAng
	ensure(initedApexDist_ > 0.);
	if (initedApexDist_ > parentShape_->RRR()/tan(halfAng)*2) 	{ cout<<" xkg ";	}
	creationPc=pc;


	advancingPc_ = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(initedApexDist_*sin(halfAng));
	receedingPc_ = intfacTen*cos( min(PI, conAngRec)+halfAng )/(initedApexDist_*sin(halfAng));
	inited_ = true;


	if(pc > initOrMaxPcHist_)  {
		initOrMinApexDistHist_ = initedApexDist_;
		initOrMaxPcHist_ = pc;
	}

}


//! sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
void CornerApex::finitCornerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj, bool overwriteTrapping)  {
	if(!exists_) return;

	int nTraps=parentShape_->eleman()->isTrappedWat(filmBlob) + (parentShape_->eleman()->isTrappedOil() || outerLayerApex_->trappingCL().first>-1);
	if( inited_/*(overwriteTrapping && nTraps < 2)*/ || nTraps < 1)  {
		//double apxDistOld=initedApexDist_;
		double apxDist;
		double conAng(oilInj ? conAngRec : conAngAdv);
		getCApexDistConAng(apxDist, conAng, pc, halfAng, intfacTen, true); ///. warnign don't send initedApexDist_ directly
		if (apxDist <=0. || apxDist > parentShape_->RRR()/tan(halfAng)*2) 	{ cout<<" xkf ";	}

		initedApexDist_ = apxDist;

		advancingPc_ = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(initedApexDist_*sin(halfAng));
		receedingPc_ = intfacTen*cos( min(PI, conAngRec)+halfAng )/(initedApexDist_*sin(halfAng));

		if(pc > initOrMaxPcHist_)  {
			initOrMinApexDistHist_ = initedApexDist_;
			initOrMaxPcHist_ = pc;
		}
		inited_ = false;
	}



}

//! sets m _pinnedApexDist m _advancingPc m _receedingPc based on pc....
void CornerApex::initCornerApex(double pc, double conAngRec, double conAngAdv, double halfAng, double intfacTen, bool oilInj)  {
	if(!exists_ || trappedCL_.first > -1 || parentShape_->eleman()->isTrappedOil() || outerLayerApex_->trappingCL().first>-1) return;
	inited_ = true;


		double  recPc = intfacTen*cos( min(PI, conAngRec+halfAng) )/(initedApexDist_*sin(halfAng));
		if (receedingPc_ < recPc) receedingPc_ = recPc;  ///. to handle strange wettability change

		advancingPc_ = intfacTen*cos( min(PI, conAngAdv)+halfAng )/(initedApexDist_*sin(halfAng));


}





void CornerApex::getCApexDistConAng(double & apexDist, double & conAng, double pc, double halfAng, double intfacTen, bool overidetrapping, bool accurat, bool debug) const
{
	double delta = accurat ? 0.: SMALL_NUM;
		///bool Warning_ensure_deactivated;
		//ensure(inited_ || trappedCL_.first>-1 || parentShape_->eleman()->isTrappedOil() || outerLayerApex_->trappingCL().first>-1);
		//ensure(exists_);

	if(!exists_)	{ apexDist = MOLECULAR_LENGTH;  return; }

	if(!overidetrapping)  {

		apexDist=initedApexDist_;
		if( trappedCL_.first>-1)  {
			pc=trappedCL_.second;
			double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));
			if (debug) cout<< "d8sAA";
			return;
		}
		else if (parentShape_->eleman()->trappingOil().first>-1)  {
			pc=parentShape_->eleman()->trappingOil().second;
			double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));
			if (debug) cout<< "d8sBB";
			return;
		}
		else if (outerLayerApex_->trappingCL().first>-1)  {
			pc=outerLayerApex_->trappingCL().second;
			double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));
			if (debug) cout<< "d8sCC";

			return;
		}
	}


	if(advancingPc_-delta <= pc && pc <= receedingPc_+delta)  {
		double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
		part = std::min(part, 0.999999);
		part = std::max(part, -0.999999);
		double hingAng(acos(part)-halfAng);
		ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
		conAng = (std::min(std::max(hingAng, 0.), PI));
		apexDist=initedApexDist_;
	}
	else if(pc < advancingPc_)  {

		conAng = parentShape_->conAngleAdv();
			if (debug) 	cout<<"   efs1 ";

		apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
		if (apexDist < initedApexDist_)  {

			double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));
			apexDist = initedApexDist_;
			if (debug) 	cout<<"   efs2 ";

		}
		else if (apexDist<0) cout<<"ErrApL<0";

	}
	else if(pc > initOrMaxPcHist_)  {   ///. recover pinning receeding contact angle

		conAng=parentShape_->minInitRecCntAng();
		apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
			if (debug) 	cout<<"   efs3 ";

	}
	else  if(pc > receedingPc_)  {
		conAng = parentShape_->conAngleRec();

		apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
		if (apexDist > initedApexDist_)  {


			double part(pc*initedApexDist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng >= -SMALL_NUM && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));

			if (debug) 	cout<<"   pcrecPcusds "<<pc<<" "<<receedingPc_<<" "<<apexDist<<" <= "<<initedApexDist_<<" "<<conAng<<"         ";

			apexDist = initedApexDist_;
			if (debug) 	cout<< apexDist / parentShape_->RRR()*tan(halfAng) <<"    ";
			if (debug) 	cout<<"   efs4 ";

		}
		else if(apexDist<initOrMinApexDistHist_)  {
			double part(pc*initOrMinApexDistHist_*sin(halfAng)/intfacTen);
			part = std::min(part, 0.999999);
			part = std::max(part, -0.999999);
			double hingAng(acos(part)-halfAng);
			ensure(hingAng+SMALL_NUM >= parentShape_->minInitRecCntAng() && hingAng <=  PI+SMALL_NUM);
			conAng = (std::min(std::max(hingAng, 0.), PI));
			if (debug) cout<<"  pcrecPcuadsdsds "<<pc<<" "<<receedingPc_<<"   ";
			apexDist=initOrMinApexDistHist_;
			if (debug) 	cout<<"   efs5 ";

		}
	}
	else
	{
		conAng = conAng;
		apexDist = (intfacTen/pc)*cos(conAng+halfAng)/sin(halfAng);
			if (debug) 	cout<<"   efs6 ";

		ensure(conAng >= 0. && conAng <=  PI);
	}

	if (debug) if (apexDist > parentShape_->RRR()/tan(halfAng)*2)
	{
		cout<<" syp ";
	}
	 if (debug) cout<<" pc"<<pc<<",conAng"<<conAng<<" "<<endl;

}
