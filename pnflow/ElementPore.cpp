#ifdef WIN32
#pragma warning(disable:4786)
#endif

//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include <cmath>
//#include <cstdlib>
//#include <algorithm>
//#include <functional>
//#include <typeinfo>
#include <vector>
//#include <ctime>
//#include <stack>
//#include <set>
#include <cassert>
//#include <utility>
//#include <map>
using namespace std;

//#include "f2c.h"
//#include "sortedEvents.h"
//#include "threeSome.h"
#include "fluid.h"
//#include "elem_Model.h"
#include "polygon.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "compareFuncs.h"
#include "inputData.h"

#include "Element.h"


Pore::Pore(const CommonData& common, int indx, dbl3 nod, double radius, double volume,
		   double volumeClay, double shapeFactor, bool insideSlvBox, bool insideSatBox,
		   double initSolverPrs, vector<Elem*>& connThroats, int typ  ) 
   : Elem(common, indx, nod, radius, volume, volumeClay, shapeFactor, connThroats.size(), true,typ)  {
	isInCalcBox_ = insideSatBox;
	isInsideSolverBox_ = insideSlvBox;
	cnctions_ = connThroats;

	/*model_->*/setGravityCorrection(node_);

	//double minRad(1e+30), maxRad(1e-30), radSum(0.);
	//for(size_t i = 0; i < cnctions_.size(); ++i)
	//{
		//double rad = cnctions_[i]->model()->RRR();
		//minRad = min(minRad, rad);
		//maxRad = max(minRad, rad);
		//radSum += rad;
	//}

	//averageAspectRatio_ = model_->RRR()*cnctions_.size()/radSum;
	//maxAspectRatio_ = model_->RRR()/maxRad;
	//minAspectRatio_ = model_->RRR()/minRad;
}


double Pore::snapOfLongitCurvature() const
{ ///. to be implemented later
	//double length = 0.;
	//double length = 0.;
	//for(i = 0; i < cnctions_.size(); ++i)
	//{
		//const Throat* throat = dynamic_cast< Throat* >(cnctions_[i]);
		//assert(throat);
		//const Elem* prj = throat->neighbouringPore(this);
		//if(prj->isEntryRes()) connToIn = true;
		//if(prj->isExitRes()) connToOut = true;
		//outOne << setw(7) << prj->indexOren();                                 // Connecting nodes
	//}

		return 0.;
}


/// In order to reach the outlet quickly during the trapping routine, we sort the
/// connecting elements according to the distance to outlet
void Pore::sortConnectingElems_DistToExit()  {
	sort(cnctions_.begin(), cnctions_.end(), DistToExitCompareThroats());
}

/// These functions checks the integrity on the network, ie that pores and throats
/// agree that they point to each other and that there are no NULL pointers lurking
/// about. The total volume and number of connections contained within the network
/// is also counted here.
void Pore::prepare2()  {
	if(index()<10)  {
		 //outD<<"p:"<<(neib(0))->index()<<"-"<<neib(1)->index()<<":"<<index()<<"\n";
		outD<<"R_p"<<index()<<":"<<RRR()<<" "<<"  X_p"<<index()<<":"<<node_.x<<" "<<node_.y<<" "<<node_.z<<" \n";
	 }
	checkConnections();
}




/// Endpores are assumed to be triangular and have zero volume and contact angles.
/// Radius and interfacial tension is set to some arbitrary value to prevent any
/// divison by zero that may occur. They are further assumed to be triangular so that
/// they wont suddenly not contain water.
InOutBoundary::InOutBoundary(const CommonData& common, int indx, dbl3 node, vector<Elem*>& connThroats)
	 : Pore(common, indx, node, 1e-12, 0., 0., sqrt(3.)/36.+(indx>1), false, false, 0., connThroats,0)///.  warning hardcode pore type, last entry
{  	 // Warning radius here has a big impact if calc_box include boundary throats, 
	         // R<1e-13 leads to Error message
	Polygon* shyp = dynamic_cast<Polygon*>(model_);
	if(shyp)  {
	  for(int i = 0; i < shyp->numCorners(); ++i)  {
		shyp->oilLayerCh()[i].setInWatFloodVec(true);
		shyp->oilLayerCh()[i].setInOilFloodVec(true);
	  }
		//polySh ape->calc R(0.);
	}

	waterSaturation_ = 0.5;

	isExitRes_ = (indx==1); //TODO source of Error
	isEntryRes_ = (indx==0);

	isConnectedToExit_ = isExitRes_;
	isConnectedToEntry_ = isEntryRes_;
	connectedToNetwork_ = false;// will be reset

	model_->setContactAngle(0.,2,180.);

}


void InOutBoundary::fillElemCentreWithOilRemoveLayersIO(double pc)  {
	if ( model()->bulkFluid() != &(comn_.oil()) )
		fillElemCentreWithOilRemoveLayers();
	model_->SetInOutletPc_pistonTypeRec(pc);
	model_->SetInOutletPc_pistonTypeAdv(pc);
	connectedToNetwork_ = false;
}

void InOutBoundary::fillElemCentreWithWaterCreateLayersIO(double pc)  {
	if ( model()->bulkFluid() != &(comn_.water()) )
		fillElemCentreWithWaterCreateLayers();
	model_->SetInOutletPc_pistonTypeRec(pc);
	model_->SetInOutletPc_pistonTypeAdv(pc);
	connectedToNetwork_ = false;

}




inline dbl3   projectAx(dbl3 nod, int iAx, const dbl3& nod2) { nod[iAx>>1]=nod2[iAx>>1]; return nod; }



void InOutBoundary::prepare2()  {
	checkConnections();
	miroredPores_.resize(nCncts());
	//flowVolumePpT_+=vol;

	for(int i = 0; i < nCncts(); ++i)  {	 Throat* tr = neiT(i);
		if(tr->neiP(0)->index()==index()) // make sure index_ is compatible in mirrored pores
		{	vector<Elem*> neis(1,tr); // dumb!!
			Pore* mirr = new Pore(comn_, index(), projectAx(tr->neiP(1)->node(),index(),node()), tr->RRR()*1.1,
				0., 0., tr->model()->shapeFactor(), false, false, 0., neis,0);
			miroredPores_[i]=mirr; 
			 tr->neiSet(0, mirr);//if(index()>1)
		}
		else if(tr->neiP(1)->index()==index())  {	vector<Elem*> neis(1,tr); // dumb!!
			Pore* mirr = new Pore(comn_, index(), projectAx(tr->neiP(0)->node(),index(),node()), tr->RRR()*1.1,
				0., 0., tr->model()->shapeFactor(), false, false, 0., neis,0);
			miroredPores_[i]=mirr; 
			tr->neiSet(1, mirr);//if(index()>1) 
		}
	  else cout<<" Errorisdj ";
	  if(index()>1) 
	  { //ChModel()->bulkFluid(&comn_.noFluid());
		  	ChModel()->waterConnection_ = false;  ChModel()->oilConnection_ = false;  ChModel()->bulkFluid_ = &comn_.noFluid();
	  }
	}
}
InOutBoundary::~InOutBoundary()  {
	for(auto& mirr:miroredPores_)
		if(mirr) { delete mirr; mirr=nullptr; }
}

