
#include "polygon.h"
#include "apex.h"
#include "cornerApex.h"
#include "layerApex.h"

#include "sortedEvents.h"
#include "compareFuncs.h"

/**
// Polygon Class Constructors
*/
VoidElem::VoidElem(Element& parent, const CommonData& common, 
	double radius, double shapeFactor, int connNum) 
	: ElemModel(parent, common, radius, connNum), m_shapeFactor(shapeFactor)
                 
{
    m_area = pow(m_R, 2.0) / (4.0 * m_shapeFactor);     //From Oren; This is ~correct for all shape
    m_SatWater = 1.0;
    m_minInitRecCntAng = 180.0;

    ensure(m_area > 0.0);

}


void Polygon::insertWatSnapEvent_IfSnapPcHgPc(SortedEvents< Apex*, PceImbCmp >& watEvents,  double globalPc)
{
    if(m_bulkFluid == &m_comn.oil())  // Water snap off when the bulk is filled by
    {                                                   // oil. Make sure there is no trapping in
		if (m_waterInCorner[0].cornerExists() && (m_waterInCorner[0].trappedCorner().first<0) && !(m_elem.isTrappedOil()))
		{
			if(m_elem.entryPc() > globalPc) 
			{
				watEvents.quickInsert(&m_elem);
			}
		}
    }
    else if( m_oilLayer[0].exists() && (m_oilLayer[0].trappedOLayer().first<0) &&
        (m_waterInCorner[0].trappedCorner().first<0 || m_elem.trappingWatBulk().first<0))   // Collapsing oil layers. Again make sure that
    {                                                                                           // trapping allows us to collapse it
        for(int i = 0; i < m_numCorners; ++i)
        {	
            if(m_oilLayer[i].LayerApex::exists(/*st ab le*/))
            {
                if(m_oilLayer[i].layerCollPc() > globalPc)  ///. Why such check
                {
                    pair<Polygon*, int> newEvent(this, i);//, collPc+m_elem.gravityCorrection()
                    watEvents.quickInsert(&m_oilLayer[i]);
                }
            }
        }
    }
}


void Polygon::insertOilSnapEvent_IfSnapPcLgPc(SortedEvents< Apex*, PceDrainCmp >& oilEvents, double globalPc)
{
    if(m_bulkFluid == &m_comn.water() && m_oilLayer[0].exists(/*st ab le*/) && m_oilLayer[0].trappedOLayer().first<0)
    {   
		ensure(m_numLayers==0);

        double snapOffPrs = calcSnapOffPressureDrain();
		
		 cout<<"  coel Pcs:"<< m_elem.entryPc() << "  " <<snapOffPrs<<endl;
        if(snapOffPrs < globalPc)
        {
            //pair<Polygon*, int> newEvent(this, -1);///, /*m_*/snapOffPrs+m_elem.gravityCorrection()
            oilEvents.quickInsert(&m_elem);
        }
    }
}




/**
 * set radius and update area and waterInCentre conductivity accordingly
 */
void VoidElem::setRadius(double newRad)
{
    m_R = newRad;
    m_area = pow(newRad, 2.0) / (4.0 * m_shapeFactor);
    ensure(m_area > 0.0);
}




///  Contact angles are assigned based on the principle of equilibrium contacr angles
/// DAng  is (modelTwo/5) Seperation angle,
void VoidElem::setContactAngle(double equilConAng, int wettClass, double DAng)
{
	

    m_conAngEquil = equilConAng;
	if(wettClass == 1)
	{
		m_cntAngRec = m_conAngEquil;
		m_cntAngAdv = m_conAngEquil;
	}
	else if(wettClass == 2)
	{
        double growthExp((PI+DAng)/PI);

        m_cntAngRec = max(0.0, growthExp*m_conAngEquil-DAng);
        m_cntAngAdv = min(PI, growthExp*m_conAngEquil);
	}
	else  if(wettClass == 3)
	{
        if(m_conAngEquil < 0.38349) m_cntAngRec = 0.0;
        else if(m_conAngEquil < 1.5289) m_cntAngRec = (0.5*exp(0.05*m_conAngEquil*180.0/PI)-1.5)*PI/180.0;
        else if(m_conAngEquil < 2.7646) m_cntAngRec = 2.0*(m_conAngEquil-1.19680);
        else m_cntAngRec = PI;

		if(m_conAngEquil < 0.38349) m_cntAngAdv = 0.0;
        else if(m_conAngEquil < 1.61268) m_cntAngAdv = 2.0*(m_conAngEquil-0.38349);
        else if(m_conAngEquil < 2.75805) m_cntAngAdv = (181.5 - 4051.0*exp(-0.05*m_conAngEquil*180.0/PI))*PI/180.0;
        else m_cntAngAdv = PI;
	}
	else  if(wettClass == 4)
	{
		m_cntAngAdv = m_conAngEquil;
		m_cntAngRec = pow(PI-1.3834263-pow(PI-m_cntAngAdv+0.004,0.45), 1.0/0.45)-0.004;
	}
	else
	{
		double PlusCoef = PI-(0.1171859 *DAng*DAng*DAng -0.6614868 *DAng*DAng + 1.632065 *DAng) ;
		double exponentCoef = 1.0-(0.01502745 *DAng*DAng*DAng -0.1015349 *DAng*DAng + 0.4734059 *DAng) ;

		m_cntAngAdv = m_conAngEquil;
		m_cntAngRec = pow(PlusCoef-pow(PI-m_cntAngAdv+0.004,exponentCoef), 1.0/exponentCoef)-0.004;
	}

    Polygon *polyShape = dynamic_cast< Polygon* >(this);
    if(polyShape)
    {
        for(int i = 0; i < polyShape->numCorners(); ++i)
        {
            polyShape->oilLayerCh()[i].advConAng(m_cntAngAdv);
        }
    }
    
    m_minInitRecCntAng=min(m_cntAngRec,m_minInitRecCntAng);

	//m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;     // TOBE initialised properly later
	//m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R;     // TOBE initialised properly later

}




/**
// Polygon Class Constructors
*/
Polygon::Polygon(Element& parent, const CommonData& common, double radius, 
                 double shapeFactor, int numCorners, int connNum)
        : VoidElem(parent, common, radius, shapeFactor, connNum),
    m_numCorners(numCorners),m_displacementType('X')
{
    m_numLayers = 0;
    m_maxConAngSpont = PI/2;

	///. Warning should be constructed as array to allow comparison 
	m_waterInCorner = new CornerApex[numCorners];
	m_oilLayer = new LayerApex[numCorners];
	//m_waterInCorner.reserve(m_numCorners);
	//m_oilLayer.reserve(m_numCorners);    
	for(int i = 0; i < numCorners; ++i)
    {
        m_waterInCorner[i].setCornerConnections(m_oilLayer+i,this,i);
        m_oilLayer[i].setLayerConnections(m_waterInCorner+i,this,i);
    }
}






/**
// Triangle Class Constructors
*/
Triangle::Triangle(Element& parent, const CommonData& common,  double radius,
                   double shapeFactor, int connNum, istringstream& inputData) 
      : Polygon(parent, common, radius, shapeFactor, 3, connNum)
{
    m_crnHafAngs.resize(m_numCorners);
    setHalfAngles();
}

void Triangle::setShapeFactor(double shapeFact)
{
    ensure(shapeFact <=  sqrt(3.0)/36.0);
    m_shapeFactor = max(shapeFact,1.0e-6);
    m_area = pow(m_R, 2.0) / (4.0 * m_shapeFactor);
    ensure(m_area > 0.0);
    //m_areaWater = m_area;
    setHalfAngles();
    //m_conductanceWater.second = SPConductance(m_areaWater, m_comn.water().viscosity());
}

/**
// This function evaluates the three half angles that make up the triangular pore. The routine
// follows the outline that was described by Patzek.
*/
void Triangle::setHalfAngles()
{
    double beta_2_min = atan((2.0/sqrt(3.0))*cos(acos(-12.0*sqrt(3.0)*m_shapeFactor)/3.0+4.0*PI/3.0));
    double beta_2_max = atan((2.0/sqrt(3.0))*cos(acos(-12.0*sqrt(3.0)*m_shapeFactor)/3.0));
    double randNum = double(rand()) / double(RAND_MAX);
	randNum = 0.5*(randNum+0.5);///.25-.75
	
    m_crnHafAngs[1] = beta_2_min + (beta_2_max - beta_2_min)*randNum;
    m_crnHafAngs[0] = 0.5*(asin((tan(m_crnHafAngs[1])+4.0*m_shapeFactor)
        * sin(m_crnHafAngs[1]) / (tan(m_crnHafAngs[1])-4.0*m_shapeFactor))-m_crnHafAngs[1]);
    m_crnHafAngs[2] = PI/2.0 - m_crnHafAngs[1] - m_crnHafAngs[0];
}


/**
// Square Class Constructors
*/
Square::Square(Element& parent, const CommonData& common, 
               double radius, int connNum, istringstream& inputData) 
   : Polygon(parent, common, radius, 0.0625, 4, connNum)
{
    m_crnHafAngs.resize(m_numCorners, PI/4.0);
}


/**
// Circle Class Constructors
*/
Circle::Circle(Element& parent, const CommonData& common, 
               double radius, int connNum, istringstream& inputData)
      : VoidElem(parent, common, radius, 1.0/(4.0*PI), connNum)
{
}




/**
// Polygon Class Destructor
*/
Polygon::~Polygon()
{
    //for(int i = 0; i < m_numCorners; ++i)
    //{
        delete[] m_waterInCorner;
        delete[] m_oilLayer;
    //}
}







///. WARNING: does not check from which throat water injection is taking place
void Polygon::fillCentreWithWaterCreateLayers(bool snapOffOverRide)
{
    m_bulkFluid = &m_comn.water(); ///K
    m_waterConnection = true;
    m_oilConnection = false;
    //ensure(m_waterInCorner[0].trappedCorner().first<0," flxpq ");
    //double snapOffPrs = ;
    double pc = m_comn.GuessCappPress()-m_elem.gravityCorrection();

    for(int i = 0; i < m_numCorners; ++i)
    {
		double snapOffPrs = m_waterInCorner[0].trappedCorner().first<0 ? calcSnapOffPressureImb() : -1.0e27;
        if( !snapOffOverRide && m_elem.entryPc() >  snapOffPrs)
        {	
            m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], m_comn.oil().interfacialTen(),false);
			//LayerApex* oillayeri = ;
            if ( m_oilLayer[i].createOLayer(pc, m_cntAngRec, m_cntAngAdv, m_maxConAngSpont, m_crnHafAngs[i],  m_comn.oil().interfacialTen(),false))
            {
				ensure(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
                ++m_numLayers;
                //this->addFillingEventToHistory(5)
            }
            else
            {
				m_waterInCorner[i].removeCorner();int Warning; ///. inconsistent film trapping and connectivity handling
			}
        }

    }

    if(m_numLayers)
    {
        m_oilConnection = true;
    }
    m_hasDisConectedCentreWCornerW = m_numLayers == m_numCorners;

    m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;

}


void Polygon::fillCentreWithOilRemoveLayers()
{
     m_oilConnection = true; m_bulkFluid = &m_comn.oil();
   

    for(int i = 0; i < m_numCorners; ++i)
    {
      
        m_oilLayer[i].LayerApex::removeLayer();//m_ exists
        

        if (!m_waterInCorner[i].cornerExists())
			m_waterInCorner[i].createFilm(m_elem.entryPc(), m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i], m_comn.oil().interfacialTen(), true);
    }
    m_numLayers = 0;
    m_hasDisConectedCentreWCornerW = false;

    m_waterConnection = m_waterInCorner[0].cornerExists();

    m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R;

}







void Circle::fillCentreWithWaterCreateLayers(bool snapOffOverRide)
{
    m_bulkFluid = &m_comn.water();///K
    m_oilConnection = false;
    m_waterConnection = true;
    m_Pc_pistonTypeAdv = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngAdv)) / m_R;

}


void Circle::fillCentreWithOilRemoveLayers()
{
    
    m_oilConnection = true; m_bulkFluid = &m_comn.oil();///K
    m_waterConnection = false;
    m_Pc_pistonTypeRec = m_comn.oil().interfacialTen()*(2.0*cos(m_cntAngRec)) / m_R;

}


void Polygon::calcOilLayerPc_syncTrappings(double pc)
{///. pc is used for error checking only

    const pair<int, double>& oilTrp = m_elem.trappingOil();
	const pair<int, double>& watTrpBulk = m_elem.trappingWatBulk();
	const pair<int, double>& watTrpFilm = m_elem.trappingWatFilm();

	//ensure (m_elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1)

	double ComWarningToFix = 1.0;

	bool oilInj =  m_comn.injectant() == &m_comn.oil();
	double tension =  m_comn.oil().interfacialTen();
	if(oilTrp.first > -1)
	{
        if (containCOil()&&watTrpBulk.first>0) cout<<"\n * Wrong trapping flags * \n";


			
		for(int i = 0; i < m_numCorners; ++i)
		{	            double oldLPc=m_oilLayer[i].layerCollPc();
		


			m_waterInCorner[i].finitCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj, true);
			if (!m_oilLayer[i].finitLayerApex(oilTrp.second, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj, true))
			{
				if (m_comn.debugMode>100) 
				{	vector<int> addCrns;
					cout<<endl<<" q3jd "; 
					cout
					<<m_elem.isInWatFloodVec()
					<<m_elem.canBeAddedToEventVec(&(m_comn.water()))
					<<m_elem.addToLayerVec(&m_comn.water(), this, addCrns)
					<<(m_elem.trappingWatFilm().first>-1)
					<<(m_elem.trappingWatBulk().first>-1)
					<<(m_oilLayer[i].trappedOLayer().first>-1) 
					<<(m_oilLayer[i].isInWatFloodVec())
					<<canNOTReconfigure(&(m_comn.water()))
					<<(m_elem.trappingOil().first>-1)

					<<"\n   "<<oilTrp.second<<" <?> "<<m_oilLayer[i].layerCollPc()<<" <?> "<<oldLPc<<" <?> "<<pc<< "    "
					<<"   "<<m_waterInCorner[i].advancingPc()<<" <?> "<<m_oilLayer[i].advancingPc()<<" <?> "<<0<<" <?> "<<pc<< "    "
					<<endl;
					//m_crnHafAngs[100000]=0;
					//exit(-1);
				}

				Pc_pin_disconnectOilLayer(i);

				
			}
			else  if (m_oilLayer[i].exists())
			{
				ensure(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
			}
			
			m_oilLayer[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
			m_comn.oil().interfacialTen(), oilInj);	
			
			m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj);
		}

        m_elem.virgin = false;
	}

	

	if(watTrpFilm.first > -1)
	{

		for(int i = 0; i < m_numCorners; ++i)
		{
			m_waterInCorner[i].finitCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj, true);
			if (!m_oilLayer[i].finitLayerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj,false))
			{ 
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" q4jd ";

			}
  			m_waterInCorner[i].CornerApex::markTrappingCorner(watTrpFilm, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
               tension, oilInj); ///. 

            m_oilLayer[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
               tension, oilInj); ///. 


			if (!m_oilLayer[i].initLayerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj))
					if (m_comn.debugMode>0) cout<<" q4jd2";


        }
	}
	
	if(watTrpBulk.first > -1)
	{

		for(int i = 0; i < m_numCorners; ++i)
		{
			m_waterInCorner[i].finitCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj,false);
			if (!m_oilLayer[i].finitLayerApex(watTrpBulk.second, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj, true))
			{
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" q8jd ";

			}
			else if (m_oilLayer[i].exists())
			{
				ensure(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
				if(!(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1])))))
				{
					//m_elem.unTrapWat(bulkBlob);
					cout<<endl
					<<m_elem.canBeAddedToEventVec(&(m_comn.oil()))
					<<m_elem.isInOilFloodVec()
					<<(m_elem.trappingWatFilm().first>-1)
					<<canNOTReconfigure(&(m_comn.oil()))
					<<(m_oilLayer[i].trappedOLayer().first>-1) 
					<< (m_elem.trappingOil().first>-1)
					<< (m_elem.trappingWatBulk().first>-1)<<endl;
					
					cout<<snapOffPrsDrainSpontNR()<<" <?> "
					<<"   "<<watTrpBulk.second<<" <?> "<<m_elem.entryPc()<< "    ";
					cout<<m_oilLayer[i].pinnedApexDist() <<"  <  "<< m_R<<" * (1.0/tan( "<<m_crnHafAngs[0]<<" ) + 1.0/tan( "<<m_crnHafAngs[1]<<endl;
					cout<<m_oilLayer[i].pinnedApexDist() <<"  <  "<< (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1])))<<endl;
					//((double*)&m_R)[100000]=0.0;

					//exit(-1);
				}
			}

            //m_oilLayer[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
               //tension, oilInj); ///. 
               
			m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj);
			if (!m_oilLayer[i].initLayerApex(watTrpBulk.second, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj))
					if (m_comn.debugMode>0) cout<<" q4jd3";

		}
        m_elem.virgin = false;
	}
	
	
	
}




void Polygon::calcOilLayerPc_markUntrappedFilms(double pc)
{///. pc is used for error checking only

    const pair<int, double>& oilTrp = m_elem.trappingOil();
	const pair<int, double>& watTrpBulk = m_elem.trappingWatBulk();
	const pair<int, double>& watTrpFilm = m_elem.trappingWatFilm();

	//ensure (m_elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1 || elem->isTrappedOil.first>-1)

	double ComWarningToFix = 1.0;

	bool oilInj =  m_comn.injectant() == &m_comn.oil();
	double tension =  m_comn.oil().interfacialTen();
	
	//if(oilTrp.first > -1)
	{

		for(int i = 0; i < m_numCorners; ++i)
		{	            double oldLPc=m_oilLayer[i].layerCollPc();
		
  			m_waterInCorner[i].CornerApex::markTrappingCorner(watTrpFilm, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
               tension, oilInj); ///. 
			m_oilLayer[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
			m_comn.oil().interfacialTen(), oilInj);	

			m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj);
			if (!m_oilLayer[i].initLayerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj, true))
			{
				if (m_comn.debugMode>100) 
				{	vector<int> addCrns;
					cout<<" q3jdu "; 
					cout
					<<m_elem.isInWatFloodVec()
					<<m_elem.canBeAddedToEventVec(&(m_comn.water()))
					<<m_elem.addToLayerVec(&m_comn.water(), this, addCrns)
					<<(m_elem.trappingWatFilm().first>-1)
					<<(m_elem.trappingWatBulk().first>-1)
					<<(m_oilLayer[i].trappedOLayer().first>-1) 
					<<(m_oilLayer[i].isInWatFloodVec())
					<<canNOTReconfigure(&(m_comn.water()))
					<<(m_elem.trappingOil().first>-1)

					<<"\n   "<<oilTrp.second<<" <?> "<<m_oilLayer[i].layerCollPc()<<" <?> "<<oldLPc<<" <?> "<<pc<< "    "
					<<"   "<<m_waterInCorner[i].advancingPc()<<" <?> "<<m_oilLayer[i].advancingPc()<<" <?> "<<0<<" <?> "<<pc<< "    "
					<<endl;
					//m_crnHafAngs[100000]=0;
					//exit(-1);
				}

				Pc_pin_disconnectOilLayer(i);
			}
			else  if (m_oilLayer[i].exists())
			{
				ensure(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
			}
			
			//m_waterInCorner[i].finitCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj,false);
			//if (!m_oilLayer[i].finitLayerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj,false))
			//{
					//if (m_comn.debugMode>0) cout<<" q4jd3";
				//Pc_pin_disconnectOilLayer(i);
//
			//}
		}
	}

	return;

	//if(watTrpFilm.first > -1)
	{

		for(int i = 0; i < m_numCorners; ++i)
		{
 			m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj);
			if (!m_oilLayer[i].initLayerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj))
			{ 
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" q4jdu ";

			}

        }

	}
	
	//if(watTrpBulk.first > -1)
	{

		for(int i = 0; i < m_numCorners; ++i)
		{
            //m_oilLayer[i].LayerApex::markCentreLayerTrappings(oilTrp, pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],
               //tension, oilInj); ///. 			

            m_waterInCorner[i].initCornerApex(pc, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj);
			if (!m_oilLayer[i].initLayerApex(watTrpBulk.second, m_cntAngRec, m_cntAngAdv, m_crnHafAngs[i],tension, oilInj))
			{
				Pc_pin_disconnectOilLayer(i);
				if (m_comn.debugMode>0) cout<<" q8jdu ";
			}
			else if (m_oilLayer[i].exists())
			{
				ensure(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1]))));
				if(!(m_oilLayer[i].pinnedApexDist() < (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1])))))
				{
					//m_elem.unTrapWat(bulkBlob);
					cout<<endl
					<<m_elem.canBeAddedToEventVec(&(m_comn.oil()))
					<<m_elem.isInOilFloodVec()
					<<(m_elem.trappingWatFilm().first>-1)
					<<canNOTReconfigure(&(m_comn.oil()))
					<<(m_oilLayer[i].trappedOLayer().first>-1) 
					<< (m_elem.trappingOil().first>-1)
					<< (m_elem.trappingWatBulk().first>-1)<<endl;
					
					cout<<snapOffPrsDrainSpontNR()<<" <?> "
					<<"   "<<watTrpBulk.second<<" <?> "<<m_elem.entryPc()<< "    ";
					cout<<m_oilLayer[i].pinnedApexDist() <<"  <  "<< m_R<<" * (1.0/tan( "<<m_crnHafAngs[0]<<" ) + 1.0/tan( "<<m_crnHafAngs[1]<<endl;
					cout<<m_oilLayer[i].pinnedApexDist() <<"  <  "<< (m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[1])))<<endl;
					//((double*)&m_R)[100000]=0.0;

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
	int * sasasas;
	if (m_waterInCorner[0].trappedCorner().first>-1 || eleman()->isTrappedOil() || m_oilLayer->trappingCL().first>-1) 
	{
		if (m_comn.debugMode>100) cout<<" nkdps ";
		return -1.0e28;
	} 

    if(m_cntAngAdv < PI/2.0 - m_crnHafAngs[0])              // Spontaneous
        snapOffPrs = snapOffPrsImbNR();
    else                                                   // Forced
		snapOffPrs = m_waterInCorner[0].cornerExists() ? m_waterInCorner[0].advancingPc() : -1.0e64;
   
    return snapOffPrs-m_elem.snapOfLongitCurvature()*m_comn.oil().interfacialTen();

}


double Polygon::calcSnapOffPressureDrain() const
{
    double snapOffPrs;
    if(m_cntAngRec > PI/2.0 + m_crnHafAngs[0])      // Spontaneous
    {
        snapOffPrs = snapOffPrsDrainSpontNR();
		ensure(snapOffPrs!=0.0);
    }
    else
    {                                                       // Forced
        snapOffPrs = m_oilLayer[0].receedingPc();
        ensure(snapOffPrs>0.0);
        if (snapOffPrs <=  0.0)
        {
			cout<<m_oilLayer[0].receedingPc()<<endl;
		}
    }


    return snapOffPrs-m_elem.snapOfLongitCurvature()*m_comn.oil().interfacialTen();
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

    double tension(m_comn.oil().interfacialTen());
	//- Tom Bultreys bug-report & fix R-> 1/R
    double snapOffPrsOne = tension * cos(m_cntAngAdv+m_crnHafAngs[0] )/((m_R/tan(m_crnHafAngs[0]) + m_R/tan(m_crnHafAngs[2]) - m_waterInCorner[2].pinnedApexDist())*sin(m_crnHafAngs[0]));

    double oldPc, errorVal(1.0);
    if(m_waterInCorner[0].trappedCorner().first>-1)
        oldPc = m_waterInCorner[0].trappedCorner().second;
    else if(m_comn.isDrainageCycle())
        oldPc = m_comn.GuessCappPress()-m_elem.gravityCorrection();
    else
        oldPc = m_comn.maxPcLastDrainCycle()-m_elem.gravityCorrection();

		const double L0pL2 = (m_R/tan(m_crnHafAngs[0]) + m_R/tan(m_crnHafAngs[1]))   ;
    for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
    {

		double apexDist2, teta2(m_cntAngAdv); 
		m_waterInCorner[1].getCApexDistConAng(apexDist2, teta2, oldPc, m_crnHafAngs[1], tension, true, true);

        double rL2 = -apexDist2*sin(m_crnHafAngs[1])/(tension*sin(teta2+m_crnHafAngs[1])); 

        double func = oldPc - tension*(cos(m_cntAngAdv)/tan(m_crnHafAngs[0]) - sin(m_cntAngAdv) + cos(teta2)/tan(m_crnHafAngs[1]) - sin(teta2)) / L0pL2;

        double funcDpc = 1 + tension*(rL2*sin(teta2)/tan(m_crnHafAngs[1]) +  rL2*cos(teta2)) /L0pL2;

        double newPc = (abs(funcDpc)>1.0e-32) ? oldPc - func / funcDpc  : oldPc - func;
        if (i>MAX_NEWT_ITR/2 ) newPc = 0.5*(newPc+oldPc);

        errorVal = fabs(newPc-oldPc)/(abs(oldPc)+1.0e-8);
        if(errorVal < EPSILON)
        {
            if(teta2 <= m_cntAngAdv+0.000001)
                return max(newPc, snapOffPrsOne)+0.0001;     // Return the pressure occuring at the highest Pc
            else
                return snapOffPrsOne+0.0001;
        }

        oldPc = newPc;
    }

/*    double tension(m_comn.oil().interfacialTen());
    double snapOffPrsOne = (tension/m_R) * (cos(m_cntAngAdv) -
            (2.0*sin(m_cntAngAdv))/(1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[1])));    // Small Contact Ang

    double oldPc, errorVal(1.0);
    if(m_waterInCorner[0].trappedCorner().first>-1)
        oldPc = m_waterInCorner[0].trappedCorner().second;
    else if(m_comn.isDrainageCycle())
        oldPc = m_comn.GuessCappPress()-m_elem.gravityCorrection();
    else
        oldPc = m_comn.maxPcLastDrainCycle()-m_elem.gravityCorrection();

    for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
    {

		double pinnedApexDistance, hingAng(m_cntAngAdv); 
		m_waterInCorner[2].getCApexDistConAng(pinnedApexDistance, hingAng, oldPc, m_crnHafAngs[2], m_comn.water().interfacialTen(), true, true);


        double hingAngDpc = -pinnedApexDistance*sin(m_crnHafAngs[2])/(tension*sin(hingAng+m_crnHafAngs[2]));

        double func = oldPc - (m_comn.oil().interfacialTen()/m_R)*((cos(m_cntAngAdv)/tan(m_crnHafAngs[0]) -
            sin(m_cntAngAdv) + cos(hingAng)/tan(m_crnHafAngs[2]) - sin(hingAng)) /
            (1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[2])));

        double funcDpc = 1 + (m_comn.oil().interfacialTen()/m_R)*((hingAngDpc*sin(hingAng)/tan(m_crnHafAngs[2]) +
            hingAngDpc*cos(hingAng)) / (1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[2])));

        double newPc = oldPc - func / funcDpc;

        errorVal = fabs((newPc-oldPc)/oldPc);
        if(errorVal < EPSILON)
        {
            if(hingAng < m_cntAngAdv)
                return max(newPc, snapOffPrsOne);     // Return the pressure occuring at the highest Pc
            else
                return snapOffPrsOne;
        }

        oldPc = newPc;
    }*/


    cerr <<"\n============  snap Off Prs Imb NR ================ " << endl
        << "Error: Failed to obtain valid value for threshold" << endl
        << "capillary pressure during snap off.              " << endl
        << "m_cntAngAdv: " << m_cntAngAdv                               << endl
        << "trapping ind" << m_elem.trappingOil().first << endl
        << "error " << errorVal<<oldPc << endl
        << "=================================================" << endl;   
        
         //exit(-1);

    return oldPc;
}

double Square::snapOffPrsImbNR() const///. to check
{
        return (m_comn.oil().interfacialTen()/m_R) * (cos(m_cntAngAdv) - sin(m_cntAngAdv));
}

/**
// The problem with these NR approaches are that they're a bit sensitive
// to the first guess.
*/
double Triangle::snapOffPrsDrainFromCorner(bool& casePossible, int cor) const
{
    if(!m_oilLayer[cor].exists()) return 0.0;

	if (m_cntAngRec <= PI/2.0 + m_crnHafAngs[cor])
	{ 
		casePossible=true;
		return m_oilLayer[cor].receedingPc();
	}
	
    double minLocalPcLastCycle(m_comn.minPcLastImbCycle()-m_elem.gravityCorrection());
    double minLocalPc(m_comn.minEverCappPress()-m_elem.gravityCorrection());

    double min_Pc = m_oilLayer[0].stablePinnedInLastCycle(minLocalPcLastCycle) ? minLocalPcLastCycle: minLocalPc;
    double tension(m_comn.oil().interfacialTen()), errorVal(1.0);
    vector< double > initGuesses(3);

    initGuesses[2] = min_Pc ;                                                 // First attempt
    initGuesses[1] = m_comn.GuessCappPress()-m_elem.gravityCorrection();    // Second attempt
    if(m_oilLayer[cor].trappedOLayer().first>-1)                          // Getting desparete now
        initGuesses[0] = m_oilLayer[cor].trappedOLayer().second;
    else
        initGuesses[0] = m_elem.trappingOil().second;

    while(!initGuesses.empty())
    {
        double oldPc = initGuesses.back();
        initGuesses.pop_back();

        for(int i = 0; i < MAX_NEWT_ITR; ++i)               // Larger contact angles  => Need to use a NR technique
        {
            double tetaCor(m_cntAngRec);// = m_oilLayer[cor].hingingConAng(oldPc,m_cntAngRec,m_crnHafAngs[cor],tension,true);// acos(part)+m_crnHafAngs[cor];
            double pinnApexThree,halfAng(m_crnHafAngs[cor]);// = m_oilLayer[cor].getApexDistance(oldPc,tetaCor,m_crnHafAngs[cor],tension);// acos(part)+m_crnHafAngs[cor];
            m_oilLayer[cor].getCAApexDist(pinnApexThree,tetaCor,halfAng,oldPc,tension);// acos(part)+m_crnHafAngs[cor];

			ensure(pinnApexThree>0);

            double tetaCorDpc = -pinnApexThree*sin(m_crnHafAngs[cor])/(tension*sin(tetaCor-m_crnHafAngs[cor]));

            double func = oldPc - (tension/m_R)*(
				(cos(m_cntAngRec)/tan(m_crnHafAngs[0]) + sin(m_cntAngRec) +  cos(tetaCor)/tan(m_crnHafAngs[cor]) + sin(tetaCor)) /
                (1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[cor]))
                );

            double funcDpc = 1 - (tension/m_R)*((tetaCorDpc*sin(tetaCor)/tan(m_crnHafAngs[cor]) -
                tetaCorDpc*cos(tetaCor)) / (1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[cor])));

            double newPc = oldPc - func / funcDpc;
            errorVal = fabs((newPc - oldPc)/oldPc);

            if(errorVal < EPSILON)
            {
                ensure(tetaCor >= m_cntAngRec || m_comn.debugMode < 1000);
                casePossible = tetaCor >= m_cntAngRec;//  && m_oilLayer[cor].freeAtPrs(newPc);
                
                return newPc;
            }

            oldPc = newPc;
        }
    }

	if ( m_comn.debugMode > 0 )
    cerr <<"==============oil Apex Meet Corner =============" << endl
        << "Error: failed to obtain valid value for threshold" << endl
        << "       capillary pressure during snap off.       " << endl
        << "Error: " << errorVal                               << endl
        << "m_cntAngRec: " << m_cntAngRec                               << endl
        << "corner: " << cor                               << endl
        << "minLocalPcLastCycle: " << minLocalPcLastCycle                               << endl
        << "minLocalPc: " << minLocalPc                               << endl
        << "min_Pc: " << min_Pc                               << endl
        << "halfAng: " << m_crnHafAngs[cor]                               << endl
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
    double tension(m_comn.oil().interfacialTen());
    double caseOneSnapOff = tension*cos(m_cntAngRec-m_crnHafAngs[0])
    /(sin(m_crnHafAngs[0])* m_R*(1.0/tan(m_crnHafAngs[0])+1.0/tan(m_crnHafAngs[2])));

    double caseTwoSnapOff = (tension/m_R) * (cos(m_cntAngRec) +
            (2.0*sin(m_cntAngRec))/(1.0/tan(m_crnHafAngs[0]) + 1.0/tan(m_crnHafAngs[1])));

    bool case3Possible(false), case4Possible(false);

    double caseThreeSnapOff = snapOffPrsDrainFromCorner(case3Possible, 2);
    double caseFourSnapOff = snapOffPrsDrainFromCorner(case4Possible, 1);

    //if(case3Possible)
        //return min(caseThreeSnapOff, caseTwoSnapOff);
    //else if(case4Possible)
        //return min(caseFourSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
    //else if(m_oilLayer[1].freeAtPrs(caseTwoSnapOff))
        //return min(caseTwoSnapOff, caseOneSnapOff);
    //else
        //return caseOneSnapOff;
    if(case3Possible && case4Possible)
        return min(caseFourSnapOff, min(caseThreeSnapOff, min(caseTwoSnapOff, caseOneSnapOff)));
    else if(case3Possible)
        return min(caseThreeSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
    else if(case4Possible)
        return min(caseFourSnapOff, min(caseTwoSnapOff, caseOneSnapOff));
    else if(m_oilLayer[1].exists() && m_cntAngRec > PI/2.0 + m_crnHafAngs[0] )  //|| m_trappedCL.first < 0;
        return min(caseTwoSnapOff, caseOneSnapOff);
    else
        return caseOneSnapOff;
}

double Square::snapOffPrsDrainSpontNR() const
{
        return (m_comn.oil().interfacialTen()/m_R) * (cos(m_cntAngRec) + sin(m_cntAngRec));
}



