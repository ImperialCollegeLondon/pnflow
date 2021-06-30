
/*---------------------------------------------------------------------------*\
2015:  Developed by Ali Q Raeini  email: a.q.raeini@imperial.ac.uk
\*---------------------------------------------------------------------------*/
 

#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <cassert>
#include <array>

#include "typses.h"
#include "results3D.h"
#include "FlowDomain.h"



#define _nl_  ((i&31)==31 ? '\n' : ' ')
#define pbak push_back









using namespace std;


results3D::results3D(const InputFile& input, const GNMData * comn, string outputFolder, const std::vector<Elem const *>*  elems_, size_t nBP2, size_t n6pPors)
  : comn_(comn), elems_(*elems_), nBSs_(nBP2), nBpPors_(n6pPors), format_(XMF_FORMAT), rkw_(1), vLinAC21I_(0), v3D_AC21I_(0),
	iWrite_(0), prefix_(outputFolder), rScaleFactor_(1.)
{

	inform=input.informative;

	istringstream ss;                  //!### Keyword **visualize** writes 3D visualization files
	if(input.giv("visualize", ss)) //!##### Keyword **visualize:**  	rScaleFactor 	thetaResulution	WriteInit	WriteOInj	WriteWInj	WriteAllSteps	WriteCorners;
	{ 
		if(inform) cout<<"Initialising full-3D vtu visualization";

		ss >> rScaleFactor_ ;


		ss>>nTheta_;
		ensure(nTheta_>2 , "not sufficient angular resolution:"+_s(nTheta_));
		ensure(nTheta_<19, "too high angular resolutione, "+_s(nTheta_)+", visualization files will be too larg");
		for (int i = 0;i<5;i++)
		{
			string tmp("F");  ss>>tmp;  char xr=tmp[0];
			if (xr=='T' || xr=='t')		v3D_AC21I_ |= (3<<(i*4)); //_visualise[i] = true;
			else if (isdigit(xr) && stoul(tmp)<15) v3D_AC21I_ |= (stoul(tmp)<<(i*4u));
			else ensure((xr=='F' || xr=='f'), "wrong data ("+tmp+").  Use visualize: rScale nTheta T(rue)/F(alse)/0-3...(5-times) ");
		}

		if(inform) cout<<" visualize: "<<(v3D_AC21I_);

		if ((rScaleFactor_<0.99 || rScaleFactor_>1.01) && v3D_AC21I_)
			cout<< "\n\n  Info: pore and throat radii will be multiplied by "<<rScaleFactor_<<" in visualization files \n"<<endl;
	}

	if(input.giv("visuaLight", ss)) //!##### Keyword **visuaLight** is outdated
	{
		if(inform) cout<<"visuaLight: initialising vtu & xmdf quasi-3D visualization";
		format_ = BOTH_FORMATS;
		for (int i = 0;i<3;++i)
		{
			string tmp("F"); ss>>tmp;  char xr=tmp[0];
			if (tmp[0]=='T' || tmp[0]=='t')		  vLinAC21I_ |=  (3<<(i*4));  //_visualiseL[i] = true;
			else if (tmp[0]=='F' || tmp[0]=='f')  vLinAC21I_ &= ~(15u<<(i*4));
			else if (isdigit(xr) && stoul(tmp)<15) { vLinAC21I_ &= ~(15u<<(i*4u)); vLinAC21I_ |= (stoul(tmp)<<(i*4u)); }
			else ensure(0, "wrong data ("+tmp+").  Use visuaLight: T(rue)/F(alse)/0-3...(3-times) ");
		}

		ss>>rkw_;
		if(inform) cout<<" vLinAC21I: "<<std::hex<<vLinAC21I_<<" rkw: "<<std::dec<<rkw_;
	}


	if(inform) cout<<endl;

}





//! write driver controlling when to write and what format to write
void results3D::write3D(double pc, double tension, bool endCycle)
{
	ensure(nBSs_>0,"internal data not initialized",2);
	int icycle=comn_->dispCycle();
	string suffix( (endCycle) ? "End" : _s(100*icycle+iWrite_));
	string prefix(prefix_+_s(icycle));
	if	  (!icycle)		            prefix += "_Init";
	else if (comn_->isDrainage() )	prefix += "_OInj";
	else                            prefix += "_WInj";


	///. ***** quasi-3D visualization  Vtu format *****
	if (format_ & VTU_FORMAT)
	if ( _1At(vLinAC21I_,icycle*4) || (endCycle  && _1At(vLinAC21I_,icycle*4+1)) )
	{ /// quasi-3D visualization
		if(inform) (cout<<" vis").flush();
			writeThroatLines(prefix+suffix+".vtu", pc, tension,icycle,iWrite_, endCycle);
		if(_1At(vLinAC21I_,icycle*4+2))
		{
			//vtuWriteThroatCornLines(prefix+"Layers"+suffix+".vtu", pc, tension,icycle,iWrite_, endCycle);
			if(inform) (cout<<"light ").flush();
		}
	}


	///. *****. Full 3D visualization  Vtu format *****
	if ( _1At(v3D_AC21I_,icycle*4) || (endCycle  && _1At(v3D_AC21I_,icycle*4+1)) )
	{
		(cout<<" visua").flush();
		vtuWritePores( prefix+"Pore"+suffix+".vtu", pc, tension);
		(cout<<"liza").flush();
		vtuWriteThroats(prefix+"Throat"+suffix+".vtu", pc, tension);
		(cout<<"tion  ").flush();
	}

	if (endCycle) iWrite_=1; //Warning double writing when all steps not requested //sync_pointxsdsdsdherei
	else 	++iWrite_;


}







float ffloatof(const Elem* slm){  float ff=slm->ffaz();
	if(ff>3.) {ff=(4.*float(int(ff)%4)+1.5)/5; }
	return max(0.6f,min(ff,2.4f)); }






int findOrInsertPoint(vector<dbl3>&  points, dbl3& point)
{
	for (vector<dbl3>::const_reverse_iterator  rip = points.rbegin(); rip != points.rbegin()+200 && rip != points.rend(); ++rip)
		if (*rip==point) return int(points.rend()-rip)-1;

	points.pbak(point);
	return points.size()-1;
}

dbl3 rotate( dbl3 n, dbl3 x, dbl3 y, double gamma)
{///. rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	return dbl3
	  (	( x[0]*(n[1]*n[1]+n[2]*n[2]) - n[0]*( x[1]*n[1]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (-x[2]*n[1]+x[1]*n[2]-n[2]*y[1]+n[1]*y[2] )*sin(gamma),
		( x[1]*(n[0]*n[0]+n[2]*n[2]) - n[1]*( x[0]*n[0]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + ( x[2]*n[0]-x[0]*n[2]+n[2]*y[0]-n[0]*y[2] )*sin(gamma),
		( x[2]*(n[0]*n[0]+n[1]*n[1]) - n[2]*( x[0]*n[0]+x[1]*n[1]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (-x[1]*n[0]+x[0]*n[1]-n[1]*y[0]+n[0]*y[1] )*sin(gamma)
	  );
}
dbl3 rotateAroundVec( dbl3 y, double gamma, dbl3 n)
{///. rotate y around n (line passing through centre, in the direction of n) http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	return dbl3
	  (	(  - n[0]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (n[1]*y[2]-n[2]*y[1])*sin(gamma),
		(  - n[1]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + (n[2]*y[0]-n[0]*y[2])*sin(gamma),
		(  - n[2]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (n[0]*y[1]-n[1]*y[0])*sin(gamma)
	  );
}


void insertHalfCorneroints(vector<dbl3>&  points, vector<int>& cellPoints, dbl3 c1, dbl3 c2, dbl3 nE1, dbl3 nE2, 
                          double lp_inrR, double lp_outR, ///<- layer location of corner premiter
							double teta1, double tetaW, double rOuter, double hafAng, double CARelax = 1.)
{
	//const double convertToRad = PI/180.;
	dbl3 dd = c2-c1;
	dbl3 normal = dd/(mag(dd)+1e-32);

	vector<dbl3> hcPs(8); // corner points

	///. radii of curvature
	double gama = abs(hafAng)*CARelax;
	double lc_iInnerR = lp_inrR*(cos(gama)+(sin(gama)/cos(gama+teta1))*(sin(gama+teta1)-1));
	double lc_iOuterR = lp_outR*(cos(gama)+(sin(gama)/cos(gama+tetaW))*(sin(gama+tetaW)-1));

	{
	 dbl3 e1 = c1+nE1*rOuter;///. edgePoint
	 dbl3 nE11 = rotateAroundVec(nE1,hafAng,normal);

	 if(hafAng>0)
	 {
		hcPs[0] = e1-lc_iInnerR*nE1;
		hcPs[1] = e1-lc_iOuterR*nE1;
		hcPs[2] = e1-lp_outR*nE11;
		hcPs[3] = e1-lp_inrR*nE11;
	 }
	 else
	 {
		hcPs[3] = e1-lc_iInnerR*nE1;
		hcPs[2] = e1-lc_iOuterR*nE1;
		hcPs[1] = e1-lp_outR*nE11;
		hcPs[0] = e1-lp_inrR*nE11;
	 }
	}


	{
	 dbl3 e1 = c2+nE2*rOuter;///. edgePoint
	 dbl3 nE11 = rotateAroundVec(nE2,hafAng,normal);
	 if(hafAng>0)
	 {
		hcPs[4] = e1-lc_iInnerR*nE2;
		hcPs[5] = e1-lc_iOuterR*nE2;
		hcPs[6] = e1-lp_outR*nE11;
		hcPs[7] = e1-lp_inrR*nE11;
	 }
	 else
	 {
		hcPs[7] = e1-lc_iInnerR*nE2;
		hcPs[6] = e1-lc_iOuterR*nE2;
		hcPs[5] = e1-lp_outR*nE11;
		hcPs[4] = e1-lp_inrR*nE11;
	 }
	}

	///. 8 points each elem 
	for (int i=0;i<8;++i)	cellPoints.pbak(findOrInsertPoint(points, hcPs[i]));

}



 

void getSolverPoreResults
(
	const GNMData *  netsim,
	const vector<Elem const *>&  elems_,
	const vector<int> & elmInds,
	vector<float> & radius,
	vector<float> & Sw,
	vector<float> & pcPiston,
	vector<float> & p,
	vector<float> & p_o,
	vector<float> & p_w,
	vector<float> & volt,
	vector<float> & elem_type,
	size_t nBpPors_
)
{

	const Fluid& c__water = netsim->water();
	const Fluid& c__oil = netsim->oil();
	const Fluid& c__elec = netsim->elec();
	vector<float> radiusTmp(elems_.size(),0.);
	vector<float> SwTmp(elems_.size(),-1.);
	vector<float> pcPistonTmp(elems_.size(),-1.);
	vector<float> pTmp(elems_.size(),-1.);
	vector<float> p_oTmp(elems_.size(),-1.);
	vector<float> p_wTmp(elems_.size(),-1.);
	vector<float> voltTmp(elems_.size(),-1.);
	vector<float> elem_typeTmp(elems_.size(),-1.);
	for(size_t i = 2; i<nBpPors_; ++i)//i < elems_.size()
	{
		const Elem* pore = elems_[i];


		double pressure(-10.);
		double volttt(-10.);
		double flowRate(-10.);
		// if(   (pore->isInsideSolverBox())
		  // )
		{
			if (pore->prevSolvrRes(c__water, 0., pressure, flowRate))
			{///. copied from amg_solver::writeVelocityEntry
				p_wTmp[i] = pressure;
			}
			if (pore->prevSolvrRes(c__elec, 0., volttt, flowRate))
			{///. copied from amg_solver::writeVelocityEntry
				///cout<< volttt<<endl;
				voltTmp[i] = volttt;
			}
			elem_typeTmp[i] = pore->rockIndex();

			if(pore->prevSolvrRes(c__oil, 0., pressure, flowRate))
			{
				p_oTmp[i] = pressure;
			}
		}
		pcPistonTmp[i] = pore->model()->Pc_pistonTypeAdv();
		pTmp[i] = pore->gravCorrectedEntryPress();

		double sat_p(1.-pore->waterSaturation());
		if(pore->isEntryOrExitRes())
		{										// Reservoirs and throats connected to reservoirs are always assumed single phase (?)
			sat_p = 1.;
		}
		SwTmp[i] = sat_p;
		radiusTmp[i] = pore->model()->RRR();
	}


	for(size_t i = 0; i < elmInds.size(); ++i)
	{
		 Sw[i] = SwTmp[elmInds[i]];
		 pcPiston[i] = pcPistonTmp[elmInds[i]];
		 p[i] = pTmp[elmInds[i]];
		 p_o[i] = p_oTmp[elmInds[i]];
		 p_w[i] = p_wTmp[elmInds[i]];
		 volt[i] = voltTmp[elmInds[i]];
		 elem_type[i] = elem_typeTmp[elmInds[i]];
		 radius[i] = radiusTmp[elmInds[i]];
	}
}



void getThroatSolverResults
(
	const GNMData *  netsim,
	const vector<Elem const *>&  elems_,
	const vector<size_t> & elmInds,
	vector<float> & radius,
	vector<float> & Sw,
	vector<float> & v_o,
	vector<float> & v_w,
	vector<float> & I,
	vector<float> & pcPiston,
	vector<float> & p,
	vector<float> & p_o,
	vector<float> & p_w,
	vector<float> & Volt,
	vector<float> & elem_type,
	int nBpPors_
)
{

	const Fluid& c__water = netsim->water();
	const Fluid& c__oil = netsim->oil();
	const Fluid& c__elec = netsim->elec();
	vector<float> radiusTmp(elems_.size(),0.);
	vector<float> SwTmp(elems_.size(),-1.);
	vector<float> v_oTmp(elems_.size(),-1.);
	vector<float> v_wTmp(elems_.size(),-1.);
	vector<float> ITmp(elems_.size(),-1.);
	vector<float> pcPistonTmp(elems_.size(),-1.);
	vector<float> pTmp(elems_.size(),-1.);
	vector<float> p_oTmp(elems_.size(),-1.);
	vector<float> p_wTmp(elems_.size(),-1.);
	vector<float> VoltTmp(elems_.size(),-0.5);
	vector<float> elem_typeTmp(elems_.size(),-1.);

	for(size_t i = nBpPors_; i < elems_.size(); ++i)
	{
		const Elem* throat = elems_[i];
		const Elem* pore1 = throat->neib(0);
		const Elem* pore2 = throat->neib(1);

		double pressure(-10.);
		double flowRate(-10.);
		//if(   (pore1->isInsideSolverBox())
			//||(pore2->isInsideSolverBox())
		  //)
		{
			throat->prevSolvrRes(c__elec, 0.5*(pore1->node().x+ pore2->node().x), pressure, flowRate);
			;//{///. copied from amg_solver::writeVelocityEntry
				ITmp[i] = flowRate;
				VoltTmp[i] = pressure;
			//}

			throat->prevSolvrRes(c__water, 0.5*(pore1->node().x+ pore2->node().x), pressure, flowRate);
				v_wTmp[i] = flowRate;
				p_wTmp[i] = pressure;

			elem_typeTmp[i] = throat->rockIndex();

			throat->prevSolvrRes(c__oil, 0.5*(pore1->node().x+ pore2->node().x), pressure, flowRate);
				v_oTmp[i] = flowRate;
				p_oTmp[i] = pressure;
		}
		pcPistonTmp[i] = throat->model()->Pc_pistonTypeAdv();
		pTmp[i] = throat->gravCorrectedEntryPress();
		double  sat_t(1.-throat->waterSaturation());
		if(pore1->isEntryOrExitRes())
		{										// Reservoirs and throats connected to phase
			sat_t = 1.;
		} else if(pore2->isEntryOrExitRes()) {
			sat_t = 1.;
		}

		SwTmp[i] = sat_t;
		radiusTmp[i] = throat->model()->RRR();
	}

	for(size_t i = 0; i < elmInds.size(); ++i)
	{
		Sw[i] = SwTmp[elmInds[i]];
		v_o[i] = v_oTmp[elmInds[i]];
		v_w[i] = v_wTmp[elmInds[i]];
		pcPiston[i] = pcPistonTmp[elmInds[i]];
		p[i] = pTmp[elmInds[i]];
		p_o[i] = p_oTmp[elmInds[i]];
		p_w[i] = p_wTmp[elmInds[i]];
		Volt[i] = VoltTmp[elmInds[i]];
		elem_type[i] = elem_typeTmp[elmInds[i]];
		radius[i] = radiusTmp[elmInds[i]];
	}
}






 string results3D::start(size_t nPoints, size_t nCells)
	{
		stringstream  str;
		str<<"<?xml version = \"1.\"?>\n"
		   <<"<VTKFile type = \"UnstructuredGrid\" version = \"0.1\" byte_order = \"LittleEndian\">\n"
		   <<"  <UnstructuredGrid>"
		   <<"	<Piece NumberOfPoints = \""<<nPoints<<"\" NumberOfCells = \""<<nCells<<"\" >\n";
		return str.str();
	}
	string  results3D::finish()
	{
		stringstream  str;
		str<<"	</Piece>\n"
		   <<"  </UnstructuredGrid>\n"
		   <<"</VTKFile>\n";
		return str.str();
	}




template<typename Type>
void writeCellData(ofstream& outp, string name, const vector<Type> & vals, string typeStr="Float32")
{
	outp<<"\t\t<DataArray type = \""<<typeStr<<"\" Name = \""<<name<<"\" format = \"ascii\">\n";
	for(size_t i = 0; i < vals.size(); ++i) outp<<vals[i]<<_nl_;
	outp<<"\t\t</DataArray>"<<endl;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addSpherePoreMesh
(
	 const vector<Elem const *>&  elems_,
	 size_t poreIndx,
	 vector<dbl3>& points,
	 vector<int>& tags,
	 vector<int>& elmInds,
	 vector<float>& ffaz,
	 vector<int>& cellPoints,
	 //vector<FacePoints>& facePoints, 		 //vector<CellFaces>& cellFaces, 		 //vector<FaceCells>& faceCells,
	 double scaleFactor,
	 unsigned int thetaResulution
)
{
	const Elem * elem = elems_[poreIndx];

	dbl3 c1(elem->node());
	dbl3 c2(elem->node());
	c2[1] += 1e-12;
	double r = elem->model()->RRR()*sqrt(scaleFactor);
	if (elem->isEntryOrExitRes())
		return ; ///. no visualization of inlet/outlet res


	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-32);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (ncc[i]<0.6){nCE[i] = 1.;break; } }

	nCE = ncc^nCE;
	nCE = nCE/(mag(nCE)+1e-32);


	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < thetaResulutionp2; ++i)
		{
			double hafAng = PI/thetaResulutionp2;
			double hafAngleAzim = 0.5*PI/thetaResulutionp2;
			nCE = rotateAroundVec(nCE,2*hafAng,ncc);
			dbl3 lAzimuth1 = ncc^nCE; ///. normal to ncc and nE1
			dbl3 nCE2 = rotateAroundVec(nCE,thetaResulutionp2*hafAngleAzim,lAzimuth1);///. edge-centre ncc vector
			for(int j = 0; j < thetaResulutionp2; ++j)
			{
				dbl3 nCE1 = nCE2;///. edge-centre ncc vector
				nCE2 = rotateAroundVec(nCE2,hafAngleAzim*2.,lAzimuth1);///. edge-centre ncc vector

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r, 0.,0.,   0,  hafAng,0.); ///. Warning: CA is not implemented for spheres
					tags.pbak(0);   elmInds.pbak(poreIndx);   ffaz.pbak(elem->model()->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1e-18,  r,  0.,0.,   0,  -hafAng,0.);///. Warning: CA is not implemented for spheres
				tags.pbak(0);   elmInds.pbak(poreIndx);   ffaz.pbak(elem->model()->containCOil());
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void results3D::vtuWritePores(string fnam, double Pc, double tension)
{

	vector<dbl3> points;
	vector<int> tags;//subTypes
	vector<int> elmInds;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
	points.reserve(elems_.size()*150);
	cellPoints.reserve(elems_.size()*300);
	tags.reserve(elems_.size()*50);
	elmInds.reserve(elems_.size()*50);
	vector<float> ffaz;
	ffaz.reserve(elems_.size()*50);


	for(int i = 0; i <  nBpPors_; ++i)
	{
		addSpherePoreMesh(elems_,i,points,tags,elmInds,ffaz,cellPoints,/*facePoints,cellFaces,faceCells,*/rScaleFactor_, nTheta_);
	}


	ofstream outp(fnam);

	outp<<results3D::start(points.size(),tags.size());



	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for(size_t i = 0; i < points.size(); ++i)
	{
		outp<<points[i][0]<< " "<<points[i][1]<< " "<<points[i][2]<< " "<<_nl_;
	}
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
	for(size_t i = 0; i < cellPoints.size(); ++i)
	{
		outp<<cellPoints[i]<<_nl_;
	}
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
	for(size_t i = 0; i < tags.size(); ++i)
	{
		outp<<8*i+8<<_nl_;
	}
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
	for(size_t i = 0; i < tags.size(); ++i)
	{
		outp<<12<<_nl_;
	}
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"      <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////

	vector<float> Sw(elmInds.size(),0.);
	vector<float> radius(elmInds.size(),0.);
	vector<float> pcPiston(elmInds.size(),0.);
	vector<float> p(elmInds.size(),0.);
	vector<float> p_o(elmInds.size(),0.);
	vector<float> p_w(elmInds.size(), 0.);
	vector<float> volt(elmInds.size(), 0.);
	vector<float> elem_type(elmInds.size(), 0.);
	getSolverPoreResults(comn_, elems_, elmInds, radius, Sw, pcPiston, p, p_o, p_w, volt, elem_type, nBpPors_);
	writeCellData( outp, "radius",  radius);
	writeCellData( outp, "subType",  tags);
	writeCellData( outp, "ffaz",  ffaz);
	writeCellData( outp, "Sw",  Sw);
	writeCellData( outp, "pcPiston",  pcPiston);
	writeCellData( outp, "p",  p);
	writeCellData( outp, "p_o",  p_o);
	writeCellData(outp, "p_w", p_w);
	writeCellData(outp, "volt", volt);
	writeCellData(outp, "type", elem_type);
	writeCellData( outp, "index",  elmInds, "Int32");

	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




	outp<<results3D::finish();
	outp.close();
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addThroatMesh(
	 const vector<Elem const *>&  elems_,
	 size_t trIndx,
	 vector<dbl3>& points,
	 vector<int>& tags,
	 vector<size_t>& elmInds,
	 vector<float>& ffaz,
	 vector<int>& cellPoints,
	 double scaleFactor,
	 bool visualizeCorners,
	 unsigned int thetaResulution
	 , double pc, double tension
	 )
{
	const Elem * elem = elems_[trIndx];

	const dbl3&  node1 = elem->neib(0)->node();
	dbl3 c1(node1.x,node1.y,node1.z);
	const dbl3&  node2 = elem->neib(1)->node();
	dbl3 c2(node2.x,node2.y,node2.z);
	double r = elem->model()->RRR()*scaleFactor;


	if (elem->neib(0)->isEntryOrExitRes())
	{
		double throatIncLength=max(((Throat*)elem)->throatLength(),1.1*elem->neib(1)->model()->RRR()+2e-9);
		c1[1] = c2[1];//y is wrong
		c1[2] = c2[2];
		c2[1] += 1e-9;	c2[2] += 1e-9;
		if (c1[0]<c2[0])			c1[0] = c2[0] - throatIncLength;
		else		c1[0] = c2[0] + throatIncLength;
	}
	if (elem->neib(1)->isEntryOrExitRes())
	{
		double throatIncLength=max(((Throat*)elem)->throatLength(),1.1*elem->neib(0)->model()->RRR()+2e-9);

		c2[1] = c1[1];//y is wrong
		c2[2] = c1[2];

		c2[1] += 1e-9;	c2[2] += 1e-9;
		if (c2[0]<c1[0])			c2[0] = c1[0]-  throatIncLength;
		else			c2[0] = c1[0]+ throatIncLength;

	}

	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1e-33);
	dbl3 nCE(0.,0.,0.);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (abs(ncc[i])<0.6){nCE[i] = 1.;break; } ; }

	nCE = ncc^nCE;
	nCE=rotateAroundVec(nCE, PI/4., ncc);
	nCE = nCE/(mag(nCE)+1e-33);

	const Polygon* shyp = dynamic_cast< const Polygon* >(elem->model());
	if(shyp && visualizeCorners)
	{

		for(int i = 0; i < shyp->numCorners(); ++i)
		{

			double hafAng = shyp->cornerHalfAngles(i);
			double rCornerOut = r/sin(hafAng);///. only an approximate
			double maxApexDist = r/tan(hafAng);///. only an approximate
			double lp_inrR = 0.;
			double lp_outR = 0.;			///. cornerAppex
			double hAng1 = PI/2.-hafAng;			///. cornerAppex
			double hAng2 = PI/2.-hafAng;			///. cornerAppex
			if (shyp->waterInCorner()[i].cornerExists())
			{
				shyp->waterInCorner()[i].getCApexDistConAng(lp_outR, hAng2, pc, hafAng, tension);		lp_outR *= scaleFactor;
				if (lp_outR>maxApexDist*scaleFactor) {lp_outR = maxApexDist*scaleFactor; if (lp_outR>maxApexDist*scaleFactor*1.5) cout<<"e"; }
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2, rCornerOut,  hafAng);
					tags.pbak(3);
				elmInds.pbak(trIndx);
				ffaz.pbak(0.);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2,  rCornerOut,  -hafAng);
					tags.pbak(3);
				elmInds.pbak(trIndx);
				ffaz.pbak(0.);
			}

			if (shyp->oilLayerConst()[i].exists(/*st ab le*/))
			{
				shyp->waterInCorner()[i].getCApexDistConAng(lp_inrR, hAng1, pc, hafAng, tension); 		lp_inrR *= scaleFactor;
				//lp_outR = shyp->oilLayerConst()[i].getApexDistance(pc, hAng2, hafAng, tension)*scaleFactor;
				//hAng2 = PI-shyp->oilLayerConst()[i].hingingConAng(pc, hAng2, hafAng, tension);
				shyp->oilLayerConst()[i].getCAApexDist(lp_outR, hAng2, hafAng, pc, tension);lp_outR*=scaleFactor;

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2,  rCornerOut,  hafAng);
				tags.pbak(10);
				elmInds.pbak(trIndx);
				ffaz.pbak(1.);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2,  rCornerOut,  -hafAng);
					tags.pbak(10);
				elmInds.pbak(trIndx);
				ffaz.pbak(1.);
			}

			if (lp_outR<maxApexDist-1e-16)///. corner has distance from inscribed circle
			{

				double lp_inrR = lp_outR;
				double lp_outR = maxApexDist; ///. ERROR, bad aproximation
				hAng1 = hAng2;
				hAng2 = 0.;
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2,  rCornerOut,  hafAng);
					tags.pbak(1);
				elmInds.pbak(trIndx);
				ffaz.pbak(shyp->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  lp_inrR,  lp_outR, hAng1,hAng2,  rCornerOut,  -hafAng);
					tags.pbak(1);
				elmInds.pbak(trIndx);
				ffaz.pbak(shyp->containCOil());
			}

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, PI/2.-hAng2,PI/2.-hAng2,  0,  PI/2.-hafAng);
			tags.pbak(2);
			elmInds.pbak(trIndx);
			ffaz.pbak(shyp->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, PI/2.-hAng2,PI/2.-hAng2,  0,  -PI/2.+hafAng);
			tags.pbak(2);
			elmInds.pbak(trIndx);
			ffaz.pbak(shyp->containCOil());

			nCE = rotateAroundVec(nCE,PI-shyp->cornerHalfAngles(i)-shyp->cornerHalfAngles((i+1)%(shyp->numCorners())),ncc);

		}

	}
	else
	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < 2*thetaResulutionp2; ++i)
		{
			double hafAng = 0.5*PI/thetaResulutionp2;
			nCE = rotateAroundVec(nCE,2*hafAng,ncc);
			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0., 0.6*(PI-hafAng),   0,  hafAng);
			tags.pbak(5);
			elmInds.pbak(trIndx);
			ffaz.pbak(elem->model()->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0., 0.6*(PI-hafAng),  0,  -hafAng);
			tags.pbak(5);
			elmInds.pbak(trIndx);
			ffaz.pbak(elem->model()->containCOil());
		}
	}


	  ///. TODO
	 // vector<CellPoints> cellPoints,
	 // vector<FacePoints> facePoints,
	 // vector<CellFaces> cellFaces,
	 // vector<FaceCells> faceCells

}




void results3D::vtuWriteThroats(string fnam, double pc, double tension)
{

	vector<dbl3> points;
	vector<int> tags;
	vector<size_t> elmInds;

	vector<int> cellPoints;
	points.reserve(elems_.size()*150);
	cellPoints.reserve(elems_.size()*300);
	tags.reserve(elems_.size()*50);
	elmInds.reserve(elems_.size()*50);
	vector<float> ffaz;
	ffaz.reserve(elems_.size()*50);



	for(size_t i = nBpPors_; i < elems_.size(); ++i)
	{
		addThroatMesh(elems_,i,points,tags,elmInds,ffaz,cellPoints,rScaleFactor_,_1At(v3D_AC21I_,9),nTheta_/2, pc, tension);
	}


	ofstream outp(fnam);


	outp<<results3D::start(points.size(),tags.size());

	outp<<"	  <Points>\n";
	outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
	for(size_t i = 0; i < points.size(); ++i)
	{
		outp<<points[i][0]<< " "<<points[i][1]<< " "<<points[i][2]<< " "<<_nl_;
	}
	outp<<"\n		</DataArray>\n";
	outp<<"	  </Points>\n";

	outp<<"	  <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"		<DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
	for(size_t i = 0; i < cellPoints.size(); ++i)
	{
		outp<<cellPoints[i]<<_nl_;
	}
	outp<<"\n		</DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
	for(size_t i = 0; i < tags.size(); ++i)
	{
		outp<<8*i+8<<_nl_;
	}
	outp<<"\n		</DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
	for(size_t i = 0; i < tags.size(); ++i)
	{
		outp<<12<<_nl_;
	}
	outp<<"\n		</DataArray>\n";
	outp<<"	  </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"	  <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////

	vector<float> radius(elmInds.size());
	vector<float> Sw(elmInds.size());
	vector<float> v_o(elmInds.size());
	vector<float> v_w(elmInds.size());
	vector<float> pcPiston(elmInds.size());
	vector<float> p_o(elmInds.size());
	vector<float> p(elmInds.size());
	vector<float> Ie(elmInds.size());
	vector<float> volt(elmInds.size());
	vector<float> p_w(elmInds.size());
	vector<float> elem_type(elmInds.size());
	getThroatSolverResults(comn_,elems_, elmInds, radius, Sw, v_o, v_w, Ie, pcPiston, p, p_o, p_w, volt, elem_type, nBpPors_);
	writeCellData( outp, "radius",  radius);
	writeCellData( outp, "tag",  tags);
	writeCellData( outp, "ffaz",  ffaz);
	writeCellData( outp, "Sw",  Sw);
	writeCellData( outp, "v_o",  v_o);
	writeCellData( outp, "v_w",  v_w);
	writeCellData( outp, "Ie",  Ie);
	writeCellData( outp, "pcPiston",  pcPiston);
	writeCellData( outp, "p",  p);
	writeCellData( outp, "p_o",  p_o);
	writeCellData( outp, "p_w",  p_w);
	writeCellData( outp, "volt",  volt);
	writeCellData( outp, "elem_type",  elem_type);
	writeCellData( outp, "index",  elmInds);

	outp<<"	  </CellData>\n";				   /////////////////////////////////////////


	outp<<results3D::finish();
	outp.close();
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////








/// ////////////////////////////////////   throat lines   ///////////////////////////////////////////////////////
void results3D::writeThroatLines(string fName, double pc, double tension, int icycl, double istp, bool endCycle)
{

	#define bcncs elems_[ib]->connections()

	ofstream outp(fName);
	int i, ib;
	const int nElms=elems_.size();
	vector<array<short,2> > btrotcpis(nElms,{{-1,-1}});

	outp<<results3D::start(nElms+elems_[0]->connections().size()+elems_[1]->connections().size(),nElms-(nBpPors_));


	outp<<"	  <Points>\n";
	{	outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
		//double Dx[2] = {0.05*(elems_[0]->node().x - elems_[1]->node().x), 0.05*(elems_[1]->node().x - elems_[0]->node().x)};
		for(i=0; i<nElms; ++i)
		{	dbl3 p(elems_[i]->node().x,elems_[i]->node().y,elems_[i]->node().z); outp<<p<< " ";		if (!((i+1)%20))	 outp<<'\n';	
			btrotcpis[i][0]=-1;btrotcpis[i][1]=-1; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i)
		 {	dbl3 p(bcncs[i]->node().x,bcncs[i]->node().y,bcncs[i]->node().z); p.x=elems_[ib]->node().x; outp<<p<< " ";		if (!((i+1)%20)) outp<<'\n';	btrotcpis[bcncs[i]->index()][ib]=i; }
		outp<<"\n		</DataArray>\n";
	}outp<<"	  </Points>\n";


	outp<<"	  <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	{	outp<<"		<DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
		for (i=nBpPors_; i<nElms; ++i)
		{
			int itrt=elems_[i]->index();
			int ip0=elems_[i]->neib(0)->index();
			int ip1=elems_[i]->neib(1)->index();

			if (ip1==1) ip1 = btrotcpis[itrt][1]+nElms+elems_[0]->connections().size();
			if (ip1==0) ip1 = btrotcpis[itrt][0]+nElms;
			if (ip0==0) ip0 = btrotcpis[itrt][0]+nElms;
			if (ip0==1) ip0 = btrotcpis[itrt][1]+nElms+elems_[0]->connections().size();

			outp<<ip0<<" "<<ip1<<" "<<itrt<<_nl_;
		}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
		for (i=0; i<nElms-(nBpPors_); ++i)
		{	outp<<3*i+3<<_nl_;	}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
		for(i=0; i<nElms-(nBpPors_); ++i)		{	outp<<21<<_nl_;	}
		outp<<"\n		</DataArray>\n";
	}outp<<"	  </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


	outp<<"	  <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"RRR\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->RRR()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<ffloatof(elems_[i])<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->saturation()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Pc\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<pc<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->index()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i) { outp<<0<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condW\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i)
			{  outp<<elems_[i]->poreToPoreCond(WTR)*1e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condO\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i)
			{  outp<<elems_[i]->poreToPoreCond(OIL)*1e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condE\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i)
			{  outp<<elems_[i]->poreToPoreCond(ELEC)*1e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qo\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->flowRate(OIL)*1e18<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qw\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->flowRate(WTR)*1e18<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i)
			//{  const auto& ec=elems_[i]->slvrCnct(WTR); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i)
			//{  auto ec=elems_[i]->slvrCnct(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"insideSatBox\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->isInCalcBox()<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->trappingCL().first<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->AmotIndxOil()<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </CellData>\n"; /////////////////////////////////////////////////////////


	outp<<"\n	  <PointData>\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"radius\" format = \"ascii\">\n";
		for(i=0; i<nElms; ++i) { outp<<elems_[i]->RRR()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i){ outp<<bcncs[i]->RRR()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L1\" format = \"ascii\">\n";
//		for(i=0; i<nBpPors_; ++i) { outp<<elems_[i]->node()*0.<<_nl_; }
//		for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->neib(0)->node()-elems_[i]->node()<<"  "; _nl_; }
//		outp<<"\n		</DataArray>\n";
//
//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L2\" format = \"ascii\">\n";
//		for(i=0; i<nBpPors_; ++i) { outp<<elems_[i]->node()*0.<<_nl_; }
//		for(i=nBpPors_; i<nElms; ++i) { outp<<elems_[i]->neib(1)->node()-elems_[i]->node()<<"  "; _nl_; }
//		outp<<"\n		</DataArray>\n";

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=0; i<nElms; ++i) { outp<<ffloatof(elems_[i])<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i){ outp<<ffloatof(elems_[ib])<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"iFPEntry\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i) { outp<<elems_[i]->inFrontOf(invf)*elems_[i]->iFPEntry()*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->inFrontOf(invf)*elems_[ib]->iFPEntry()*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"PcEntry\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i) { outp<<elems_[i]->inFrontOf(invf)*elems_[i]->PcEntry()*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->inFrontOf(invf)*elems_[ib]->PcEntry()*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=0; i<nElms; ++i) { outp<<elems_[i]->index()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->index()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i) { outp<<elems_[i]->rockIndex()+0.1*double(size_t(elems_[i]->index())>=nBpPors_)+0.1*(double(elems_[i]->index()<2))<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<0.1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_o\" format = \"ascii\">\n";
		//for(i=0; i<nBpPors_; ++i)       { outp<<dynamic_cast<const Pore*>(elems_[i ])->solverPrs(OIL)*tension<<_nl_; }
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->avgPressure(OIL)*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i) { outp<<dynamic_cast<const Pore*>(elems_[ib])->solverPrs(OIL)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_w\" format = \"ascii\">\n";
		//for(i=0; i<nBpPors_; ++i)      { outp<<dynamic_cast<const Pore*>(elems_[i])->solverPrs(WTR)*tension<<_nl_; }
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->avgPressure(WTR)*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<dynamic_cast<const Pore*>(elems_[ib])->solverPrs(WTR)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_o\" format = \"ascii\">\n";
		//for(i=0; i<nBpPors_; ++i)       { outp<<dynamic_cast<const Pore*>(elems_[i ])->solverErr(OIL)*tension<<_nl_; }
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->avgPressure(OIL)*0<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i) { outp<<dynamic_cast<const Pore*>(elems_[ib])->solverErr(OIL)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_w\" format = \"ascii\">\n";
		//for(i=0; i<nBpPors_; ++i)      { outp<<dynamic_cast<const Pore*>(elems_[i])->solverErr(WTR)*tension<<_nl_; }
		//for(i=nBpPors_; i<nElms; ++i) { outp<<dynamic_cast<const Throat*>(elems_[i])->avgPressure(WTR)*0<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<dynamic_cast<const Pore*>(elems_[ib])->solverErr(WTR)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i)  { outp<<float(elems_[i]->AmotIndxOil())<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i) { outp<<float(elems_[ib]->AmotIndxOil())<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i)	{auto& ec=elems_[i]->slvrCnct(WTR); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){auto& ec=elems_[ib]->slvrCnct(WTR); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i)	{ auto& ec=elems_[i]->slvrCnct(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ auto& ec=elems_[ib]->slvrCnct(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagE\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i)	{ auto& ec=elems_[i]->slvrCnct(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ auto& ec=elems_[ib]->slvrCnct(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=0; i<nElms; ++i) { outp<<elems_[i]->trappingCL().first<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->trappingCL().first<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=0; i<nElms; ++i) { outp<<elems_[i]->waterSaturation()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->waterSaturation()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"cachInd\" format = \"ascii\">\n";
		//for(i=0; i<nElms; ++i) { outp<<elems_[i]->cachInd<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<int(bcncs.size()); ++i){ outp<<elems_[ib]->cachInd<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </PointData>\n"; /////////////////////////////////////////////////////////


	outp<<results3D::finish();
	outp.close();
}









#undef bcncs




//---------------------------------------------------------------------'
