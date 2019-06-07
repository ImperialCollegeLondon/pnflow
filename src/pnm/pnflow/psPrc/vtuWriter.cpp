
/*---------------------------------------------------------------------------*\
2015:  Developed by Ali Q Raeini  email: a.qaseminejad-raeini09@imperial.ac.uk
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

using namespace std;


#include "typses.h"
#include "vtuWriter.h"
#include "netsim.h"

#define _pi   3.14159265359



#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/types.h>
#include <sys/stat.h>
#endif





vtuWriter::vtuWriter(const string& KeywordData,const Netsim * netsim, const string title_res): m_comn(netsim), m_fileNamePrefix("res3D"),m_rScaleFactor(1.0),iWrite(0)
{
	for (int i = 0;i<3;i++) m_visualise[i] = false;

	if (!KeywordData.empty())
	{ 
		m_FullOrLight ="Light";
		cout<<"Initialising VTK visualization";
		stringstream data;
		data<<KeywordData;
		data>>m_FullOrLight>>m_rScaleFactor ;
		m_fileNamePrefix=title_res;
		if (*m_fileNamePrefix.rbegin()=='/' || *m_fileNamePrefix.rbegin()=='\\')
		{
				cout<<"creating folder: "<<m_fileNamePrefix
				#if defined(_WIN32)
					<< mkdir(m_fileNamePrefix.c_str()) // check also _mkdir
				#else 
					<< mkdir(m_fileNamePrefix.c_str(), 0733) // notice that 777 is different than 0777
				#endif
				<<endl;
		}
		data>>m_thetaResulution;
		if (m_thetaResulution<3 )
			cout<< "\n\n   Error: not sufficient angular resolution"<<m_thetaResulution<<"\n\n"<<endl;
		else if (m_thetaResulution>18 )
			cout<< "\n\n   Warning: too high resolution, visualization files will be too large\n\n"<<endl;

		for (int i = 0;i<4;i++)
		{
			string tmp;
			data>>tmp;
			if (tmp[0] == 'T' || tmp[0] == 't')        m_visualise[i] = true;
			else if (tmp[0] == 'F' || tmp[0] == 'f')	m_visualise[i] = false;
			else cout<< "\n\nError :  wrong data in input file \"" << tmp
			<<"\".  Only T(rue) or F(alse) is accepted\n\n"<<endl;
		}

		if ((m_rScaleFactor<0.99 || m_rScaleFactor>1.01) && m_visualise[0])
			cout<< "\n\n  Info: pore and throat radii will be multiplied by "<<m_rScaleFactor<<" in visualization files \n"<<endl;


		iWrite = 0;
	}


	cout<<endl;

}


void vtuWriter::vtuWrite( const vector<Element const *> *  elems, size_t nPors, double pc, double tension)
{

	string suffix;
	if      (!iWrite)
	{	if (!m_visualise[0] )return ;
		else suffix = toStr(m_comn->floodingCycle())+"_Init";
	}else if (m_comn->oilInjection() )
	{	if (!m_visualise[1] )return ;
		else suffix = toStr(m_comn->floodingCycle())+"_OInj";
	}else if ( ! m_comn->oilInjection() )
	{	if (!m_visualise[2] )return ;
		else suffix = toStr(m_comn->floodingCycle())+"_WInj";
	}else suffix = "Cycle";


	cout<<" visua";cout.flush();
	if(m_FullOrLight[0]=='T')
	{
	vtuWritePores( suffix+"Pore",  elems, nPors);
	cout<<"liza";cout.flush();
	vtuWriteThroats(suffix+"Throat", elems, nPors, pc, tension);
	cout<<"tion "<<endl;
	}

	vtuWriteThroatLines(suffix, *elems, nPors, pc, tension);
	iWrite++;
}

#define _nl_  ((i&31)==31 ? '\n' : ' ')

#define pbak push_back

float ffloatof(const Element* slm){  float ff=slm->ffaz();
    if(ff>3.0) {ff=(4.0*float(int(ff)%4)+1.5)/5;}
    return max(0.6f,min(ff,2.4f));    }






int findOrInsertPoint(vector<dbl3>&  points, dbl3& point)
{

	for (vector<dbl3>::const_reverse_iterator  rip = points.rbegin(); rip != points.rbegin()+200 && rip != points.rend(); ++rip)
		if (*rip == point) return int(points.rend()-rip)-1;

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


void insertHalfCorneroints(vector<dbl3>&  points, vector<int>& cellPoints, dbl3 c1, dbl3 c2, dbl3 nE1, dbl3 nE2, double apex_inrR, double apex_outR,
							double conAng1, double conAng2, double rOuter, double hafAng, double CARelax = 1.0)
{
	//const double convertToRad = _pi/180.0;
	dbl3 dd = c2-c1;
	dbl3 normal = dd/(mag(dd)+1.0e-32);

	vector<dbl3> hcPoints(8);

	///. radii of curvature
	double gama = abs(hafAng)*CARelax;
	double lc_iInnerR = apex_inrR*(cos(gama)+(sin(gama)/cos(gama+conAng1))*(sin(gama+conAng1)-1));
	double lc_iOuterR = apex_outR*(cos(gama)+(sin(gama)/cos(gama+conAng2))*(sin(gama+conAng2)-1));

	dbl3 e1 = c1+nE1*rOuter;///. edgePoint
	dbl3 nE11 = rotateAroundVec(nE1,hafAng,normal);

	if(hafAng>0)
	{
		hcPoints[0] = e1-lc_iInnerR*nE1;
		hcPoints[1] = e1-lc_iOuterR*nE1;
		hcPoints[2] = e1-apex_outR*nE11;
		hcPoints[3] = e1-apex_inrR*nE11;
	}
	else
	{
		hcPoints[3] = e1-lc_iInnerR*nE1;
		hcPoints[2] = e1-lc_iOuterR*nE1;
		hcPoints[1] = e1-apex_outR*nE11;
		hcPoints[0] = e1-apex_inrR*nE11;
	}



	 e1 = c2+nE2*rOuter;///. edgePoint
	 nE11 = rotateAroundVec(nE2,hafAng,normal);
	if(hafAng>0)
	{
		hcPoints[4] = e1-lc_iInnerR*nE2;
		hcPoints[5] = e1-lc_iOuterR*nE2;
		hcPoints[6] = e1-apex_outR*nE11;
		hcPoints[7] = e1-apex_inrR*nE11;
	}
	else
	{
		hcPoints[7] = e1-lc_iInnerR*nE2;
		hcPoints[6] = e1-lc_iOuterR*nE2;
		hcPoints[5] = e1-apex_outR*nE11;
		hcPoints[4] = e1-apex_inrR*nE11;
	}

	///. 8 points each elem 
	for (int i=0;i<8;++i)	cellPoints.pbak(findOrInsertPoint(points, hcPoints[i]));

}





void getSolverPoreResults
(
	const Netsim *  netsim,
	const vector<Element const *> *  elems,
	const vector<int> & elmInds,
	vector<float> & radius,
	vector<float> & Sw,
	vector<float> & pcPiston,
	vector<float> & p,
	vector<float> & p_o,
	vector<float> & p_w,
	vector<float> & volt,
	vector<float> & elem_type,
	size_t nPors
)
{

	const Water* m_c_water = &netsim->water();
	const Oil* m_c_oil = &netsim->oil();
	vector<float> RTmp((*elems).size(),0.0);
	vector<float> SwTmp((*elems).size(),-1.0);
	vector<float> pcPistonTmp((*elems).size(),-1.0);
	vector<float> pTmp((*elems).size(),-1.0);
	vector<float> p_oTmp((*elems).size(),-1.0);
	vector<float> p_wTmp((*elems).size(),-1.0);
	vector<float> voltTmp((*elems).size(),-1.0);
	vector<float> elem_typeTmp((*elems).size(),-1.0);
    for(size_t i = 1; i <= nPors; ++i)//i < (*elems).size()
	{
		const Element* pore = (*elems)[i];


		double pressure(-10.0);
		double volttt(-10.0);
		double flowRate(-10.0);
        // if(   (pore->isInsideSolverBox())
          // )
		{
			int resistSolve = 0;
			if (pore->prevSolvrRes(m_c_water, resistSolve, 0.0, pressure, flowRate))
			{///. copied from amg_solver::writeVelocityEntry
				p_wTmp[i] = pressure;
			}
			if (pore->prevSolvrRes(m_c_water, 2, 0.0, volttt, flowRate))
			{///. copied from amg_solver::writeVelocityEntry
				///cout<< volttt << endl;
				voltTmp[i] = volttt;
			}
			elem_typeTmp[i] = pore->iRockType();

			if(pore->prevSolvrRes(m_c_oil, resistSolve, 0.0, pressure, flowRate))
			{
				p_oTmp[i] = pressure;
			}
		}
		pcPistonTmp[i] = pore->model()->Pc_pistonTypeAdv();
		pTmp[i] = pore->model()->gravCorrectedEntryPress();

		double sat_p(1.0-pore->waterSaturation());
		if(pore->isEntryOrExitRes())
		{			                            // Reservoirs and throats connected to reservoirs are always assumed single phase (?)
			sat_p = 1.0;
		}
		SwTmp[i] = sat_p;
		RTmp[i] = pore->model()->radius();
	}


    for(size_t i = 0; i < elmInds.size(); ++i)
    {
		radius[i] = RTmp[elmInds[i]];
		 Sw[i] = SwTmp[elmInds[i]];
		 pcPiston[i] = pcPistonTmp[elmInds[i]];
		 p[i] = pTmp[elmInds[i]];
		 p_o[i] = p_oTmp[elmInds[i]];
		 p_w[i] = p_wTmp[elmInds[i]];
		 volt[i] = voltTmp[elmInds[i]];
		 elem_type[i] = elem_typeTmp[elmInds[i]];
	}
}



void getThroatSolverResults
(
	const Netsim *  netsim,
	const vector<Element const *> *  elems,
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
	int nPors
)
{

	const Water * m_c_water = &netsim->water();
	const Oil * m_c_oil = &netsim->oil();
	vector<float> RTmp((*elems).size(),0.0);
	vector<float> SwTmp((*elems).size(),-1.0);
	vector<float> v_oTmp((*elems).size(),-1.0);
	vector<float> v_wTmp((*elems).size(),-1.0);
	vector<float> ITmp((*elems).size(),-1.0);
	vector<float> pcPistonTmp((*elems).size(),-1.0);
	vector<float> pTmp((*elems).size(),-1.0);
	vector<float> p_oTmp((*elems).size(),-1.0);
	vector<float> p_wTmp((*elems).size(),-1.0);
	vector<float> VoltTmp((*elems).size(),-0.5);
	vector<float> elem_typeTmp((*elems).size(),-1.0);

    for(size_t i = nPors+2; i < (*elems).size(); ++i)
	{
		const Element* throat = (*elems)[i];
		const Element* pore1 = throat->connection(0);
		const Element* pore2 = throat->connection(1);

		double pressure(-10.0);
		double flowRate(-10.0);
        //if(   (pore1->isInsideSolverBox())
            //||(pore2->isInsideSolverBox())
          //)
		{
			bool resistSolve = 2;
			throat->prevSolvrRes(m_c_water, resistSolve, 0.5*(pore1->node()->xPos()+ pore2->node()->xPos()), pressure, flowRate);
			;//{///. copied from amg_solver::writeVelocityEntry
				ITmp[i] = flowRate;
				VoltTmp[i] = pressure;
			//}

			resistSolve = 0;
			throat->prevSolvrRes(m_c_water, resistSolve, 0.5*(pore1->node()->xPos()+ pore2->node()->xPos()), pressure, flowRate);
				v_wTmp[i] = flowRate;
				p_wTmp[i] = pressure;

			elem_typeTmp[i] = throat->iRockType();

			throat->prevSolvrRes(m_c_oil, resistSolve, 0.5*(pore1->node()->xPos()+ pore2->node()->xPos()), pressure, flowRate);
				v_oTmp[i] = flowRate;
				p_oTmp[i] = pressure;
		}
		pcPistonTmp[i] = throat->model()->Pc_pistonTypeAdv();
		pTmp[i] = throat->model()->gravCorrectedEntryPress();
		double  sat_t(1.0-throat->waterSaturation());
		if(pore1->isEntryOrExitRes())
		{			                            // Reservoirs and throats connected to phase
			sat_t = 1.0;
		} else if(pore2->isEntryOrExitRes()) {
			sat_t = 1.0;
		}

		SwTmp[i] = sat_t;
		RTmp[i] = throat->model()->radius();
	}

	for(size_t i = 0; i < elmInds.size(); ++i)
	{
		radius[i] = RTmp[elmInds[i]];
		Sw[i] = SwTmp[elmInds[i]];
		v_o[i] = v_oTmp[elmInds[i]];
		v_w[i] = v_wTmp[elmInds[i]];
		pcPiston[i] = pcPistonTmp[elmInds[i]];
		p[i] = pTmp[elmInds[i]];
		p_o[i] = p_oTmp[elmInds[i]];
		p_w[i] = p_wTmp[elmInds[i]];
		Volt[i] = VoltTmp[elmInds[i]];
		elem_type[i] = elem_typeTmp[elmInds[i]];
	}
}






 string vtuWriter::start(size_t nPoints, size_t nCells)
	{
		stringstream  str;
		str<<"<?xml version = \"1.0\"?>\n"
		   <<"<VTKFile type = \"UnstructuredGrid\" version = \"0.1\" byte_order = \"LittleEndian\">\n"
		   <<"  <UnstructuredGrid>"
		   <<"	<Piece NumberOfPoints = \""<<nPoints<<"\" NumberOfCells = \""<<nCells<<"\" >\n";
		return str.str();
	}
	 string  vtuWriter::finish()
	{
		stringstream  str;
		str<<"	</Piece>\n"
		   <<"  </UnstructuredGrid>\n"
		   <<"</VTKFile>\n";
		return str.str();
	}




template<typename Type>
void writeCellData(ofstream& outp, string name, const vector<Type> & data, string typeStr="Float32")
{
	outp<<"        <DataArray type = \""<<typeStr<<"\" Name = \""<<name<<"\" format = \"ascii\">\n";
	for(size_t i = 0; i < data.size(); ++i)
	{
		outp << data[i] << " ";
		if ((i+1)%20 == 0)     outp<<'\n';
	}
	outp<<"        </DataArray>"<<endl;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addSpherePoreMesh
(
	 const vector<Element const *> *  elems,
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
	const Element * elem = (*elems)[poreIndx];

	const Node *  node1 = elem->node();
	dbl3 c1(node1->xPos(),node1->yPos(),node1->zPos());
	const Node *  node2 = elem->node();
	dbl3 c2(node2->xPos(),node2->yPos(),node2->zPos());
	c2[1] += 1.0e-12;
	double r = elem->model()->radius()*sqrt(scaleFactor);
	if (elem->isEntryOrExitRes())
		return ; ///. no visualization of inlet/outlet res


	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1.0e-32);
	dbl3 nCE(0.0,0.0,0.0);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (ncc[i]<0.6){nCE[i] = 1.0;break;} }

	nCE = ncc^nCE;
	nCE = nCE/(mag(nCE)+1.0e-32);

    const Polygon* polyShape = dynamic_cast< const Polygon* >(elem->model());
    if(0 && polyShape) ///. for now ignore corners
    { }
	else
	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < thetaResulutionp2; ++i)
		{
			double hafAng = _pi/thetaResulutionp2;
			double hafAngleAzim = 0.5*_pi/thetaResulutionp2;
			nCE = rotateAroundVec(nCE,2*hafAng,ncc);
			dbl3 lAzimuth1 = ncc^nCE; ///. normal to ncc and nE1
			dbl3 nCE2 = rotateAroundVec(nCE,thetaResulutionp2*hafAngleAzim,lAzimuth1);///. edge-centre ncc vector
			for(int j = 0; j < thetaResulutionp2; ++j)
			{
				dbl3 nCE1 = nCE2;///. edge-centre ncc vector
				nCE2 = rotateAroundVec(nCE2,hafAngleAzim*2.0,lAzimuth1);///. edge-centre ncc vector

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1.0e-18,  r, 0.0,0.0,   0,  hafAng,0.0); ///. Warning: CA is not implemented for spheres
					tags.pbak(0);
				elmInds.pbak(poreIndx);
				ffaz.pbak(elem->model()->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1.0e-18,  r,  0.0,0.0,   0,  -hafAng,0.0);///. Warning: CA is not implemented for spheres
				tags.pbak(0);
				elmInds.pbak(poreIndx);
				ffaz.pbak(elem->model()->containCOil());
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void vtuWriter::vtuWritePores(string suffix, const vector<Element const *> *  elems, size_t nPors)
{

	vector<dbl3> points;
	vector<int> tags;
	vector<int> elmInds;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
	points.reserve((*elems).size()*150);
	cellPoints.reserve((*elems).size()*300);
	tags.reserve((*elems).size()*50);
	elmInds.reserve((*elems).size()*50);
	vector<float> ffaz;
	ffaz.reserve((*elems).size()*50);


    for(size_t i = 0; i <  nPors+2; ++i)
    {
        addSpherePoreMesh(elems,i,points,tags,elmInds,ffaz,cellPoints,/*facePoints,cellFaces,faceCells,*/m_rScaleFactor, m_thetaResulution);
    }


	stringstream fileNamepp;
    fileNamepp<< m_fileNamePrefix<<suffix<<"_"<<100+vtuWriter::iWrite<<".vtu";
    ofstream outp(fileNamepp.str().c_str());

	outp<<vtuWriter::start(points.size(),tags.size());



	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)
    {
        outp << points[i][0]<< " " << points[i][1]<< " " << points[i][2]<< " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)
    {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < tags.size(); ++i)
    {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < tags.size(); ++i)
    {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"      <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////

	vector<float> Sw(elmInds.size(),0.0);
	vector<float> radius(elmInds.size(),0.0);
	vector<float> pcPiston(elmInds.size(),0.0);
	vector<float> p(elmInds.size(),0.0);
	vector<float> p_o(elmInds.size(),0.0);
	vector<float> p_w(elmInds.size(), 0.0);
	vector<float> volt(elmInds.size(), 0.0);
	vector<float> elem_type(elmInds.size(), 0.0);
	getSolverPoreResults(m_comn, elems, elmInds, radius, Sw, pcPiston, p, p_o, p_w, volt, elem_type, nPors);
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




    outp<<vtuWriter::finish();
    outp.close();
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addThroatMesh(
	 const vector<Element const *> *  elems,
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
	const Element * elem = (*elems)[trIndx];

	const Node *  node1 = elem->connection(0)->node();
	dbl3 c1(node1->xPos(),node1->yPos(),node1->zPos());
	const Node *  node2 = elem->connection(1)->node();
	dbl3 c2(node2->xPos(),node2->yPos(),node2->zPos());
	double r = elem->model()->radius()*scaleFactor;


	if (elem->connection(0)->isEntryOrExitRes())
	{
		double throatIncLength=max(((Throat*)elem)->length(),1.1*elem->connection(1)->model()->radius()+2e-9);
		c1[1] = c2[1];//y is wrong
		c1[2] = c2[2];
		c2[1] += 1e-9;	c2[2] += 1e-9;
		if (c1[0]<c2[0])			c1[0] = c2[0] - throatIncLength;
		else		c1[0] = c2[0] + throatIncLength;
	}
	if (elem->connection(1)->isEntryOrExitRes())
	{
		double throatIncLength=max(((Throat*)elem)->length(),1.1*elem->connection(0)->model()->radius()+2e-9);

		c2[1] = c1[1];//y is wrong
		c2[2] = c1[2];

		c2[1] += 1e-9;	c2[2] += 1e-9;
		if (c2[0]<c1[0])			c2[0] = c1[0]-  throatIncLength;
		else			c2[0] = c1[0]+ throatIncLength;

	}

	dbl3 c1c2 = c2-c1;
	dbl3 ncc = c1c2/(mag(c1c2)+1.0e-33);
	dbl3 nCE(0.0,0.0,0.0);	///. pick a corner point for the first subElem
	for(size_t i = 0;i<3;i++){ if (abs(ncc[i])<0.6){nCE[i] = 1.0;break;} ;}

	nCE = ncc^nCE;
	nCE=rotateAroundVec(nCE, _pi/4.0, ncc);
	nCE = nCE/(mag(nCE)+1.0e-33);

    const Polygon* polyShape = dynamic_cast< const Polygon* >(elem->model());
    if(polyShape && visualizeCorners)
    {

		for(int i = 0; i < polyShape->numCorners(); ++i)
        {

			double hafAng = polyShape->cornerHalfAngles(i);
			double rCornerOut = r/sin(hafAng);///. only an approximate
			double maxApexDist = r/tan(hafAng);///. only an approximate
			double apex_inrR = 0.0;
			double apex_outR = 0.0;			///. cornerAppex
			double conAng1 = _pi/2.0-hafAng;			///. cornerAppex
			double conAng2 = _pi/2.0-hafAng;			///. cornerAppex
			if (polyShape->waterInCorner()[i].cornerExists())
			{
				polyShape->waterInCorner()[i].getCApexDistConAng(apex_outR, conAng2, pc, hafAng, tension);		apex_outR *= scaleFactor;
				if (apex_outR>maxApexDist*scaleFactor) {apex_outR = maxApexDist*scaleFactor; if (apex_outR>maxApexDist*scaleFactor*1.5) cout<<"e";}
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2, rCornerOut,  hafAng);
					tags.pbak(3);
				elmInds.pbak(trIndx);
				ffaz.pbak(0.0);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					tags.pbak(3);
				elmInds.pbak(trIndx);
				ffaz.pbak(0.0);
			}

			if (polyShape->oilLayerConst()[i].exists(/*st ab le*/))
			{
				polyShape->waterInCorner()[i].getCApexDistConAng(apex_inrR, conAng1, pc, hafAng, tension); 		apex_inrR *= scaleFactor;
				//apex_outR = polyShape->oilLayerConst()[i].getApexDistance(pc, conAng2, hafAng, tension)*scaleFactor;
				//conAng2 = _pi-polyShape->oilLayerConst()[i].hingingConAng(pc, conAng2, hafAng, tension);
				polyShape->oilLayerConst()[i].getCAApexDist(apex_outR, conAng2, hafAng, pc, tension);apex_outR*=scaleFactor;

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2,  rCornerOut,  hafAng);
				tags.pbak(10);
				elmInds.pbak(trIndx);
				ffaz.pbak(1.0);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					tags.pbak(10);
				elmInds.pbak(trIndx);
				ffaz.pbak(1.0);
			}

			if (apex_outR<maxApexDist-1.0e-16)///. corner has distance from inscribed circle
			{

				double apex_inrR = apex_outR;
				double apex_outR = maxApexDist; ///. ERROR, bad aproximation
				conAng1 = conAng2;
				conAng2 = 0.0;
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2,  rCornerOut,  hafAng);
					tags.pbak(1);
				elmInds.pbak(trIndx);
				ffaz.pbak(polyShape->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  apex_inrR,  apex_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					tags.pbak(1);
				elmInds.pbak(trIndx);
				ffaz.pbak(polyShape->containCOil());
			}

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, _pi/2.0-conAng2,_pi/2.0-conAng2,  0,  _pi/2.0-hafAng);
			tags.pbak(2);
			elmInds.pbak(trIndx);
			ffaz.pbak(polyShape->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, _pi/2.0-conAng2,_pi/2.0-conAng2,  0,  -_pi/2.0+hafAng);
			tags.pbak(2);
			elmInds.pbak(trIndx);
			ffaz.pbak(polyShape->containCOil());

			nCE = rotateAroundVec(nCE,_pi-polyShape->cornerHalfAngles(i)-polyShape->cornerHalfAngles((i+1)%(polyShape->numCorners())),ncc);

		}

	}
	else
	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < 2*thetaResulutionp2; ++i)
        {
			double hafAng = 0.5*_pi/thetaResulutionp2;
			nCE = rotateAroundVec(nCE,2*hafAng,ncc);
			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0.0, 0.6*(_pi-hafAng),   0,  hafAng);
			tags.pbak(5);
			elmInds.pbak(trIndx);
			ffaz.pbak(elem->model()->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0.0, 0.6*(_pi-hafAng),  0,  -hafAng);
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




void vtuWriter::vtuWriteThroats(string suffix, const vector<Element const *> *  elems, size_t nPors, double pc, double tension)
{

	vector<dbl3> points;
	vector<int> tags;
	vector<size_t> elmInds;

	vector<int> cellPoints;
    points.reserve((*elems).size()*150);
    cellPoints.reserve((*elems).size()*300);
    tags.reserve((*elems).size()*50);
    elmInds.reserve((*elems).size()*50);
	vector<float> ffaz;
    ffaz.reserve((*elems).size()*50);



    for(size_t i = nPors+2; i < (*elems).size(); ++i)
    {
			addThroatMesh(elems,i,points,tags,elmInds,ffaz,cellPoints,m_rScaleFactor,m_visualise[3],m_thetaResulution/2, pc, tension);
    }


	stringstream fileNamepp;
	fileNamepp<< m_fileNamePrefix<<suffix<<"_"<<100+vtuWriter::iWrite<<".vtu";
	ofstream outp(fileNamepp.str().c_str());

	outp<<vtuWriter::start(points.size(),tags.size());

	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)
    {
        outp << points[i][0]<< " " << points[i][1]<< " " << points[i][2]<< " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)
    {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < tags.size(); ++i)
    {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < tags.size(); ++i)
    {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp<<'\n';
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



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
	getThroatSolverResults(m_comn,elems, elmInds, radius, Sw, v_o, v_w, Ie, pcPiston, p, p_o, p_w, volt, elem_type, nPors);
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

	outp<<"	  </CellData>\n"; /////////////////////////////////////////////////////////


    outp<<vtuWriter::finish();
    outp.close();
}







////////////////////////////////////////////////////////////////////////////////////////////////////////////////





/// ////////////////////////////////////   throat lines   ///////////////////////////////////////////////////////
void vtuWriter::vtuWriteThroatLines(string fName, const vector<Element const *> & elems, size_t nPors, double pc, double tension)
{

	#define bcncs elems[ib]->connections()

	stringstream fileNamepp;
	fileNamepp<< m_fileNamePrefix<<fName<<100+vtuWriter::iWrite<<".vtu";
	ofstream outp(fileNamepp.str().c_str());


	size_t i, ib;
	vector<array<short,2> > btrotcpis(elems.size(),{{-1,-1}});

	outp<<vtuWriter::start(elems.size()+elems[0]->connections().size()+elems[1]->connections().size(),elems.size()-(nPors+2));


	outp<<"	  <Points>\n";
	{	outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
		double Dx[2] = {0.05*(elems[0]->node()->xPos() - elems[1]->node()->xPos()), 0.05*(elems[1]->node()->xPos() - elems[0]->node()->xPos())};
		for(i=0; i<elems.size(); ++i)
		{	dbl3 p(elems[i]->node()->xPos(),elems[i]->node()->yPos(),elems[i]->node()->zPos()); outp<<p<< " ";		if (!((i+1)%20))	 outp<<'\n';	btrotcpis[i][0]=-1;btrotcpis[i][1]=-1;}
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i)
		 {	dbl3 p(bcncs[i]->node()->xPos(),bcncs[i]->node()->yPos(),bcncs[i]->node()->zPos()); p.x=elems[ib]->node()->xPos()+Dx[ib]; outp<<p<< " ";		if (!((i+1)%20)) outp<<'\n';	btrotcpis[bcncs[i]->index()][ib]=i;}
		outp<<"\n		</DataArray>\n";
	}outp<<"	  </Points>\n";


	outp<<"	  <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	{	outp<<"		<DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
		for (i=nPors+2; i<elems.size(); ++i)
		{
			int itrt=elems[i]->index();
			int ip0=elems[i]->connection(0)->index();
			int ip1=elems[i]->connection(1)->index();

			if (ip1==1) ip1 = btrotcpis[itrt][1]+elems.size()+elems[0]->connections().size();
			if (ip1==0) ip1 = btrotcpis[itrt][0]+elems.size();
			if (ip0==0) ip0 = btrotcpis[itrt][0]+elems.size();
			if (ip0==1) ip0 = btrotcpis[itrt][1]+elems.size()+elems[0]->connections().size();

			outp<<ip0<<" "<<ip1<<" "<<itrt<<_nl_;	
		}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
		for (i=0; i<elems.size()-(nPors+2); ++i)
		{	outp<<3*i+3<<_nl_;	}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
		for(i=0; i<elems.size()-(nPors+2); ++i)		{	outp<<21<<_nl_;	}
		outp<<"\n		</DataArray>\n";
	}outp<<"	  </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


	outp<<"	  <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"RRR\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->RRR()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<ffloatof(elems[i])<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->saturation()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Pc\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<pc<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->index()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<0<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condW\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[0]*1.0e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condO\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[1]*1.0e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condE\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[2]*1.0e18<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qo\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->flowRate(OIL)*1.0e18<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qw\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->flowRate(WATER)*1.0e18<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i)
			//{  const auto& ec=elems[i]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";   _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i)
			//{  auto ec=elems[i]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"insideSatBox\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->isInCalcBox()<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->trappingCL().first<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->AmotIndxOil()<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </CellData>\n"; /////////////////////////////////////////////////////////


	outp<<"\n	  <PointData>\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"radius\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->RRR()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<bcncs[i]->RRR()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L1\" format = \"ascii\">\n";
//		for(i=0; i<nPors+2; ++i) { outp<<elems[i]->node()*0.0<<_nl_; }
//		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->connection(0)->node()-elems[i]->node()<<"  ";  _nl_; }
//		outp<<"\n		</DataArray>\n";
//
//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L2\" format = \"ascii\">\n";
//		for(i=0; i<nPors+2; ++i) { outp<<elems[i]->node()*0.0<<_nl_; }
//		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->connection(1)->node()-elems[i]->node()<<"  ";  _nl_; }
//		outp<<"\n		</DataArray>\n";

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<ffloatof(elems[i])<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<ffloatof(elems[ib])<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"iFPEntry\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->inFrontOf(invf)*elems[i]->iFPEntry()*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->inFrontOf(invf)*elems[ib]->iFPEntry()*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"PcEntry\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->inFrontOf(invf)*elems[i]->PcEntry()*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->inFrontOf(invf)*elems[ib]->PcEntry()*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->index()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->index()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->rockIndex()+0.2*double(size_t(elems[i]->index())>=nPors+2)+0.1*(double(elems[i]->index()<2))<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<0.1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_o\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)       { outp<<dynamic_cast<const Pore*>(elems[i ])->solverPrs(OIL)*tension<<_nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(OIL)*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<dynamic_cast<const Pore*>(elems[ib])->solverPrs(OIL)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_w\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)      { outp<<dynamic_cast<const Pore*>(elems[i])->solverPrs(WATER)*tension<<_nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(WATER)*tension<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<dynamic_cast<const Pore*>(elems[ib])->solverPrs(WATER)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_o\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)       { outp<<dynamic_cast<const Pore*>(elems[i ])->solverResidual(OIL)*tension<<_nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(OIL)*0<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<dynamic_cast<const Pore*>(elems[ib])->solverResidual(OIL)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_w\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)      { outp<<dynamic_cast<const Pore*>(elems[i])->solverResidual(WATER)*tension<<_nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(WATER)*0<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<dynamic_cast<const Pore*>(elems[ib])->solverResidual(WATER)*tension<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)  { outp<<float(elems[i]->AmotIndxOil())<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<float(elems[ib]->AmotIndxOil())<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{auto& ec=elems[i]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){auto& ec=elems[ib]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{ auto& ec=elems[i]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ auto& ec=elems[ib]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagE\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{ auto& ec=elems[i]->solverConect(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ auto& ec=elems[ib]->solverConect(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->trappingCL().first<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->trappingCL().first<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->waterSaturation()<<_nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->waterSaturation()<<_nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"cachInd\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->cachInd<<_nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->cachInd<<_nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </PointData>\n"; /////////////////////////////////////////////////////////



	outp<<vtuWriter::finish();
	outp.close();
}










// ************************************************************************* //
