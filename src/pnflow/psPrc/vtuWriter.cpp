
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
#define _nl_  if(!((i+1)%20)) outp<<"\n"


#if defined _MSC_VER
#include <direct.h>
#elif defined __GNUC__
#include <sys/types.h>
#include <sys/stat.h>
#endif




void vtuWriter::vtuWrite( const vector<Element const *> *  elems, size_t m_numPores, double pc, double intfacTen)
{

	string suffix;
	if      (!iWrite)
	{	if (!m_visualise[0] )return ;
		else suffix = myto_string(m_comn->floodingCycle())+"_Init";
	}else if (m_comn->oilInjection() )
	{	if (!m_visualise[1] )return ;
		else suffix = myto_string(m_comn->floodingCycle())+"_OInj";
	}else if ( ! m_comn->oilInjection() )
	{	if (!m_visualise[2] )return ;
		else suffix = myto_string(m_comn->floodingCycle())+"_WInj";
	}else suffix = "Cycle";


	cout<<" visua";cout.flush();
	if(m_FullOrLight[0]=='T')
	{
	vtuWritePores( suffix+"Pore",  elems, m_numPores);
	cout<<"liza";cout.flush();
	vtuWriteThroats(suffix+"Throat", elems, m_numPores, pc, intfacTen);
	cout<<"tion "<<endl;
	}

	vtuWriteThroatLines(suffix, *elems, m_numPores, pc, intfacTen);
	iWrite++;
}


float ffTofloat(const Element* slm){  float ff=slm->ffaz();
    if(ff>3.0) {ff=(4.0*float(int(ff)%4)+1.5)/5;}
    //if(!slm->azad()) ff=(4.0*ff-1.5)/3;  
    return max(0.6f,min(ff,2.4f));    }
	// {	fluidf ff(slm->ffaz());if(ff&OIL) return 1+0.025*(ff-OIL);   else if(ff&WATER)  return 0.025*(ff-WATER);   else  return -0.1; }


dbl3 rotate( dbl3 n, dbl3 x, dbl3 y, double gamma)
{///. rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	return dbl3
	  (	( x[0]*(n[1]*n[1]+n[2]*n[2]) - n[0]*( x[1]*n[1]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (-x[2]*n[1]+x[1]*n[2]-n[2]*y[1]+n[1]*y[2] )*sin(gamma),
		( x[1]*(n[0]*n[0]+n[2]*n[2]) - n[1]*( x[0]*n[0]+x[2]*n[2]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + ( x[2]*n[0]-x[0]*n[2]+n[2]*y[0]-n[0]*y[2] )*sin(gamma),
		( x[2]*(n[0]*n[0]+n[1]*n[1]) - n[2]*( x[0]*n[0]+x[1]*n[1]-n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (-x[1]*n[0]+x[0]*n[1]-n[1]*y[0]+n[0]*y[1] )*sin(gamma)
	  );
}
dbl3 rotateAroundVec( dbl3 n, dbl3 y, double gamma)
{///. rotate y around n (line passing through centre, in the direction of n) http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	return dbl3
	  (	(  - n[0]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[0]*cos(gamma) + (n[1]*y[2]-n[2]*y[1])*sin(gamma),
		(  - n[1]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[1]*cos(gamma) + (n[2]*y[0]-n[0]*y[2])*sin(gamma),
		(  - n[2]*( -n[0]*y[0]- n[1]*y[1]-n[2]*y[2] ) )*(1-cos(gamma)) + y[2]*cos(gamma) + (n[0]*y[1]-n[1]*y[0])*sin(gamma)
	  );
	//-n*(n&y)*(1-cos(gamma)) + y*cos(gamma) + n^y*sin(gamma);
}




int findOrInsertPoint(vector<dbl3>&  points, dbl3& point)
{

	for (vector<dbl3>::const_reverse_iterator  rip = points.rbegin(); rip != points.rbegin()+200 && rip != points.rend(); ++rip)
		if (*rip == point) return int(points.rend()-rip)-1;

	points.push_back(point);
	return points.size()-1;
}



void insertHalfCorneroints(vector<dbl3>&  points, vector<int>& cellPoints, dbl3 c1, dbl3 c2, dbl3 nE1, dbl3 nE2, double appexDist_inrR, double appexDist_outR,
							double conAng1, double conAng2, double rOuter, double hafAng, double CARelax = 1.0)
{
	//const double convertToRad = _pi/180.0;
	dbl3 dd = c2-c1;
	dbl3 normal = dd/(mag(dd)+1.0e-32);

	vector<dbl3> hcPoints(8);

	dbl3 e1 = c1+nE1*rOuter;///. edgePoint
	dbl3 nE11 = rotateAroundVec(normal,nE1,hafAng);

	///. radii of curvature

	double gama = abs(hafAng)*CARelax;
	double lc_InnerR = appexDist_inrR*(cos(gama)+(sin(gama)/cos(gama+conAng1))*(sin(gama+conAng1)-1));
	double lc_OuterR = appexDist_outR*(cos(gama)+(sin(gama)/cos(gama+conAng2))*(sin(gama+conAng2)-1));

	if(hafAng>0)
	{
		hcPoints[0] = e1-lc_InnerR*nE1;
		hcPoints[1] = e1-lc_OuterR*nE1;
		hcPoints[2] = e1-appexDist_outR*nE11;
		hcPoints[3] = e1-appexDist_inrR*nE11;
	}
	else
	{
		hcPoints[3] = e1-lc_InnerR*nE1;
		hcPoints[2] = e1-lc_OuterR*nE1;
		hcPoints[1] = e1-appexDist_outR*nE11;
		hcPoints[0] = e1-appexDist_inrR*nE11;
	}



	 e1 = c2+nE2*rOuter;///. edgePoint
	 nE11 = rotateAroundVec(normal,nE2,hafAng);
	if(hafAng>0)
	{
		hcPoints[4] = e1-lc_InnerR*nE2;
		hcPoints[5] = e1-lc_OuterR*nE2;
		hcPoints[6] = e1-appexDist_outR*nE11;
		hcPoints[7] = e1-appexDist_inrR*nE11;
	}
	else
	{
		hcPoints[7] = e1-lc_InnerR*nE2;
		hcPoints[6] = e1-lc_OuterR*nE2;
		hcPoints[5] = e1-appexDist_outR*nE11;
		hcPoints[4] = e1-appexDist_inrR*nE11;
	}


	for (int i=0;i<8;++i)
	{///. 8 points each elem
		int id=(findOrInsertPoint(points, hcPoints[i]));
		cellPoints.push_back(id);
	}

}





void getSolverPoreResults
(
	const Netsim *  netsim,
	const vector<Element const *> *  elems,
	const vector<int> & cellPores,
	vector<float> & saturation,
	vector<float> & pcPiston,
	vector<float> & p,
	vector<float> & p_o,
	vector<float> & p_w,
	vector<float> & volt,
	vector<float> & elem_type,
	size_t m_numPores
)
{

	const Water* m_c_water = &netsim->water();
	const Oil* m_c_oil = &netsim->oil();
	vector<float> saturationTmp((*elems).size(),-1.0);
	vector<float> pcPistonTmp((*elems).size(),-1.0);
	vector<float> pTmp((*elems).size(),-1.0);
	vector<float> p_oTmp((*elems).size(),-1.0);
	vector<float> p_wTmp((*elems).size(),-1.0);
	vector<float> voltTmp((*elems).size(),-1.0);
	vector<float> elem_typeTmp((*elems).size(),-1.0);
    for(size_t i = 1; i <= m_numPores; ++i)//i < (*elems).size()
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
		saturationTmp[i] = sat_p;
	}


    for(size_t i = 0; i < cellPores.size(); ++i)
    {
		 saturation[i] = saturationTmp[cellPores[i]];
		 pcPiston[i] = pcPistonTmp[cellPores[i]];
		 p[i] = pTmp[cellPores[i]];
		 p_o[i] = p_oTmp[cellPores[i]];
		 p_w[i] = p_wTmp[cellPores[i]];
		 volt[i] = voltTmp[cellPores[i]];
		 elem_type[i] = elem_typeTmp[cellPores[i]];
	}
}



void getThroatSolverResults
(
	const Netsim *  netsim,
	const vector<Element const *> *  elems,
	const vector<size_t> & cellPores,
	vector<float> & saturation,
	vector<float> & v_o,
	vector<float> & v_w,
	vector<float> & I,
	vector<float> & pcPiston,
	vector<float> & p,
	vector<float> & p_o,
	vector<float> & p_w,
	vector<float> & Volt,
	vector<float> & elem_type,
	int m_numPores
)
{

	const Water * m_c_water = &netsim->water();
	const Oil * m_c_oil = &netsim->oil();
	vector<float> saturationTmp((*elems).size(),-1.0);
	vector<float> v_oTmp((*elems).size(),-1.0);
	vector<float> v_wTmp((*elems).size(),-1.0);
	vector<float> ITmp((*elems).size(),-1.0);
	vector<float> pcPistonTmp((*elems).size(),-1.0);
	vector<float> pTmp((*elems).size(),-1.0);
	vector<float> p_oTmp((*elems).size(),-1.0);
	vector<float> p_wTmp((*elems).size(),-1.0);
	vector<float> VoltTmp((*elems).size(),-0.5);
	vector<float> elem_typeTmp((*elems).size(),-1.0);

    for(size_t i = m_numPores+2; i < (*elems).size(); ++i)
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

		saturationTmp[i] = sat_t;
	}

	for(size_t i = 0; i < cellPores.size(); ++i)
	{
		saturation[i] = saturationTmp[cellPores[i]];
		v_o[i] = v_oTmp[cellPores[i]];
		v_w[i] = v_wTmp[cellPores[i]];
		pcPiston[i] = pcPistonTmp[cellPores[i]];
		p[i] = pTmp[cellPores[i]];
		p_o[i] = p_oTmp[cellPores[i]];
		p_w[i] = p_wTmp[cellPores[i]];
		Volt[i] = VoltTmp[cellPores[i]];
		elem_type[i] = elem_typeTmp[cellPores[i]];
	}
}





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
		if ((i+1)%20 == 0)     outp << "\n";
	}
	outp<<"        </DataArray>"<<endl;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addSpherePoreMesh
(
	 const vector<Element const *> *  elems,
	 size_t poreIndx,
	 vector<dbl3>& points,
	 vector<int>& subTypes,
	 vector<int>& cellPores,
	 vector<float>& alpha,
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
			nCE = rotateAroundVec(ncc,nCE,2*hafAng);
			dbl3 lAzimuth1 = ncc^nCE; ///. normal to ncc and nE1
			dbl3 nCE2 = rotateAroundVec(lAzimuth1,nCE,thetaResulutionp2*hafAngleAzim);///. edge-centre ncc vector
			for(int j = 0; j < thetaResulutionp2; ++j)
			{
				dbl3 nCE1 = nCE2;///. edge-centre ncc vector
				nCE2 = rotateAroundVec(lAzimuth1,nCE2,hafAngleAzim*2.0);///. edge-centre ncc vector

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1.0e-18,  r, 0.0,0.0,   0,  hafAng,0.0); ///. Warning: CA is not implemented for spheres
					subTypes.push_back(0);
				cellPores.push_back(poreIndx);
				alpha.push_back(elem->model()->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE1,  -nCE2,  1.0e-18,  r,  0.0,0.0,   0,  -hafAng,0.0);///. Warning: CA is not implemented for spheres
				subTypes.push_back(0);
				cellPores.push_back(poreIndx);
				alpha.push_back(elem->model()->containCOil());
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void vtuWriter::vtuWritePores(string suffix, const vector<Element const *> *  elems, size_t m_numPores)
{

	vector<dbl3> points;
	vector<int> subTypes;
	vector<int> cellPores;

	vector<int> cellPoints;
	//vector<FacePoints> facePoints;
	//vector<CellFaces> cellFaces;
	//vector<FaceCells> faceCells;
	points.reserve((*elems).size()*150);
	cellPoints.reserve((*elems).size()*300);
	subTypes.reserve((*elems).size()*50);
	cellPores.reserve((*elems).size()*50);
	vector<float> alpha;
	alpha.reserve((*elems).size()*50);


    for(size_t i = 0; i <  m_numPores+2; ++i)
    {
        addSpherePoreMesh(elems,i,points,subTypes,cellPores,alpha,cellPoints,/*facePoints,cellFaces,faceCells,*/m_rScaleFactor, m_thetaResulution);
    }


	stringstream fileNamepp;
    fileNamepp<< m_fileNamePrefix<<suffix<<"_"<<100+vtuWriter::iWrite<<".vtu";
    ofstream outp(fileNamepp.str().c_str());

	outp<<vtuWriter::start(points.size(),subTypes.size());



	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)
    {
        outp << points[i][0]<< " " << points[i][1]<< " " << points[i][2]<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)
    {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)
    {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)
    {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"      <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////

	vector<float> saturation(cellPores.size(),0.0);
	vector<float> pcPiston(cellPores.size(),0.0);
	vector<float> p(cellPores.size(),0.0);
	vector<float> p_o(cellPores.size(),0.0);
	vector<float> p_w(cellPores.size(), 0.0);
	vector<float> volt(cellPores.size(), 0.0);
	vector<float> elem_type(cellPores.size(), 0.0);
	getSolverPoreResults(m_comn, elems, cellPores, saturation, pcPiston, p, p_o, p_w, volt, elem_type, m_numPores);
	writeCellData( outp, "subType",  subTypes);
	writeCellData( outp, "alpha",  alpha);
	writeCellData( outp, "saturation",  saturation);
	writeCellData( outp, "pcPiston",  pcPiston);
	writeCellData( outp, "p",  p);
	writeCellData( outp, "p_o",  p_o);
	writeCellData(outp, "p_w", p_w);
	writeCellData(outp, "volt", volt);
	writeCellData(outp, "type", elem_type);
	writeCellData( outp, "index",  cellPores, "Int32");

	outp<<"      </CellData>\n"; /////////////////////////////////////////////////////////




    outp<<vtuWriter::finish();
    outp.close();
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void addThroatMesh(
	 const vector<Element const *> *  elems,
	 size_t trIndx,
	 vector<dbl3>& points,
	 vector<int>& subTypes,
	 vector<size_t>& cellPores,
	 vector<float>& alpha,
	 vector<int>& cellPoints,
	 double scaleFactor,
	 bool visualizeCorners,
	 unsigned int thetaResulution
	 , double pc, double intfacTen
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
	nCE=rotateAroundVec( ncc, nCE, _pi/4.0);
	nCE = nCE/(mag(nCE)+1.0e-33);

    const Polygon* polyShape = dynamic_cast< const Polygon* >(elem->model());
    if(polyShape && visualizeCorners)
    {

		for(int i = 0; i < polyShape->numCorners(); ++i)
        {

			double hafAng = polyShape->cornerHalfAngles(i);
			double rCornerOut = r/sin(hafAng);///. only an approximate
			double maxApexDist = r/tan(hafAng);///. only an approximate
			double appexDist_inrR = 0.0;
			double appexDist_outR = 0.0;			///. cornerAppex
			double conAng1 = _pi/2.0-hafAng;			///. cornerAppex
			double conAng2 = _pi/2.0-hafAng;			///. cornerAppex
			if (polyShape->waterInCorner()[i].cornerExists())
			{
				polyShape->waterInCorner()[i].getCApexDistConAng(appexDist_outR, conAng2, pc, hafAng, intfacTen);		appexDist_outR *= scaleFactor;
				if (appexDist_outR>maxApexDist*scaleFactor) {appexDist_outR = maxApexDist*scaleFactor; if (appexDist_outR>maxApexDist*scaleFactor*1.5) cout<<"e";}
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2, rCornerOut,  hafAng);
					subTypes.push_back(3);
				cellPores.push_back(trIndx);
				alpha.push_back(0.0);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					subTypes.push_back(3);
				cellPores.push_back(trIndx);
				alpha.push_back(0.0);
			}

			if (polyShape->oilLayerConst()[i].exists(/*st ab le*/))
			{
				polyShape->waterInCorner()[i].getCApexDistConAng(appexDist_inrR, conAng1, pc, hafAng, intfacTen); 		appexDist_inrR *= scaleFactor;
				//appexDist_outR = polyShape->oilLayerConst()[i].getApexDistance(pc, conAng2, hafAng, intfacTen)*scaleFactor;
				//conAng2 = _pi-polyShape->oilLayerConst()[i].hingingConAng(pc, conAng2, hafAng, intfacTen);
				polyShape->oilLayerConst()[i].getCAApexDist(appexDist_outR, conAng2, hafAng, pc, intfacTen);appexDist_outR*=scaleFactor;

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2,  rCornerOut,  hafAng);
				subTypes.push_back(10);
				cellPores.push_back(trIndx);
				alpha.push_back(1.0);

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					subTypes.push_back(10);
				cellPores.push_back(trIndx);
				alpha.push_back(1.0);
			}

			if (appexDist_outR<maxApexDist-1.0e-16)///. corner has distance from inscribed circle
			{

				double appexDist_inrR = appexDist_outR;
				double appexDist_outR = maxApexDist; ///. ERROR, bad aproximation
				conAng1 = conAng2;
				conAng2 = 0.0;
				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2,  rCornerOut,  hafAng);
					subTypes.push_back(1);
				cellPores.push_back(trIndx);
				alpha.push_back(polyShape->containCOil());

				insertHalfCorneroints( points,cellPoints,c1,c2,  nCE, nCE,  appexDist_inrR,  appexDist_outR, conAng1,conAng2,  rCornerOut,  -hafAng);
					subTypes.push_back(1);
				cellPores.push_back(trIndx);
				alpha.push_back(polyShape->containCOil());
			}

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, _pi/2.0-conAng2,_pi/2.0-conAng2,  0,  _pi/2.0-hafAng);
			subTypes.push_back(2);
			cellPores.push_back(trIndx);
			alpha.push_back(polyShape->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r, _pi/2.0-conAng2,_pi/2.0-conAng2,  0,  -_pi/2.0+hafAng);
			subTypes.push_back(2);
			cellPores.push_back(trIndx);
			alpha.push_back(polyShape->containCOil());

			nCE = rotateAroundVec(ncc,nCE,_pi-polyShape->cornerHalfAngles(i)-polyShape->cornerHalfAngles((i+1)%(polyShape->numCorners())));

		}

	}
	else
	{
		int thetaResulutionp2 = (thetaResulution+1)/2;

		for(int i = 0; i < 2*thetaResulutionp2; ++i)
        {
			double hafAng = 0.5*_pi/thetaResulutionp2;
			nCE = rotateAroundVec(ncc,nCE,2*hafAng);
			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0.0, 0.6*(_pi-hafAng),   0,  hafAng);
			subTypes.push_back(5);
			cellPores.push_back(trIndx);
			alpha.push_back(elem->model()->containCOil());

			insertHalfCorneroints( points,cellPoints,c1,c2,  -nCE, -nCE,  0,  r,  0.0, 0.6*(_pi-hafAng),  0,  -hafAng);
			subTypes.push_back(5);
			cellPores.push_back(trIndx);
			alpha.push_back(elem->model()->containCOil());
		}
	}


	  ///. TODO
	 // vector<CellPoints> cellPoints,
	 // vector<FacePoints> facePoints,
	 // vector<CellFaces> cellFaces,
	 // vector<FaceCells> faceCells

}




void vtuWriter::vtuWriteThroats(string suffix, const vector<Element const *> *  elems, size_t m_numPores, double pc, double intfacTen)
{

	vector<dbl3> points;
	vector<int> subTypes;
	vector<size_t> cellPores;

	vector<int> cellPoints;
    points.reserve((*elems).size()*150);
    cellPoints.reserve((*elems).size()*300);
    subTypes.reserve((*elems).size()*50);
    cellPores.reserve((*elems).size()*50);
	vector<float> alpha;
    alpha.reserve((*elems).size()*50);



    for(size_t i = m_numPores+2; i < (*elems).size(); ++i)
    {
			addThroatMesh(elems,i,points,subTypes,cellPores,alpha,cellPoints,/*facePoints,cellFaces,faceCells,*/m_rScaleFactor,m_visualise[3],m_thetaResulution/2, pc, intfacTen);
    }


	stringstream fileNamepp;
	fileNamepp<< m_fileNamePrefix<<suffix<<"_"<<100+vtuWriter::iWrite<<".vtu";
	ofstream outp(fileNamepp.str().c_str());

	outp<<vtuWriter::start(points.size(),subTypes.size());



	outp<<"      <Points>\n";
	outp<<"        <DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"ascii\">\n";
    for(size_t i = 0; i < points.size(); ++i)
    {
        outp << points[i][0]<< " " << points[i][1]<< " " << points[i][2]<< " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Points>\n";

	outp<<"      <Cells>\n";///// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/
	outp<<"        <DataArray type = \"Int32\" Name = \"connectivity\" format = \"ascii\">\n";
    for(size_t i = 0; i < cellPoints.size(); ++i)
    {
        outp << cellPoints[i] << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }

	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)
    {
        outp << 8*i+8 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";

	outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
    for(size_t i = 0; i < subTypes.size(); ++i)
    {
        outp << 12 << " ";
        if ((i+1)%20 == 0)     outp << "\n";
    }
	outp<<"\n        </DataArray>\n";
	outp<<"      </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//



	outp<<"	  <CellData Scalars = \"alpha\">\n"; //////////////////////////////////////

	vector<float> saturation(cellPores.size());
	vector<float> v_o(cellPores.size());
	vector<float> v_w(cellPores.size());
	vector<float> pcPiston(cellPores.size());
	vector<float> p(cellPores.size());
	vector<float> p_o(cellPores.size());
	vector<float> Ie(cellPores.size());
	vector<float> volt(cellPores.size());
	vector<float> p_w(cellPores.size());
	vector<float> elem_type(cellPores.size());
	getThroatSolverResults(m_comn,elems, cellPores, saturation, v_o, v_w, Ie, pcPiston, p, p_o, p_w, volt, elem_type, m_numPores);
	writeCellData( outp, "subType",  subTypes);
	writeCellData( outp, "alpha",  alpha);
	writeCellData( outp, "saturation",  saturation);
	writeCellData( outp, "v_o",  v_o);
	writeCellData( outp, "v_w",  v_w);
	writeCellData( outp, "Ie",  Ie);
	writeCellData( outp, "pcPiston",  pcPiston);
	writeCellData( outp, "p",  p);
	writeCellData( outp, "p_o",  p_o);
	writeCellData( outp, "p_w",  p_w);
	writeCellData( outp, "volt",  volt);
	writeCellData( outp, "elem_type",  elem_type);
	writeCellData( outp, "index",  cellPores);

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
		{	dbl3 p(elems[i]->node()->xPos(),elems[i]->node()->yPos(),elems[i]->node()->zPos()); outp<<p<< " ";		if (!((i+1)%20))	 outp<<"\n";	btrotcpis[i][0]=-1;btrotcpis[i][1]=-1;}
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i)
		 {	dbl3 p(bcncs[i]->node()->xPos(),bcncs[i]->node()->yPos(),bcncs[i]->node()->zPos()); p.x=elems[ib]->node()->xPos()+Dx[ib]; outp<<p<< " ";		if (!((i+1)%20)) outp<<"\n";	btrotcpis[bcncs[i]->index()][ib]=i;}
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

			outp<<ip0<<" "<<ip1<<" "<<itrt<<" ";	_nl_;	
		}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"Int32\" Name = \"offsets\" format = \"ascii\">\n";
		for (i=0; i<elems.size()-(nPors+2); ++i)
		{	outp<<3*i+3<<" ";	_nl_;	}
		outp<<"\n		</DataArray>\n";

		outp<<"	<DataArray type = \"UInt8\" Name = \"types\" format = \"ascii\">\n";
		for(i=0; i<elems.size()-(nPors+2); ++i)		{	outp<<21<<" ";	_nl_;	}
		outp<<"\n		</DataArray>\n";
	}outp<<"	  </Cells>\n";// //  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//


	outp<<"	  <CellData Scalars = \"ffaz\">\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"RRR\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->RRR()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<ffTofloat(elems[i])<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->saturation()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Pc\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<pc<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->index()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<0<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condW\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[0]*1.0e18<<" ";   _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condO\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[1]*1.0e18<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"condE\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i)
			{  outp<<elems[i]->m_conductance[2]*1.0e18<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qo\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->flowRate(OIL)*1.0e18<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"qw\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->flowRate(WATER)*1.0e18<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i)
			//{  const auto& ec=elems[i]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";   _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i)
			//{  auto ec=elems[i]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"insideSatBox\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->isInCalcBox()<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->trappingCL().first<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->AmotIndxOil()<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </CellData>\n"; /////////////////////////////////////////////////////////


	outp<<"\n	  <PointData>\n"; //////////////////////////////////////
	{
		outp<<"		<DataArray type = \"Float32\" Name = \"radius\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->RRR()<<" ";  _nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<bcncs[i]->RRR()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L1\" format = \"ascii\">\n";
//		for(i=0; i<nPors+2; ++i) { outp<<elems[i]->node()*0.0<<" ";  _nl_; }
//		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->connection(0)->node()-elems[i]->node()<<"  ";  _nl_; }
//		outp<<"\n		</DataArray>\n";
//
//		outp<<"		<DataArray type = \"Float32\" NumberOfComponents = \"3\"  Name = \"L2\" format = \"ascii\">\n";
//		for(i=0; i<nPors+2; ++i) { outp<<elems[i]->node()*0.0<<" ";  _nl_; }
//		for(i=nPors+2; i<elems.size(); ++i) { outp<<elems[i]->connection(1)->node()-elems[i]->node()<<"  ";  _nl_; }
//		outp<<"\n		</DataArray>\n";

		outp<<"		<DataArray type = \"Float32\" Name = \"ffaz\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<ffTofloat(elems[i])<<" ";  _nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<ffTofloat(elems[ib])<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"iFPEntry\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->inFrontOf(invf)*elems[i]->iFPEntry()*tension<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->inFrontOf(invf)*elems[ib]->iFPEntry()*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"PcEntry\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->inFrontOf(invf)*elems[i]->PcEntry()*tension<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->inFrontOf(invf)*elems[ib]->PcEntry()*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"index\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->index()<<" ";  _nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->index()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"type\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->rockIndex()+0.2*double(size_t(elems[i]->index())>=nPors+2)+0.1*(double(elems[i]->index()<2))<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<0.1<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_o\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)       { outp<<dynamic_cast<const Pore*>(elems[i ])->solverPrs(OIL)*tension<<" ";  _nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(OIL)*tension<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<dynamic_cast<const Pore*>(elems[ib])->solverPrs(OIL)*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"p_w\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)      { outp<<dynamic_cast<const Pore*>(elems[i])->solverPrs(WATER)*tension<<" ";  _nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(WATER)*tension<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<dynamic_cast<const Pore*>(elems[ib])->solverPrs(WATER)*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_o\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)       { outp<<dynamic_cast<const Pore*>(elems[i ])->solverResidual(OIL)*tension<<" ";  _nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(OIL)*0<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<dynamic_cast<const Pore*>(elems[ib])->solverResidual(OIL)*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"Err_w\" format = \"ascii\">\n";
		//for(i=0; i<nPors+2; ++i)      { outp<<dynamic_cast<const Pore*>(elems[i])->solverResidual(WATER)*tension<<" ";  _nl_; }
		//for(i=nPors+2; i<elems.size(); ++i) { outp<<dynamic_cast<const Throat*>(elems[i])->avgPressure(WATER)*0<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<dynamic_cast<const Pore*>(elems[ib])->solverResidual(WATER)*tension<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"oilWetability\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)  { outp<<float(elems[i]->AmotIndxOil())<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i) { outp<<float(elems[ib]->AmotIndxOil())<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;


		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagW\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{auto& ec=elems[i]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){auto& ec=elems[ib]->solverConect(WATER); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagO\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{ auto& ec=elems[i]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ auto& ec=elems[ib]->solverConect(OIL); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"solverflagE\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i)	{ auto& ec=elems[i]->solverConect(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ auto& ec=elems[ib]->solverConect(ELEC); outp<<ec.isPassed()*4+ec.isSearched()*2+ec.isInSlvrBox()*1<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"trpIndx\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->trappingCL().first<<" ";  _nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->trappingCL().first<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		outp<<"		<DataArray type = \"Float32\" Name = \"Sw\" format = \"ascii\">\n";
		for(i=0; i<elems.size(); ++i) { outp<<elems[i]->waterSaturation()<<" ";  _nl_; }
		for(ib=0; ib<2; ++ib)
		 for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->waterSaturation()<<" ";  _nl_; }
		outp<<"\n		</DataArray>"<<endl;

		//outp<<"		<DataArray type = \"Float32\" Name = \"cachInd\" format = \"ascii\">\n";
		//for(i=0; i<elems.size(); ++i) { outp<<elems[i]->cachInd<<" ";  _nl_; }
		//for(ib=0; ib<2; ++ib)
		 //for(i=0; i<bcncs.size(); ++i){ outp<<elems[ib]->cachInd<<" ";  _nl_; }
		//outp<<"\n		</DataArray>"<<endl;

	}outp<<"	  </PointData>\n"; /////////////////////////////////////////////////////////



	outp<<vtuWriter::finish();
	outp.close();
}










// ************************************************************************* //
