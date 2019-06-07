#ifndef RESULTS3D_H
#define RESULTS3D_H

/*---------------------------------------------------------------------------*\
Developed by (2015): Ali Q Raeini  email: a.qaseminejad-raeini09@imperial.ac.uk
\*---------------------------------------------------------------------------*/


#include <sstream>
#include "Element.h"

typedef   std::vector<int> CellPoints;
typedef std::vector<int> FacePoints;
typedef std::vector<int> CellFaces;
typedef std::vector<int> FaceCells;

class Netsim;


class results3D
{

public:
	results3D(const std::string& KeywordData,const Netsim * netsim, std::string outputfolder);

	void vtuWrite(const std::vector<Element const *> *  rockLattices, size_t m_numPores, double pc, double intfacTen);

	std::string fileNamePrefix;
	std::string start(size_t nPoints, size_t nCells);
	std::string  finish();
	void vtuWritePores(std::string suffix, const std::vector<Element const *> *  rockLattices, size_t m_numPores);
	void vtuWriteThroats(std::string suffix, const std::vector<Element const *> *  rockLattices, size_t m_numPores, double pc, double intfacTen);
	void vtuWriteThroatLines(string fName, const vector<Element const *> & elems, size_t nPors, double pc, double tension);

private:

	const Netsim * m_comn;
	
	bool visualise_[4];
	std::string FullOrLight_;
	std::string fileNamePrefix_;
	double rScaleFactor_;
	unsigned int thetaResulution_;
	unsigned int iWrite;
};





#endif

