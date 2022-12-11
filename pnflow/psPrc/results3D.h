#ifndef RESULTS3D_CnM_H // possible conflict with gnflow
#define RESULTS3D_CnM_H

/*---------------------------------------------------------------------------*\
Developed by: Ali Q Raeini (2015-2019)
\*---------------------------------------------------------------------------*/

#include <memory>
#include <sstream>
#include "Element.h"



class GNMData;





class results3D
{

public:
	results3D(const InputFile& input,const GNMData* comn, ststr outputfolder, const std::vector<Elem const *>*  elems, size_t nBP2=0, size_t n6pPors=0);
	void init(size_t nBP2, size_t n6pPors) { nBSs_=nBP2 ; nBpPors_=n6pPors; };

	void write3D(double pc, double intfacTen, bool endCycle = false);

			private: ///. depricated

				ststr  start(size_t nPoints, size_t nCells);
				ststr  finish();
				void vtuWritePores(ststr suffix, double pc, double intfacTen);
				void vtuWriteThroats(ststr suffix, double pc, double intfacTen);
				void writeThroatLines(ststr fName, double pc, double tension, int icycl, double tstp, bool endCycle);

	void writeThroatLinesXmf(ststr suffix, double pc, double tension, int icycl, double tstp, bool endCycle);


private:

	const GNMData*  comn_;
	const std::vector<Elem const*>&  elems_;
	int nBSs_;
	int nBpPors_;

	enum Format { 	 VTU_FORMAT=1, 	 XMF_FORMAT=2, 	 BOTH_FORMATS=3  };
	Format format_;
	unsigned int      rkw_;  //. rx4+keepx2+writex1
	unsigned int      vLinAC21I_;
	unsigned int      v3D_AC21I_;
	unsigned int   iWrite_;
	ststr          prefix_;
	double         rScaleFactor_;
	unsigned int   nTheta_;


public :

	bool        inform;

};





#endif
