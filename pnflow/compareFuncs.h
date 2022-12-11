#ifndef COMPARE_FUNC_H
#define COMPARE_FUNC_H

#ifdef WIN32
#pragma warning(disable:4786)
#endif




#include "Element.h"


class ElemRadCmpInc
{ public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{   return (elem1->RRR() < elem2->RRR());    }

	bool operator() (const VoidElem* elem1, const VoidElem* elem2) const
	{   return (elem1->RRR() < elem2->RRR());    }
};

class ElemRadCmpRed
{
public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{  return (elem1->RRR() > elem2->RRR());  }

	bool operator() (const VoidElem* elem1, const VoidElem* elem2) const
	{  return (elem1->RRR() > elem2->RRR());  }
};

class ElemGCmpInc
{ public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{   return (pow(elem1->RRR(),3)/elem1->flowVolumeX() < pow(elem2->RRR(),3)/elem2->flowVolumeX());    }
};


class ElemVolCmpRed
{
public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{
		return (elem1->flowVolume() > elem2->flowVolume());
	}
};

class ElemGCmpRed
{
public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{
		return (elem1->model()->shapeFactor() > elem2->model()->shapeFactor());
	}
};


class PceImbCmp
{
public:
	bool operator() (const Apex* elem1, const Apex* elem2) const
	{
		double pc1(elem1->gravCorrectedEntryPress());
		double pc2(elem2->gravCorrectedEntryPress());
		if(pc1 == pc2)
			return elem1 < elem2;
		else
			return pc1 < pc2;
	}
};

class PceDrainCmp
{
public:
	bool operator() (const Apex* elem1, const Apex* elem2) const
	{
		double pc1(elem1->gravCorrectedEntryPress());
		double pc2(elem2->gravCorrectedEntryPress());
		if(pc1 == pc2)
			return elem1 < elem2;
		else
			return pc1 > pc2;
	}
};





class FracWettInc
{
public:
	bool operator() (std::pair<double, Elem*> prop1, std::pair<double, Elem*> prop2) const
	{
		return (prop1.first < prop2.first);
	}
};

class FracWettDec
{
public:
	bool operator() (std::pair<double, Elem*> prop1, std::pair<double, Elem*> prop2) const
	{
		return (prop1.first > prop2.first);
	}
};

class TrappingWatStorageCmp
{
public:
	bool operator() (std::pair<Elem*, FluidBlob> elm1, std::pair<Elem*, FluidBlob> elm2) const
	{
		return (elm1.first->index() > elm2.first->index());
	}
};

/**
// In trapping routine we want to advance to the exits as fast as possible. In
// order to achieve this we want to always first retrieve the element cloasest to
// the exit first. Connecting rock elements are theerfor sorted according to
// distance to exit. Throats do not have a position associated. We must therefor
// identify the common pore and select the throat where the other pore is closest
// to the exit.
*/
class DistToExitCompareThroats
{
public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{
		const Elem* pOneOne = elem1->neib(0);
		const Elem* pOneTwo = elem1->neib(1);
		const Elem* pTwoOne = elem2->neib(0);
		const Elem* pTwoTwo = elem2->neib(1);

		if(pOneOne == pTwoOne)
			return pOneTwo->node().x > pTwoTwo->node().x;
		else if(pOneOne == pTwoTwo)
			return pOneTwo->node().x > pTwoOne->node().x;
		else if(pOneTwo == pTwoOne)
			return pOneOne->node().x > pTwoTwo->node().x;
		else
			return pOneOne->node().x > pTwoOne->node().x;
	}
};

/**
// Comparing pores we directly obtain a position from their node
*/
class DistToExitComparePores
{
public:
	bool operator() (const Elem* elem1, const Elem* elem2) const
	{
		return elem1->node().x > elem2->node().x;
	}
};
/**
// How are the elements stored in the compressed row format sparse matrix
// that is passed to the solver. Petsc requires the elements to be in the
// same order as they occur in the matrix. AMG requires the diagonal to be
// first.
*/
class poreIndexCompare
{
public:
	bool operator() (std::pair<int,double> colOne, std::pair<int,double> colTwo) const
	{
		return (colOne.first < colTwo.first);
	}
};

/**
// The diagonal (positive) will always be larger than the off diagonal (that are
// negative)
*/
class poreDiagonalFirst
{
public:
	bool operator() (std::pair<int,double> colOne, std::pair<int,double> colTwo) const
	{
		return (colOne.second > colTwo.second);
	}
};

class throatIndexCompare
{
public:
	bool operator() (std::pair<const Elem*, double> thOne, std::pair<const Elem*, double> thTwo) const
	{
		return (thOne.first->indexOren() < thTwo.first->indexOren());
	}
};





#endif
