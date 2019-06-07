#ifndef COMPARE_FUNC_H
#define COMPARE_FUNC_H

#ifdef WIN32
#pragma warning(disable:4786)
#endif




#include "Element.h"

class ElemRadCmpInc
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->model()->radius() < elem2->model()->radius());
    }
};

class ElemRadCmpRed
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->model()->radius() > elem2->model()->radius());
    }
};

//class ElemEquivRadCmpRed
//{
//public:
    //bool operator() (const Element* elem1, const Element* elem2) const
    //{
        //return (elem1->model()->equivalentRadius() > elem2->model()->equivalentRadius());
    //}
//};
 
class ElemVolCmpRed
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->flowVolume() > elem2->flowVolume());
    }
};

class ElemGCmpRed
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->model()->shapeFactor() > elem2->model()->shapeFactor());
    }
};

class ElemRandCmp
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->randomAssInt() < elem2->randomAssInt());
    }
};

class lenToRadRatioCmp
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        return (elem1->lenToRadRatio() > elem2->lenToRadRatio());
    }
};


class PceImbCmp
{
public:
    bool operator() (const Apex* elem1, const Apex* elem2) const
    {
        register double pc1(elem1->gravCorrectedEntryPress());
        register double pc2(elem2->gravCorrectedEntryPress());
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
        register double pc1(elem1->gravCorrectedEntryPress());
        register double pc2(elem2->gravCorrectedEntryPress());
        if(pc1 == pc2)
            return elem1 < elem2;
        else
            return pc1 > pc2;
    }
};

/*


class CoalesceWatFillCmp  ///. PceImbCmp
{
public:
    bool operator() (const Apex* elem1, const Apex* elem2) const
    {
        // if(elem1.third() == elem2.third())
            // return elem1.second() < elem2.second();
        // else
            // return elem1.third() < elem2.third();
            if (elem1->gravCorrectedEntryPress() == elem2->gravCorrectedEntryPress()) return elem1->subIndex() < elem2->subIndex();
            else return elem1->gravCorrectedEntryPress() < elem2->gravCorrectedEntryPress();
    }
};

class CoalesceOilFillCmp ///. PceDrainCmp
{
public:
    bool operator() (const Apex* elem1, const Apex* elem2) const
    {
        // if(elem1.third() == elem2.third())
            // return elem1.second() > elem2.second();
        // else
            // return elem1.third() > elem2.third();
            if (elem1->gravCorrectedEntryPress() == elem2->gravCorrectedEntryPress()) return elem1->subIndex() > elem2->subIndex();
            else return elem1->gravCorrectedEntryPress() > elem2->gravCorrectedEntryPress();
    }
};

 class PcColCmp ///. PceImbCmp
{
public:
    bool operator() (const  Apex* elem1, const Apex* elem2) const
    {
        register double pc1(elem1->gravCorrectedEntryPress());
        register double pc2(elem2->gravCorrectedEntryPress());
        if(pc1 == pc2)
            return elem1->subIndex() < elem2->subIndex();
        else
            return pc1 < pc2;
    }
};

class PcRefCmp ///. PceDrainCmp
{
public:
    bool operator() (const Apex* elem1, const Apex* elem2) const
    {
        register double pc1(elem1->gravCorrectedEntryPress());
        register double pc2(elem2->gravCorrectedEntryPress());
        if(pc1 == pc2)
            return elem1->subIndex() > elem2->subIndex();
        else
            return pc1 > pc2;
    }
};


*/







class FracWettInc
{
public:
    bool operator() (pair<double, Element*> prop1, pair<double, Element*> prop2) const
    {
        return (prop1.first < prop2.first);
    }
};

class FracWettDec
{
public:
    bool operator() (pair<double, Element*> prop1, pair<double, Element*> prop2) const
    {
        return (prop1.first > prop2.first);
    }
};

class TrappingWatStorageCmp
{
public:
    bool operator() (pair<Element*, FluidBlob> elm1, pair<Element*, FluidBlob> elm2) const
    {
        return (elm1.first->latticeIndex() > elm2.first->latticeIndex());
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
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        const Pore* pOneOne = dynamic_cast< const Pore* >(elem1->connection(0));
        const Pore* pOneTwo = dynamic_cast< const Pore* >(elem1->connection(1));
        const Pore* pTwoOne = dynamic_cast< const Pore* >(elem2->connection(0));
        const Pore* pTwoTwo = dynamic_cast< const Pore* >(elem2->connection(1));

        if(pOneOne == pTwoOne)
            return pOneTwo->node()->distToExit() < pTwoTwo->node()->distToExit();
        else if(pOneOne == pTwoTwo)
            return pOneTwo->node()->distToExit() < pTwoOne->node()->distToExit();
        else if(pOneTwo == pTwoOne)
            return pOneOne->node()->distToExit() < pTwoTwo->node()->distToExit();
        else
            return pOneOne->node()->distToExit() < pTwoOne->node()->distToExit();
    }
};

/**
// Comparing pores we directly obtain a position from their node
*/
class DistToExitComparePores
{
public:
    bool operator() (const Element* elem1, const Element* elem2) const
    {
        const Pore* pOne = dynamic_cast< const Pore* >(elem1);
        const Pore* pTwo = dynamic_cast< const Pore* >(elem2);

        return pOne->node()->distToExit() < pTwo->node()->distToExit();
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
    bool operator() (pair< int, double > colOne, pair< int, double > colTwo) const
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
    bool operator() (pair< int, double > colOne, pair< int, double > colTwo) const
    {
        return (colOne.second > colTwo.second);
    }
};

class throatIndexCompare
{
public:
    bool operator() (pair<const Element*, double> thOne, pair<const Element*, double> thTwo) const
    {
        return (thOne.first->orenIndex() < thTwo.first->orenIndex());
    }
};





#endif
