#ifndef NETSTATISTICS_H
#define NETSTATISTICS_H
/*---------------------------------------------------------------------------*\
Developed by (2015): Ali Q Raeini  email: a.qaseminejad-raeini09@imperial.ac.uk
\*---------------------------------------------------------------------------*/


#include <vector>
#include "../Element.h"

void printCornerAngStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printCornerNumStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printShapeFactorStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printRadiusStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printAspectRatioStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printCoordinaNumStatistics( const std::vector<Element const *> *  rockLattices, int m_numPores);
void printDistanceMapStatistics( const std::vector<Element const *> &  rockLattices, int m_numPores);




//void writeOilGangliaDistribution(string fileName, const std::vector< std::vector<Element*> >	& trappedRegionsOil );
//void writeWatGangliaDistribution( const std::vector< std::vector< std::pair<Element*,FluidBlob> > > & trappedRegionsWat );


#endif

