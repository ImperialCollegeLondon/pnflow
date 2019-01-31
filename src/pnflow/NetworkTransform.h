#include <random>
#include "polygon.h"

double sqr(double a) {return a*a;}

template<typename T, T> class MemberGetter;
template <class T, typename R, R (T::*mf)() const>
class MemberGetter <R (T::*)() const, mf> { public:
  int operator()(const T* elm) const {return (elm->*mf)();}
};

template<typename T, T> class UnitScalar;
template <class T, typename R, R (T::*mf)() const>
class UnitScalar <R (T::*)() const, mf> { public:
  int operator()(const T* elm) const  {return 1.0;}
};

template<typename T, T> class ComparerInc;
template <class T, typename R, R (T::*mf)() const>
class ComparerInc <R (T::*)() const, mf> { public:
    bool operator() (const T* e1, const T* e2) const  {   return (e1->*mf)() < (e2->*mf)();    }
};

template<typename T, T> class ComparerDec;
template <class T, typename R, R (T::*mf)() const>
class ComparerDec <R (T::*)() const, mf> { public:
    bool operator() (const T* e1, const T* e2) const  {   return (e1->*mf)() > (e2->*mf)();    }
};





template<typename T, size_t N>
std::array<T,N>& operator *= (std::array<T,N>& vs, T sf)
{
	for(auto& v : vs) v*=sf;
	return vs;
}
template<typename T, size_t N>
T sum(std::array<T,N>& vs)
{
	T s=0;
	for(auto v : vs) s+=v;
	return s;
}


//template <typename T>
class WeibalParam
{
 public:
 
	WeibalParam() : minV(0.0), maxV(0.0), deltaExp(-0.2), etaExp(-3.0), correlation("rand") {}

	inline void read(std::istream& data)
	{
		data>>minV;
		if(data.good()) data >>maxV;    else maxV=minV;

		if(data.good()) data >>deltaExp;
		if(data.good()) data >>etaExp;
		if(data.good()) data >>correlation;
		cout<<"WeibalParam: "<<minV<<"  "<<maxV<<"  "<<deltaExp<<"  "<<etaExp<<"  "<<correlation<<"  "<<endl;
	}
 public:

	double minV, maxV;
	double deltaExp, etaExp;
	string correlation;
};

class NetworkTransform
{  

template<class ET, typename T, T> class Correlate;

	template<class ET, class ET2, typename R, R (ET2::*mfcor)() const> //typename MFScale, typename MFWeight, 
	class Correlate <ET, R (ET2::*)() const, mfcor> { public:
		 vector<double> randCorrelate(vector<ET*> & elms, string correlate,  double minv, double maxv, double deltaExp, double etaExp, NetworkTransform& randGen)
		{
			if(correlate[2] == 'a' || correlate[2] == 'A')      // rMax
			  sort(elms.begin(), elms.end(), ComparerInc<double(ET2::*)() const, mfcor>());
			else if(correlate[2] == 'i' || correlate[2] == 'I') // rMin
			  sort(elms.begin(), elms.end(), ComparerDec<double(ET2::*)() const, mfcor>());
			else  shuffle(elms.begin(), elms.end(), randGen.randomGenerator());            // Random

			vector<double> randomfield(elms.size());
			for(size_t i = 0; i < randomfield.size(); ++i)
			{
			  randomfield[i] = randGen.weibull(minv, maxv, deltaExp, etaExp);
			}
			sort(randomfield.begin(), randomfield.end(), greater<double>());

			return randomfield;
		}
		void operator()(vector<ET*> & elms, string order, NetworkTransform& randGen)
			{
				if(order[2] == 'a' || order[2] == 'A')      // rMax
				  sort(elms.begin(), elms.end(), ComparerInc<double(ET2::*)() const, mfcor>());
				else if(order[2] == 'i' || order[2] == 'I') // rMin
				  sort(elms.begin(), elms.end(), ComparerDec<double(ET2::*)() const, mfcor>());
				else  shuffle(elms.begin(), elms.end(), randGen.randomGenerator());            // Random
			}
	};
	template<typename T>
	vector<T> randfield(size_t n, WeibalParam& wb, NetworkTransform& randGen)
	{
		vector<T> rf(n);
		for(size_t i = 0; i < rf.size(); ++i)
		  rf[i] = randGen.weibull(wb.minV, wb.maxV, wb.deltaExp, wb.etaExp);
		sort(rf.begin(), rf.end(), greater<double>());
		return rf;
	}




	void  scale(Throat* trot, double Rtppsf[3], double Gtppsf[3])
	{
		trot->ChModel()->setRadius(Rtppsf[0]*trot->RRR());

		//trot->m_sagittalKc = ((trot->m_Rtpp[0]-trot->m_Rtpp[1])/trot->m_LhTroppt[0]+  (trot->m_Rtpp[0]-trot->m_Rtpp[2])/trot->m_LhTroppt[1])*PI/(2.0*trot->m_LhTroppt[2]);

	}

	void  transform(Throat* tshap, double Rsf, double Gsf=1.0)
	{
		double	Rtppsf[3]={1.0,1.0,1.0};		Rtppsf[0]=Rsf;
		double	Gtppsf[3]={1.0,1.0,1.0};		Gtppsf[0]=Gsf;
		scale(tshap, Rtppsf, Gtppsf);
	}

	void  fixRadius(Throat* tshap)
	{ 
		double	Rtppsf[3]={1.0,1.0,1.0};
		double	Gtppsf[3]={1.0,1.0,1.0};
		
		Rtppsf[0]=min(min(1.0, 0.999*tshap->connection(0)->RRR()/tshap->RRR() ) , 0.999*tshap->connection(1)->RRR()/tshap->RRR() );
		scale(tshap, Rtppsf, Gtppsf);
	}

	void  transform(Pore* por, double Rsf, double Gsf=1.0)
	{
		por->ChModel()->setRadius(Rsf*por->RRR());
		for(int i=0; i<por->connectionNum();++i)
		{
			Throat* tshap = dynamic_cast<Throat*>(por->connection(i));
			if(tshap)
			{
				int neiI = (tshap->connection(0) == por) ? 0 : 1;
				double	Rtppsf[3]={1.0,1.0,1.0};		Rtppsf[neiI+1]=0.5*Rsf+0.5;
				double	Gtppsf[3]={1.0,1.0,1.0};		Gtppsf[neiI+1]=0.5*Gsf+0.5;
				scale(tshap, Rtppsf, Gtppsf);
			}
		}
	}

 public:
	NetworkTransform(unsigned int randSeed, mstream& out)
	: m_randomGenerator(randSeed), m_uniformRand01(0.0,1.0), m_out(out) {}
	

	void getModif(	vector<Pore*>& pores2BAltered, vector<Throat*>& throats2BAltered, const InputFile & input, istringstream& data, const vector<Element*> & elemans, size_t nPors)
	{

		bool volBased(false), totalFraction(false), ffClust(true);
		double fraction(1.0);
		//int clustDiam1(0), clustDiam2;,          deltaClustD(-0.2), etaClustD(-3.0)
		string  spatialDistrib("rand");
		WeibalParam  clustDs; clustDs.correlation = "rand";


		if(data.good()) 
		{
			
			if(data.good()) data >> fraction;
			if(data.good())
			{	char vB('x');    data >>vB;  	volBased = (vB=='T'||vB=='t'||vB=='V'||vB=='v');	
				if (!(vB=='T'||vB=='t'||vB=='V'||vB=='v'||vB=='F'||vB=='f'||vB=='N'||vB=='n'))	{ m_out<<" Error Wrong choice for \"volume-/number- based\" fraction, expected T/t/V/v or F/f/N/n"<<endl; exit(-1); }
			}
			if(data.good()) 
			{ 	char tF('x');  data >> tF;  totalFraction = (tF=='T'||tF=='t'); 
			   if (!(tF=='T'||tF=='t'||tF=='F'||tF=='f'||tF=='O'||tF=='o'))	{ m_out<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}
			}
			
			if(data.good()) data >> spatialDistrib;
			if (spatialDistrib[1] == 'o' || spatialDistrib[1] == 'O')
			{
				clustDs.read(data);
				clustDs.maxV+=0.999999;
				if (data.fail()) m_out <<"Error: cant read cluster size "<<endl;
				input.Assert(clustDs.minV>1 && clustDs.maxV>=clustDs.minV, "Error: wrong cluster length");
				if (data.good()) 
				{	char oInW('T'); data >> oInW;		ffClust = (oInW == 'Y' || oInW == 'y' || oInW == 'T' || oInW == 't');
					if (!(oInW=='T'||oInW=='t'||oInW=='F'||oInW=='f'||oInW=='O'||oInW=='o'))	{ m_out<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}
				}
			}
			else clustDs.correlation = spatialDistrib;
		}


		vector<Pore*> poresToSeed;poresToSeed.reserve(nPors+2);
		vector<int> clustDiams;
		{
		 for(size_t i = 1; i < 1+nPors; ++i)
			if (/*elemans[i]->exists(Ff) &&*/ dynamic_cast<Pore*>(elemans[i]))  poresToSeed.push_back(dynamic_cast<Pore*>(elemans[i]));


			Correlate<Pore, double(Element::*)()const, &Element::RRR>()(poresToSeed, clustDs.correlation, *this);
			clustDiams = randfield<int>(poresToSeed.size(),clustDs, *this);

		}

		pores2BAltered.reserve(poresToSeed.size());
		throats2BAltered.reserve(poresToSeed.size()*3);


		double ffInvadedVol(0.0), TotalVol(0.0), slctedVol(0.0);
		int totElem(0);

		for(size_t i = 1; i < elemans.size(); ++i)
		{
			//if (dynamic_cast<VoidElem*>(elemans[i]))
			if (i!=size_t(nPors+1))
			{
				//if(elemans[i]->exists(Ff))
				{
					if(i < size_t(nPors+2)) ++totElem;                        // Only count pores
					ffInvadedVol += elemans[i]->flowVolume();
				}
				TotalVol += elemans[i]->flowVolume();
			}
		}

		double targetFrac = fraction;
		if (totalFraction)  targetFrac = min(fraction*TotalVol/ffInvadedVol, 1.0);
		if (ffInvadedVol/TotalVol>1.0e-12)
		{

			vector<int> clusterIndx(elemans.size()+2, 0);

			int numFracWetted(0), clusterIdx(1);
			if(spatialDistrib[1] == 'o' || spatialDistrib[1] == 'O')    // Spatial correlation approach
			{
				bool clct_more(true);


				targetFrac = ffClust ? targetFrac: 1.0-targetFrac;

				while(clct_more)
				{
					set< Element * > frontier;
					Pore* elem = poresToSeed.back(); poresToSeed.pop_back();
					int   clustDiam = clustDiams.back(); clustDiams.pop_back();
						//int clustDiam = (*this).weibull(clustDiam1,clustDiam2, deltaClustD, etaClustD) +0.5;

					if (clusterIndx[elem->latticeIndex()]==0)
					{

						frontier.insert(elem);
						typedef set<Element*>::iterator ItrSet;
						int frontExpansions(0);
						while(!frontier.empty() && clct_more && frontExpansions < clustDiam)  // Front expansions will go over thraots as well as pores, wheras
						{                                                                           // cluster diam is wrt pores  => Don't divide by 2
							set<Element*> oldFrontier = frontier;
							frontier.clear();
							++frontExpansions;
							for(ItrSet itrElem = oldFrontier.begin(); itrElem != oldFrontier.end(); ++itrElem)
							{

								if( clusterIndx[(*itrElem)->latticeIndex()] != clusterIdx && clct_more)
								{
									if(clusterIndx[(*itrElem)->latticeIndex()] == 0)// && (*itrElem)->exists(Ff))
									{
										clusterIndx[(*itrElem)->latticeIndex()] = clusterIdx;                    /// Set flag
										slctedVol += (*itrElem)->flowVolume();
										if(dynamic_cast<Pore*>(*itrElem)) ++numFracWetted;   /// Only count pores

										clct_more = volBased ? slctedVol/ffInvadedVol < targetFrac :  double(numFracWetted)/totElem < targetFrac;

										if(!clct_more)			break;
									}

									for(int i = 0; i < (*itrElem)->connectionNum(); ++i)
										if ( clusterIndx[(*itrElem)->connection(i)->latticeIndex()] != clusterIdx
											  && (dynamic_cast<Throat*>((*itrElem)->connection(i)) || dynamic_cast<Pore*>((*itrElem)->connection(i))) )
											frontier.insert((*itrElem)->connection(i));
								}
							}
						}
						++clusterIdx;
					}
				}
				numFracWetted = 0;
				slctedVol = 0.0;
				if(ffClust)
				{
					for(size_t i = 0; i < elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->latticeIndex()] > 0)// && elemans[i]->exists(Ff))
						{
							if( i < nPors+2 && dynamic_cast<Pore*>(elemans[i]))
							{
								++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
								pores2BAltered.push_back(dynamic_cast<Pore*>(elemans[i]));
							}
							else if( dynamic_cast<Throat*>(elemans[i]))
							{
								++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
								throats2BAltered.push_back(dynamic_cast<Throat*>(elemans[i]));
							}
						}
					}
				}
				else
				{
					for(size_t i = 0; i < elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->latticeIndex()] == 0)// && elemans[i]->exists(Ff))
						{
							if( i < nPors+2 && dynamic_cast<Pore*>(elemans[i]))
							{
								++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
								pores2BAltered.push_back(dynamic_cast<Pore*>(elemans[i]));
							}
							else if( dynamic_cast<Throat*>(elemans[i]))
							{
								++numFracWetted;  slctedVol+=elemans[i]->flowVolume();
								throats2BAltered.push_back(dynamic_cast<Throat*>(elemans[i]));
							}
						}
					}
				}
				m_out<< "Number of correlated regions: " << clusterIdx-1                         << endl;

			}
			else ///. spatialDistrib == rand
			{
				//vector< Pore*> > toBeAltered;

				//for(int i = 2; i < nPors+2; ++i)
				//{
					//if(elemans[i]->exists(Ff) && dynamic_cast<Pore*>(elemans[i]))
					//{
						//pair<double, Pore*> rockEntry;
						//rockEntry.second = dynamic_cast<Pore*>(elemans[i]);
						//if(spatialDistrib[1] == 'a' || spatialDistrib[1] == 'A')        ///. rand // Random
						//{
							//rockEntry.first = (*this).rand01();
						//}
						//else if(spatialDistrib[0] == 'r' || spatialDistrib[0] == 'R')   ///. rMax/rMin // Radius
							//rockEntry.first = rockEntry.second->RRR();
						//else
						//{
							//m_out << endl << endl << "Error: Did not recognize fractional wetting model: " << spatialDistrib  << endl  << endl << endl;
							//exit(-1);
						//}
						//toBeAltered.push_back(rockEntry);
					//}
				//}
				//if(spatialDistrib[2] == 'a' || spatialDistrib[2] == 'A') ///. rM'a'x
					//sort(toBeAltered.begin(), toBeAltered.end(), FracWettInc());
				//else
					//sort(toBeAltered.begin(), toBeAltered.end(), FracWettDec()); ///. rMin/rand

				int targetNum(fraction*poresToSeed.size());

				while(!poresToSeed.empty() &&
					((volBased && slctedVol/ffInvadedVol < targetFrac) || (!volBased && numFracWetted < targetNum)))
				{
					++numFracWetted;
					Pore* elem = poresToSeed.back(); poresToSeed.pop_back();
					slctedVol += elem->flowVolume();
					pores2BAltered.push_back(elem);
					int pind = elem->latticeIndex();
					clusterIndx[pind]=1;
					for(int t = 0; t < elem->connectionNum(); ++t)
					{
						Throat* throat = dynamic_cast<Throat*>(elem->connection(t));
						if( throat && clusterIndx[throat->latticeIndex()] == 0) // throat->exists(Ff) &&
						{
							int tind = throat->latticeIndex();
							double randNum = (*this).rand01();
							if(randNum < fraction)
							{
								throats2BAltered.push_back(throat);
								clusterIndx[tind]=1;
								slctedVol += throat->flowVolume();
							}
						}
					}
				}
			}
			m_out
				<< "selected volume (fraction): " << slctedVol/ffInvadedVol   << endl
				<< "selected volume (of total net pore volume): " << slctedVol/TotalVol   << endl
				<< "Number of pores selected: " << numFracWetted                            << endl
				<< "================================================================="  << endl;

		}
		else m_out
			<< "\nWarning:  nothing selected, too low volume fraction: " << ffInvadedVol/TotalVol             << endl
			<< "=================================================================\n"      << endl;




		//setContactAngles(pores2BAltered, throats2BAltered, minVal, maxVal, deltaExp, etaExp, cntctAngModel, modelTwoSepAng,  globalCorrelation, nPors, (*this));


	}


	
	void modify(const InputFile & input, vector<Element*> & elemans, size_t nPors, double boxVolume)//, dbl3 boxSize
	{
		///. modify network
		//IncreaseClay     0.3   0.3  -0.2  -3.0   rand  ;
		//ConvertToClay    0.3   0.3  -0.2  -3.0   rand  T  ;
		//ScaleRadius     0.9   1.1  -0.2  -3.0   rand;
		//ScalePoreRadius     0.9   1.1  -0.2  -3.0   rand;
		//ScaleThroatRadius     0.9   1.1  -0.2  -3.0   rand;
		//MatchMICP     MICPCurve.txt  0.485  45;

		istringstream data;
		if (input.getData(data,"AddClay"))
		{  m_out<<"Adding clay porosity: "<<endl;

			WeibalParam  wbdist;
			wbdist.read(data);
			input.Assert(!data.fail(), "AddClay", "wrong data",true);

			vector<Throat*> throats;throats.reserve(elemans.size()-nPors+2);
			for(size_t i = 2+nPors; i < elemans.size(); ++i)
				if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));

			vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
			getModif( pores2BAltered,  throats2BAltered,  input,  data, elemans, nPors);

			//vector<double> claypores = randfield<double>(pores2BAltered.size(),minVal, maxVal, deltaExp, etaExp, *this);
			//Correlate<Pore, double(Element::*)()const, &Element::RRR>()(pores2BAltered, wbdist.correlation, *this);
			
			vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
			Correlate<Throat, double(Element::*)()const, &Element::RRR>()(throats2BAltered, wbdist.correlation, *this);

			for(size_t i = 0; i < throats2BAltered.size(); ++i)
			{
				//Throat* trot = throats2BAltered[i];
				//double clayFrac = rands[i];
				//trot->setClayVolume(clayFrac*trot->flowVolume()/(1.0-clayFrac));
			}
		}

		if (input.getData(data,"FillWithClay"))
		{
			m_out<<"Filling void space by clay : "<<endl;
			WeibalParam  wbdist;
			wbdist.read(data);
			input.Assert(!data.fail(), "FillWithClay", "wrong data",true);

			vector<Throat*> throats;throats.reserve(elemans.size()-nPors+2);
			for(size_t i = 2+nPors; i < elemans.size(); ++i)
				if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));

			vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
			getModif( pores2BAltered,  throats2BAltered,  input,  data, elemans, nPors);

			vector<double> clayporesRad = randfield<double>(pores2BAltered.size(),wbdist, *this);
			Correlate<Pore, double(Element::*)()const, &Element::RRR>()(pores2BAltered, wbdist.correlation, *this);
			
			vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
			Correlate<Throat, double(Element::*)()const, &Element::RRR>()(throats2BAltered, wbdist.correlation, *this);

			double totalFlowVolume(0.0), totalClayVolSum(0.0);
			for(size_t i = 0; i < elemans.size(); ++i)
			{
				totalFlowVolume -= elemans[i]->flowVolume();
				totalClayVolSum -= elemans[i]->clayVolume();
			}
			//double selectedFlowVolume(0.0), selectedClayVolSum(0.0);
			//for(size_t i = 0; i < throats2BAltered.size(); ++i)
			//{
				//double totalFlowVolume = throats2BAltered[i]->flowVolume();
				//double totalClayVolSum = throats2BAltered[i]->clayVolume();
			//}

			for(size_t i = 0; i < throats2BAltered.size(); ++i)
			{
				Throat* trot = throats2BAltered[i];
				double netVolSum, clayVolSum;
				double clayFrac = rands[i];
				double newNetVol = trot->flowVolume()*(1.0-clayFrac);
				double newClayVol = trot->flowVolume()*(clayFrac) +trot->clayVolume();
				trot->adjustVolume(newNetVol, newClayVol, netVolSum, clayVolSum);
				transform(trot, sqrt(1.0-clayFrac));
				//for (int i=0;i<2;++i)
				//{
					//Pore* por = dynamic_cast<Pore*>(trot->connection(i));
					//double newNetVol = por->flowVolume()*(1.0-clayFrac);
					//double newClayVol = por->flowVolume()*(clayFrac) +por->clayVolume();
					//por->adjustVolume(newNetVol, newClayVol, netVolSum, clayVolSum);
					//transform(por, 0.5+0.5*sqrt(1.0-clayFrac));
				//}
			}

			for(size_t i = 0; i < pores2BAltered.size(); ++i)
			{
				Pore* por = pores2BAltered[i];
				double netVolSum, clayVolSum;
				double clayFrac = clayporesRad[i];
				double newNetVol = por->flowVolume()*(1.0-clayFrac);
				double newClayVol = por->flowVolume()*(clayFrac) +por->clayVolume();
				por->adjustVolume(newNetVol, newClayVol, netVolSum, clayVolSum);
				transform(por, sqrt(1.0-clayFrac));
			}

			for(size_t i = 2+nPors; i < elemans.size(); ++i)
				if (dynamic_cast<Throat*>(elemans[i]))	  fixRadius(dynamic_cast<Throat*>(elemans[i]));

				//exit(-1);
		}

		double radiusScale;
		if (input.getType(radiusScale,"scaleRadius"))
		{
			m_out<<" scalling radius "<<endl;
			WeibalParam  wbdist;
			wbdist.read(data);
			input.Assert(!data.fail(), "scaleRadius", "wrong data",true);

			vector<Throat*> throats;throats.reserve(elemans.size()-nPors+2);
			for(size_t i = 2+nPors; i < elemans.size(); ++i)
				if (dynamic_cast<Throat*>(elemans[i]))	  throats.push_back(dynamic_cast<Throat*>(elemans[i]));

			vector<Pore*> pores2BAltered; vector<Throat*> throats2BAltered;
			getModif( pores2BAltered,  throats2BAltered,  input,  data, elemans, nPors);

			//vector<double> claypores = randfield<double>(pores2BAltered.size(),minVal, maxVal, deltaExp, etaExp, *this);
			//Correlate<Pore, double(Element::*)()const, &Element::RRR>()(pores2BAltered, wbdist.correlation, *this);

			vector<double> rands = randfield<double>(throats2BAltered.size(),wbdist, *this);
			Correlate<Throat, double(Element::*)()const, &Element::RRR>()(throats2BAltered, wbdist.correlation, *this);

			for(size_t i = 0; i < throats2BAltered.size(); ++i)
			{
				Throat* trot = throats2BAltered[i];
				double scaleFact = rands[i];
				transform(trot, scaleFact);
			}
		}
		if (input.getType(radiusScale,"scaleRadiusMICP"))
		{
			m_out<<"scaleMICP: ToBeImplemented"<<endl;
		}
	}
 private:

	std::default_random_engine&  randomGenerator() const  { return m_randomGenerator; }
	double  rand01() const   { return m_uniformRand01(m_randomGenerator); }
	inline double weibull(double minv, double maxv, double deltaExp, double etaExp) const // deltaExp = pow(beta,etaExp)
	{	if(deltaExp < 0.0 || etaExp < 0.0)    return minv + (maxv-minv)*rand01(); // Uniform Distribution
		else  return (maxv-minv) * pow(-deltaExp*log(1.0-rand01()*(1.0-exp(-1.0/deltaExp))), 1.0/etaExp) + minv; // Weibull truncated up to 1.0, scaled between min and max  	// return (maxv-minv) * pow(-deltaExp*log(    rand01()*(1.0-exp(-1.0/deltaExp))+exp(-1.0/deltaExp)), 1.0/etaExp) + minv;   
	}
	mutable std::default_random_engine		         m_randomGenerator;
	mutable std::uniform_real_distribution<double>  m_uniformRand01;
	mstream& m_out;
};














