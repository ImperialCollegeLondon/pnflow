#include <random>
#include "polygon.h"


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



class NetworkTransform
{  

template<class ET, typename T, T> class Correlate;

	template<class ET, class ET2, typename R, R (ET2::*mfcor)() const> //typename MFScale, typename MFWeight, 
	class Correlate <ET, R (ET2::*)() const, mfcor>
	{ public:
		 vector<double> randCorrelate(vector<ET*> & elms, string correlate,  double minv, double maxv, double deltaExp, double etaExp, NetworkTransform& randGen)
		{
			if(correlate[2] == 'a' || correlate[2] == 'A')      // rMax
			  sort(elms.begin(), elms.end(), ComparerInc<double(ET2::*)() const, mfcor>());
			else if(correlate[2] == 'i' || correlate[2] == 'I') // rMin
			  sort(elms.begin(), elms.end(), ComparerDec<double(ET2::*)() const, mfcor>());
			else  shuffle(elms.begin(), elms.end(), randGen.randomGenerator());            // Random

			vector<double> rnds(elms.size());
			for(size_t i = 0; i < rnds.size(); ++i)   rnds[i] = randGen.weibull(minv, maxv, deltaExp, etaExp);
			sort(rnds.begin(), rnds.end(), greater<double>());
			return rnds;
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
	vector<T> randfield(size_t n, WeibulParam& wb, NetworkTransform& randGen)
	{
		vector<T> rf(n);
		for(auto& rfi:rf) rfi = randGen.weibull(wb.minV, wb.maxV, wb.deltaExp, wb.etaExp);
		sort(rf.begin(), rf.end(), greater<double>());
		return rf;
	}




	void  scale(Throat* trot, double Rtppsf[3], double Gtppsf[3])
	{
		trot->ChModel()->setRadius(Rtppsf[0]*trot->RRR());

		//trot->sagittalKc_ = ((trot->Rtpp_[0]-trot->Rtpp_[1])/trot->LhTroppt_[0]+  (trot->Rtpp_[0]-trot->Rtpp_[2])/trot->LhTroppt_[1])*PI/(2.0*trot->LhTroppt_[2]);

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
	: randomGenerator_(randSeed), uniformRand01_(0.0,1.0), out_(out) {}


	void getModif(	vector<Pore*>& pores2BAltered, vector<Throat*>& throats2BAltered, const InputFile & input, istringstream& ins, const vector<Element*> & elemans, size_t nBpPors);


	void modify(const InputFile & input, vector<Element*> & elemans, size_t nBpPors, double boxVolume);
 private:

	std::default_random_engine&  randomGenerator() const  { return randomGenerator_; }
	double  rand01() const   { return uniformRand01_(randomGenerator_); }
	double weibull(double minv, double maxv, double deltaExp, double etaExp) const // deltaExp = pow(beta,etaExp)
	{	if(deltaExp < 0.0 || etaExp < 0.0)    return minv + (maxv-minv)*rand01(); // Uniform Distribution
		else  return (maxv-minv) * pow(-deltaExp*log(1.0-rand01()*(1.0-exp(-1.0/deltaExp))), 1.0/etaExp) + minv; // Weibull truncated up to 1.0, scaled between min and max  	// return (maxv-minv) * pow(-deltaExp*log(    rand01()*(1.0-exp(-1.0/deltaExp))+exp(-1.0/deltaExp)), 1.0/etaExp) + minv;   
	}
	mutable std::default_random_engine				 randomGenerator_;
	mutable std::uniform_real_distribution<double>  uniformRand01_;
	mstream& out_;
};



void setShapeFactFromFile(vector<Element*>& elems, const string& fileName, double lowCO, double hiCO);

void modifyInscribedRadii(int model, const string& options, vector<Element*>& elems, const string& fileName, bool writeToFile, int numPts, bool forPores = false);
void modifyShapeFactor(int model, const string& options, vector<Element*>& elems, const string& fileName, bool writeToFile, int numPts);
void setRadiiFromFile(vector<Element*>& elems, const string& fileName, double lowCO, double hiCO, bool forPores);


inline double weibull(double minv, double maxv, double delta, double eta) 
{
	double randm = double(rand())*(1.0/double(RAND_MAX));
	if(delta < 0.0 && eta < 0.0) return minv + (maxv-minv)*randm;                   // Uniform Distribution
	else    return (maxv-minv) * pow(-delta*log(randm*(1.0-exp(-1.0/delta))+exp(-1.0/delta)), 1.0/eta) + minv;
}





