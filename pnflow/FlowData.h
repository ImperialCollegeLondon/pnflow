#ifndef FLOWDATA_H
#define FLOWDATA_H

#include <random>
#include <array>
#include <stack>

#include "InputFile.h"
#include "typses.h"
#include "fluid.h"


#define ISw 0
#define IPc 1



class RockType // shall not be constrcuted in parLoop necause of RockTypeCounter
{
 public:
	//RockType(int index)
	//{
		//index_ = index; //++RockTypeCounter;
		//maxRcRandnessF_ = 0.;
		//maxRcRandness_  = 0.; //between 0 and 1 (voxels)
		//RcRoughBias_    = 0.; //between -1 and 1  (voxels), positive values decrease the radius
	//}
	RockType(const InputFile& inp, int index): index_(index), volume_(0.)  {
		index_ = index; //++RockTypeCounter;
		maxRcRandnessF_ = inp.getOr("maxPcRandness", 0.01);
		maxRcRandness_ = inp.getOr("maxRcRandness", 0.50); //between 0 and 1 (voxels)
		RcRoughBias_ = inp.getOr("RcRoughBias", 0.00); //between -1 and 1  (voxels), positive values decrease the radius
	}
	virtual double Sw(double Pc, int icycle) const = 0;
	virtual double Pc(double Sw, int icycle) const = 0;
	virtual ~RockType() {};
	int index() const { return index_;}
//protected :
	int index_;
	double maxRcRandnessF_;
	double maxRcRandness_;
	double RcRoughBias_;
	mutable double volume_; // used by satTracker
	//static thread_local int     RockTypeCounter;

};

class SolidWall : public RockType {
 public:
	SolidWall(const InputFile& inp) : RockType(inp,0) {}
	double Sw(double Pc, int icycle) const { return 0.; } //used by SaturationTracker
	double Pc(double Sw, int icycle) const { return 1e29*(2*(icycle%2)-1.); }
};


//template <typename T>
struct WeibulParam {
	WeibulParam() : minV(0.), maxV(0.), deltaExp(-0.2), etaExp(-3.), correlation("rand") {}
	WeibulParam(std::istream& ins) : minV(0.), maxV(0.), deltaExp(-0.2), etaExp(-3.), correlation("rand")  { read(ins); }

	inline void read(std::istream& ins)  {
		ins>>minV;
		if(ins.good()) ins >>maxV;    else maxV=minV;

		if(ins.good()) ins >>deltaExp;
		if(ins.good()) ins >>etaExp;
		if(ins.good()) ins >>correlation;
		std::cout<<"WeibulParam: "<<minV<<"  "<<maxV<<"  "<<deltaExp<<"  "<<etaExp<<"  "<<correlation<<"  "<<std::endl;
	}

	double minV, maxV;
	double deltaExp, etaExp;
	std::string correlation;
};

struct Weibul1  {
	Weibul1() : minV(0.), delV(0.), lambDV(-0.6), kk(-3.), kInv(1/kk), cor("rand") {}
	Weibul1(std::istream& ins) : minV(0.), delV(0.), lambDV(-0.6), kk(-3.), kInv(1/kk), cor("rand") { read(ins); }

	void read(std::istream& ins)  {
		ins>>minV;
		if(ins.good()) { ins >>delV;  delV-=minV;  } // maxV

		if(ins.good()) ins >>lambDV; // == lambda here
		if(ins.good()) { ins >>kk; }  kInv=1./kk;
		cdf1=1.-exp(-pow(lambDV,-kk)); //! ==CDF(1.) = CDF(minV+x/delV), before scale by delV
		if(ins.good()) ins >>cor;
		std::cout<<"Weibul1: "<<minV<<"  "<<minV+delV<<"  "<<lambDV<<"  "<<kk<<"  "<<cor<<"  "<<std::endl;
		lambDV*=delV; // == lambda*delV here
	}
	void scale(double s)        { delV*=s;  lambDV*=s; }
	void scaleLambda(double s)  { lambDV*=s; cdf1=1-exp(-pow(lambDV/delV,-kk));} //
	double cdf(double xx) const { return (1.-exp(-pow((xx-minV)/lambDV,kk)))/cdf1; } // todo add unif
	double quantile(double pp) const {
		return (lambDV<0. || kInv<0.)  ? minV + delV*pp : // Uniform Distribution
		       minV+lambDV*pow(-log(1.-pp*cdf1), kInv);  /// Weibull truncated up to 1., shifted by minV, scaled by delV (lambDV=lambda*delV here)
	}

	double minV;
//private:
	double delV;
	double lambDV, kk, kInv, cdf1; //!<   cdf1: CDF(1)=1-exp(-1/lambda^k)
	std::string cor;//!<  correlation
};
inline std::istream&  operator >> (std::istream& in, Weibul1& wb) { wb.read(in); return in; }



class GNMData  {
 public:

	GNMData(const InputFile& inp, std::string prtFileName)//, const FlowDomain* domain
		: input_(inp), watoil_("WaterOilInterface",inp,1), //_watrok("WaterSolidInterface",inp,0), _oilrok("OilSolidInterface",inp,0),
		  water_(WTR,inp, "Water", 1./watoil_.sigma(),1), oil_(OIL,inp, "Oil", 1./watoil_.sigma(),1),
		  elec_(ELEC,inp, "Electricity", 1./watoil_.sigma(),0), clay_(CLAY,inp, "Clay", 1./watoil_.sigma(),0),
		  randomGenerator_(inp.getOr("RAND_SEED", (unsigned)time(nullptr))), uniformRand01_(0.,1.),
		  out_(prtFileName,inp.getOr("outPrtStd", 3))
	{
		rockTypes_.push_back(new SolidWall(inp));
		floodingCycle_ = 0;
		if(inp.informative) cout<<"ift:"<<watoil_.sigma()<<endl;
		if(inp.informative) cout<<"Water,  mu:"<<water_.viscosity()*watoil_.sigma()<<" rho:"<<water_.density()*watoil_.sigma()<<" R:"<<water_.resistivity()<<" "<<endl;
		if(inp.informative) cout<<"Oil,    mu:"<<  oil_.viscosity()*watoil_.sigma()<<" rho:"<<  oil_.density()*watoil_.sigma()<<" R:"<<  oil_.resistivity()<<" "<<endl;
		if(inp.informative) cout<<"Clay,   R:"<<  clay_.resistivity()<<" "<<endl;
	}

	int  dispCycle() const {return floodingCycle_;}
	void incrementFloodCycle() {++floodingCycle_;}
	bool isDrainage() const {return floodingCycle_%2==1;}

	template<class T> void lookupOpt(const string& ky, T &val) { if(input_.lookup(ky,val))   out_<<"   "+ky+" "<<val<<endl; }

	default_random_engine& randomGenerator() const  { return randomGenerator_; }
	double  rand01() const   { return uniformRand01_(randomGenerator_); }
	double weibullObsolete(double minv, double maxv, double delta, double eta) const {// delta = pow(beta,eta)
		if(delta < 0. || eta < 0.)    return minv + (maxv-minv)*rand01(); // Uniform Distribution
		else  return minv+(maxv-minv) * pow(-delta*log(1.-rand01()*(1.-exp(-1./delta))), 1./eta);  /// Weibull truncated up to 1., scaled between min and max
		  /// wiki:  eta -> k ,  delta -> lambda^k
	}

	double weibul1(const Weibul1& wb) const { return wb.quantile(rand01()); }

	const Fluid&  oil()           const { return oil_; }
	const Fluid&  water()         const { return water_; }
	const Fluid&  elec()          const { return elec_; }
	const Fluid&  clay()          const { return clay_; }
	const Fluid& fluid(fluidf ff) const { return (&water_)[ff-1]; }
	const Interface& watoil()     const { return watoil_; }
	double sigmaOW()              const { return watoil_.sigma(); }
	mstream& outstream()          const { return out_; }
	const InputFile&	input()     const { return input_; }
	const stvec<RockType*>& rockTypes() const { return rockTypes_; }
	stvec<RockType*> &	   rockTypesCh()     { return rockTypes_; }
	bool	informative()           const { return input_.informative; }

	void addKcSwKrsQs(const std::array<double,8>& KcSwKrsQs)  {
		int icy = dispCycle();
		while(len(KrQsss_)<=icy) KrQsss_.push_back(stvec<std::array<double,8>>());
		if(KrQsss_[icy].size()) {
			if(isDrainage())  {
				ensure(KrQsss_[icy].back()[2]+1e-8>=KcSwKrsQs[2]);
				ensure(KrQsss_[icy].back()[3]<=KcSwKrsQs[3]+1e-8);  }
			else  {
				ensure(KrQsss_[icy].back()[2]<=KcSwKrsQs[2]+1e-8);
				ensure(KrQsss_[icy].back()[3]+1e-8>=KcSwKrsQs[3]);  }
		}
		KrQsss_[icy].push_back(KcSwKrsQs);
	}

	void writeResultData(bool wantRelPerm, bool wantResIdx, const std::string& title) const;

	stvec<stvec<std::array<double,8>>>  KrQsss_;

 protected:
	const InputFile&                    input_;

	Interface                           watoil_; //, _watrok, _oilrok
	Fluid                               water_;
	Fluid                               oil_; ///. TODO make const
	Fluid                               elec_;
	Fluid                               clay_;

	int                                 floodingCycle_;

	vector<RockType*>                   rockTypes_;


	mutable default_random_engine              randomGenerator_;
	mutable uniform_real_distribution<double>  uniformRand01_;

	mutable mstream                    out_;
};

#endif //FLOWDATA_H
