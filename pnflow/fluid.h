#ifndef FLUID_H
#define FLUID_H
#define FLUIDF_SkipH // used in XdmfI.h to check if fluidf is defined

#include <string>
#include "InputFile.h"

using namespace std;

//ff&(AP0|AP1|ATR)==(AP0|AP1|ATR)  -> duplicate vis
//OTHER=ff^(OIL|WTR)
typedef  unsigned short ffT;
enum fluidf  : ffT
{
	NONE     = 0,
	WTR      = 1, /// water
	OIL      = 1<<1, ///. 2 -> oil
	ELEC     = WTR|OIL,
	GAS      = 1<<2,
	ANF      = WTR|OIL|GAS,
	CLAY     = 1<<4,
	mj_ELEC     = 4,

// GNM_1:
	SOLID    = 1<<5,
	NEIF0    = 1<<6, //p0  neiL
	NEIF1    = 1<<7, //p1  neiL
	NEIF     = 1<<8, //both neiL

	WTRNEI  = NEIF|WTR,
	WTRNEI0 = NEIF0|WTRNEI,
	WTRNEI1 = NEIF1|WTRNEI,
	WTRNEI01 = WTRNEI0|WTRNEI1,

	OILNEI  = NEIF|OIL,
	OILNEI0 = NEIF0|OILNEI,
	OILNEI1 = NEIF1|OILNEI,
	OILNEI01 = OILNEI0|OILNEI1,


//! interf conf
	TRAP = 1<<9,
	AtPL = 1<<10, // poer layer connected to throat
	AtTC = 1<<11,
	//AtPT = AtPC|AtTC,
	AtLL = 1<<12, // Layer interface, "separating from centre"
	AtNP = 1<<13, // piston interface separating from pore centre
	//AtPL = AtLL|AtNP,
	//AtNT = 1<<14,
	//AtLT = AtLL|AtNT,
	//AtLB = AtLL|AtNP|AtNT,
//! half T index
	IHT = 1<<15,
	LEFTSIDE = 1<<15,	//!< closer to an onOutlet_ at the LeFT SIDe
// GNM_1;
};
inline fluidf  othrOW(fluidf ff)             {return fluidf(ff^(OIL|WTR));}
inline bool    isAll(fluidf ff, fluidf FFs)  {return (ff&FFs)==FFs;}
inline bool    hasAny(fluidf ff, fluidf FFs) {return (ff&FFs);}
inline ffT     iof(fluidf ff)                {return (ff>>1)&3;}//works for water oil and gas only

inline fluidf  operator | (fluidf f1, fluidf f2) {return fluidf(ffT(f1) | f2);}
inline fluidf  operator & (fluidf f1, fluidf f2) {return fluidf(ffT(f1) & f2);}

inline fluidf& operator |= (fluidf& f1, fluidf f2) {(ffT&)(f1)|=ffT(f2); return f1;}
inline fluidf& operator &= (fluidf& f1, fluidf f2) {(ffT&)(f1)&=ffT(f2); return f1;}
inline fluidf& operator ^= (fluidf& f1, fluidf f2) {(ffT&)(f1)^=ffT(f2); return f1;}

class Fluid
{
  public:

	Fluid(): ff_(NONE) {};
	Fluid(fluidf ff): ff_(ff) {};
	Fluid(fluidf ff, const InputFile& input, std::string nam, double scaleFactor, int importance)
	: ff_(ff), name(nam)
	{
		std::istringstream data;
		if(input.giv(name, data,importance))
		{
			data >> viscosity_ >> resistivity_ >> density_;

			input.checkEndOfData(data, name, "\nExpected:\n  "+nam+":  viscosity resistivity density",importance==2);
		 }else{
			resistivity_=input.getOr(name+"Resistivity",ff==OIL?1000.:1.);
			viscosity_=input.getOr(name+"Viscosity", 0.001);
			density_=input.getOr(name+"Density", 1000.);
		}
		viscosity_*=scaleFactor;
		density_*=scaleFactor;
	};

	void setFluidProps(double viscosity, double interfacialTen, double resistivity, double density, double epsilon, double zeta)
	{
		viscosity_ = viscosity;
		//interfacialTen_ = interfacialTen;
		resistivity_ = resistivity;
		density_ = density;
		epsilon_ = epsilon;//mj_
		zeta_ = zeta;//mj_
	}///. TODO enable auto initialization from input data


	double viscosity() const {return viscosity_;}
	double resistivity() const {return resistivity_;}
	double density() const {return density_;}
	double epsilon() const { return epsilon_; }//mj_
	double zeta() const { return zeta_; }//mj_

	bool isOil()  const { return ff_==OIL; };
	fluidf ff()  const { return ff_; };

	friend class GNMData; ///. scales viscosity and density

	const fluidf ff_;
	const std::string name;

//protected:

	double                  viscosity_;
	double                  resistivity_;
	double                  density_;
	double                  epsilon_;//mj_
	double                  zeta_; //mj_

};




class Interface
{
  public:

	Interface(std::string nam, const InputFile& input, int importance = 0)
	: name(nam), interfacialTen_(0.03)
	, epsilon_(0.), zeta_(1e37)//mj_
	{
		std::istringstream data;
		if(input.giv(name, data, importance))
		{
			data >> interfacialTen_ ;
			data >> epsilon_ >> zeta_;//mj_
		}
	};


	double sigma() const {return interfacialTen_;}
	double epsilon() const { return epsilon_; }//mj_
	double zeta() const { return zeta_; }//mj_

	const std::string 				name;

protected:

	double                  interfacialTen_;
	double                  epsilon_;//mj_
	double                  zeta_; //mj_

};


inline double dirInj(fluidf ff) {return 2.*(double(ff)-1.)-1.;}
inline char name(fluidf ff) {return ff==OIL ? 'O' : ff==WTR ? 'W' : 'E';}


#endif
