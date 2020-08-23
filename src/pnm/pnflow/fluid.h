#ifndef FLUID_H
#define FLUID_H
#define FLUIDF_SkipH // used in XdmfI.h to check if fluidf is defined

#include "globals.h" // defines _begins_
#include "InputFile.h"

//keep compatible with gnflow
typedef  unsigned short ffT;
enum fluidf  : ffT
{
	NONE     = 0,
	WTR      = 1, /// water
	OIL      = 1<<1, ///. 2 -> oil
	ELEC     = WTR|OIL,
	CLAY     = 4
};




class Fluid
{
public:

	//Fluid(double viscosity, double interfacialTen, double resistivity, double density) : viscosity_(viscosity),
		//interfacialTen_(interfacialTen), resistivity_(resistivity), density_(density) {}
	//Fluid(const Fluid& fluid) : viscosity_(fluid.viscosity_), interfacialTen_(fluid.interfacialTen_),
		//resistivity_(fluid.resistivity_) {}

	Fluid(fluidf ff): ff_(ff) {};
	Fluid(): ff_(NONE) {};
	Fluid(fluidf ff, const InputFile& input, std::string nam, int importance=0)
		: ff_(ff), name(nam)
	{
		std::istringstream data;
		if(input.getData(data, name,importance))
		{
			data >> viscosity_ >> resistivity_ >> density_;

			input.checkEndOfData(data, name, "\nExpected:\n  "+nam+":  viscosity resistivity density",importance==2);
		 }else{
			resistivity_=input.getOr(ff==WTR?1.0:1000.0,name+"Resistivity");
			viscosity_=input.getOr(0.001,name+"Viscosity");
			density_=input.getOr(1000.0,name+"Density");
		}
		//viscosity_*=scaleFactor;
		//density_*=scaleFactor;
	};

	void setFluidProps(double viscosity, double interfacialTen, double resistivity, double density, double epsilon, double zeta)
	{
		viscosity_ = viscosity;
		//interfacialTen_ = interfacialTen; 
		resistivity_ = resistivity; 
		density_ = density;
	}///. TODO enable auto initialization from input data


	double viscosity() const {return viscosity_;}
	//double interfacialTen() const {return interfacialTen_;}
	double resistivity() const {return resistivity_;}
	double density() const { return density_; }

	bool isOil()  const { return ff_==OIL; };
	fluidf ff()  const { return ff_; };

	friend class CommonData; ///. scales viscosity and density

	fluidf                  ff_;
	const std::string name;

//protected:

	double                  viscosity_;
	//double                  interfacialTen_;
	double                  resistivity_;
	double                  density_;

};




class Interface
{
  public:

	Interface(std::string nam, const InputFile& input, int importance = 0)
	: name(nam), interfacialTen_(0.03)
	{
		std::istringstream data;
		if(input.getData(data, name, importance)) 
		{
			data >> interfacialTen_ ;
		}
	};


	double sigma() const {return interfacialTen_;}

	const std::string 				name;

protected:

	double                  interfacialTen_;

};

#endif

