#ifndef FLUID_H
#define FLUID_H

class Fluid
{
public:

    //Fluid(double viscosity, double interfacialTen, double resistivity, double density) : m_viscosity(viscosity),
        //m_interfacialTen(interfacialTen), m_resistivity(resistivity), m_density(density) {}
//
    //Fluid(const Fluid& fluid) : m_viscosity(fluid.m_viscosity), m_interfacialTen(fluid.m_interfacialTen),
        //m_resistivity(fluid.m_resistivity) {}
        
    Fluid(){};
    void setFluidProps(double viscosity, double interfacialTen, double resistivity, double density, double epsilon, double zeta)
    {
		m_viscosity = viscosity;
        m_interfacialTen = interfacialTen; 
        m_resistivity = resistivity; 
		m_density = density;
		m_epsilon = epsilon;
		m_zeta = zeta;
	}///. TODO enable auto initialization from input data
	
    virtual ~Fluid() {}

    double viscosity() const {return m_viscosity;}
    double interfacialTen() const {return m_interfacialTen;}
    double resistivity() const {return m_resistivity;}
	double density() const { return m_density; }
	double epsilon() const { return m_epsilon; }
	double zeta() const { return m_zeta; }

    virtual bool isOil()  const = 0;

protected:

    double                  m_viscosity;
    double                  m_interfacialTen;
    double                  m_resistivity;
	double                  m_density;
	double                  m_epsilon;
	double                  m_zeta;

};

class Oil : public Fluid
{
public:
	
	bool isOil() const {return true;};
    static int ind()  {return 1;};

    //Oil(double viscosity, double interfacialTen, double resistivity, double density) : Fluid(viscosity,
        //interfacialTen, resistivity, density) {}
    //Oil(const Oil& oil) : Fluid(oil) {}

private:

};

class Water : public Fluid
{
public:
	bool isOil() const { return false; };
	static int ind() { return 0; };

	//Water(double viscosity, double interfacialTen, double resistivity, double density) : Fluid(viscosity,
	//interfacialTen, resistivity, density) {}
	//Water(const Water& water) : Fluid(water) {}

private:

};

class Clay : public Fluid
{
public:
	bool isOil() const { return false; };
	static int ind() { return 3; };

	//Water(double viscosity, double interfacialTen, double resistivity, double density) : Fluid(viscosity,
	//interfacialTen, resistivity, density) {}
	//Water(const Water& water) : Fluid(water) {}

private:

};


class Surface : public Fluid
{
public:
	bool isOil() const { return false; };
	static int ind() { return 4; };

	//Water(double viscosity, double interfacialTen, double resistivity, double density) : Fluid(viscosity,
	//interfacialTen, resistivity, density) {}
	//Water(const Water& water) : Fluid(water) {}

private:

};

#endif

