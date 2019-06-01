#ifndef TYPSES_H
#define TYPSES_H
/*---------------------------------------------------------------------------*\
Developed by (2017): Ali Q Raeini  email: a.qaseminejad-raeini09@imperial.ac.uk
\*---------------------------------------------------------------------------*/


#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <map>

#ifndef verySmall
	#define verySmall  1.0e-31
#endif


//typedef  int Int;
typedef  std::array<int,3> int3;
typedef  std::array<int3,3> int3x3;


class dbl3
{
public:

	double	x;
	double	y;
	double	z;

	dbl3() {}
	dbl3(double r, double s, double t)  { x = r;  y = s;  z = t; }
	dbl3(int3 n)  { x = n[0];  y = n[1];  z = n[2]; }
	dbl3(const double* d)  { x = d[0];  y = d[1];  z = d[2]; }

	double&       operator [](long k) {	return ((&x)[k]); }
	const double& operator [](long k) const  { return ((&x)[k]); }

	dbl3& operator +=(const dbl3& v)  { x += v.x;  y += v.y;  z += v.z;  return (*this); }
	dbl3& operator -=(const dbl3& v)  { x -= v.x;  y -= v.y;  z -= v.z;  return (*this); }
	dbl3& operator *=(double t)  { x *= t;  y *= t;  z *= t;  return (*this); }
	dbl3& operator /=(double t)  {  double f = 1.0 / t;  x *= f;  y *= f;  z *= f;  return (*this); }
	dbl3& operator ^=(const dbl3& v)  { double r, s;  r=y*v.z-z*v.y;  s=z*v.x-x*v.z;  z=x*v.y-y*v.x;  x=r; y=s; 	return (*this); }
	dbl3& operator *=(const dbl3& v)  { x *= v.x;  y *= v.y;  z *= v.z;  return (*this); }
	dbl3  operator -(void) const  { return (dbl3(-x, -y, -z)); }
	dbl3  operator +(const dbl3& v) const  { return (dbl3(x+v.x, y+v.y, z+v.z)); }
	dbl3  operator -(const dbl3& v) const  { return (dbl3(x-v.x, y-v.y, z-v.z)); }
	dbl3  operator *(double t) const  { return (dbl3(x*t, y*t, z*t)); }
	dbl3  operator /(double t) const  { double f = 1.0 / t;  return (dbl3(x*f, y*f, z*f)); }
	double operator &(const dbl3& v) const  { return (x*v.x+y*v.y+z*v.z); }
	dbl3  operator ^(const dbl3& v) const  { return (dbl3(y*v.z-z*v.y,  z*v.x-x*v.z,  x*v.y-y*v.x)); }
	dbl3  operator *(const dbl3& v) const  { return (dbl3(x*v.x, y*v.y, z*v.z)); }
	dbl3  operator /(const dbl3& v) const  { return (dbl3(x/v.x, y/v.y, z/v.z)); }
	bool operator ==(const dbl3& v) const  { return ((x-v.x)*(x-v.x) < verySmall) && ((y-v.y)*(y-v.y) < verySmall) && ((z-v.z)*(z-v.z) < verySmall); }
	bool operator !=(const dbl3& v) const  { return ((x-v.x)*(x-v.x) >= verySmall) || ((y-v.y)*(y-v.y) >= verySmall) || ((z-v.z)*(z-v.z) >= verySmall); }

};

inline dbl3 rotateAroundLine(dbl3 y, double gamma,  dbl3 n, dbl3 x)
{///. rotate y around line passing through x, in the direction of n, http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return dbl3(
	 	( x.x*(n.y*n.y+n.z*n.z) - n.x*( x.y*n.y+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (-x.z*n.y+x.y*n.z-n.z*y.y+n.y*y.z )*s,
		( x.y*(n.x*n.x+n.z*n.z) - n.y*( x.x*n.x+x.z*n.z-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + ( x.z*n.x-x.x*n.z+n.z*y.x-n.x*y.z )*s,
		( x.z*(n.x*n.x+n.y*n.y) - n.z*( x.x*n.x+x.y*n.y-n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (-x.y*n.x+x.x*n.y-n.y*y.x+n.x*y.y )*s );
}
inline dbl3 rotateAroundVec(const dbl3 y, double gamma, dbl3 n)
{///. rotate y around n (line passing through centre, in the direction of n) http://inside.mines.edu/~gmurray/ArbitraryAxisRotation
	double s = sinf(gamma),   c = cosf(gamma);
	double k = 1.0 - c;
	return dbl3(
		(  - n.x*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.x*c + (n.y*y.z-n.z*y.y)*s,
		(  - n.y*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.y*c + (n.z*y.x-n.x*y.z)*s,
		(  - n.z*( -n.x*y.x- n.y*y.y-n.z*y.z ) )*k + y.z*c + (n.x*y.y-n.y*y.x)*s );
}



inline dbl3  operator*(double t, const dbl3& v) { return dbl3(t*v.x, t*v.y, t*v.z); }
inline double  mag(const dbl3& v) { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
inline double  magSqr(const dbl3& v) { return (v.x*v.x+v.y*v.y+v.z*v.z); }
inline dbl3  norm(const dbl3& v)  { return  v/(mag(v)+1.0e-300); }

inline dbl3  operator*( int3 n, dbl3 x) { return dbl3(n[0]*x[0], n[1]*x[1], n[2]*x[2]);}
inline int3  operator-( int3 n, int3 x) { return int3{{n[0]-x[0], n[1]-x[1], n[2]-x[2]}};}
inline int3  operator*( int s, int3 n) { return int3{{s*n[0], s*n[1], s*n[2]}};}
inline int3  operator*( double s, int3 n) { return int3{{int(s*n[0]), int(s*n[1]), int(s*n[2])}};}
inline int3  operator/( int3 n, int s) { return int3{{n[0]/s, n[1]/s, n[2]/s}};}
inline int3& operator+=( int3& n, int3 x) { n[0]+=x[0]; n[1]+=x[1]; n[2]+=x[2]; return n;}



template<class T>
class piece
{
  protected:
	piece(): d(0), dn(0) {};
  public:

	piece(T* dd, int n): d(dd), dn(d+n) {};
	piece(const piece& p): d(p.d), dn(p.dn) {};//! note data hold by piece are not const unless piece is const itself
	piece(const std::vector<T>& vs): d(&vs[0]), dn(d+vs.size()) {};

	T* begin() const {return d;};
	T* end()   const {return dn;};
	const T* cbegin() const {return d;};
	const T* cend()   const {return dn;};
	T* operator()()   const {return d;};
	//const T& operator[](int i) const {return d[i];};
	T& operator[](int i) const {return d[i];};
	//const T* operator()() const {return d;};
	size_t size() const {return dn-d;};




	piece& operator +=(const piece& v)  { for(auto& a:*this){ a += v[&a-d];};  return (*this); }
	piece& operator -=(const piece& v)  { for(auto& a:*this){ a -= v[&a-d];};  return (*this); }
	piece& operator *=(const piece& v)  { for(auto& a:*this){ a *= v[&a-d];};  return (*this); }
	piece& operator /=(const piece& v)  { for(auto& a:*this){ a /= v[&a-d];};  return (*this); }
	piece& operator +=(T v)  { for(auto& a:*this){ a += v;};  return (*this); }
	piece& operator -=(T v)  { for(auto& a:*this){ a -= v;};  return (*this); }
	piece& operator *=(double t)             { for(auto& a:*this){ a *= t;};  return (*this); }
	piece& operator *=(int t)             { for(auto& a:*this){ a *= t;};  return (*this); }
	piece& operator /=(double t)  { return (*this)*=(1.0/t); }
	T sum() const  { T sm=0; for(auto a:*this){ sm += a;}  return sm; }


  //protected:
	T*  d;
	T*  dn;
};

template<class T>
class lazyvec: public piece<T>
{
	using piece<T>::d;
	using piece<T>::dn;

 public:

	lazyvec(): piece<T>(0,0) {};
	lazyvec(int siz): piece<T>(new T[siz],siz) {};  
	lazyvec(size_t siz, T val): piece<T>(new T[siz],siz) {  std::fill(d, dn, val);}
	lazyvec(const lazyvec& v): piece<T>(new T[v.size()],v.size())  {  std::memcpy(d, v.d, v.size()*sizeof(T)); }
	lazyvec(lazyvec&& v): piece<T>(v.d,v.size())  {  v.d=0; }
	lazyvec(const T* dd, int nn): piece<T>(new T[nn],nn)            { std::memcpy(d, dd, nn*sizeof(T)); }
	lazyvec(const T* dd, const T* de): piece<T>(new T[de-dd],de-dd)    { std::memcpy(d, dd, (de-dd)*sizeof(T));}
	~lazyvec(){ if(d) delete[] d; };

	void operator=(const lazyvec& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::memcpy(d, v.d, v.size()*sizeof(T));  dn=d+v.size(); };
	void operator=(const std::vector<T>& v){ { if(d) delete[] d; }  d = new T[v.size()]; std::memcpy(d, v.data(), v.size()*sizeof(T));  dn=d+v.size(); };

	lazyvec  operator -(void) const  { lazyvec res(*this);  for(auto& a:*this){ a = -a;};  return (*this);  }
	lazyvec  operator +(const piece<T>& v) const { return lazyvec(*this)+=v;  }
	lazyvec  operator -(const piece<T>& v) const { return lazyvec(*this)-=v;  }
	lazyvec  operator *(const piece<T>& v) const { return lazyvec(*this)*=v;  }
	lazyvec  operator /(const piece<T>& v) const { return lazyvec(*this)/=v;  }
	lazyvec  operator *(double t) const { lazyvec tmp(*this); tmp*=t;     return tmp; }
	lazyvec  operator *(int t) const { lazyvec tmp(*this); tmp*=t;     return tmp; }
	lazyvec  operator /(double t) const { lazyvec tmp(*this); tmp*=1.0/t; return tmp; }

   void resize(int nn)  
   { { if(d) delete[] d; }   if(nn) {d=new T[nn]; dn=d+nn;} else {d=0; dn=0;} }
   void pbak(T& vj)
   {	if(d) 
		{ T* od=d;  d=new T[dn+1-od];  
			std::memcpy(d, od, (dn-od)*sizeof(T));
			dn=dn+1-od; *(dn-1)=vj;  delete[] od;  }
		else    { d=new T[1];   *d=vj;   dn=d+1; }
	}
};

inline 	dbl3 round(const dbl3& vec) 	{ 	dbl3 res(vec); 	res.x=round(vec.x); res.y=round(vec.y); res.z=round(vec.z); 	return res; 	}
template<class T> 	lazyvec<T> round(const lazyvec<T>& vec) 	{ 	lazyvec<T> res(vec); 	for(auto& v:res) v=round(v); 	return res; 	}


#if __cplusplus > 201103L
 using dbl2 = std::pair<double,double>;
#else
 #define dbl2  std::pair<double,double>
#endif



#ifdef VMMLIB__VECTOR__HPP
 template< size_t M, typename T >
 Vctr<M,T>::Vctr(const dbl3& v3) 	{ 	array[ 0 ] = v3.x;	array[ 1 ] = v3.y;	array[ 2 ] = v3.z; 	}
#endif


#ifndef TOSTR
 #define TOSTR
 #if __cplusplus > 201103L
 #define toStr std::to_string
 #else
 template<typename T> std::string toStr(const T& n){  std::ostringstream stm;  stm<<n;  return stm.str();  }
 #endif
#endif



//! static vector, ideally lazyvec, but testing needed
#if __cplusplus >= 201103L
 template<class T> using  vars = lazyvec<T>; // or std::vector<T>;
 using  dbls = vars<double>;
 using  dbl3s = vars<dbl3>;
 using  ints = vars<int>;
#else
 #define vars  std::vector
 #define dbls  std::vector<double>
 #define dbl3s  std::vector<dbl3>
 #define ints  std::vector<int>
#endif


template<class T>
double sumdbl(const piece<T>& ps)  { double sm=1.0e-300; for(auto a:ps){ sm += a;}  return sm; }
template<class T>
double sumdbl(const piece<T>& ps, const piece<T>& ws)  { double sm=1.0e-300; const T* p = ps.d-1, *w=ws.d-1; while(++p<ps.dn){ sm += *w * *(++p);}  return sm; }
inline double sumdbl(const dbl3& ps, const dbl3& ws)  { return ps.x+ws.x + ps.y+ws.y + ps.z+ws.z; }
inline double sumdbl(const dbl3& ps)  { return ps.x+ ps.y+ ps.z; }
template<class T>//! Note: can be made more efficient if assuming non-zero vec
T sumvars(const piece<T>& ps, const piece<T>& ws)  { T sm; sm*=0.0; const T* p = ps.d-1, *w=ws.d-1; while(++p<ps.dn){ sm += *w * *(++p);}  return sm; }
//inline double sumdbl(const piece<double> ps, const piece<double> ws)  { double sm=1.0e-300; double* p = ps.d-1,  *w=ws.d-1; 
	 //while(++p<ps.dn){ sm += *w * *(++p);}  return sm; }









inline std::ostream& operator<< (std::ostream& out, const dbl3& node)
{
	std::ios_base::fmtflags flgs=out.flags();
	out.flags(std::ios::showpoint | std::ios::scientific);	
	out << std::setprecision(5)  << node.x   <<" " << node.y   <<" " << node.z;
	out.flags(flgs);
	return out;
}
//inline std::ofstream& operator<< (std::ofstream& out, const dbl3& node)
//{
//	std::ios_base::fmtflags flgs=out.flags();
//	out.flags(std::ios::showpoint | std::ios::scientific);	
//	out << std::setprecision(8)  << node.x   <<" " << node.y   <<" " << node.z;
//	out.flags(flgs);
//	return out;
//}

inline std::ostream& operator<< (std::ostream& out, const int3& ijk)
{
    out << ijk[0] <<" "<< ijk[1]<<" " << ijk[2];
	return out;
}
inline std::istream& operator>> (std::istream& in, dbl3& node)
{
	in >> node.x >> node.y >> node.z;
	return in;
}
inline std::istream& operator>> (std::istream& in, int3& ijk)
{
    in >> ijk[0] >> ijk[1] >> ijk[2];
    return in;
}


template<typename T1, typename T2>
inline std::istream& operator>> (std::istream& in, std::pair<T1,T2>& interval)
{
	in >> interval.first >> interval.second;
	return in;
}

template<typename T1, typename T2>
inline std::ostream& operator<< (std::ostream& out, std::pair<T1,T2>& interval)
{
	out  << interval.first<<" "<< interval.second;
	return out;
}

template<typename Type>
inline std::vector<Type>& operator*= (std::vector<Type>& vs, double scale)
{
	for(Type& v : vs) v*=scale;
	return vs;
}

template<typename T>
std::ostream & operator << (std::ostream & out, const std::vector<T>& vec)
{
	if(vec.size() && vec.size()<10)  for (auto v : vec) out << v << '\t';
	else                             for (auto v : vec) out << v << '\n';
	return out;
}

template<typename T1, typename T2>
std::ostream & operator << (std::ostream & out, const std::map<T1,T2>& vec)
{
	for (auto v : vec) out << v << '\n';
	return out;
}


template<typename T>
std::ostream & operator << (std::ostream & out, const piece<T>& vec)
{
	if(vec.size() && vec.size()<10 )  for (auto v : vec) out << v << '\t';
	else                             for (auto v : vec) out << v << '\n';
	return out;
}




#endif


