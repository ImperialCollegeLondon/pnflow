#ifndef GLOBALS_H
#define GLOBALS_H

// definition of global class and macros used to debugging and testing  

const static double PI = 3.14159265358979;

class globals
{
    static int    debugLevel_;
    //static size_t typeCounter_;//dynamic type index, not used yet, see https://codereview.stackexchange.com/questions/44936/unique-type-id-in-c
public:
    //template<typename T>
    //static size_t getTypeCounter()
    //{
        //static size_t id = typeCounter_++;
        //return id;
    //}
    static int getDebugLevel()
    {
        return debugLevel_;
    }
    static void setDebugLevel(int dbgLvl)
    {
        debugLevel_=dbgLvl;
    }
};
//template<typename T> const auto  TypeID = globals::getTypeCounter<T>;
#define debugLevel    globals::getDebugLevel()





//using toStr = std::to_string  is bad in decimal notation
template<typename T> std::string toStr(const T& n){  std::ostringstream ss;  ss<<n;  return ss.str();  }
template<typename T> T strTo(const std::string &s){  std::istringstream ss(s);  T t;  ss>>t;  return t; }



//- Debugging  testing macros

inline void _cerr_(std::string msg="", bool xit=false) // for debugger breakpoints: don't optimize out please !!!
	{ 	 if(xit) throw std::runtime_error(msg); 	 else  std::cerr<< msg <<std::endl; 	 }

// Variable argument macro trick
 #define ERR_HDR(isOk) " Error  "+std::string("in ")+ std::string(__FUNCTION__)+", " \
                       +std::string(__FILE__)+":"+toStr(__LINE__) +std::string(":  ")+std::string(#isOk)
 #define ensure1(isOk)           if(!(isOk)) _cerr_(ERR_HDR(isOk))
 #define ensure2(isOk, msg)      if(!(isOk)) _cerr_(ERR_HDR(isOk)+" \\\n *** "+msg+" *** \n")
 #define ensure3(isOk, msg, xit) if(!(isOk)) _cerr_(ERR_HDR(isOk)+" \\\n *** "+msg+" *** \n", xit)
 #define GET_MACRO3(_1,_2,_3,NAME,...) NAME

//! Validation/production phase ensure/assert. Usage:
//!   \code{.cpp} ensure(condition, "message", throw_on_error=false); \endcode
 #define ensure(...)  GET_MACRO3(__VA_ARGS__, ensure3, ensure2,ensure1, "Only 1 to 3 args please")(__VA_ARGS__)


//#define _debugCompile_
//! Debug error message, similar to ensure, but printed only if _debugCompile_ is defined and debugLevel is not zero
#ifndef _debugCompile_
  #define d_assert(...)
#else
  //#define d_assert(...) ensure((!debugLevel)||__VA_ARGS__)
  #define d_assert(...) if(debugLevel) ensure(__VA_ARGS__)
#endif



#endif
