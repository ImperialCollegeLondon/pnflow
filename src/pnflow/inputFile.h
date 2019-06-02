#ifndef INPUTFILE_H
#define INPUTFILE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include "typses.h"

using namespace std;


// For compatibility with old codes, use the InputFile instead from include directory if you can


template < typename T > string myto_string( const T& n )
{
	std::ostringstream stm ;
	stm << n ;
	return stm.str() ;
}


#define outD OnDemandStream::dbgFile
///. class to write .dbg  file for debugging
class OnDemandStream
{
	public:
	ofstream prtFile;
	bool opened;
	OnDemandStream(): opened(false){};
	OnDemandStream(string fileName,int dbgMode) : opened(dbgMode) { if (dbgMode && !fileName.empty()) { prtFile.open(fileName.c_str());  opened=true;} };
	~OnDemandStream();
	OnDemandStream& operator<<(ostream& (*pfun)(ostream&)) {
		pfun(prtFile); return *this;}
	void open(string fileName) {if (!fileName.empty()) {prtFile.open(fileName.c_str()); opened=true;}};
	void close() {if (opened) {prtFile.close(); opened=false;}};
	void flush() {prtFile.flush();};

	static  OnDemandStream    dbgFile;//thread_local

};


template <class T>
OnDemandStream& operator<< (OnDemandStream& st, T val)
{
	#ifdef _debugCompile_
	if (st.opened)	{  (st.prtFile) << val;  st.prtFile.flush(); }
	#endif
   return st;
}
template <>
inline OnDemandStream& operator<< <double> (OnDemandStream& st, double val)
{
	#ifdef _debugCompile_
	if (st.opened)	{  if (val>1.0e10)	(st.prtFile)<<"+~"; else if (val<-1.0e10)	(st.prtFile)<<"-~"; else (st.prtFile)<<val;  	st.prtFile.flush();  }
	#endif
   return st;
}
inline OnDemandStream::~OnDemandStream() { outD<<"\nend"<<endl;}

//inline bool replace(const std::string& from, const std::string& to, std::string& str) {
//    size_t start_pos = str.find(from);    if(start_pos == std::string::npos)        return false;
//    str.replace(start_pos, from.length(), to);    return true;
//}

class InputFile
{

	static bool safeGetline(std::istream& is, std::string& sr)
	{
		sr.clear();
		bool continueReading=true;
		std::istream::sentry se(is, true);
		std::streambuf* sb = is.rdbuf();

		for(;;) {
		    int c = sb->sbumpc();
		    switch (c) {
		    case '\n':
		        return continueReading;
		    case '\r':
		        if(sb->sgetc() == '\n')
		            sb->sbumpc();
		        return continueReading;
		    case EOF:
		        // Also handle the case when the last line has no line ending
		        if(sr.empty())
		            is.setstate(std::ios::eofbit);
		        return false;
			case '#':
		        sr += (char)c;
		        return false;
			case ';':
		        return false;
			case ',':
			case ':':
		        sr += '\t';
		    default:
		        sr += (char)c;
		    }
		}
	}


public:

    InputFile(const std::string& fileName)
	{
		verbose = false;
		std::string prevKeyword("NO_KEYWORD_READ");
		cout<< "Loading file: " << fileName; cout.flush();
		std::ifstream in(fileName.c_str());
		if (!in) {	std::cerr << "\n\nError: Unable to open input file, " << fileName << endl;	exit(-1); }





		while(in.good())
		{
			std::string keyword, dataString;
			std::string bufferStr;
			bool oktocontinue=safeGetline(in,bufferStr);
			if(bufferStr.empty() || bufferStr[0] == '%') continue;


			{
				size_t beginKey=std::min(bufferStr.find_first_not_of(" \n\r\t"), bufferStr.size());
				size_t keyEnding=std::min(bufferStr.find_first_of("%"), bufferStr.size());
				size_t endKey=std::min(bufferStr.find_last_of(":="), keyEnding);
				if (endKey < keyEnding) endKey = std::min(bufferStr.find_last_not_of(" :=", endKey)+1, endKey);
				else endKey = std::min(bufferStr.find_first_of(" \n\r\t", beginKey+1), endKey);

				keyword = bufferStr.substr(beginKey, endKey-beginKey);
				if(keyword.empty() || keyword[0] == '%' || keyword[0] == '#' || keyword[0] == ';') continue ;

				beginKey=std::min(bufferStr.find_first_not_of(" :=",endKey), keyEnding);
				if (bufferStr[beginKey]!='%' && bufferStr[beginKey]!='#')
				{
					dataString = bufferStr.substr(beginKey, keyEnding-beginKey);
					if (!dataString.empty()) dataString += '\n';
				}

			}


			if(keyword.size() < 1 || keyword.size() > 100)
			{	std::cerr << "\n\nError: Data file contains errors after keyword: " << prevKeyword << endl << bufferStr;	exit(-1);
			}
			else
			{
				while(oktocontinue)
				{
					oktocontinue=safeGetline(in,bufferStr);
					size_t beginKey=std::min(bufferStr.find_first_not_of(" \n\r\t"), bufferStr.size());
					if (bufferStr[beginKey]!='%' && bufferStr[beginKey]!='#')
					{
						size_t keyEnding=std::min(keyword.find_first_of("%#"), bufferStr.size());
						dataString += bufferStr.substr(beginKey, keyEnding-beginKey) + '\n';
					}
				}


				m_parsedData.push_back(
				std::pair< std::string, std::string >(keyword, dataString.substr(0, std::min(dataString.find_last_not_of(" \n\r\t#")+1, dataString.size())-0))
				);
				prevKeyword = keyword;
				if (keyword=="verbose")
				{
					verbose = true;
					if (m_parsedData.rbegin()->second[0]=='f' || m_parsedData.rbegin()->second[0]=='F' || m_parsedData.rbegin()->second[0]=='n' || m_parsedData.rbegin()->second[0]=='N') verbose = false;
				}
				if(verbose) cout<< keyword << ":\n" ;
				if(verbose) cout<< m_parsedData.rbegin()->second << endl ;
				
			}
		}
		std::pair< std::string, std::string > dataEntry("", "");///. add empty data BC
		m_parsedData.push_back(dataEntry);///. add empty data BC

		in.close();
		m_baseFileName = fileName.substr(0, fileName.rfind('.'));
		
		cout<< endl;

	}





	inline void echoKeywords(std::ostream& out) const
	{
		for(unsigned i = 0; i < m_parsedData.size(); ++i)
			out << m_parsedData[i].first    << endl
				<< m_parsedData[i].second   << endl
				<< endl;
	};
   
    inline const std::string& keywordData(const std::string& keyword) const;
   
    template<typename T>
    bool getVar(T& var, const string& keyword) const
    {
		istringstream data;
		if(getData(data, keyword))	{	data>>var;	if (verbose) std::cout<<" "<<keyword<<" = "<<var<<";"<<endl; return true;	}
		else										return false;
	}
	
    bool getVar(bool& var, const string& keyword) const
    {
		istringstream data;
		if(getData(data, keyword))	{ char varC;	data>>varC;	var = (varC == 'T' || varC == 't' || varC == 'Y' || varC == 'y' || varC == '1'); cout<<keyword<<" = "<<var<<";"<<endl; return true; 	}
		else										return false;
	}

    template<typename Type>
	Type getOr(Type var, const std::string& keyword)	 const {  getVar(var, keyword);  return var;  }

    inline void errorMsg(const std::string& keyword) const;
    inline void missingDataErr(const std::string& keyword) const;
    inline void errorInDataCheck(std::istringstream& data, const std::string& keyword) const;

    inline bool getData(std::istringstream& data, const std::string& keyword) const;

	inline void Assert(bool isOK, const std::string& keyword, const std::string extraMessage="", bool severe = true) const
	{
	  if (!isOK)
	  { std::cout<< endl
			  << "========================================================"   << endl
			  << "Error: " << extraMessage                                    << endl
			  << "  File: " << m_baseFileName                                     << endl
			  << "  Keyword: "<< keyword <<"  "<< keywordData(keyword)<<";" << endl
			  << "========================================================"   << endl;
		 if (severe) exit(-1);
	  }
	}


protected:

    void removeComments(std::string& data) const{
		std::string::iterator itr;
		itr = std::find(data.begin(), data.end(), '%');
		if(itr != data.end()) data.erase(itr,data.end());
	};


    std::vector< std::pair< std::string, std::string > >         m_parsedData;
    std::string                                                  m_baseFileName;
    bool verbose;

};


/**
// Returns error message is error occurs during reading of data string
*/
inline void InputFile::errorMsg(const std::string& keyword) const
{
    std::cerr << endl
        << "================================"   << endl
        << "Error while reading input file. "   << endl
        << "         " << m_baseFileName        << endl
        << "Keyword: " << keyword               << endl
        << "================================"   << endl;    exit(-1);
}

/**
// Returns error message when required keyword is missing
*/
inline void InputFile::missingDataErr(const std::string& keyword) const
{
    std::cerr << endl
        << "============================================== "   << endl
        << "Error: could not find keyword " << keyword         << endl
        << "       in file " << m_baseFileName                 << endl
        << "============================================== "   << endl;    exit(-1);
}

/**
// Retrieves data sting based on supplied keyword, and connects a
// string stream to it. Data is removed from storage after having been
// retrived
*/
inline bool InputFile::getData(std::istringstream& data, const std::string& keyword) const
{
    for(unsigned i = 0; i < m_parsedData.size(); ++i)
    {
        if(m_parsedData[i].first == keyword)
        {
            data.str(m_parsedData[i].second);
            return true;
        }
    }
    return false;
}

inline const std::string& InputFile::keywordData(const std::string& keyword) const
{
    for(unsigned i = 0; i < m_parsedData.size(); ++i)
    {
        if(m_parsedData[i].first == keyword)
        {
			if(verbose) cout<<"Reading "<<keyword<<endl;
            return (m_parsedData[i].second);
        }
    }
    return m_parsedData[m_parsedData.size()-1].second;///. empty string
}

inline void InputFile::errorInDataCheck(std::istringstream& data, const std::string& keyword) const
{
    char leftOvers;
    data >> leftOvers;
    if(data)
    {
        cout<< endl
            << "==========================================================="   << endl
            << "Error: Too much data read for keyword: " << keyword            << endl
            << "Data read: " << data.str()                                     << endl
            << "==========================================================="   << endl;
    }
}







#endif

