#ifndef INPUTFILE_H
#define INPUTFILE_H

// Input data file used by 3D image processing, network extraction, 
// flow simulation  and other codes
// Developed by Ali Qaseminejad Raeini. See the documentation of the 
// relevant codes for user guids and contact details, 





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

using std::string;

#if defined _MSC_VER
 #include <direct.h>
#elif defined __GNUC__
 #include <sys/types.h>
 #include <sys/stat.h>
 #include <unistd.h>
#endif

using namespace std;







inline int readInt(std::istringstream& in) {int n = 0; in>>n; return n;}
inline bool readBoolOr(std::string str, istringstream& in) { in>>str; return str[0] == 'T' || str[0] == 't';}



/// #define MULTILINEKEYDATA 1
/// #define SINGLELINEKEYDATA 0
class InputFile //! InputFile is a general input file reader, with some flexibility to chose the keyword endings etc
{
public:

	InputFile(bool multiline=true)
	:	verbose(false), informative(true), multiline_(multiline)
	{
		_parsedData.push_back(std::pair< std::string, std::string >("end", ""));//. add empty data at the end
	};

   InputFile(const InputFile& input, const std::string& title)
	:	_parsedData(input._parsedData), _fileName(input._fileName), _folder(input._folder), _name(title)
	  , verbose(input.verbose), informative(input.informative), multiline_(input.multiline_)
	{
		setKeyword("name",title);
		renameKeys("stage1","processed stage1 , ...");
		renameKeys("stage1_1","processed stage1_1 , ...");
		globals::setDebugLevel(input.getOr(0,"debugLevel"));

		if(verbose) std::cout<< " created " << _name <<std::endl;
	}

	InputFile(const std::string& fileName, bool multiline=true)
	:	verbose(false), informative(true), multiline_(multiline)
	{
		read(fileName);

		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if( (_parsedData[i].first=="include" || _parsedData[i].first=="append" ) && !_parsedData[i].second.empty() )
			{
				_parsedData[i].first = "included";
				read(_parsedData[i].second);
			}

		std::string wdir=getOr(std::string(),"workDir");
		if(!wdir.empty() && wdir!="PWD" &&  wdir!="pwd")
		{
			if(wdir=="inputDir")  {
				size_t eslshp=fileName.find_last_of("\\/")+1;
				if(eslshp<fileName.size())  wdir=fileName.substr(0,eslshp);  }
			cout<<"Changing working directory: "<<wdir<<": "<<
			chdir(wdir.c_str())             << endl;
		}

		setTitle();

		if(verbose) std::cout<< " read " << _fileName <<std::endl;
	};

	void setTitle()
	{	// call this afters etting name and/or prefix
		std::string prf=getOr(std::string(""),"prefix");
		if(prf.size())
		{
			_folder.resize(0);
			size_t slashloc=prf.find_first_of("\\/");
			if (slashloc<prf.size()-1)
			{
				_folder=prf.substr(0,slashloc+1);
				cout<<"Creating folder: "<<_folder<<"  "
				#if defined(_WIN32)
					<< mkdir(_folder.c_str()) //. check also _mkdir
				#else
					<< mkdir(_folder.c_str(), 0733) // notice that 777 is different than 0777
				#endif
					<<endl;
				prf=prf.substr(slashloc+1);
			}

			if (prf.size()>1 && (*prf.rbegin()=='/' || *prf.rbegin()=='\\') )
			{
				_folder+=prf;
				cout<<"Creating folder: "<<_folder<<"  "
					#if defined(_WIN32)
						<< mkdir(_folder.c_str()) //. check also _mkdir
					#else
						<< mkdir(_folder.c_str(), 0733) // notice that 777 is different than 0777
					#endif
					<<endl;
				prf="";
			}
		}

		if( getVar(_name, "name") || getVar(_name, "TITLE") || getVar(_name, "title") )
			                _name = prf+_name;
		else if(prf.size()) _name = prf;
		else
		{  prf = getOr(getOr(std::string(""),"network"),"networkFile");
			if(prf.empty()) prf = getOr(getOr(prf,"ElementDataFile"),"imageFile");
			if (prf.size()>7 && prf.substr(prf.size()-3,3)==".gz") prf = prf.substr(0,prf.size()-3);
			size_t dotloc=prf.find_last_of("."); if (dotloc<prf.size()) prf.erase(prf.find_last_of("."), std::string::npos);
			size_t slashloc=prf.find_last_of("\\/"); if (slashloc<prf.size()) prf=prf.substr(slashloc+1);
			_name = prf;
		}
		std::cout<< " output base name: " << _name <<std::endl;
	}

	void addKeyword(std::string key, std::string data)
	{

		if (_parsedData.rbegin()->first!="end") std::cout<<"Error in input data handelling"<<std::endl;
		_parsedData.rbegin()->first=key;   _parsedData.rbegin()->second=data;
		_parsedData.push_back(std::pair< std::string, std::string >("end", ""));//. add empty data at the end
	}
	void setKeyword(std::string key, std::string data)
	{
		if(verbose) cout<<"setting keyword "<<key<<" "<<data<<endl;

		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == key)
			{	_parsedData[i].second = data; 	return;	}
		addKeyword(key,data);
	}
	void renameKeys(std::string key, std::string newkey)
	{
		if(verbose) cout<<"replacing keys "<<key<<" with "<<newkey<<endl;
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == key)
				_parsedData[i].first = newkey;
		return;
	}

	bool safeGetline(std::istream& is, std::string& sr, bool noeq=false)
	{
		sr.clear();
		auto begl = is.tellg();
		std::istream::sentry se(is, true);
		std::streambuf* sbr = is.rdbuf();
		for(;;) {
			int cr = sbr->sbumpc();
			switch (cr) {

				case '\\': sr+=(char) sbr->sbumpc(); break; //. read next

				case '/':	if(sbr->sgetc()!='/') {  sr += '/';  break;  }
				case '%':
				case '#':
					while (sbr->sbumpc()!='\n');
					sr += '\t';
					return multiline_;

				case '=':
				case ':':
					if(noeq)
					{	sr.clear(); is.seekg(begl);
						return false;
					}
					sr += '\t';	return true;;
				case ',':  sr += '\t';	break;

				case EOF: // Also handle the case when the last line has no line ending
					if(sr.empty())   is.setstate(std::ios::eofbit);
					return false;

				case '{': case '}': cr='\t';  break;
				//case '{':  cr='\t'; do{ sr+=cr; cr=sbr->sbumpc(); }while(cr!='}' && cr!=EOF); cr='\t';  break;
				case '\'': cr='\t'; do{ sr+=cr; cr=sbr->sbumpc(); }while( cr!='\''&& cr!=EOF );   cr='\t';  break;
				case '"':  cr='\t'; do{ sr+=cr; cr=sbr->sbumpc(); }while( cr!='"' && cr!=EOF );   cr='\t';  break;
				case '[':  cr='\t'; do{ sr+=cr; cr=sbr->sbumpc(); }while( cr!=']' && cr!=EOF );   cr='\t';  break;
				case '(':  cr='\t'; do{ sr+=cr; cr=sbr->sbumpc(); }while( cr!=')' && cr!=EOF );   cr='\t';  break;

				case '\r':
					if(sbr->sgetc()=='\n') sbr->sbumpc();
				case '\n':
					cr=sbr->sgetc();
					return (cr=='\n' || cr=='\r') ? false //! double new lines are treated as end of keyword
					                              : multiline_;

				case ';':				return false;
				default:   sr += (char)cr;
			}
		}
	}


	int read(std::string fileName, int importance=2)
	{
		verbose = false;
		_fileName = fileName;

		cout<< " Loading file: " << fileName; cout.flush();
		std::ifstream in(fileName.c_str());
		if (!in) {	std::cerr << "\n\n"<<(importance>1 ? "Error" : "Warning")<<": Unable to open input file, " << fileName << endl;	if (importance>1) exit(-1);  return 0;}
		return read(in);
	}
	int read(std::istream& in)
	{
		if(_parsedData.size() && _parsedData.back().first=="end") _parsedData.pop_back();

		std::string prevKeyword("NO_KEYWORD_READ");
		while(in.good())
		{
			std::string keyword, dataString;
			std::string bufferStr;
			bool oktocontinue=safeGetline(in,bufferStr,false);
			if(bufferStr.empty()) continue;

			{
				size_t beginKey=bufferStr.find_first_not_of(" \t");
				if (beginKey == std::string::npos) continue ;
				size_t keyEnding= bufferStr.find_last_not_of(" \t")+1;
				size_t endKey=bufferStr.find_first_of(":=");
				if (endKey == std::string::npos) endKey = std::min(bufferStr.find_first_of(" \t", beginKey+1), keyEnding);
				if (endKey == std::string::npos) continue ;

				keyword = bufferStr.substr(beginKey, endKey-beginKey);

				beginKey=bufferStr.find_first_not_of(" \t",endKey);
		
				if (beginKey != std::string::npos) dataString = bufferStr.substr(beginKey, keyEnding-beginKey);

				if (keyword=="end")					break;

				if(keyword.size() < 1 || keyword.size() > 100)	{
					std::cerr << "\n\nError: Data file contains errors after keyword: " << prevKeyword << endl << bufferStr<< endl;
					std::cerr << "beginKey: " << beginKey << endl;
					std::cerr << "endKey: "   << beginKey << endl;
					exit(-1);}
			}

			{
				while(oktocontinue)
				{
					oktocontinue=safeGetline(in,bufferStr, true);
					size_t beginKey=bufferStr.find_first_not_of(" \t");
					if (beginKey==std::string::npos) continue;
					size_t keyEnding=bufferStr.find_last_not_of(" \t");
					dataString += '\n' + bufferStr.substr(beginKey, keyEnding-beginKey+1) ;
				}

				prevKeyword = keyword;
				if (keyword=="debugLevel")
				{
					globals::setDebugLevel(atoi(dataString.c_str()));
					verbose = debugLevel%2;
				}

				_parsedData.push_back(std::pair<std::string, std::string>(keyword, dataString));
				if(verbose) cout<< keyword << ":\n" << _parsedData.rbegin()->second << endl ;
			}
		}

		_parsedData.push_back(std::pair< std::string, std::string >("end", ""));//. add empty data at the end

		//if(_parsedData.back().first!="end") _parsedData.push_back({"end",""});
		cout<< endl;
		return 1;
	}



	inline void echoKeywords(std::ostream& out) const
	{
		out<<"//"<<"!-*- C -*-! input keywords:\n{";
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			out <<" "<< _parsedData[i].first  << ":\t"
				 << _parsedData[i].second  << "\n\n";
		out<< "}"<< endl;
	}

	inline const std::string& keywordData(const std::string& keyword, int importance=0) const
	{
		if(verbose) cout<<"Reading "<<keyword<<endl;
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == keyword)   return (_parsedData[i].second);

		Assert(importance<1, keyword, "missing keyword", importance>1);
		return _parsedData[_parsedData.size()-1].second;//. empty string
	}

	bool getData(std::istringstream& data, const std::string& keyword, int importance=0) const
	{
		data.clear();
		for(unsigned i = 0; i < _parsedData.size(); ++i)
			if(_parsedData[i].first == keyword) {
				data.str(_parsedData[i].second);	return true; }

		Assert(importance<1, keyword, "missing keyword", importance>1);
		return false;
	}

	bool getVar(bool& var, const std::string& keyword) const
	{
		istringstream data;
		if(getData(data, keyword))
		{
			char varC;  data>>varC;  var = (varC=='T' || varC=='t' || varC=='Y' || varC=='y' || varC=='1');
			if(verbose) cout<<keyword<<" = "<<var<<";"<<endl;
			return true;
		}
		else
			return false;
	}

	template<typename Type>
    bool getVar(Type& var, const std::string& keyword) const
    {
		istringstream data;
		if(getData(data, keyword))
		{	data>>var;	if(verbose) std::cout<<" "<<keyword<<" = "<<var<<";"<<endl; return true; }
		else						return false;
	}

    template<typename Type>
	Type getOr(Type var, const std::string& keyword)	 const {  getVar(var, keyword);  return var;  }


	inline void Assert(bool isOK, const std::string& keyword, const std::string extraMessage="", bool severe = true) const
	{
		if (!isOK)
		{	std::cout << "\n\nError in file: " << _fileName
			      << "\n  in keyword "<<keyword<<": "<< keywordData(keyword,0)<<";"
			      << "\n  "<< extraMessage<<"\n"<< endl;
			if (severe) exit(-1);
		}
	}

	inline void checkEndOfData(std::istringstream& data, const std::string& keyword, bool severe = true) const
	{
		 Assert(!data.fail(), keyword,"Incomplete/bad data", severe);
		 char leftOvers;  data >> leftOvers;
		 Assert(!data, keyword,"Too much data", severe);
	}

	std::string outputName() const { return _folder+_name; }
	std::string prefix() const { return _folder; }
	std::string name() const { return _name; }
	std::string fileName() const { return _fileName; }


	const std::vector< std::pair< std::string, std::string > >&  data() const {return _parsedData;};

protected:

	std::vector< std::pair< std::string, std::string > >        _parsedData;
	std::string                                                 _fileName;
	std::string                                                 _folder;
	std::string                                                 _name;
	bool                                                        verbose;

public:
	//. extra:
	bool                           informative;
	bool                           multiline_;
};



#endif





//TODO: move to separate file
#ifndef ONDEMANDSTREAM_H
#define ONDEMANDSTREAM_H

/// output stream to write both to std::cout, and .prt file
class mstream
{
 public:
	enum: unsigned char {PRTF = 1, STDO = 2};
	ofstream prtFile;
	unsigned char outps;
	mstream(std::string fileName, unsigned char po=PRTF|STDO)
	: outps(po) { if (!fileName.empty()) prtFile.open(fileName.c_str()); if(!prtFile) outps&=~PRTF;};
	~mstream(void){};
	mstream& operator<<(ostream& (*pfun)(ostream&)) {if(outps&PRTF) pfun(prtFile); if(outps&STDO) pfun(std::cout); return *this;}
	ofstream& fileStream() {return prtFile;};
};

template <class T>
mstream& operator<< (mstream& st, T val)
{
  if(st.outps & mstream::PRTF) st.prtFile << val;
  if(st.outps & mstream::STDO) std::cout<< val;
  return st;
}

/// output stream to write .dbg  file for debugging
class OnDemandStream
{
	public:
	ofstream prtFile;
	bool opened;
	OnDemandStream(): opened(false){};
	OnDemandStream(std::string fileName,int dbgMode) : opened(dbgMode) { if (dbgMode && !fileName.empty()) { prtFile.open(fileName.c_str());  opened=true;} };
	~OnDemandStream(void){};
	OnDemandStream& operator<<(ostream& (*pfun)(ostream&)) {
		pfun(prtFile); return *this;}
	void open(std::string fileName) {if (!fileName.empty()) {prtFile.open(fileName.c_str()); opened=true;}};
	void close() {if (opened) {prtFile.close(); opened=false;}};
	void flush() {if (opened) {prtFile.flush();} cout.flush(); };

	 static OnDemandStream    dbgFile; //thread_local

};

#define outD OnDemandStream::dbgFile


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


#endif
