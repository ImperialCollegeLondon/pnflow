#ifndef CONFIG_H
#define CONFIG_H


#include "inputFile.h"
#include "voxelImage.h"
#include "ElementGNE.h"



#define MICROCTDATA 1
#define BINARYDATA 2
#define ASCIIDATA 3
#define MHDDATA 4

#ifndef PORORANGE_H
#define PORORANGE_H
class poroRange :
public std::pair<unsigned char,unsigned char>
{
public:
	poroRange(unsigned char lower, unsigned char upper) : std::pair<unsigned char,unsigned char>(lower,upper){};
	poroRange(std::string nam, unsigned char lower,unsigned char upper)
	  : std::pair<unsigned char,unsigned char>(lower,upper),name(nam){};
	poroRange(const poroRange& copyt): std::pair<unsigned char,unsigned char>(copyt),name(copyt.name) {};
	bool outside(unsigned char value){return this->first > value  ||  value > this->second; };
	std::string name;
};
#endif //PORORANGE_H


inline int generate_input_nextract(string fileName, string opts)
{
	if (fileName.empty()) fileName = "input_image.mhd";

	std::ifstream in(fileName.c_str());
	if (in) {
		cout << "\n\nError: file "<<fileName<<" exists," << endl
		     <<      "  to run simulation: rerun with "<< fileName<<" as the only argument "<< endl
		     <<      "  to regenerate: delete it and try again" << endl;
		return 0; }

	ofstream of(fileName);
	if(opts=="-g")
	{
	 of	<<"ObjectType =  Image\n"
		<<"NDims =       3\n"
		<<"ElementType = MET_UCHAR\n"
		<<"ElementByteOrderMSB = False\n"
		<<"ElementNumberOfChannels = 1\n"
		<<"CompressedData = True\n\n"
		<<"HeaderSize = 0\n"
		<<"DimSize =    	1000	1000	1000\n"
		<<"ElementSize = 	1.6 	1.6 	1.6\n"
		<<"Offset =      	0   	0   	0\n"
		<<"\n"
		<<"ElementDataFile = input_image.raw.gz\n"
		<<"\n\n"
		<<"//! The above keywords are compatible with mhd format and \n"
		<<"//! can be used to open the Image in Fiji/ImageJ or Paraview\n"
		<<"//! The following commands are optional, remove the \"//\" to activate them\n"
		<<"\n"
		<<"//DefaultImageFormat = tif\n"
		<<"\n"
		<<"//!______________  image processing  commands _________________\n"
		<<"\n"
		<<"//! crop image to  [ Nxyz_begin  Nxyz_end )"
		<<"//cropD                0 0 0    300 300 300 \n"
		<<"\n"
		<<"//! flip x direction with y or z \n"
		<<"//direction z\n"
		<<"\n"
		<<"//! manipulate voxel values:\n"
		<<"//!    range   [start...end] -> value \n"
		<<"//replaceRange   0     127       0 \n"
		<<"//replaceRange   128   255       1 \n"
		<<"\n"
		<<"//! threshold image: range -> 0 (void-space), rest->1 (solid) \n"
		<<"//threshold   0  128 \n"
		<<"\n\n\n";
	 of	<<"//!_______________  network extraction keywords __________________\n"
		<<"//!______(should be after image processing  commands above) ______\n"
		<<"\n"
		<<"//title:   output_network"
		<<"\n"
		<<"//write_all:	true; // use `write_all` is a memorable alternative to all visualization keywords "
		<<"// write_radius:	true\n"
		<<"// write_statistics:	true\n"
		<<"// write_elements:	true\n"
		<<"// write_poreMaxBalls:	true\n"
		<<"// write_throatMaxBalls:	true\n"
		<<"// write_throats:	true\n"
		<<"// write_poroats:	true\n"
		<<"// write_hierarchy:	true\n"
		<<"// write_medialSurface:	true\n"
		<<"// write_throatHierarchy:	true\n"
		<<"// write_vtkNetwork:	true\n"
		<<endl;
	}
	else
		cout<<"Error unknown option (first argument)"<<endl;

		of.close();

		cout <<" file "<<fileName<<" generated, edit: set Image size, name etc, and rerun\n";	
	return 0;
}





class inputDataNE : public InputFile
{
public:
 inputDataNE(std::string file)
	: InputFile(file,false), precision(1), X0(0.0,0.0,0.0),  datatype(BINARYDATA), invalidSeg{-10000, 255}
	{}

 inputDataNE(std::string file, int minvoid, int maxvoid, int outsiderange, bool readCfg)
	: InputFile(false),  precision(1), X0(0.0,0.0,0.0),  datatype(MHDDATA), invalidSeg{-10000, 255}
	{
		if(!readCfg)
		{
			setKeyword("ElementDataFile",file);
			setKeyword("void_range",toStr(minvoid)+" "+toStr(maxvoid));
			if(outsiderange<256) setKeyword("outside_range", toStr(outsiderange));
		}
		else	read(file);
		setTitle();
		//echoKeywords(cout);
	}

 void init(bool verbos=true)
 {

	cout<< "Reading inputDataNE data:"<<endl;

	std::istringstream inputKeyData;


	Assert(getVar(imgfileName,"ElementDataFile") || getVar(imgfileName,"imageFile") || !verbos, "ElementDataFile or imageFile", "keyword not found");
	if(verbos)
		cout<<" image file: "<<imgfileName<<endl;



	Assert(getData(inputKeyData, "DimSize") ||  getData(inputKeyData, "imageSize") || !verbos, "DimSize or imageSize", "keyword not found");
		inputKeyData >> nx >> ny >> nz;
	if(verbos)
		cout<< " image size: " << nx << " " << ny << " " << nz <<endl;


	Assert(getVar(precision, "ElementSpacing") || getVar(precision, "voxelSize") || getVar(precision, "ElementSize") || !verbos, "ElementSpacing or voxelSize", "keyword not found");
	if(verbos)
		cout<< " voxel size: " <<precision <<endl;


	std::string dataType("binary char");
	if (getVar(dataType,"ElementType") || getVar(dataType,"fileType"))
	{
		if (dataType == "binary")	datatype = BINARYDATA;
		if (dataType == "MET_UCHAR")	datatype = MHDDATA;
		else if	(dataType == "microct")	                        datatype = MICROCTDATA;
		else if	(dataType == "ascii")	                        datatype = ASCIIDATA;
		else cout<<"wrong data type, going with binary"<<endl;
		cout<<"  file type: "<<dataType <<endl;
	}

	if (!getVar(imgfrmt,"DefaultImageFormat")) imgfrmt=".tif";
	if(imgfrmt[0]!='.') imgfrmt="."+imgfrmt;
	suffix(imgfrmt);
	cout<<"DefaultImageFormat: "<<imgfrmt<<endl;
	

	cout<<" voxel indices:"<<endl;
	if (getData(inputKeyData,"voidSpace"))
	{
		std::string rockTypeName;
		inputKeyData>>rockTypeName;
		cout<<"  "<<0<<": rockTypeName = "<<rockTypeName<<endl;
		poroRange ithRockType(rockTypeName,0,0);

		_rockTypes.push_back(ithRockType);
	}
	else
	{
		poroRange ithRockType("void",0,0);
		_rockTypes.push_back(ithRockType);
		cout<<"  "<<0<<": void voxels "<<endl;
	}

	int nRTypes(0);
	if(getData(inputKeyData,"porousRocks"))
	{
		inputKeyData >> nRTypes;
		for(int i = 1; i <=  nRTypes; ++i)
		{
				std::string rockTypeName;
				inputKeyData>>rockTypeName;
				cout<<"  "<<i<<": porousRocks = "<<rockTypeName<<endl;
				poroRange ithRockType(rockTypeName,i,i);
				_rockTypes.push_back(ithRockType);
		}
		cout<< "  number of rock types: "<<_rockTypes.size()<<endl;
	}

	segValues.resize(256, _rockTypes.size());

	for(int i = 0; i <=  nRTypes; ++i)
	{
		if(getData(inputKeyData, _rockTypes[i].name+"_range") || getData(inputKeyData, _rockTypes[i].name+"_thresholds"))
		{
			int lower,upper;
			inputKeyData >> lower >> upper;
			 _rockTypes[i].first = lower;
			 _rockTypes[i].second = upper;
			if (_rockTypes[i].first > _rockTypes[i].second)
				cout<<"  Wrong entries for keyword \""<<_rockTypes[i].name+"_range"<<"\":\n"<<keywordData(_rockTypes[i].name+"_thresholds")<<"\n lower value is higher than upper value"<<endl;
		}

		for(size_t j = _rockTypes[i].first; j <=  _rockTypes[i].second; ++j)
			segValues[j] = i;

		cout<<"  "<< _rockTypes[i].name<<" voxel values: ["<<int(_rockTypes[i].first)<< " "<<int(_rockTypes[i].second)<<"]"<<endl;
	}

	cout<<"  Voxel value indices:";
	for(size_t i = 0; i <=  12; ++i)		cout<<" "<<segValues[i];
	cout<<" ... "<<endl;




 }


 void readImage()
 {
	cout<< "\nLoading voxel data, format:"<<datatype <<" fileName:"<<fileName()<<endl;

	VImage.reset(nx,ny,nz,0);
	if (datatype == MICROCTDATA)	VImage.readMicroCT(imgfileName);
	else if (datatype == BINARYDATA)	{if (!VImage.readBin(imgfileName))   {cout<<"\nError: didn't read binary image!\n"<<endl; exit(-1);}}
	else if (datatype == ASCIIDATA)	{if (!VImage.readAscii(imgfileName)) {cout<<"\nError: didn't read ascii image!\n"<<endl; exit(-1);}}
	else
	{
		VImage.reset(0,0,0,255);
		VImage.readFromHeader(fileName());
		int3 siz=VImage.size3();
		nx= siz[0];  ny= siz[1];  nz= siz[2];
		precision = VImage.dx()[0];
		X0=VImage.X0();
	}
	if(!VImage.size3()[0]) {  cout<<"\nError: no image read!\n"<<endl; exit(-1);}
	VImage.printInfo();

	nInside= (long long)(nx)*ny*nz;
	std::istringstream inputKeyData;
	if(getData(inputKeyData, "outside_range") || getData(inputKeyData,"outside_thresholds"))
	{
		int lower,upper;
		inputKeyData >> lower >> upper;
		forAllcp(VImage) if(lower<=(*cp) && (*cp)<=upper) --nInside;

	}
 }

 void createSegments()
 {




	nVxlVs.resize(_rockTypes.size()+1,0);

	std::vector<segment> segTmp(nx+1);
 	segs_.resize(nz,std::vector<segments>(ny));

 	for (int iz = 0; iz<nz; iz++)
 	 for (int iy = 0; iy<ny; iy++)
 		if(segs_[iz][iy].s != NULL) cout<<"ERROR"<<endl;


	for (int iz = 0; iz<nz; ++iz)
	{
 		for (int iy = 0; iy<ny; ++iy)
 		{
			int cnt = 0;
			int  currentSegValue = 257;
			for (int ix = 0; ix<nx; ++ix)
			{
				unsigned char vV = VImage(ix,iy,iz);
				if (segValues[vV] != currentSegValue)
				{
					currentSegValue = segValues[vV];
					segTmp[cnt].start = ix;
					segTmp[cnt].value = currentSegValue;
					++cnt;
				}
				++(nVxlVs[segValues[vV]]);
			}

			segments & ss = segs_[iz][iy];
			ss.cnt = cnt;
			ss.reSize(cnt+1);
			if (segs_[iz][iy].s == NULL)		{cout<<"\n   iz "<<iz<<" iy "<<iy<<" cnt"<<cnt<<" ERROR XX XX"<<endl;	}

			for (int i = 0; i<cnt; ++i)
			{
				ss.s[i].start = segTmp[i].start;
				ss.s[i].value = segTmp[i].value;
				if (i>0 && ss.s[i].value == ss.s[i-1].value) 	cout<<"\n   ERROR XXX"<<i<<" "<<int(ss.s[i-1].value)<<" "<<int(ss.s[i].value)<<"    "<<int(segValues[3])<<"    "<<int(VImage(ss.s[i-1].start,iy,iz))<<" "<<int(VImage(ss.s[i].start,iy,iz))<<endl;
			}
			ss.s[cnt].start = nx;
			ss.s[cnt].value = 254;
 		}
 	}


	cout<< endl;
	for (int i = 0;i<int(_rockTypes.size());i++)
		cout<<" "<<i<< ". " << _rockTypes[i].name<< ": " << nVxlVs[i] << " voxels, " <<  (nVxlVs[i]*(100.0/nx/ny/nz))<< "%"<<endl;
	cout<< endl;

 }


 const segment* segptr(int i, int j, int k) const
 {
 	if (i<0 || j<0 || k<0 || i>= nx || j>= ny || k>= nz)  return &invalidSeg;

 	const segments& s = segs_[k][j];
 	for (int p = 0; p<s.cnt; ++p)
 		if (i >= s.s[p].start && i < s.s[p+1].start)	  return s.s+p;

	cout<<"Error can not find segment at "<<i<<" "<<j<<" "<<k<<" nSegs: "<<s.cnt<<endl;
 	return (s.s+s.cnt);
 }


	std::string netsufix() const { return (flowBaseDir.empty() ? "DS0" : (flowBaseDir.back()=='/' ? "DS1": "DS4")); }
public:

	int nx, ny, nz;
	double precision;
	dbl3 X0;
	int datatype;
	std::string imgfrmt;
	long long nInside;

	std::string imgfileName;
	std::string flowBaseDir;

	std::vector< std::vector<segments> > segs_;
	segment invalidSeg;

	std::vector<size_t> nVxlVs;
	std::vector<int> segValues;
	std::vector<poroRange> _rockTypes;
	voxelImage VImage;

};

#endif
