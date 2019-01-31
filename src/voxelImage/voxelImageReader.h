/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.


This file is part of voxelImage library, a small c++ template library  
developed by Ali Qaseminejad Raeini for handelling 3D raw images.


Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:	a.qaseminejad-raeini09@imperial.ac.uk

\*-------------------------------------------------------------------------*/

#include <sstream>
//~ #include <streambuf>


namespace MCTProcessing
{
template<typename T> bool ignore( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	if (ins.good() && vxlImage.size3()[0]==-1) std::cout<<" ";
	
	return true;
}

template<typename T> bool fillHoles( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	unsigned int maxHoleSize;
	ins>>maxHoleSize;

		std::cout<<"fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<std::endl;
		vxlImage.fillHoles(maxHoleSize);

		vxlImage.FaceMedian06(1,5);
		//~ vxlImage.FaceMedian07(2,5);
		//~ vxlImage.FaceMedian07(2,5);
		return true;
}

template<typename T> bool selectPore( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
		std::cout<<"  converting to binary (0 and 1):"<<std::endl
			 <<"  selecting pore (->0) with values between:";
		unsigned int  thresholdMin=0,thresholdMax=0;
		ins>>thresholdMin;
		ins>>thresholdMax;

		std::cout<<" "<<int(thresholdMin)<<"  and "<<int(thresholdMax)<<"  inclusive."<<std::endl;
		vxlImage.threshold101(thresholdMin,thresholdMax);
		return true;
}

template<typename T> bool rescale( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	//~ thresholdImage=true;
	(std::cout<<"  rescaling voxel values to [ ").flush();
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(std::cout<<thresholdMin<<", "<<thresholdMax<<" ]    ").flush();
	rescale(vxlImage,T(thresholdMin),T(thresholdMax));
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool growPore( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
		std::cout<<"  growing voxels:"<<std::endl;
		int voxelValueTogrow; ins>>voxelValueTogrow;
 		char growingAlgorithm; ins>>growingAlgorithm;

		while (ins.good())		  // loop while extraction from file is possible
		{
			if (growingAlgorithm!='f')
			{
				if(voxelValueTogrow==0)
					vxlImage.growPore();
				else if (voxelValueTogrow==0)
					vxlImage.shrinkPore();
				else
				{
					std::cerr<<"growing is only implemented for binary images: "<<
					"selected voxel value to grow is "<<voxelValueTogrow << ", which is not acceptable"<<std::endl;
					return false;//error occurred
				}
			}
			else
			{
				std::cerr<<"selected growing algorithm: "<<growingAlgorithm<<
				" the only implemented algorithm is f which stands for faceGrowing"<<std::endl;
				return false;//error occurred
			}

			ins>>voxelValueTogrow;
			ins>>growingAlgorithm;
		}
		std::cout<<" done"<<std::endl;
		return true;
}


template<typename T> bool resample( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	double nResample=1;
		ins>>nResample, std::cout<<"  resampling factor: "<<nResample<<std::endl;
		vxlImage = resample(vxlImage,nResample);
		return true;
}


template<typename T> bool resampleMax( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	double nResample=1;
		ins>>nResample, std::cout<<" resampleMax factor: "<<nResample<<std::endl;
		vxlImage = resampleMax(vxlImage,nResample);
		return true;
}

template<typename T> bool redirect( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
		char direction;
		ins>>direction;
		(std::cout<<direction<<", swapping x and "<<direction<<" directions").flush();

		vxlImage.rotate(direction);
		std::cout<<std::endl;
		return true;
}

template<typename T> bool replaceRange( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	int  thresholdMin(0),thresholdMax(0); ///. Warning don't use T, uchar wont work
	ins >> thresholdMin >> thresholdMax;

	int  value=(thresholdMin+thresholdMax)/2; ///. Warning don't use T, uchar wont work
	ins >> value;

	std::cout<<" Replacing range  ["<<thresholdMin<<"  "<<thresholdMax<<"] with "<<value<<";   ";
	replaceRange(vxlImage,T(thresholdMin),T(thresholdMax),T(value));
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool crop( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	int cropBegin[3], cropEnd[3];

	std::cout<<"Crop:   ";
	for (int i=0; i<3;++i)   ins>>cropBegin[i] >>cropEnd[i],  std::cout<<cropBegin[i]<<' '<<cropEnd[i]<<"	";
	std::cout<<' '<<std::endl;

	//cropEnd[0]+=1; cropEnd[1]+=1; cropEnd[2]+=1;
	vxlImage.crop(cropBegin,cropEnd);
	return true;
}


template<typename T> bool cropD( std::stringstream & ins, voxelImageT<T>& vxlImage)
{
	int3 cropBegin(0,0,0), cropEnd=vxlImage.sizeu3();
	int nLayers(0); int value(1);
	std::cout<<"cropD:   ";
	ins>>cropBegin[0] >>cropBegin[1] >>cropBegin[2];  std::cout<<" "<<cropBegin[0] <<" "<<cropBegin[1] <<" "<<cropBegin[2]<<" -- ";  
	ins>>cropEnd[0] >>cropEnd[1] >>cropEnd[2];		std::cout<<cropEnd[0] <<" "<<cropEnd[1] <<" "<<cropEnd[2]<<"  +  ";;  
	ins >> nLayers >> value;
	std::cout<<nLayers<<" layers of "<<value<<std::endl;
	vxlImage.cropD(cropBegin,cropEnd,nLayers,value);
	return true;
}

template<typename T> bool write( std::stringstream & ins, voxelImageT<T>& voximage)
{
	std::string outName("dump.tif");
	ins >> outName;
	voximage.write(outName);
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool writeUchar( std::stringstream & ins, voxelImageT<T>& voximage)
{
	std::string outName("dump.tif");
	ins >> outName;
	double minv=-0.5, maxv=255.0;
	ins>>minv>>maxv;
	double delv=255.499999999/(maxv-minv);
	(std::cout<<minv<<" "<<maxv).flush();
	voxelImageT<unsigned char> voxels(voximage.size3(),voximage.dx(),voximage.X0(),255);
	forAlliii_(voxels) voxels(iii)=std::max(0,std::min(255,int(delv*(voximage(iii)-minv))));
	voxels.write(outName);
	(std::cout<<".").flush();
	return true;
}

template<typename T> bool read( std::stringstream & ins, voxelImageT<T>& voximage)
{
	int3 nnn = voximage.size3();
	//vxlImage.reset(int3(0,0,0),0);

	std::string fname;
	ins>>fname;
	std::cout<<"  reading from  image "<<fname<<std::endl;
	if(fname.size()>4)
	{ 
		if (  (fname.size()>=4 && fname.compare(fname.size()-4,4,".tif") == 0 )
			||	(fname.size()>=7 && fname.compare(fname.size()-7,7,".raw.gz") == 0 )
			||	(fname.size()>=4 && fname.compare(fname.size()-4,4,".raw") == 0 )
			)
		{
			  voximage.reset(nnn,0);
			  voximage.readBin(fname);
		}
		else voximage.readFromHeader(fname,0);
	}
	return true;
}
template<typename T> bool medianFilter( std::stringstream & ins, voxelImageT<T> & voximage)
{
	int nIterations(1); 
	ins >> nIterations;
	(std::cout<<"  median Filter, nIterations: "<<nIterations).flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage=median(voximage);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool FaceMedian06( std::stringstream & ins, voxelImageT<T> & voximage)
{
	if(ins.peek()=='?') { ins.str("thereshold0(2), thereshold1(4),  nIterations(1)"); return true; }
	int thereshold0(2), thereshold1(4),  nIterations(1); 
	ins >> thereshold0>> thereshold1>> nIterations;
	(std::cout<<"  FaceMedian06: "<<thereshold0<<" "<<thereshold1<<" "<<nIterations<<"     ").flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage.FaceMedian06(thereshold0,thereshold1);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}



template<typename T> bool PointMedian026( std::stringstream & ins, voxelImageT<T> & voximage)
{
	int nIterations(1),  thereshold0(11), thereshold1(15); 
	ins >> thereshold0>> thereshold1>> nIterations;
	(std::cout<<"  PointMedian026: "<<thereshold0<<" "<<thereshold1<<" "<<nIterations<<"     ").flush();
	voximage.growBox(2);
	for (int i=0; i<nIterations; ++i)
	{
		voximage.PointMedian026(thereshold0,thereshold1);
	}
	voximage.shrinkBox(2);
	(std::cout<<".").flush();
	return true;
}


template<typename T> bool circleOut( std::stringstream & ins, voxelImageT<T> & voximage)
{

	char d='z';
	ins >> d;
	int i = std::max<int>(d-'x',0);
	int X0(voximage.size3()[(i+1)%3]/2), Y0(voximage.size3()[(i+2)%3]/2);
	int R((X0+Y0)/2);

	ins >> X0 >> Y0 >> R;
	(std::cout<<"  circleOut: dir="<<d<<",  X0="<<X0 <<"  Y0="<<Y0  <<"  R="<<R ).flush();

	circleOut(voximage,X0,Y0,R,d);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool maskWriteFraction( std::stringstream & ins, voxelImageT<T> & voximage)
{
	int maskvv(2); 
	T minIelm(1), maxIelm=std::numeric_limits<T>::max();
	std::string maskname, outName("maskWriteFraction.txt");
	ins >> maskname >> outName >> maskvv >> minIelm >> maxIelm;
	(std::cout<<"  maskWriteFraction:  mask:"<<maskname <<"  outName:"<<outName<<"  maskvv:"<<maskvv  <<"  minIelm:"<<minIelm<<"  maxIelm:"<<maxIelm ).flush();

	maskWriteFraction(voximage,maskname,outName,maskvv,minIelm,maxIelm);

	(std::cout<<".").flush();
	return true;
}


template<typename T> bool Offset( std::stringstream & ins, voxelImageT<T> & voximage)
{
	dbl3 offset; 
	ins >> offset;
	(std::cout<<"  Offset:"<<offset<<" " ).flush();
	voximage.X0Ch()=offset;
	(std::cout<<".").flush();
	return true;
}


 



template<typename T> std::unordered_map<std::string,bool(*)( std::stringstream&, voxelImageT<T>&)> namedProcesses()
{

	typedef bool(*ProcessP)( std::stringstream&  ins, voxelImageT<T>& vxlImage);
	return std::unordered_map<std::string,ProcessP>{
		{  ""      , ProcessP(& ignore)},// ProcessP can be removed if using g++
		{  ";"		   , ProcessP(& ignore )},
		{  "fillHoles"   , ProcessP(& fillHoles )},
		{  "rescale"		, ProcessP(& rescale )},
		{  "pore"		, ProcessP(& selectPore )},
		{  "threshold"   , ProcessP(& selectPore )},
		{  "threshold101"   , ProcessP(& selectPore )},
		{  "resample"	, ProcessP(& resample )},
		{  "Offset"   , ProcessP(& Offset )},
		{  "direction"   , ProcessP(& redirect )},
		{  "crop"		, ProcessP(& crop )},
		{  "cropD"	   , ProcessP(& cropD )},
		{  "resampleMax" , ProcessP(& resampleMax )},
		{  "replaceRange", ProcessP(& replaceRange )},
		{  "write"  , ProcessP(& write )},
		{  "writeUchar"  , ProcessP(& writeUchar )},
		{  "read"  , ProcessP(& read )},
		{  "medianFilter"  ,ProcessP(& medianFilter )},
		{  "FaceMedian06"  ,ProcessP(& FaceMedian06 )},
		{  "PointMedian026"  ,ProcessP(& PointMedian026 )},
		{  "circleOut"  ,ProcessP(& circleOut )},
		{  "maskWriteFraction"  ,ProcessP(& maskWriteFraction )}	};
}


}

inline void getAmiraHeaderSize(const std::string& inputName, int3& n, dbl3& dx_, dbl3& X0_, int& nSkipBytes)
{
	std::ifstream headerFile(inputName);

	std::cout<<" .am:"<<inputName<<": "<<std::endl;
	while (true)
	{
		std::string tmpStr;
		headerFile>>tmpStr;
		std::cout<<tmpStr<<":"<<std::endl;

		std::stringstream keywordData;
		if(headerFile.peek()!='\n') headerFile.get (*(keywordData.rdbuf()));
		if (headerFile.fail()) {std::cout<<"Error reading "<<inputName<<",  after "<<tmpStr<<std::endl; break;}
		std::string tmp;
		if (tmpStr == "define")
		{
			keywordData >> tmp;
			keywordData >>n;
			if (tmp != "Lattice") std::cout<<" Warning: define != Lattice n3, read: "<<tmp<<std::endl;
		}
		else if (tmpStr == "BoundingBox")
		{
			keywordData >> X0_[0]>>dx_[0]>> X0_[1]>>dx_[1]>> X0_[2]>>dx_[2];
			dx_=(dx_-X0_)/n;
			std::cout<<"  dx: "<<dx_<<",   X0: "<<X0_<<std::endl;
		}
		else if (tmpStr=="@1")
		{
			std::cout<<std::endl;
			break;
		}
		else
		{
			std::cout<<"       "<<keywordData.str() <<std::endl;
		}
	}
	nSkipBytes = headerFile.tellg();
	std::cout<<"nSkipBytes:"<<nSkipBytes<<std::endl;
}


template<typename T>
void voxelImageT<T>::readFromHeader
(
	std::ifstream& headerFile,
	std::string header,
	int processKeys,
	std::string inputName
)
{
	int3 n(0,0,0);
	std::string BinaryData="XXX";
	bool X0read=false, dxread=false;
	double unit_=1.0;
	int nSkipBytes(0);
	#ifdef TIFLIB
	if ((header.size()>4 && header.compare(header.size()-4,4,".tif") == 0))
	{
		(std::cout<<  " reading tif file "<<header<<" ").flush();
		readTif(*this, header);
		std::cout<<  "."<<std::endl;
		return;
	}
	else
	#endif
	if ((header.size()>4 && header.compare(header.size()-4,4,".mhd") == 0) || (header.size()>3 && header.compare(header.size()-3,3,".py") == 0))
	{
		std::cout<<" mhd:"<<header<<": "<<std::endl;
		while (true)
		{
			std::string tmpStr;
			std::streampos begLine = headerFile.tellg();
			headerFile>>tmpStr;

			

			//~ ObjectType = Image
			//~ NDims = 3
			//~ Offset = 0 0 0
			//~ ElementSize = 8 8 8
			//~ DimSize = 200 225 153
			//~ ElementType = MET_UCHAR
			//~ ElementDataFile = Ketton100.raw
			std::stringstream keywordData;
			if(headerFile.peek()!='\n') headerFile.get (*(keywordData.rdbuf()));
			if (headerFile.fail()) break;
			std::string tmp;
			if (tmpStr == "ObjectType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "Image") std::cout<<" Warning: ObjectType != Image :="<<tmp<<std::endl;
			}
			else if (tmpStr == "NDims")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "3") std::cout<<" Warning: NDims != 3 :="<<tmp<<std::endl;
			}
			else if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				if (tmp != "MET_UCHAR") std::cout<<" Warning: ElementType != MET_UCHAR :="<<tmp<<std::endl;
			}
			else if (tmpStr == "Offset")
			{
				keywordData >> tmp; keywordData>>	X0_[0]>>X0_[1]>>X0_[2] ;
				std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2]<<std::endl ;
				X0read=true;
			}
			else if (tmpStr == "ElementSize" || tmpStr == "ElementSpacing")
			{
				keywordData >> tmp; keywordData>>	dx_[0]>>dx_[1]>>dx_[2] ;
				std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; 
				dxread=true;
				if(dx_[0]>0.01)
				{
					std::cout<<"	 Warning: too large dx (="<<dx_[0]<<"), assuming unit is um. ";
					unit_ = 1.0e-6;
				}
				std::cout<<std::endl; 
			}
			else if (tmpStr == "DimSize")
			{
				keywordData >> tmp; keywordData>>	n[0]>>n[1]>>n[2];
				std::cout<<" Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "<<std::endl; 
			}
			else if (tmpStr == "ElementDataFile")
			{
				keywordData >> tmp; if (inputName.empty()) keywordData >> inputName;

				size_t islash=header.find_last_of("\\/");
				if (islash<header.size() && inputName[0]!='/' &&  inputName[1]!=':') inputName=header.substr(0,islash+1)+inputName;
				std::cout<<" ElementDataFile = "<<inputName<<"	"<<std::endl;
			}
			else if (tmpStr == "BinaryData")
			{
				keywordData >> tmp; keywordData >> BinaryData;
				std::cout<<" BinaryData = "<<BinaryData<<"	"<<std::endl;
			}
			else if (tmpStr == "DefaultImageFormat")
			{
				std::string defSuffix;
				keywordData >> tmp; keywordData >> defSuffix;
				std::cout<<" OutputFormat = "<<defSuffix<<", suffix:"<<suffix(defSuffix)<<"	"<<std::endl; ///. sets suffix+format
			}
			else if (tmpStr == "Unit")
			{
				keywordData >> tmp; keywordData >> unit_;
				std::cout<<" Unit, OneMeter = "<<unit_<<std::endl;
			}
			else if (tmpStr == "HeaderSize")
			{
				keywordData >> tmp; keywordData >> nSkipBytes;
				std::cout<<"HeaderSize, nSkipBytes = "<<nSkipBytes<<std::endl;
			}
			else if (tmpStr!="BinaryDataByteOrderMSB" && tmpStr!="ElementByteOrderMSB" && tmpStr!="CompressedData" &&  tmpStr!="CompressedDataSize" &&  tmpStr!="TransformMatrix" &&
					 tmpStr!="ElementNumberOfChannels" && tmpStr!="CenterOfRotation" && tmpStr!="AnatomicalOrientation" && tmpStr!="AnatomicalOrientation")
			{
				headerFile.clear();
				headerFile.seekg(begLine);
				(std::cout<<" ; ").flush();
				break;
			}

		}


	}
	else if ((header.size()>3 && header.compare(header.size()-3,3,".am") == 0))
	{
		getAmiraHeaderSize(header, n,dx_,X0_,nSkipBytes);
		inputName=header;
	}
	else
	{
		std::cout<<" (depricated) _header:"<<header<<","<<std::endl;

		char tmpc;
		for (int i=0; i<8;++i)   headerFile>>tmpc, std::cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (header.size()>7 && header.compare(header.size()-7,7,"_header") == 0)  inputName=header.substr(0,header.size()-7);
		headerFile>>n[0]>>n[1]>>n[2];						// number of variables (dimension of
		std::cout<<"\n Nxyz: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"   "; std::cout.flush();
		headerFile>>	dx_[0]>>dx_[1]>>dx_[2] ;
		std::cout<<" dX: "<< dx_[0]<<"  "<<dx_[1]<<"  "<<dx_[2]<<"   "; std::cout.flush();
		headerFile>>	X0_[0]>>X0_[1]>>X0_[2] ;
		std::cout<<" X0: "<<  X0_[0]<<"  "<<X0_[1]<<"   "<<X0_[2] <<" um"<< std::endl;
		if (!headerFile)	 { std::cout<<"  Incomplete/bad header, aborting"<<std::endl; exit(-1);}
		//if (!headerFile)	 { std::cout<<"  Incomplete/bad header, continuing anyway"<<std::endl; }
		if(dx_[0]>0.01)
		{
			std::cout<<"Warning: too large dx (="<<dx_[0]<<"), assuming unit is um"<<std::endl;
			unit_ = 1.0e-6;
		}
	}



	dx_*=unit_;
	X0_*=unit_;
	if(std::abs(unit_-1.0)<epsT(float)) std::cout<<"unit= "<<unit_<<" => dx= "<<dx_<<", X0= "<<X0_<<std::endl;
	this->reset(n,1);
	if( !inputName.empty() && inputName!="NO_READ" && processKeys!=2 )
	{
	  if (inputName.compare(inputName.size()-4,4,".tif") == 0)
	  {
			dbl3 dx=dx_, X0=X0_;
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
			if(X0read) X0_=X0;
			if(dxread) dx_=dx;
	  }
	  else if ((inputName.compare(inputName.size()-4,4,".raw") == 0 && BinaryData!="False") || BinaryData=="True")
	  {
			bool readingImage = this->readBin(inputName, nSkipBytes);
			assert(readingImage);
	  }
	  else if (inputName.compare(inputName.size()-3,3,".am") == 0)
	  {
			dbl3 dx=dx_, X0=X0_;
			bool readingImage = this->readBin(inputName, nSkipBytes);
			assert(readingImage);
			if(X0read) X0_=X0;
			if(dxread) dx_=dx;
	  }
	  else if (inputName.size()>7 && inputName.compare(inputName.size()-7,7,".raw.gz") == 0)
	  {
			bool readingImage = this->readBin(inputName);
			assert(readingImage);
	  }
	  else
	  {
		std::ifstream in(inputName.c_str());
		assert(in);
		if(nSkipBytes) in.ignore(nSkipBytes);
		readAscii(in);
	  }
	}

	typedef bool(*ProcessP)( std::stringstream&  ins, voxelImageT<T>& vxlImage);


	std::unordered_map<std::string,ProcessP> name_Processes = MCTProcessing::namedProcesses<T>();

	if (processKeys)
	{


	while (true)
	{
		std::streampos begLine = headerFile.tellg();
		std::string tmpStr;
		headerFile>>tmpStr;
		//bool validKey=false;
		//cout<<tmpStr<<endl;///. keep me
		if (headerFile.fail())
		{std::cout<<" Finished reading "<<header<<":/  "<<headerFile.tellg()<<std::endl;  break; }
		else if (tmpStr[0]=='#' || tmpStr[0]=='\'' || tmpStr[0]=='/' || tmpStr[0]=='%')
		{
			headerFile.ignore(10000,'\n');
			//validKey=true;
		}
		else
		{
			auto paer = name_Processes.find(tmpStr);
			if (paer!=name_Processes.end())	
			{
				(std::cout<<" "<<tmpStr<<": ").flush();
				std::stringstream keywordData;
				if(headerFile.peek()!='\n') headerFile.get (*(keywordData.rdbuf()));
				(*(paer->second))(keywordData,*this);
				std::cout<<std::endl;
				//validKey=true;
			}
			else
			{	std::cout<<"  read "<<header<<" util entry \""<<tmpStr<<"\":/ "<<std::endl;
				headerFile.clear();
				headerFile.seekg(begLine);
				break;
			}
		}
	}


	}

}










inline std::unique_ptr<voxelImageTBase> readImage
(
	std::string headerName,
	int processKeys = 1
)
{

	std::cout<<"Openning header file: "<<headerName<<std::endl;
	std::ifstream headerFile(headerName.c_str());
	if(!headerFile)  {std::cout<<"\n\n\nError: can not open header file, "<<headerName<<std::endl<<std::endl; }
	else
	{
	#ifdef TIFLIB
	if (headerName.size()>4 && headerName.compare(headerName.size()-4,4,".tif") == 0)
		{ headerFile.close(); return readTif(headerName); }
	#endif
	if (headerName.size()>3 && headerName.compare(headerName.size()-3,3,".am") == 0)
	{
	std::cout<<"reading unsigned short .am file: "<<headerName<<std::endl;
		headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(headerName,processKeys));
	}
	if (headerName.size()>4 && headerName.compare(headerName.size()-4,4,".mhd") == 0)
	{
		while (true)
		{
			std::string tmpStr;
			headerFile>>tmpStr;



			std::stringstream keywordData;
			if(headerFile.peek()!='\n') headerFile.get (*(keywordData.rdbuf()));
			if (headerFile.fail()) break;
			std::string tmp;
			if (tmpStr == "ElementType")
			{
				keywordData >> tmp; keywordData >> tmp;
				headerFile.close();
				if (tmp=="MET_UCHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned char>(headerName, processKeys)); }
				if (tmp=="MET_CHAR")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<char>(headerName, processKeys)); }
				if (tmp=="MET_USHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned short>(headerName, processKeys)); }
				if (tmp=="MET_SHORT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<short>(headerName, processKeys)); }
				if (tmp=="MET_UINT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<unsigned int>(headerName, processKeys)); }
				if (tmp=="MET_INT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<int>(headerName, processKeys)); }
				if (tmp=="MET_FLOAT")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<float>(headerName, processKeys)); }
				if (tmp=="MET_DOUBLE")
				 { headerFile.close();	return std::unique_ptr<voxelImageTBase>(new voxelImageT<double>(headerName, processKeys)); }
				  
			}

		}
	 }
	}
	
	headerFile.close();
	return std::unique_ptr<voxelImageTBase>(new voxelImage(headerName, processKeys));

}




template<typename T>
void readConvertFromHeader
(	voxelImageT<T>& vxlImg,
	std::string headerName,
	int processKeys = 1
)
{
	std::unique_ptr<voxelImageTBase> vxlImgTup = readImage(headerName,processKeys);
	voxelImageTBase* vxlImgT = vxlImgTup.get();
	
	bool red = false;
	{auto vxlImage = dynamic_cast<voxelImageT<char>* >(vxlImgT);			if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" chars "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned char>* >(vxlImgT);   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" ucars "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<short>* >(vxlImgT);		   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" shrts "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned short>* >(vxlImgT);  if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" usrts "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<int>* >(vxlImgT);			 if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" intgs "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<unsigned int>* >(vxlImgT); 	if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" uints "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<float>* >(vxlImgT);		   if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" flots "; } }
	{auto vxlImage = dynamic_cast<voxelImageT<double>* >(vxlImgT);		  if(vxlImage) { vxlImg.resetFrom(*vxlImage); red=true; std::cout<<"read into "<<vxlImg.size3()<<" dobls "; } }

	if(!red) std::cout<<"\n\ncan not convert image\n\n"<<std::endl;
	if(!red) std::cerr<<"\n\ncan not convert image\n\n"<<std::endl;
}
