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
Ali Q Raeini: a.qaseminejad-raeini09@imperial.ac.uk
\*-------------------------------------------------------------------------*/

#include "voxelImage.h"
#include <map>
#include <limits>   // std::numeric_limits

#ifdef ZLIB
#include "zfstream.h"
#endif

#ifdef TIFLIB
#include "voxelTiff.h"
#endif

//using namespace std;

#ifndef voxelImageMacros_H
#define voxelImageMacros_H



#ifdef OpenMP
 #define OMPragma(_args_) _Pragma(_args_) 
#else
 #define OMPragma(_args_)
#endif

 //#define forAllitr12(datas3s_M1,datas3s_M2)   for (const auto itr1 datas3s_M1.data_.begin(), itr2 datas3s_M2.data_.begin(); itr2<datas3s_M2.end() ; ++itr1, ++itr2)
 #define 	forAlliii_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for ( size_t iii=0; iii<datas3s_M.data_.size() ; ++iii )
 
 #define forAllvp_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for(auto vp=datas3s_M.data_.data(); vp<&(*datas3s_M.data_.cend()); ++vp )
 #define 	forAll_vp_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for(auto vp=datas3s_M.data_.data()+datas3s_M.data_.size()-1; vp>&(*datas3s_M.data_.cbegin())-1; --vp )
 
 
 #define 	forAllkvp_(datas3s_M)   \
	OMPragma("omp parallel for") 		\
 	for ( int k=0; k<(datas3s_M).size3()[2] ; ++k )   \
 	for(register auto vp=datas3s_M.data_.data()+k*datas3s_M.nij_, _ve_=datas3s_M.data_.data()+(k+1)*datas3s_M.nij_; vp<_ve_; ++vp )
 #define 	forAllk_vp_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for ( int k=0; k<(datas3s_M).size3()[2] ; ++k )   \
 	for(register auto vp=datas3s_M.data_.data()+(k+1)*datas3s_M.nij_-1, _ve_=datas3s_M.data_.data()+k*datas3s_M.nij_-1; vp>_ve_; --vp )

 #define 	forAllkji(datas3s_M)   \
 	for ( int k=0; k<(datas3s_M).size3()[2] ; ++k )   \
 	for ( int j=0; j<(datas3s_M).size3()[1] ; ++j )   \
 	for ( int i=0; i<(datas3s_M).size3()[0] ; ++i )
 
 #define 	forAllkji_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for ( int k=0; k<(datas3s_M).size3()[2] ; ++k )   \
 	for ( int j=0; j<(datas3s_M).size3()[1] ; ++j )   \
 	for ( int i=0; i<(datas3s_M).size3()[0] ; ++i )
 #define 	forAllkji_1_(datas3s_M)   \
 	OMPragma("omp parallel for") 		\
 	for ( int k=1; k<(datas3s_M).size3()[2]-1 ; ++k )   \
 	for ( int j=1; j+1<(datas3s_M).size3()[1] ; ++j )   \
 	for ( int i=1; i+1<(datas3s_M).size3()[0] ; ++i )
 

//#define 	forAllkji_1(datas3s_M)   
	//for ( int k=1; k+1<(datas3s_M).size3()[2] ; ++k )   
	//for ( int j=1; j+1<(datas3s_M).size3()[1] ; ++j )   
	//for ( int i=1; i+1<(datas3s_M).size3()[0] ; ++i )
#define forAllkji_1(datas3s_M)   \
	for (register int k=1,_ke_=(datas3s_M).size3()[2]-1; k<_ke_; ++k )   \
	for (register int j=1,_je_=(datas3s_M).size3()[1]-1; j<_je_; ++j )   \
	for (register int i=1,_ie_=(datas3s_M).size3()[0]-1; i<_ie_; ++i )
 
 #define forAllvp_seq(datas3s_M)   for(auto vp=datas3s_M.data_.data(); vp<&(*datas3s_M.data_.cend()); ++vp )
 #define forAllcp(datas3s_M)   for(auto cp=datas3s_M.data_.cbegin(); cp<datas3s_M.data_.cend(); ++cp)
 #define forAllcp_seq_(datas3s_M)   for(const auto* cp=datas3s_M.data_.data(), *_ve_=&(*datas3s_M.data_.cend()); cp<_ve_; ++cp )
 #define forAllvv_seq(datas3s_M)   for (auto vv : datas3s_M.data_ )
 #define forAllvr_seq(datas3s_M)   for (auto& vr : datas3s_M.data_ )
 




 
#define forAllNei(r1,r2) \
for (int k_nei_rel=r1;k_nei_rel<=r2;++k_nei_rel) \
for (int j_nei_rel=r1;j_nei_rel<=r2;++j_nei_rel) \
for (int i_nei_rel=r1;i_nei_rel<=r2;++i_nei_rel)
#define nei(datas3s_M,i_M,j_M,k_M) (datas3s_M)(i_M+i_nei_rel, j_M+j_nei_rel, k_M+k_nei_rel)

#define forAllNeiInt(ri1,ri2,rj1,rj2,rk1,rk2) \
for (int k_nei_abs = rk1;k_nei_abs <= rk2;++k_nei_abs) \
for (int j_nei_abs = rj1;j_nei_abs <= rj2;++j_nei_abs) \
for (int i_nei_abs = ri1;i_nei_abs <= ri2;++i_nei_abs)
#define neib(datas3s_M) (datas3s_M)(i_nei_abs, j_nei_abs, k_nei_abs)

#define sum6Nei(datas3s_M, vp_) \
		(datas3s_M.v_i(-1,vp_)+datas3s_M.v_i( 1,vp_)+\
		 datas3s_M.v_j(-1,vp_)+datas3s_M.v_j( 1,vp_)+\
		 datas3s_M.v_k(-1,vp_)+datas3s_M.v_k( 1,vp_))

#define min7Nei(datas3s_M, vp_) \
		min(*vp_,min(min( \
		min(datas3s_M.v_i(-1,vp_),datas3s_M.v_i( 1,vp_)),\
		min(datas3s_M.v_j(-1,vp_),datas3s_M.v_j( 1,vp_))),\
		min(datas3s_M.v_k(-1,vp_),datas3s_M.v_k( 1,vp_))))
#define min5Nei(datas3s_M, vp_) \
		min(*vp_, min( \
		min(datas3s_M.v_i(-1,vp_),datas3s_M.v_i( 1,vp_)),\
		min(datas3s_M.v_j(-1,vp_),datas3s_M.v_j( 1,vp_))))


//#define OPERATOR(M_func_,Type_) static_cast<Type_ const & (*) (Type_ const &, Type_ const &)>(M_func_<Type_>)

#endif  /// voxelImageMacros_H

template<typename  T> inline T maxNei(const voxelImageT<T>& image, int i, int j, int k, int r11, int r22)
{
	T maxx=std::numeric_limits<T>::min();
	forAllNei(r11,r22)
	{
		maxx=std::max(maxx,nei(image,i,j,k));
	}
	return maxx;
}

#define  accumulatedbl  accumulate
template<typename T,typename TRes>
TRes accumulate(const voxelField<T>& vf, TRes const & (*operatorFunc)(TRes const &,  TRes const&) , TRes result=0)
{	
	//result = std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);   return result; 

    OMPragma("omp parallel")
    {
        TRes res_private = result;
        OMPragma("omp for nowait")
        for(auto vp=vf.data_.begin(); vp<vf.data_.end(); ++vp )
        {
				res_private=operatorFunc(res_private,*vp);
		}
        OMPragma("omp critical")
        result = operatorFunc(res_private,result);
    }
    return result;
}

template<typename T,typename TRes, typename TFunc>
TRes accumulate(const voxelField<T>& vf, TFunc operatorFunc , TRes result=0)
{	
	//result = std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);   return result; 

	OMPragma("omp parallel")
	{
		TRes res_private = result;
		OMPragma("omp for nowait")
		for(auto vp=vf.data_.begin(); vp<vf.data_.end(); ++vp )
		{
			res_private=operatorFunc(res_private,*vp);
		}
		OMPragma("omp critical")
		result = operatorFunc(res_private,result);
	}
	return result;
}

template<typename T>
int accumulateT(const voxelImage& vf, std::function<T(T, T)> operatorFunc, T result=0)
{
	return std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);
}





template<typename Type>   void voxelField<Type>::reset(int3 n)
{
	nij_=size_t(n[0])*n[1];
	this->data_.resize(n[2]*nij_);
	nnn_=n;
}

template<typename Type>   void voxelField<Type>::reset(int3 n, Type value)
{
	this->data_.resize(0);
	nij_=size_t(n[0])*n[1];
	this->data_.resize(size_t(n[2])*nij_,value);
	nnn_=n;
}




template<typename Type>   void voxelField<Type>::getSize(int& n1, int& n2, int& n3) const
{
	n3 = (*this).size3()[2];
	if (n3>0)
	{
	  n1 = (*this).size3()[0];
	  n2 = (*this).size3()[1];
	}
	else
	{
	  n1 = 0;
	  n2 = 0;
	}
}



template<typename Type>
int3 voxelField<Type>::sizeu3() const
{
   return nnn_;
}

template<typename Type>
int3 voxelField<Type>::size3() const
{
   return nnn_;
}





//read order sensitive
template<typename Type>   void voxelField<Type>::readAscii(std::ifstream& in)
{
	forAllvr_seq((*this))
				in>>vr;

}

//read order sensitive
template<typename Type>   bool voxelField<Type>::readAscii(std::string fileName)
{
	std::cout<<  " reading ascii file "<<fileName<<std::endl;
	std::ifstream in(fileName.c_str());
	if(!in)  {std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl;
			   return false;}


	readAscii(in);

	in.close();
	return !in.fail();
}




//read order sensitive
template<typename Type>   void voxelField<Type>::readMicroCT(std::string fileName)
{
	std::cout<<  " reading micro-CT file "<<fileName<<std::endl;
	std::ifstream in(fileName.c_str());
	if(!in)  {std::cout<<"\n\nError: can not open image file, "<<fileName<<std::endl<<std::endl; exit(-1);}

	char tmpc;
	for ( int i = 0; i<8;i++)   in>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)
	int n[3];
	double  xmin[3];
	double  xmax[3];

	in>>n[2]>>n[1]>>n[0];						// number of variables (dimension of
	//in>>	dx[0]>>dx[1]>>dx[2] ;
	in>>	xmin[0]>>xmin[1]>>xmin[2] ;
	in>>	xmax[0]>>xmax[1]>>xmax[2] ;

	readAscii(in);

	in.close();
}


template<typename Type>   bool voxelField<Type>::readBin(std::string fileName, int nSkipBytes)
{
	int3 n = size3();
	size_t siz = size_t(n[0])*n[1]*n[2];

	#ifdef TIFLIB
	if(fileName.compare(fileName.size()-4,4,".tif")==0)
	{ 
		(std::cout<<  " reading tif file "<<fileName<<" size:"<<siz<<"*"<<sizeof(Type)).flush();
		readTif(*this, fileName);
		std::cout<<  "."<<std::endl;
		return true;
	}
	#endif 


	if(fileName.compare(fileName.size()-3,3,".gz")==0)
	{

	 #ifdef ZLIB
		if(std::ifstream(fileName.c_str()).good())
		{
			(std::cout<<  " reading compressed file "<<fileName<<" size:"<<siz<<"*"<<sizeof(Type)).flush();
			gzifstream  in(fileName.c_str());
			//in << setcompression(Z_NO_COMPRESSION);
			in.read(reinterpret_cast<char*>(&((*this)(0,0,0))), siz*sizeof(Type) );
			in.close();

			std::cout<<"."<<std::endl;
			return true;
		}else std::cout<<"Error: can not be read "<<fileName<<std::endl;	
	 #endif

		fileName=fileName.substr(0,fileName.size()-3);
		std::cout<<".gz not read or not supported, trying "<<fileName<<" instead"<<std::endl;	

	}

	(std::cout<<  " reading binary file "<<fileName<<" size:"<<siz<<"*"<<sizeof(Type)).flush();
	std::ifstream in (fileName.c_str(), std::ios::in | std::ios::binary);
	if(!in)  {std::cout<<"\n\n  ***** Error: could not open image file, "<<fileName<<"  *****\n"<<std::endl;		return false;}
	if(nSkipBytes) in.ignore(nSkipBytes);
	in.read(reinterpret_cast<char*>(&((*this)(0,0,0))), siz*sizeof(Type) );

	if (!in) 	 std::cout<<  "\n\n ***** Error in reading "<<fileName<<", only  read "<<in.gcount()<<"B ***** \n"<<std::endl;
	in.close();
	std::cout<<  "."<<std::endl;
	return true;

}


template<typename Type>   bool voxelField<Type>::readBin(std::string fileName,
				int iS,int iE ,	int jS,int jE ,	int kS,int kE, int nSkipBytes
)
{
	if(fileName.compare(fileName.size()-4,4,".tif")==0)
	{ 
		voxelImageT<Type> vxls(fileName);
		if (vxls.size3()[0]!=iE-iS)	std::cout<<"Error in reading "<<fileName<<", unexpected size: Nx="<<vxls.size3()[0]<<"  !="<<iE-iS<<std::endl;
		setBlock(iS,jS,kS, vxls);
		return true;
	}
	if(fileName.compare(fileName.size()-3,3,".gz")==0)
	{
		voxelField<Type> vxls(int3(iE-iS, jE-jS, kE-kS));
		vxls.readBin(fileName);
		setBlock(iS,jS,kS, vxls);
		return true;
	}
	
	(std::cout<<  " reading binary file "<<fileName).flush();
	std::cout<<"\n"<<"@@=  "; std::cout.flush();
	std::ifstream in (fileName.c_str(), std::ios::in | std::ios::binary);
	if(!in)  {std::cout<<"\n\n  Error: can not open image file, "<<fileName<<std::endl<<std::endl;		return false;}
	if(nSkipBytes) in.ignore(nSkipBytes);
	int k = kS;
	for ( ;k < kE;k++)
	{
		for ( int j = jS;j < jE;j++)
		{
			if(in)  in.read(reinterpret_cast<char*>(&((*this)(iS,j,k))), (iE-iS)*sizeof(Type) );
		}
		if (!in) break;
	}

	if (!in)	std::cout<<  "\n\n ***** Error in reading "<<fileName<<" ***** \n"<<"only "<<k<<" layers read"<<std::endl;

	in.close();
	std::cout<<  "."<<std::endl;
	return true;
}


template<typename Type>   void voxelField<Type>::writeBin(std::string fileName) const
{
	int3 n = size3();

	if(fileName.compare(fileName.size()-4,4,".tif")==0)
	{
	 #ifdef TIFLIB
		std::cout<<  "\n writting tif file "<<fileName<<";  size: "<<n<<" "; std::cout.flush();
		writeTif(*this, fileName);
		std::cout<<  "."<<std::endl;
		return;
	 #else
		fileName = fileName.substr(0,fileName.size()-4)+suffix();
	 #endif //TIFLIB
	}

	#ifdef ZLIB
	if(fileName.compare(fileName.size()-3,3,".gz")==0)
	{
		std::cout<<  "\n writting compressed file "<<fileName<<";  size: "<<n; std::cout.flush();
		gzofstream  of((fileName).c_str());
		of << setcompression(Z_DEFAULT_COMPRESSION);
		assert(of);
		if(data_.size())
		of.write(reinterpret_cast<const char*>(&((*this)(0,0,0))), (size_t(n[0])*n[1]*n[2]) * sizeof(Type));
		of.flush();
		of.close();
		std::cout<<  "."<<std::endl;
		return;
	}
	#else
	if(fileName.compare(fileName.size()-3,3,".gz")==0) fileName=fileName.substr(0,fileName.size()-3);
	#endif


	std::cout<<  " writting binary file "<<fileName<<";  size: "<<n; std::cout.flush();
	auto mode = std::ios::out | std::ios::binary;
	if(fileName.compare(fileName.size()-3,3,".am")==0) //-Warning  remember to call Writeheader before this if suffix is .amam
	{
		char cs[]="xxx";
		std::ifstream is(fileName.c_str());		if(is)	{ is.seekg (3, is.end);	is.get(cs,3); }		is.close();
		if(cs[0]!='@' || cs[1]!='1' || cs[2]!='\n')	writeHeader(fileName,{0,0,0},n);
		mode = mode | std::ios::app;
	}
	std::ofstream of (fileName.c_str(), mode);
	assert(of);
	if(data_.size())
	of.write(reinterpret_cast<const char*>(&((*this)(0,0,0))), (size_t(n[0])*n[1]*n[2]) * sizeof(Type));
	of.flush();
	of.close();
	std::cout<<  "."<<std::endl;

}




template<typename Type>   void voxelField<Type>::writeBin(std::string fileName,
							   int iS,int iE , int jS,int jE , int kS,int kE ) const
{


	if(fileName.compare(fileName.size()-4,4,".tif")==0)
	{
	 #ifdef TIFLIB
		std::cout<<  " writting tif file "<<fileName<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE; std::cout.flush();
		writeTif(*this, fileName, iS, iE,  jS, jE ,  kS, kE);
		std::cout<<  "."<<std::endl;
		return;
	 #else
		fileName = fileName.substr(0,fileName.size()-4)+suffix();
	 #endif //TIFLIB
	}

	#ifdef ZLIB
	if(fileName.compare(fileName.size()-3,3,".gz")==0)
	{
		std::cout<<  " writting compressed file "<<fileName<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE; std::cout.flush();
		gzofstream  of((fileName).c_str());
		of << setcompression(Z_DEFAULT_COMPRESSION);
		assert(of);
		if(data_.size())
		 for ( int k = kS;k < kE;k++)
			for ( int j = jS;j < jE;j++)
			{
				of.write(reinterpret_cast<const char*>(&((*this)(iS,j,k))), (iE-iS) * sizeof(Type));
			}
		of.flush();
		of.close();
		std::cout<<  "."<<std::endl;
		return;
	}
	#endif

	std::cout<<  " writting binary file "<<fileName<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE; std::cout.flush();
	auto mode = std::ios::out | std::ios::binary;
	if(fileName.compare(fileName.size()-3,3,".am")==0)
	{
		char cs[]="xxx";
		std::ifstream is(fileName.c_str());		if(is)	{ is.seekg (3, is.end);	is.get(cs,3); }		is.close();
		if(cs[0]!='@' || cs[1]!='1' || cs[2]!='\n')	writeHeader(fileName,{iS,jS,kS},{iE ,jE,kE});
		mode = mode | std::ios::app;
	}
	std::ofstream of(fileName.c_str(), mode);
	assert(of);
	if(data_.size())
	 for ( int k = kS;k < kE;k++)
		for ( int j = jS;j < jE;j++)
		{
			of.write(reinterpret_cast<const char*>(&((*this)(iS,j,k))), (iE-iS) * sizeof(Type));
		}
	of.flush();
	of.close();
	std::cout<<  "."<<std::endl;

}

template<typename Type>   void voxelField<Type>::writeAscii(std::string fileName,int iS,int iE , int jS,int jE , int kS,int kE) const
{
	std::cout<<  " writting ascii file "<<fileName<<";  "; std::cout.flush();

	std::ofstream of (fileName.c_str());
	assert(of);

	for ( int k = kS;k < kE;k++)
	  for ( int j = jS;j < jE;j++)
	  {
		for ( int i = iS;i < iE;i++)
		   of<<((*this)(i,j,k))<<' ';
		of<<"\n";
	  }
	of<<std::endl;
	of.close();
	std::cout<<  "."<<std::endl;
}

template<>  inline void voxelField<unsigned char>::writeAscii(std::string fileName,int iS,int iE , int jS,int jE , int kS,int kE) const
{
	std::cout<<  " writting ascii file "<<fileName<<";  "; std::cout.flush();

	std::ofstream of (fileName.c_str());
	assert(of);

	for ( int k = kS;k < kE;k++)
	  for ( int j = jS;j < jE;j++)
	  {
		for ( int i = iS;i < iE;i++)
		   of<<int((*this)(i,j,k))<<' ';
		of<<"\n";
	  }
	of<<std::endl;
	of.close();
	std::cout<<  "."<<std::endl;
}

template<typename Type>   void voxelField<Type>::writeAscii(std::string fileName) const
{
	int3 imgsize = size3();
	writeAscii(fileName,0, imgsize[0] ,0,imgsize[1] ,0, imgsize[2]);
}

template<typename Type>   void voxelField<Type>::writeRotatedXZ(std::ofstream& of) const
{
	for (int i = 0;i<(*this).size3()[0];i++) //reversed order with k
	{
		for (int j = 0;j<(*this).size3()[1];j++)
		{
			for (int k = 0;k<(*this).size3()[2];k++) //reversed order with i
			{
				of<<double((*this)(i,j,k))<<' ';
			}
			of<<std::endl;
		}
	}
	std::cout<<  "writeRotatedXZ."<<std::endl;

}




template<typename T>
void voxelImageT<T>::cropD( int3 cropBegin,  int3 cropEnd,int emptylayers, T emptylayersValue)
{
	crop(cropBegin[0],cropEnd[0]-1,cropBegin[1],cropEnd[1]-1,cropBegin[2],cropEnd[2]-1,emptylayers,emptylayersValue);
}


template<typename T>
void voxelImageT<T>::crop( int cropBegin[3],  int cropEnd[3],int emptylayers, T emptylayersValue)
{
	crop(cropBegin[0],cropEnd[0],cropBegin[1],cropEnd[1],cropBegin[2],cropEnd[2],emptylayers,emptylayersValue);
}


template<typename T>
void voxelImageT<T>::crop(
							int iS, int iE ,
							int jS, int jE ,
							int kS, int kE  ,
							int emptylayers, T emptylayersValue
)
{
	(std::cout<<  "  cropping, "<<  "   ["<<iS<<" "<<iE+1 <<  ")  ["<<jS<<" "<<jE+1<< ")  ["<<kS<<" "<<kE+1<<") ").flush();
	if (emptylayersValue) (std::cout<<  ", adding "<<emptylayers<<" layers of "<< double(emptylayersValue)<<"  ").flush();

	X0_[0]=X0_[0]+(int(iS)-int(emptylayers))*dx_[0];   X0_[1]=X0_[1]+(int(jS)-int(emptylayers))*dx_[1];   X0_[2]=X0_[2]+(int(kS)-int(emptylayers))*dx_[2];

	voxelImageT<T> tmp=*this;


	this->reset(iE+1-iS+2*emptylayers,jE+1-jS+2*emptylayers,kE+1-kS+2*emptylayers, emptylayersValue);


	for ( int k=0; k<kE+1-kS; k++ )
		for ( int j=0; j<jE+1-jS; ++j )
			//for ( int i=0; i<iE+1-iS; ++i )
			//{
				//(*this)[k+emptylayers][j+emptylayers][i+emptylayers]=tmp[k+kS][j+jS][i+iS];
					std::copy(&tmp(iS,j+jS,k+kS), &tmp(iE+1,j+jS,k+kS), &(*this)(emptylayers,j+emptylayers,k+emptylayers));

			//}
		//}
	//}


}




template<typename T>
void voxelField<T>::setLayer(int k, const T* Values)
{
	//(*this)[k]=Values;
	std::copy(Values, Values+(*this).size3()[0]*(*this).size3()[1], &(*this)(0,0,k));
	
}

template<typename T>
void voxelField<T>::setSlice(char dir, int ijk, T vv)
{
	if(dir=='i')
	 for ( size_t iii=ijk; iii<(*this).data_.size() ; iii+=(*this).size3()[0] )
		(*this).data_[iii]=vv;
	else if(dir=='j')
	 for ( int k=0; k<(*this).size3()[2] ; k++ )
		std::fill( &(*this)(0,ijk,k), &(*this)(0,ijk+1,k),vv);
	else if(dir=='k')
		std::fill( &(*this)(0,0,ijk), &(*this)(0,0,ijk+1),vv);
	else 	std::cout<<"Error: wrong dir "<<dir<<std::endl;
}

template<typename T>
void voxelField<T>::replacexLayer(int i, int fromi)
{

	 for ( int k=0; k<(*this).size3()[2] ; k++ )
	{
		for ( int j=0; j<(*this).size3()[1] ; ++j )
		{
			//for ( int i=0; i<Values.size3()[0] ; ++i )
			//{
				(*this)(i,j,k)=(*this)(fromi,j,k);
			//}
		}
	}

}
template<typename T>
void voxelField<T>::replaceyLayer(int j, int fromj)
{

	 for ( int k=0; k<(*this).size3()[2] ; k++ )
	{
		//for ( int j=0; j<Values.size3()[1] ; ++j )
		//{
			for ( int i=0; i<(*this).size3()[0] ; ++i )
			{
				(*this)(i,j,k)=(*this)(i,fromj,k);
			}
		//}
	}

}


template<typename T>
void voxelField<T>::setBlock(int n1, int n2, int n3, const voxelField<T>& Values)
{
	forAllkji_(Values)
			(*this)(i+n1,j+n2,k+n3)=Values(i,j,k);
}
template<typename T>
void voxelField<T>::setFrom(const voxelField<T>&Values, int n1, int n2, int n3)
{
	forAllkji_(*this)	(*this)(i,j,k)=Values(i+n1,j+n2,k+n3);
}

template<typename T>
template<typename T2>
void voxelImageT<T>::resetFrom(const voxelImageT<T2>&Values)
{
	dx_= Values.dx();
	X0_ = Values.X0();
	this->reset(Values.sizeu3(),0);
	forAlliii_((*this))
			(*this)(iii)=Values(iii);
}
template<typename T>
void voxelImageT<T>::setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3)
{
	dx_= Values.dx();
	X0_[0] = Values.X0()[0]+n1*dx_[0];
	X0_[1] = Values.X0()[1]+n2*dx_[1];
	X0_[2] = Values.X0()[2]+n3*dx_[2];
	forAllkji_(*this)
			(*this)(i,j,k)=Values(i+n1,j+n2,k+n3);
}


template<typename T>
void voxelImageT<T>::growBox(int nLayers)
{

	int3 n = (*this).sizeu3();
	(*this).crop(0,n[0]-1,0,n[1]-1,0,n[2]-1, nLayers,1);//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	for (int i=0; i<nLayers ; i++ )
	{
		(*this).replaceyLayer(n[1]+nLayers+i, n[1]+nLayers-1);
		(*this).replaceyLayer(i, nLayers);
		(*this).replacexLayer(n[0]+nLayers+i, n[0]+nLayers-1);
		(*this).replacexLayer(i, nLayers);
		(*this).setLayer(n[2]+nLayers+i, &(*this)(0,0,(n[2]+nLayers-1)));
		(*this).setLayer(i, &(*this)(0,0,nLayers));
	}
}




template<typename T>
voxelImageT<T>  resample(const voxelImageT<T>& img, double nReSampleNotSafe)//  TODO to be tested
{
	voxelImageT<T> tmp;
	if (nReSampleNotSafe < .999)
	{
		int nReS=1.0/nReSampleNotSafe+0.5;
		tmp.reset(nReS*img.size3());
		forAllkji_(tmp)
			tmp(i,j,k)=img((0.5+i)/nReS,(0.5+j)/nReS,(0.5+k)/nReS);

		tmp.dxCh()=img.dx()/nReS; //fixed
		tmp.X0Ch()=img.X0()/nReS; //fixed
		return tmp;
	}
	else if (nReSampleNotSafe > 1.001)
	{
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset(img.size3()/nReS);
		forAllkji_(tmp)
		{
			int neiSum=0;
			forAllNei(0,nReS-1)
			{
				neiSum+=nei(img,i*nReS,j*nReS,k*nReS);
			}
			tmp(i,j,k)=(0.5+double(neiSum)/(nReS*nReS*nReS));
		}

		tmp.dxCh()=img.dx()*nReS; //fixed
		tmp.X0Ch()=img.X0()*nReS; //fixed
		return tmp;
	}
	else return img;
}




template<typename T>
voxelImageT<T> resampleMax(const voxelImageT<T>& img, double nReSampleNotSafe)//  TODO to be tested
{
	voxelImageT<T> tmp;
	if (nReSampleNotSafe < .999)
	{
		int nReS=1.0/nReSampleNotSafe+0.5;
		tmp.reset(nReS*img.size3());
		forAllkji_(img)
			tmp(i,j,k)=img((0.5+i)/nReS,(0.5+j)/nReS,(0.5+k)/nReS);

		tmp.dxCh()/=nReS; //fixed
	}
	else if (nReSampleNotSafe > 1.001)
	{
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset(img.size3()/nReS);
		forAllkji_(img)
		{
			T neiSum=std::numeric_limits<T>::min();
			forAllNei(0,nReS-1)
			{
				neiSum=std::max(neiSum, nei(img,i*nReS,j*nReS,k*nReS));
			}
			tmp(i,j,k)=neiSum;//(0.5+double(neiSum)/(nReS*nReS*nReS));
		}

		tmp.dxCh()*=nReS;
	}
	return tmp;
}


template<typename T>
void voxelImageT<T>::rotate(char direction)
{// wrong X0
	int n1,n2,n3;

	(std::cout<<" x<->"<<direction<<" ").flush();
	voxelField<T>::getSize(n1,n2,n3);
	if (direction=='z' || direction=='Z')
	{
		//int nMinTmp=nMin_[0];
		//nMin_[0]=nMin_[2];
		//nMin_[2]=nMinTmp;
		{
			double X0Tmp=X0_[0];
			X0_[0]=X0_[2];
			X0_[2]=X0Tmp;
			double dxTmp=dx_[0];
			dx_[0]=dx_[2];
			dx_[2]=dxTmp;
		}
		voxelImageT<T> tmp=*this;
		this->reset(n3,n2,n1,0);
		size_t nij =this->nij_;
		OMPragma("omp parallel for")
		for ( int k=0; k<n3 ; k++ )
			for ( int j=0; j<n2 ; ++j )
			{
				//for (int i=0; i<n1; ++i) (*this)(k,j,i)=tmp(i,j,k);
				for ( int i=1; i<n1 ; i+=2 )
				{
					const T& vv0 = tmp(i,j,k);
					const T  vv1 = *(&vv0-1);
					*(&((*this)(k,j,i)=vv0)-nij)=vv1;
				}
				if(n1%2) (*this)(k,j,n1-1)=tmp(n1-1,j,k);
			}
	}
	else if (direction=='y' || direction=='Y')
	{
		//int nMinTmp=nMin_[0];
		//nMin_[0]=nMin_[1];
		//nMin_[1]=nMinTmp;
		{
			double X0Tmp=X0_[0];
			X0_[0]=X0_[1];
			X0_[1]=X0Tmp;
			double dxTmp=dx_[0];
			dx_[0]=dx_[1];
			dx_[1]=dxTmp;
		}

		voxelImageT<T> tmp=*this;
		this->reset(n2,n1,n3,0);
		for ( int k=0; k<n3 ; k++ )
			for ( int j=0; j<n2 ; ++j )
			{	//for (int i=1; i<n1; ++i) (*this)(j,i,k)=tmp(i,j,k);
				for ( int i=1; i<n1 ; i+=2 )
				{
					const T& vv0 = tmp(i,j,k);
					const T  vv1 = *(&vv0-1);
					*(&((*this)(j,i,k)=vv0)-(*this).nnn_[0])=vv1;
				}
				if(n1%2) (*this)(j,n1-1,k)=tmp(n1-1,j,k);

			}
				
	}
	else if (direction=='-')
	{
		std::cout<<" -> flipping image,  x origin will be invalid "<<std::endl;
		voxelImageT<T> tmp=*this;
		for ( int k=0; k<n3 ; k++ )
			for ( int j=0; j<n2 ; ++j )
				for ( int i=0; i<n1 ; ++i )
					(*this)(n1-1-i, j, k)=tmp(i,j,k);
	}
	else
	{
		std::cout<<"\n\nSwapping "<<direction<<" and x directions(!?!), impossibile  >-( "<<std::endl;
		std::cerr<<"Swapping "<<direction<<" and x directions(!?!), impossibile >-( \n\n"<<std::endl;
	}

}


template<typename T>
void voxelImageT<T>::PointMedian026(int thereshold0,int thereshold1)
{
	unsigned long nChanged(0);
	dbl3 doubletmp(0,0,0);
	int3 n = voxelField<T>::size3();
	for (int i=0;i<3;i++) n[i]=n[i]+2;
	
	voxelImageT<T> voxls(n,doubletmp,doubletmp,1);
		voxls.setBlock(0, 0, 0, (*this));
		voxls.setBlock(2, 2, 2, (*this));
		voxls.setBlock(1, 1, 1, (*this));


		forAllkji_1_(voxls)
		{
			int neiSum=0;
			forAllNei(-1,1)
			{
				neiSum+=nei(voxls,i,j,k);
			}
			neiSum-=voxls(i,j,k);
			if (neiSum <= thereshold0  && (*this)(i-1,j-1,k-1))
			{
				(*this)(i-1,j-1,k-1)=0;
				++nChanged;
			}
			else if (neiSum >= thereshold1  && !((*this)(i-1,j-1,k-1)))
			{
				(*this)(i-1,j-1,k-1)=1;
				++nChanged;
			}
		}

	std::cout<<"PointMedian026  changed: "<<nChanged<<std::endl;

}





//~
//#define forAllFaceNei \/
//for (int k_nei_m=-1;k_nei_m<2;k_nei_m++) \/
//for (int j_nei_m=-1;j_nei_m<2;j_nei_m++) \/
//for (int i_nei_m=-1;i_nei_m<2;i_nei_m++)
//~
//#define nei(i,j,k) (*this)[k+k_nei_m][j+j_nei_m][i+i_nei_m]


template<typename T>
void voxelImageT<T>::FaceMedian06(int thereshold0,int thereshold1)
{
	unsigned long nChanged(0);
	dbl3 doubletmp(0,0,0);
	int3 n = voxelField<T>::size3();
	for (int i=0;i<3;i++) n[i]=n[i]+2;

	voxelImageT<T> voxls=*this;
	voxls.growBox(1);

	OMPragma("omp parallel for reduction(+:nChanged)") 		
	for ( int k=1; k<(voxls).size3()[2]-1 ; ++k )   
	for ( int j=1; j<(voxls).size3()[1]-1 ; ++j )   
	for ( int i=1; i<(voxls).size3()[0]-1 ; ++i )
	{		const T* vp=&voxls(i,j,k);
			int neiSum = sum6Nei(voxls, vp);

			if (neiSum <= thereshold0 && (*this)(i-1,j-1,k-1))
			{
				(*this)(i-1,j-1,k-1)=0;
				++nChanged;
			}
			else if (neiSum >= thereshold1 && !((*this)(i-1,j-1,k-1)))
			{
				(*this)(i-1,j-1,k-1)=1;
				++nChanged;
			}
		}

	std::cout<<"FaceMedian06  changed: "<<nChanged<<std::endl;
}


template<typename T>
void voxelImageT<T>::shrinkPore()
{
	voxelImageT<T> voxls=*this;


	forAllkji_1_(voxls)
	{
		T* vp = &voxls(i,j,k);
		if (*vp==0 && ( sum6Nei(voxls, vp) ) )
			(*this)(i,j,k)=1;
	}



	OMPragma("omp parallel for") 		
	for ( int k=1; k<voxls.size3()[2]-1 ; k++ )
	{
		for ( int j=1; j<voxls.size3()[1]-1 ; ++j )
		{
			//for ( int i=0; i<voxls.size3()[0]-1 ; ++i )
			{	int i=0;
				if (voxls(i,j,k)==0 && ( voxls(i,j+1,k) || voxls(i,j-1,k) || voxls(i,j,k+1) || voxls(i,j,k-1) ) )
					(*this)(i,j,k)=1;

				i=voxls.size3()[0]-1;
				if (voxls(i,j,k)==0 && ( voxls(i,j+1,k) || voxls(i,j-1,k) || voxls(i,j,k+1) || voxls(i,j,k-1) ) )
					(*this)(i,j,k)=1;
		   }
		}
	}


	OMPragma("omp parallel for") 		
	for ( int k=1; k<voxls.size3()[2]-1 ; k++ )
	{
		//for ( int j=0; j<voxls.size3()[1]-1 ; ++j )
		{
			for ( int i=1; i<voxls.size3()[0]-1 ; ++i )
			{
				int j=0;
				if (voxls(i,j,k)==0 && ( voxls(i-1,j,k) || voxls(i+1,j,k) || voxls(i,j,k-1) || voxls(i,j,k+1) ) )
					(*this)(i,j,k)=1;

				j=voxls.size3()[1]-1;
				if (voxls(i,j,k)==0 && ( voxls(i-1,j,k) || voxls(i+1,j,k) || voxls(i,j,k-1) || voxls(i,j,k+1) ) )
					(*this)(i,j,k)=1;
		   }
		}
	}

	//for ( int k=0; k<voxls.size3()[2]-1 ; k++ )
	{
		OMPragma("omp parallel for") 		
		for ( int j=1; j<voxls.size3()[1]-1 ; ++j )
		{
			for ( int i=1; i<voxls.size3()[0]-1 ; ++i )
			{
				int k=0;
				if (voxls(i,j,k)==0 && ( voxls(i-1,j,k) || voxls(i+1,j,k) || voxls(i,j-1,k) || voxls(i,j+1,k) ) )
					(*this)(i,j,k)=1;

				k=voxls.size3()[2]-1;
				if (voxls(i,j,k)==0 && ( voxls(i-1,j,k) || voxls(i+1,j,k) || voxls(i,j-1,k) || voxls(i,j+1,k) ) )
					(*this)(i,j,k)=1;
		   }
		}
	}


}

template<typename T>
class mapComparer  {  public: bool operator() (std::pair<const T,short>& i1, std::pair<const T,short> i2) {return i1.second<i2.second;}  };

template<typename T>
void voxelImageT<T>::mode(short minDif)
{
	//long long nChanges = 0;
	voxelImageT<T> voxls=*this;
   forAllkji_1_(voxls)
	{
		register T* vp = &((*this)(i,j,k));
		register const T pID = *vp;

		short nSames(0);
		std::map<T,short> neis;///.  ID-counter

		T
		neiPID = voxls.v_i(-1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls.v_i( 1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls.v_j(-1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls.v_j( 1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls.v_k(-1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;
		neiPID = voxls.v_k(1,vp);
		if (neiPID != pID  ) 	 ++(neis.insert(std::pair<T,short>(neiPID,0)).first->second); else ++nSames;

		if(nSames+minDif/2<3)
		{
			 auto neitr = std::max_element(neis.begin(), neis.end(), mapComparer<T>());
			 if (neitr->second >= nSames+minDif)
			{
				//++nChanges; 
				(*this)(i,j,k) = neitr->first;
			}
		}
	}
	//( std::cout<<"  nMedian: "<< std::left<<nChanges<<"  \n").flush();
}







template<typename T>
void voxelImageT<T>::growPore() // optimized function, should be further optimized as it is frequently used
{

	voxelImageT<T> voxls=*this;


	forAllkji_1_(voxls)
	{

		if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i+1,j,k) ||
							  !voxls(i,j-1,k) || !voxls(i,j,k+1) || !voxls(i,j+1,k) || !voxls(i,j,k-1) ) )
			(*this)(i,j,k)=0;
	}


	OMPragma("omp parallel for") 		
	for ( int k=1; k<voxls.size3()[2]-1 ; k++ )
	{
		for ( int j=1; j<voxls.size3()[1]-1 ; ++j )
		{
			// for ( int i=0; i<voxls.size3()[0]-1 ; ++i )
			{	int i=0;
				if (voxls(i,j,k) && ( !voxls(i+1,j,k) || !voxls(i,j-1,k)  || !voxls(i,j,k+1) || !voxls(i,j+1,k) || !voxls(i,j,k-1) ) )
					(*this)(i,j,k)=0;

				i=voxls.size3()[0]-1;
				if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i,j-1,k)  || !voxls(i,j,k+1) || !voxls(i,j+1,k) || !voxls(i,j,k-1) ) )
					(*this)(i,j,k)=0;
		   }
		}
	}


	OMPragma("omp parallel for") 		
	for ( int k=1; k<voxls.size3()[2]-1 ; k++ )
	{
		// for ( int j=0; j<voxls.size3()[1]-1 ; ++j )
		{
			for ( int i=1; i<voxls.size3()[0]-1 ; ++i )
			{	int j=0;

				if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i+1,j,k)  || !voxls(i,j,k+1) || !voxls(i,j+1,k) || !voxls(i,j,k-1) ) )
					(*this)(i,j,k)=0;

				j=voxls.size3()[1]-1;
				if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i+1,j,k) || !voxls(i,j-1,k)  || !voxls(i,j,k+1) || !voxls(i,j,k-1) ) )
					(*this)(i,j,k)=0;
		   }
		}
	}

	// for ( int k=0; k<voxls.size3()[2]-1 ; k++ )
	{
		OMPragma("omp parallel for") 		
		for ( int j=1; j<voxls.size3()[1]-1 ; ++j )
		{
			for ( int i=1; i<voxls.size3()[0]-1 ; ++i )
			{	int k=0;

				if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i+1,j,k) || !voxls(i,j,k+1) || !voxls(i,j+1,k) || !voxls(i,j-1,k) ) )
					(*this)(i,j,k)=0;

				k=voxls.size3()[2]-1;
				if (voxls(i,j,k) && ( !voxls(i-1,j,k) || !voxls(i+1,j,k) || !voxls(i,j-1,k) || !voxls(i,j+1,k) || !voxls(i,j,k-1) ) )
					(*this)(i,j,k)=0;
		   }
		}
	}

}








template<typename T>
void voxelImageT<T>::NOT(const voxelImageT& data2)
{
	forAlliii_((*this))	(*this)(iii)= (*this)(iii) && !data2(iii);
}
template<typename T>
void voxelImageT<T>::AND(const voxelImageT& data2)
{
	forAlliii_((*this))
				(*this)(iii)= (*this)(iii) && data2(iii);
}
template<typename T>
void voxelImageT<T>::OR(const voxelImageT& data2)
{
	forAlliii_((*this))
				(*this)(iii)= (*this)(iii) || data2(iii);
}

template<typename T>
void voxelImageT<T>::XOR(const voxelImageT& data2)
{
	forAlliii_((*this))
		(*this)(iii)= (*this)(iii) != data2(iii);
}

template<typename T>
void voxelImageT<T>::maxEq(const voxelImageT& data2)
{
	forAlliii_((*this))
				(*this)(iii)= max((*this)(iii), data2(iii));
}

template<typename T>
void voxelImageT<T>::minEq(const voxelImageT& data2)
{
	forAlliii_((*this))
				(*this)(iii)= min((*this)(iii), data2(iii));
}


template<typename T>
void voxelImageT<T>::threshold101(T theresholdMin,T  theresholdMax)
{
	forAllvp_((*this))
	{	T vv = *vp;
		*vp=	( vv < theresholdMin )  || ( vv > theresholdMax  );
	}

}

template<typename T>
void rescale(voxelImageT<T>& img, T theresholdMin,T  theresholdMax)
{
	T vmin = std::numeric_limits<T>::max();
	T vmax = std::numeric_limits<T>::min();
	int deltaT = theresholdMax - theresholdMin;

	OMPragma("omp parallel for reduction(min:vmin) reduction(max:vmax)")	
	forAllcp(img) {vmin = std::min(*cp,vmin); vmax = std::max(*cp,vmin); }
	std::cout<<"   vmin:"<<int(vmin)<<"   vmax:"<<int(vmax)<<"  ";
	vmax = std::max(T(vmin+1),vmax);
	forAllvp_(img)
	{	*vp = theresholdMin + (deltaT*(*vp - vmin))/(vmax-vmin);
	}

}




template<typename T>
void voxelImageT<T>::fillHoles(int maxHoleRadius)
{
	std::cout<<"  filling small isolated parts: "<<std::flush;
	voxelImageT<T> dataTmp=*this;
		std::cout<<"-"<<std::flush;

	dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();
	for ( int i=0 ; i < 6 ; ++i )
		{ dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush();}
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	dataTmp.growPore();
	for ( int i=0 ; i < 4 ; ++i )
		{ dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	if ( maxHoleRadius > 1)
	{
		for ( int i=0 ; i < maxHoleRadius ; ++i )
			{ dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();}
		for ( int i=0 ; i < maxHoleRadius*4 ; ++i )
			{ dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;

		for ( int i=0 ; i < maxHoleRadius ; ++i )
			{ dataTmp.growPore(); std::cout<<".";std::cout.flush();}
		for ( int i=0 ; i < maxHoleRadius*3 ; ++i )
			{ dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;
	}
	std::cout<<"."<<std::endl;

}

template<typename T>
void voxelImageT<T>::writeAConnectedPoreVoxel(std::string fileName) const
{
	std::cout<<" finding a connected pore voxel:";
	for ( int nShrink=4; nShrink>=0  ; nShrink-- )
	{
		voxelImageT<T> dataTmp=*this;
		for ( int iShrink=1; iShrink<=nShrink ; iShrink++ )
		{
			dataTmp.shrinkPore();
		}

		std::cout<<" nShrink: "<<nShrink<<std::endl;;
		dataTmp.printInfo();

		for (int k=dataTmp.size3()[2]*3/4-2; k>dataTmp.size3()[2]*1/8+1  ; k-- )
		 for (int j=dataTmp.size3()[1]*3/4-2; j>dataTmp.size3()[1]*1/8+1  ; j-- )
			for (int i=dataTmp.size3()[0]*3/4-2; i>dataTmp.size3()[0]*1/8+1  ; i-- )
			{
				if (dataTmp(i,j,k)==0 && 
					 !(   dataTmp(i+1,j+1,k+1) || dataTmp(i-1,j+1,k+1) || dataTmp(i-1,j-1,k+1) || dataTmp(i+1,j-1,k+1)
						|| dataTmp(i-1,j-1,k-1) || dataTmp(i+1,j-1,k-1) || dataTmp(i-1,j+1,k-1) || dataTmp(i+1,j+1,k-1)
						|| dataTmp(i-1,j,k) || dataTmp(i+1,j,k)
						|| dataTmp(i,j-1,k) || dataTmp(i,j+1,k)
						|| dataTmp(i,j,k-1) || dataTmp(i,j,k+1)
					  )
					)
				{
					std::ofstream of(fileName.c_str());
					assert(of);
					of<<(i+0.5)*dx_[0]+X0_[0]<<" "<<(j+0.5)*dx_[1]+X0_[1]<<" "<<(k+0.5)*dx_[2]+X0_[2]<<std::endl;
					of.close();
					std::cout<<" found  ("<<i<<" "<<j<<" "<<k<<") -> "<<double(dataTmp(i,j,k))<<std::endl;
					std::cout<<" found  xyz:  "<<i*dx_[0]+X0_[0]<<" "<<j*dx_[1]+X0_[1]<<" "<<k*dx_[2]+X0_[2]<<std::endl;
					return ;
				}
			}
	}
	std::cout<<" \n			  ----  ERRORR   -----	  \n"<<std::endl;
	std::cout<<" \n\n ----  didn't find  a connected pore voxel   -----  \n\n"<<std::endl;
}

template<typename T>
void replaceRange(voxelImageT<T>& vImage, int minvi, int  maxvi, int midvi)
{
	(std::cout<<minvi<<":"<<maxvi<<"->"<<midvi<<"    ").flush();
	const T minv=minvi,   maxv=maxvi,  midv=midvi;
	forAllvp_(vImage)
	{	T vv = *vp;
		if (minv<=vv && vv<=maxv)
			*vp=midv;
	}
}


template<typename T>
void voxelImageT<T>::printInfo() const
{
	int3 n=(*this).size3();
	if(std::is_integral<T>::value)
	{
		unsigned long long nPores=0;
		unsigned long long nTotal=0;
		std::cout<<"  Calculating image porosity: (void==0, valid!=255):"<<std::endl;
		OMPragma("omp parallel for reduction(+:nPores) reduction(+:nTotal)")
		forAllcp((*this))  { nPores += (*cp==0); nTotal += (*cp!=255);}

		std::cout << "   total porosity: " << nPores<<"/ ("<<n[0]<<"*"<<n[1]<<"*"<<n[2]<<") = "<< double(nPores)/(double(n[0])*n[1]*n[2]) << std::endl;
		std::cout << "   validPorosity: " << nPores<<"/"<<nTotal<<" = "<< double(nPores)/double(nTotal) << std::endl;
		std::cout << "   dx: " << dx_<<",  X0: "<< X0_ << std::endl;
	}
	else
	{
		double minv=1.0e64, maxv=-1.0e64, avgv=0.0;
		OMPragma("omp parallel for reduction(min:minv) reduction(max:maxv) reduction(+:avgv)")
		forAllcp((*this))  { double val = toScalar(*cp); minv=std::min(minv,val); maxv=std::max(maxv,val); avgv+=val;
			
				if(!std::isfinite(avgv))
				{
					std::cout<<"overflow or ???";
					write("dump_overflow_dump.dat");
					exit(-1);
				}
			 }
		std::cout << " min: "<<minv  << " max: "<<maxv  << " avg: "<<avgv/(double(n[0])*n[1]*n[2])<< std::endl;
	}

}

template<typename T>
double voxelImageT<T>::volFraction(T vv1,T vv2) const
{
	unsigned long long nPores=0;
	int3 n=(*this).size3();
	OMPragma("omp parallel for reduction(+:nPores)")
	forAllcp((*this))	nPores += ( vv1<=*cp && *cp<=vv2 );

	return double(nPores)/(double(n[0])*n[1]*n[2]);

}

template<typename T>
void voxelImageT<T>::writeHeader(std::string outputName) const
{
	voxelField<T>::writeHeader(outputName, int3(0,0,0), voxelField<T>::sizeu3(),dx_,X0_);
}

template<typename T>
void voxelField<T>::writeHeader(std::string outName, int3 iS, int3 iE, dbl3 dx, dbl3 X0) const
{

	if(dx[0]<0)
	{
		std::cerr<<"Error negative dx, writing abs value instead" ;
		dx[0]=std::abs(dx[0]);  dx[1]=std::abs(dx[1]);  dx[2]=std::abs(dx[2]);
	}
	dbl3 xmin(X0+iS*dx);
	dbl3 xmax(X0+iE*dx);
	int3	n(iE-iS);
	if (outName.compare(outName.size()-3,3,".am")==0)
	{
		std::string typ="uchar";
		if (typeid(T)==typeid(char)) typ="char";
		else if (typeid(T)==typeid(short)) typ="short";
		else if (typeid(T)==typeid(unsigned short)) typ="ushort";
		else if (typeid(T)==typeid(int)) typ="int";
		else if (typeid(T)==typeid(int)) typ="uint";
		else if (typeid(T)==typeid(float)) typ="float";
		else if (typeid(T)==typeid(double)) typ="double";
		else if (typeid(T)==typeid(float3)) typ="float[3]";

		std::ofstream out(outName.c_str());
		assert(out);
		out	<<"# Avizo BINARY-LITTLE-ENDIAN 2.1\n\n\n"
			<<"define Lattice "<<n[0]<<" "<<n[1]<<" "<<n[2]<<"\n\n"
			<<"Parameters {\n    Units {\n        Coordinates \"m\"\n    }\n";
			if(typ=="float[3]")
				out <<"    XLabExperiment {\n        viscosity 0.001,\n        inputPressure 1,\n        outputPressure 0,\n        flowRate 1111111\n    }\n";
			
		out	<<"    Content \""<<n[0]<<"x"<<n[1]<<"x"<<n[2]<<" "<<typ<<", uniform coordinates\",\n"
			<<"    BoundingBox "<<xmin[0]<<" "<<xmax[0]<<" "<<xmin[1]<<" "<<xmax[1]<<" "<<xmin[2]<<" "<<xmax[2]<<",\n"
			<<"    CoordType \"uniform\"\n}\n\n"
			<<"Lattice { "<<typ<<" Data } @1\n\n# Data section follows\n@1\n";

	}else
	if (outName.size()<7 || outName.compare(outName.size()-7,7,"_header")!=0)
	{
		int islash=outName.find_last_of("\\/"); if (islash>=int(outName.size())) islash=-1;
		std::string title=outName.substr(islash+1);
		if (outName.size()>=4 && outName.compare(outName.size()-4,4,".mhd")==0)
			title=title.substr(0,title.size()-4)+suffix();
		else if (outName.size()>=7 && outName.compare(outName.size()-7,7,".raw.gz")==0)
			outName=outName.substr(0,outName.size()-7)+".mhd";
		else
			outName=outName.substr(0,outName.find_last_of("."))+".mhd";


		std::string typeNmeVTK="MET_UCHAR";
		if (typeid(T)==typeid(char)) typeNmeVTK="MET_CHAR";
		else if (typeid(T)==typeid(short)) typeNmeVTK="MET_SHORT";
		else if (typeid(T)==typeid(unsigned short)) typeNmeVTK="MET_USHORT";
		else if (typeid(T)==typeid(int)) typeNmeVTK="MET_INT";
		else if (typeid(T)==typeid(int)) typeNmeVTK="MET_UINT";
		else if (typeid(T)==typeid(float)) typeNmeVTK="MET_FLOAT";
		else if (typeid(T)==typeid(double)) typeNmeVTK="MET_DOUBLE";

		std::ofstream outputHeaderFile(outName.c_str());
		assert(outputHeaderFile);
		outputHeaderFile
			 <<"ObjectType =  Image"<<std::endl
			 <<"NDims =	   3"<<std::endl
			 <<"ElementType = "<<typeNmeVTK <<std::endl
			 <<"ElementByteOrderMSB = False\n"//(typeNmeVTK=="MET_UCHAR" ? "\n":)
			 <<"ElementNumberOfChannels = 1\n"//(typeNmeVTK=="MET_UCHAR" ? "\n":)
			 <<"CompressedData = "<<((title.size()>=3 && title.compare(title.size()-3,3,".gz")==0) ? "True" : "False")<<std::endl
			 <<"\nDimSize =		"<<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl
			 <<"ElementSize =	"<<dx[0]<<" "<<"  " <<dx[1]<<" "<<"  " <<dx[2]<<std::endl
			 <<"Offset =		 "<<X0[0]<<" "<<"  " <<X0[1]<<" "<<"  " <<X0[2]<<std::endl
			 <<"ElementDataFile = "<<title<<std::endl <<std::endl;
			 if(dx[0]>=0.001)
				outputHeaderFile <<"Unit = "<<1<<std::endl;

			outputHeaderFile <<std::endl <<std::endl;
	}
	else
	{
	  std::ofstream outputHeaderFile(outName.c_str());					 // file pointer
	  assert(outputHeaderFile);
	  outputHeaderFile
		 <<"Nxyz"<<std::endl
		 <<"dxX0"<<std::endl
		 <<n[0]<<" "<<n[1]<<" "<<n[2]<<std::endl
		 <<dx[0]<<" "<<"  " <<dx[1]<<" "<<"  " <<dx[2]<<std::endl
		 <<xmin[0]<<" "<<"  " <<xmin[1]<<" "<<"  " <<xmin[2]<<std::endl
		 <<"\n\nComments:"<<std::endl
		 <<" first 9 entries above are:"<<std::endl
		 <<"	Nx Ny Nz"<<std::endl
		 <<"	dx dy dz"<<std::endl
		 <<"	Xo Yo Zo"<<std::endl
		 <<" Nx, Ny and Nz  count for the number of columns, rows and layers respectively as written in the file"<<std::endl
		 <<" Optional keywords (move above Comments to activate):"<<std::endl
		 <<"	crop		0  299   0  299   0  299 "<<std::endl
		 <<"	pore 		0 0 "<<std::endl
		 <<"	resample	1"<<std::endl
		 <<"	direction	z"<<std::endl
		 <<"	..... "<<std::endl
		 <<std::endl;
	}
}


template<typename T>
void voxelField<T>::writeNoHdr(std::string outName) const
{
	if (outName.compare(outName.size()-4,4,".mhd") == 0 )
		outName = outName.substr(0,outName.size()-4)+suffix();

	if (outName.compare(outName.size()-4,4,".dat") == 0 || outName.compare(outName.size()-4,4,".txt") == 0)
	{
		this->writeAscii(outName);
	}
	else if(outName!="NO_WRITE")
	{
		this->writeBin(outName);
	}

}

template<typename T>
void voxelImageT<T>::writeNoHdr(std::string outName) const
{
	if(outName.compare(outName.size()-3,3,".am")==0) writeHeader(outName);// Warning: remember to append to this in writeBin
	voxelField<T>::writeBin(outName);
}

template<typename T>
void voxelImageT<T>::write(std::string outName) const
{

	if (outName.compare(outName.size()-4,4,".dat") == 0 || outName.compare(outName.size()-4,4,".txt") == 0)
	{
		this->writeAscii(outName);
		writeHeader(outName+"_header");
	}
	else if(outName!="NO_WRITE")
	{
		if (outName.compare(outName.size()-4,4,".mhd") == 0)
		{	outName = outName.substr(0,outName.size()-4)+suffix();
			writeHeader(outName);
		}else if(outName.compare(outName.size()-4,4,".tif")!=0 )//|| outName.compare(outName.size()-3,3,".am")!=0
			writeHeader(outName);
		this->writeBin(outName);
	}
}

template<typename T>
voxelImageT<T> median(const voxelImageT<T>& vImage)
{
	//unsigned long nChanged(0);
	(std::cout<<"  median ").flush();
	voxelImageT<T> voxls=vImage;
	forAllkji_1_(vImage)
	{  const T* vp=&vImage(i,j,k);
		std::array<T,7> vvs={{ *vp,
								vImage.v_i(-1,vp), vImage.v_i( 1,vp),
								vImage.v_j(-1,vp), vImage.v_j( 1,vp),
								vImage.v_k(-1,vp), vImage.v_k( 1,vp)
								}};

		std::nth_element(vvs.begin(),vvs.begin()+3,vvs.end());
		//nChanged+=voxls(i,j,k) != vvs[3];	
		voxls(i,j,k) = vvs[3];
	}

	//(std::cout<<nChanged<<", ").flush();
	return voxls;
}


template<typename T>
void circleOut(voxelImageT<T>& vImage, int X0,int Y0,int R, char dir = 'z', T outVal=std::numeric_limits<T>::max())//  TODO to be tested
{
	int rlim = R*R;
	if (dir=='z')
	{
	 forAllkji_(vImage)
		  if((i-X0)*(i-X0)+(j-Y0)*(j-Y0)>rlim)
		    vImage(i,j,k)=outVal;
	}
	else if (dir=='x')
	{
	 for ( int k=0; k<vImage.size3()[2]; k++ )
	  for ( int j=0; j<vImage.size3()[1]; ++j )
		if((j-X0)*(j-X0)+(k-Y0)*(k-Y0)>rlim)
	      for ( int i=0; i<vImage.size3()[0]; ++i )
		    vImage(i,j,k)=outVal;
	} else std::cout<<"Error: bad direction "<<dir<<std::endl;
}



template<typename T>
void maskWriteFraction(voxelImageT<T>& vImage, std::string maskname, std::string outName, unsigned char maskvv, T minIelm, T maxIelm)//  TODO to be tested
{
	voxelImageT<unsigned char> mask(maskname);
	T maxvv = std::min(maxIelm, accumulate(vImage,(std::max<T>)));//(const T& (*)(const T&, const T&))(std::max<T>)
	std::cout<<"  maxvv:"<<maxvv<<std::endl;
	std::vector<int> nMasked(maxvv+3,0);
	std::vector<int> nNotmsk(maxvv+3,0);
	
	forAllkji_(vImage)
	{	T vv=vImage(i,j,k);
		if(minIelm<=vv && vv<=maxIelm)
		{
			if(mask(i,j,k)==maskvv)		++nMasked[vv];
			else                       ++nNotmsk[vv];
		}
	}
	std::cout<<" Mask Info:"<<std::endl;
	mask.printInfo();
	//mask.write("dumpMask.mhd");
	std::ofstream outf(outName);
	if(outf) 
	{
		std::cout<<"  Writting "<<outName<<std::endl;
		for(T i=minIelm; i<=maxvv ;++i) outf<<double(nMasked[i])/(nMasked[i]+nNotmsk[i]+1.0e-38)<<std::endl;
	}
	else std::cout<<"  Can not open file "<<outName<<" for writting"<<std::endl;
}

