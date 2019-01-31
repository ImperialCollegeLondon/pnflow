/**************************************************************************************


pnextract algorithm version 0.3, March. 2017

For more information, please visit Imperial College pore-scale modelling website:
http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/
or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk
***************************************************************************************/

#include <array>



#ifndef RELEASE_DATE
 #define RELEASE_DATE  __DATE__
#endif


#include "blockNet.h"
#include "writers.h"

void GNMAddCorrelations(string netName, dbl2 minRtForDNS, dbl2 scaleCrlByDNS);




#ifndef MAIN

inline void usage(int detailed=1)
{

	std::cout<<"Pore Network Extraction: pnextract version " << RELEASE_DATE << std::endl;
	if(detailed)
	{
		cout<<"Usage:\n  pnextract input_image.mhd\n";
		cout<<"Usage:\n  pnextract input_image.mhd";

		std::cout<<" For more information, please visit Imperial College pore-scale modelling website:"<<std::endl
			 <<"http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling"<<std::endl
			 <<"or contact Ali Q. Raeini by email: a.qaseminejad-raeini@imperial.ac.uk"<<endl;
		std::cout<<"\nUsage:"<< std::endl;
		std::cout<<"  pnextract vxlImage.mhd\n"<< std::endl;
	}
}

int debugLevel=0;
int main(int argc, char* argv[])
{

	usage(0);

	string inputName;
	if(argc>1)  inputName = argv[1];
	else  {     cout << "Please input data file: ";    getline(cin, inputName);    }
	if (inputName.empty()) { usage(); inputName = "vxlImage.mhd"; }
	if (inputName=="-h") { usage(); exit(0);}
	if (inputName=="-g") return generate_input_nextract(string(argc>2? argv[2] : ""), inputName);
	srand(1001);
	try
	{

		inputDataNE cfg(inputName);
			nextract(cfg,true);

	}
	catch (exception &exc) {  cerr<<"\n\n Error in pnextract: \n   "<<exc.what()<<"\n Aborting! \n" <<endl;	return 1; }
	catch (...)            {  cerr<<"\n\n Unknown Error in pnextract! \n Aborting! \n"<<endl;	return 1; }

	return 0;
}

#endif

int nextract(inputDataNE& cfg, bool verbose)
{

		if (cfg.getOr(false,"write_all")) // use `write_all` as a rememberable alternative for all other visualization keywords 
		{
			cfg.addKeyword("write_radius","true");
			cfg.addKeyword("write_statistics","true");
			cfg.addKeyword("write_elements","true");
				 cfg.addKeyword("write_poreMaxBalls","true");
				 cfg.addKeyword("write_throatMaxBalls","true");
				 cfg.addKeyword("write_throats","true");
			cfg.addKeyword("write_poroats","true");
			cfg.addKeyword("write_hierarchy","true");
			cfg.addKeyword("write_medialSurface","true");
			cfg.addKeyword("write_throatHierarchy","true");
			cfg.addKeyword("write_vtkNetwork","true");
			cfg.addKeyword("write_StatoilFormat","true");
		}
		if (cfg.getOr(false,"write_Statoil")) cfg.addKeyword("write_StatoilFormat","true");



	cfg.init(verbose);
	cfg.readImage();
	cfg.createSegments();


	medialSurface* refs;
	blockNetwork mpn(refs,cfg);
	mpn.createMedialSurface(refs,cfg,0);

			if (cfg.getOr(false,"write_radius"))	ballRadiiToVoxel(mpn).writeBin(cfg.name()+"_radius"+cfg.imgfrmt);

	mpn.CreateVElem(0);

			if (cfg.getOr(true,"write_elements"))	mpn.VElems.write(cfg.name()+"_VElems.mhd");

	mpn.createNewThroats(refs);



		mpn.writeStatoilFormat();

			if (cfg.getOr(false,"write_hierarchy"))
																vtuWriteMbMbs(cfg.name()+"_mbHierarchy"+toStr(0), mpn.allSpace,  mpn.poreIs,  mpn.VElems, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_medialSurface"))
																 vtuWriteMedialSurface(cfg.name()+"_medialSurface"+toStr(0), mpn.allSpace,  mpn.poreIs,  mpn.VElems, cfg.precision,mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_poroats"))         VElemsPlusThroats(mpn).writeBin(cfg.name()+"_poroats"+cfg.imgfrmt,1,cfg.nx,1,cfg.ny,1,cfg.nz);
			if (cfg.getOr(false,"write_throatHierarchy")) vtuWriteThroatMbMbs(cfg.name()+"_throatHierarchy", mpn.throatIs,  mpn.poreIs,  mpn.VElems, cfg.precision,mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_vtkNetwork"))      vtuWritePores(cfg.name()+"_pores",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());
			if (cfg.getOr(false,"write_vtkNetwork"))      vtuWriteTHroatSpheres(cfg.name()+"_throatsBalls",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());



		int outputBlockSize = 0; //. developed by Tom Bultreys
		if(cfg.getVar(outputBlockSize,"outputBlockSize"))
			cout << "OutputBlockSize:" << outputBlockSize << endl;

        if (!outputBlockSize) {
            if (cfg.getOr(false,"write_throats"))                  VThroats(mpn).writeBin(cfg.name()+"_throats"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr(false,"write_poreMaxBalls"))             poreMaxBalls(mpn).writeBin(cfg.name()+"_poreMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
            if (cfg.getOr(false,"write_throatMaxBalls"))           throatMaxBalls(mpn).writeBin(cfg.name()+"_throatMBs"+cfg.imgfrmt,0,cfg.nx,0,cfg.ny,0,cfg.nz);
        } else {

            int blockNumber = 1;
            int beginSlice = 0;
            int endSlice = 0;

            while (endSlice < cfg.nz-1){
                cout << " WRITING BLOCK \n";
                beginSlice = (blockNumber-1) * outputBlockSize;
                endSlice = min(blockNumber * outputBlockSize, cfg.nx-1);
                if (cfg.getOr(false,"write_throats"))                  VThroats(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_throats" + toStr(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr(false,"write_poreMaxBalls"))             poreMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_poreMBs" + toStr(blockNumber) + cfg.imgfrmt, 0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                if (cfg.getOr(false,"write_throatMaxBalls"))           throatMaxBalls(mpn, beginSlice, endSlice).writeBin(cfg.name()+"_throatMBs" + toStr(blockNumber) + cfg.imgfrmt,0 , cfg.nx, 0,cfg.ny,0, endSlice-beginSlice);
                blockNumber ++;
            }
		  }

																if (cfg.getOr(false,"write_vtkNetwork"))	vtuWriteThroats(cfg.name()+"_throats",  mpn.poreIs, mpn.throatIs, cfg.precision, mpn.VElems.X0()+mpn.VElems.dx());

	cout<<endl<<cfg.name()<<endl;
	cout<<"***  " <<mpn.poreIs.size()-2<<" pores, "<<mpn.throatIs.size()<<" throats,   ratio: "<<double(mpn.throatIs.size())/(mpn.poreIs.size()-2.0+1e-6)<<"  ***"<<endl;
	cout<<"end"<<endl;


 return 0;
}

