#ifndef READSETCAS_H
#define READSETCAS_H
void readSetCAs(istringstream& data, const vector<Elem*>& elemans, int nBpPors, mstream& out_)
{ 
	string elmTyp("pore"), fnam;
	int CAMdl(4);
	const double piy180 = acos(-1.) / 180.;
	double mdl2SepAng(90.*piy180);

	data >>fnam >> elmTyp>>CAMdl>>mdl2SepAng;
	out_<<"CA model"<<CAMdl; if(CAMdl==2||CAMdl==5) out_<<",  CA_adv - CA_rec = "<<mdl2SepAng; out_<<endl;
	out_<<"reading "<<elmTyp<<" contact angles from file "<<fnam<<endl;

	if(elmTyp[0]=='p')
	{
		dbls CAs(nBpPors,-1.);
		ifstream ContAngFile(fnam);
		ensure(ContAngFile);
		std::string line;
		std::getline(ContAngFile,line);
		for(size_t i=0; i<CAs.size(); ++i)  { ContAngFile>>CAs[i];   CAs[i]*=piy180; }
		ensure(CAs.back()>=0., "Wrong pore contact angles, expected "+_s(nBpPors)+" values of [0-180] in file "+fnam, -1);
		for(int i=0; i<nBpPors; ++i)
			if (dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
				dynamic_cast<VoidElem*>(elemans[i]->ChModel())->setContactAngle(CAs[i], CAMdl, mdl2SepAng);

		for(size_t i = nBpPors; i < elemans.size(); ++i)
		{
			VoidElem* trot=dynamic_cast<VoidElem*>(elemans[i]->ChModel());
			if (trot)
				trot->setContactAngle(0.5*(CAs[elemans[i]->neib(0)->index()]+CAs[elemans[i]->neib(1)->index()]), CAMdl, mdl2SepAng);
		}
	}
	else if(elmTyp[0]=='a')
	{
		vector<double> CAs(elemans.size(),-1.);
		ifstream ContAngFile(fnam);
		ensure(ContAngFile);
		std::string line;
		std::getline(ContAngFile,line);
		for(size_t i=0; i<CAs.size(); ++i) ContAngFile>>CAs[i];
		ensure(CAs.back()>=0., "Wrong pore+throat contact angles, expected "+_s(elemans.size())+" values of [0-180] in file "+fnam, -1);

		for(int i=0; i<len(elemans); ++i)
			if (dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
				dynamic_cast<VoidElem*>(elemans[i]->ChModel())->setContactAngle(CAs[i], CAMdl, mdl2SepAng);
	}
	else if(elmTyp[0]=='c')
	{
		out_<<"Wrong elem type for contact angles: "<<elmTyp<<endl;
		out_<<" only pore is accepted "<<elmTyp<<endl;
		exit(-1);
	}
	else
	{
		out_<<"Wrong elem type for contact angles: "<<elmTyp<<endl;
		out_<<" only pore and corner are accepted "<<elmTyp<<endl;
		exit(-1);
	}


}
#endif

