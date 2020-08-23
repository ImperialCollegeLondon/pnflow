

#include "readSetCAs.h"



void setContactAngles(vector<VoidElem*>& pores, vector<VoidElem*>& trots, double minCA, double maxCA,
							  double delta, double eta, int CAMdl, double CAMdl2SepAng, string  CACrl, int nBpPors, const CommonData& comn)
{

 if(CACrl.size()>4 && CACrl[4] == 'A') //rMaxAll rMinAll
 {
	const size_t nps = pores.size();
	const size_t nes = pores.size()+trots.size();
	vector<double> conAngles(nes);
	for(size_t i = 0; i < nes; ++i)  conAngles[i] = comn.weibull(minCA, maxCA, delta, eta);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	vector<VoidElem*> elems(nes);
	for(size_t i = 0; i < pores.size()  ; ++i) elems[i]=pores[i];
	for(size_t i = 0; i < trots.size(); ++i) elems[i+nps]=trots[i];

	if(CACrl[2] == 'a' || CACrl[2] == 'A')      // rMax
	  sort(elems.begin(), elems.end(), ElemRadCmpRed());
	else if(CACrl[2] == 'i' || CACrl[2] == 'I') // rMin
	  sort(elems.begin(), elems.end(), ElemRadCmpInc());
	else
	  shuffle(elems.begin(), elems.end(), comn.randomGenerator());            // Random


	for(size_t i = 0; i < nes; ++i)
	{
		if(auto elem = dynamic_cast<VoidElem*>(elems[i])) elem->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
		//else
		//if(auto elem = dynamic_cast<VoidElem*>(elems[i])) elem->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
	}
 }
 else if(CACrl.size()>4 && CACrl[4] == 'E')//rMaxEach rMinEach
 {
  {
	vector<double> conAngles(pores.size());
	for(size_t i = 0; i < conAngles.size(); ++i)
	  conAngles[i] = comn.weibull(minCA, maxCA, delta, eta);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(CACrl[2] == 'a' || CACrl[2] == 'A')      // rMax
	  sort(pores.begin(), pores.end(), ElemRadCmpRed());
	else if(CACrl[2] == 'i' || CACrl[2] == 'I') // rMin
	  sort(pores.begin(), pores.end(), ElemRadCmpInc());
	else
	  shuffle(pores.begin(), pores.end(), comn.randomGenerator());            // Random


	for(size_t i = 0; i < pores.size(); ++i)
		pores[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
  }

  {
	vector<double> conAngles(trots.size());
	for(size_t i = 0; i < conAngles.size(); ++i)
	  conAngles[i] = comn.weibull(minCA, maxCA, delta, eta);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(CACrl[2] == 'a' || CACrl[2] == 'A')      // rMax
	  sort(trots.begin(), trots.end(), ElemRadCmpRed());
	else if(CACrl[2] == 'i' || CACrl[2] == 'I') // rMin
	  sort(trots.begin(), trots.end(), ElemRadCmpInc());
	else
	  shuffle(trots.begin(), trots.end(), comn.randomGenerator());            // Random


	for(size_t i = 0; i < trots.size(); ++i)
		trots[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
  }
 }
 else
 {
	vector<double> conAngles(pores.size());
	for(size_t i = 0; i < conAngles.size(); ++i)
	  conAngles[i] = comn.weibull(minCA, maxCA, delta, eta);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(CACrl[2] == 'a' || CACrl[2] == 'A')      // rMax
	  sort(pores.begin(), pores.end(), ElemRadCmpRed());
	else if(CACrl[2] == 'i' || CACrl[2] == 'I') // rMin
	  sort(pores.begin(), pores.end(), ElemRadCmpInc());
	else
	  shuffle(pores.begin(), pores.end(), comn.randomGenerator());            // Random


	vector<double> conAnglesMaped(nBpPors,-1.0);
	for(size_t i = 0; i < pores.size(); ++i)
	{
		pores[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
		conAnglesMaped[pores[i]->index()]=conAngles[i];
	}

	for(size_t k = 0; k < trots.size(); ++k)
	{
		const Element* pOne = trots[k]->connection(0);
		const Element* pTwo = trots[k]->connection(1);

		if( conAnglesMaped[pOne->index()]>=0 && conAnglesMaped[pTwo->index()]>=0 /*clusterIndx[pOne->index()]>0 && clusterIndx[pTwo->index()]>0*/ )
		{
			if(comn.rand01()>0.5)
			{
				 trots[k]->setContactAngle(conAnglesMaped[pOne->index()], CAMdl, CAMdl2SepAng);
				 d_assert(conAnglesMaped[pOne->index()]>=0);
			}
			else
			{
				 trots[k]->setContactAngle(conAnglesMaped[pTwo->index()], CAMdl, CAMdl2SepAng);
				 d_assert(conAnglesMaped[pTwo->index()]>=0);
			}
		}
		else if(conAnglesMaped[pTwo->index()]>=0 /*clusterIndx[pTwo->index()] > 0 &&*/)
			trots[k]->setContactAngle(conAnglesMaped[pTwo->index()], CAMdl, CAMdl2SepAng);
		else
		{
			if (conAnglesMaped[pOne->index()]>=0)
				trots[k]->setContactAngle(conAnglesMaped[pOne->index()], CAMdl, CAMdl2SepAng);
			else
				trots[k]->setContactAngle(comn.weibull(minCA, maxCA, delta, eta), CAMdl, CAMdl2SepAng);
			//d_assert(conAnglesMaped[pOne->index()]>=0);
		}
		d_assert(trots[k]->conAngleAdv()>0.0001);
//        cout<<trots[k]->conAngleAdv()<<" ";
	}
 }
 if(comn.informative) cout<<"Contact Angles ["<<int(minCA*180/PI)<<"-"<<int(maxCA*180/PI)<<"] set for "<<pores.size()<<" pores and "<<trots.size()<<" trots. Crl:"<<CACrl<<endl;
}


#ifndef READSETCAS_H
void readSetCAs(istringstream& data, const vector<Element*>& elemans, int nBpPors, mstream& out_)
{
		out_<<"\nreadSetCAs: Not supported "<<endl;
}
#endif


///. set init contact angles
void applyInitWettability(const InputFile& inp, const vector<Element*>& elemans, int nBpPors, mstream& out_, const CommonData& comn)
{
	out_<<"\nSetting initial (oil-injection) contact angles:"<<endl;

	istringstream data;
	double minAng(0.0), maxAng(0.0),  weiDelta(-10.0), weiEta(-2.0), CAMdl2SepAng(25.2);
	int       CAMdl = 1;
	string    CACrl = "rand";


	if (inp.getData(data,"INIT_CONT_ANG",0))
	{
		data >> CAMdl >> minAng >> maxAng >> weiDelta >> weiEta;
		if (data.good()) data >> CACrl;
		if (data.good()) data >> CAMdl2SepAng;
		inp.checkEndOfData(data, "INIT_CONT_ANG");
	}
	else if(inp.getData(data,"INIT_CON_ANG",0)) /// Maintain backward compatibility
	{
		data >> minAng >> maxAng >> weiDelta >> weiEta;
		if (data.good()) data >> CACrl;
		inp.checkEndOfData(data, "INIT_CON_ANG");
	}
	else if( (!inp.giv("READ_INIT_CA",data,0)) && inp.giv("EQUIL_CON_ANG",data,0) ) /// Maintain backward compatibility
	{
		out_<<"\nINIT_CONT_ANG not found,\n -> Setting contact angles from EQUIL_CON_ANG:"<<endl;

		data >> CAMdl >> minAng >> maxAng;
		if (data.good()) data >> weiDelta >> weiEta;
		if (data.good()) data >> CACrl;
		if (data.good()) data >> CAMdl2SepAng;

		inp.checkEndOfData(data, "EQUIL_CON_ANG","EQUIL_CON_ANG:  CAModel  minAng  maxAng  weiDelta  weiEta CACrl CAMdl2SepAng ;");
	}
	minAng *= acos(-1.0) / 180.0;
	maxAng *= acos(-1.0) / 180.0;
	CAMdl2SepAng *= acos(-1.0) / 180.0;

	vector<VoidElem*> pores2Set;  pores2Set.reserve(nBpPors);
	for(int i = 2; i < nBpPors; ++i)
		if (auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(WTR))  pores2Set.push_back(elem);
	vector<VoidElem*> trots2Set;trots2Set.reserve(elemans.size()-nBpPors);
	for(size_t i = nBpPors; i < elemans.size(); ++i)
		if (auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(WTR))  trots2Set.push_back(elem);

	setContactAngles(pores2Set, trots2Set, minAng, maxAng, weiDelta, weiEta, CAMdl, CAMdl2SepAng,  CACrl, nBpPors, comn);
	if (inp.getData(data,"READ_INIT_CA",0))
	{
		readSetCAs(data, elemans, nBpPors, out_);
	}
}




/** Change the advancing contact angle depending on several criteria set out in
 * inp file. Wettability alterations will only occur in elems that have been
 * drained */
void applyFWettabilityChange(const InputFile& inp, const vector<Element*>& elemans, int nBpPors, mstream& out_, const CommonData& comn)
{

	out_<<"\nWettability change, setting water-injection contact angles:"<<endl;


	string  CACrl("rand");
	int     CAMdl = 1;
	double  CAMdl2SepAng(0.0);
	istringstream data;

	int warn=inp.giv("READ_ALTR_CA",data,0) ? 0 : 1;

	if(	inp.giv("INIT_CONT_ANG",data,0) || inp.giv("INIT_CON_ANG",data,0) || inp.giv("READ_INIT_CA",data,0) ) //! if INIT_CONT_ANG not set don't reset the EQUIL_CON_ANG, it is already set
	if(inp.giv("EQUIL_CON_ANG",data,warn))
	{		///set cycle 2 contact angles
		double minAng(0.0), maxAng(0.0),  weiDelta(-0.2), weiEta(-3.0), CAMdl2SepAng(25.2);
		warn=0;
		data >> CAMdl >> minAng >> maxAng;
		if (data.good()) data >> weiDelta >> weiEta;
		if (data.good()) data >> CACrl;
		if (data.good()) data >> CAMdl2SepAng;

		inp.checkEndOfData(data, "EQUIL_CON_ANG","EQUIL_CON_ANG  CAModel  minAng  maxAng  weiDelta  weiEta CACrl CAMdl2SepAng ;");

		minAng *= acos(-1.0) / 180.0;
		maxAng *= acos(-1.0) / 180.0;
		CAMdl2SepAng *= acos(-1.0) / 180.0;

		vector<VoidElem*> pores2Set;pores2Set.reserve(nBpPors);
		vector<VoidElem*> trots2Set;trots2Set.reserve(elemans.size()-nBpPors);

		for(int i = 2; i < nBpPors; ++i)
			if(auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(OIL))  pores2Set.push_back(elem);

		for(size_t i = nBpPors; i < elemans.size(); ++i)
			if(auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(OIL))  trots2Set.push_back(elem);

		setContactAngles(pores2Set, trots2Set, minAng, maxAng, weiDelta, weiEta, CAMdl, CAMdl2SepAng,  CACrl, nBpPors, comn);

	}


	if( inp.giv("FRAC_CONT_OPT",data) || inp.giv("FRAC_CON_ANG",data) )
	{
		out_ << "================================================================"    << endl
				<< "Applying fractional wetting"                                       << endl;


		bool volBased(false), totalFraction(false), oilClustInWat(true);
		double fraction, minAng, maxAng, weiDelta(-0.2), weiEta(-3.0), deltaClustD(-0.2), etaClustD(-3.0);
		int totElem(0), clustDiam1(0), clustDiam2;
		string spatialDistrib("rand"), clustSeedScheme("rMax");

		if( inp.getData(data,"FRAC_CONT_OPT") )
		{

			char volBase('x'), totalFrac('x');
			data >> fraction >> volBase >> totalFrac  >> spatialDistrib;
			volBased = (volBase == 'T' || volBase == 't' || volBase == 'V' || volBase == 'v');
			if (!(volBase=='T' || volBase=='t' || volBase == 'V' || volBase == 'v' || volBase=='F' || volBase=='f' || volBase=='N' || volBase=='n'))	{ out_<<" Error Wrong choice for \"volume-/number- based\" fraction, expected T/t/V/v or F/f/N/n"<<endl; exit(-1); }
			totalFraction = (totalFrac == 'T' || totalFrac == 't');
			if (!(totalFrac=='T' || totalFrac=='t' || totalFrac=='F' || totalFrac=='f' || totalFrac=='O' || totalFrac=='o' || totalFrac=='W' || totalFrac=='w'))	{ out_<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}

			char oInW('T');
			if (spatialDistrib[0]=='C' || spatialDistrib[0]=='c')
			{
				data >> oInW >> clustDiam1 >> clustDiam2;
				if (data.good()) data >> deltaClustD >> etaClustD;
				if (data.good()) data >> clustSeedScheme;                  // Maintain backward compatibility
				oilClustInWat = (oInW == 'W' || oInW == 'w' || oInW == 'T' || oInW == 't');
				inp.Assert(clustDiam1>0 && clustDiam2>=clustDiam1, "FRAC_CONT_OPT:  fraction  volBase spatialDistrib  oilInWCluster clustDiam1 clustDiam2 delta  eta clustSeedScheme");
			}
			else clustSeedScheme = spatialDistrib;
			inp.checkEndOfData(data, "FRAC_CONT_OPT","",false);


			if ( inp.getData(data,"FRAC_CONT_ANG",2) )
			{
				data >> CAMdl >> minAng >> maxAng;
				if (data.good()) data >> weiDelta >> weiEta;
				if (data.good()) data >> CACrl;
				if (data.good()) data >> CAMdl2SepAng;

				inp.checkEndOfData(data, "FRAC_CONT_ANG");

				minAng *= acos(-1.0) / 180.0;
				maxAng *= acos(-1.0) / 180.0;
				CAMdl2SepAng *= acos(-1.0) / 180.0;

			}
		}
		else if( inp.getData(data,"FRAC_CON_ANG") )     ///. Maintain backward compatibility
	   {

			char volBase('x');
			data >> fraction >> volBase  >> minAng >> maxAng >> weiDelta >> weiEta >> spatialDistrib;
			minAng *= acos(-1.0) / 180.0;
			maxAng *= acos(-1.0) / 180.0;
			volBased = (volBase == 'T' || volBase == 't');
			if (!(volBase == 'T' || volBase == 't' || volBase == 'F' || volBase == 'f'))
				out_<<"      Error Wrong value for volBase, expected T/t or F/f"<<endl;

			char oInW('T');
			if (spatialDistrib[0]=='C' || spatialDistrib[0]=='c')
			{
				data >> clustDiam1 >> oInW;
				oilClustInWat = (oInW == 'T' || oInW == 't');
				inp.Assert(clustDiam1>1, "FRAC_CON_ANG");
				clustDiam2 = clustDiam1;
				clustSeedScheme = "rMax";
			}
			else clustSeedScheme = spatialDistrib;
			inp.checkEndOfData(data, "FRAC_CON_ANG","",false);
	   }




		vector<VoidElem*> poresToSeed;poresToSeed.reserve(nBpPors);
	   {
			for(int i = 2; i < nBpPors; ++i)
				if (auto elem=dynamic_cast<VoidElem*>(elemans[i]))  if(elem->exists(OIL))  poresToSeed.push_back(elem);

			vector<int> clustDiams(poresToSeed.size());
			for(size_t i = 0; i < clustDiams.size(); ++i)
			  clustDiams[i] = comn.weibull(clustDiam1, clustDiam2, deltaClustD, etaClustD)+0.5;
			sort(clustDiams.begin(), clustDiams.end(), greater<int>());

			if(clustSeedScheme[2] == 'a' || clustSeedScheme[2] == 'A')      // rMax
				  sort(poresToSeed.begin(), poresToSeed.end(), ElemRadCmpInc());
			else if(clustSeedScheme[2] == 'i' || clustSeedScheme[2] == 'I') // rMin
				  sort(poresToSeed.begin(), poresToSeed.end(), ElemRadCmpRed());
			else  shuffle(poresToSeed.begin(), poresToSeed.end(), comn.randomGenerator());            // Random


	   }
		vector<VoidElem*> pores2BAltered; poresToSeed.reserve(poresToSeed.size());
		vector<VoidElem*> trots2BAltered; trots2BAltered.reserve(poresToSeed.size()*3);


		double oilInvadedVol(0.0), TotalVol(0.0), oilwettedVol(0.0);

		for(size_t i = 2; i < elemans.size(); ++i)
		{
			if (dynamic_cast<const VoidElem*>(elemans[i]->model()))
			{
				if(elemans[i]->model()->exists(OIL))
				{
					if(i < size_t(nBpPors)) ++totElem;                        // Only count pores
					oilInvadedVol += elemans[i]->flowVolumeX();
				}
				TotalVol += elemans[i]->flowVolumeX();
			}
		}

		double targetFrac = fraction;
		if (totalFraction)  targetFrac = min(fraction*TotalVol/oilInvadedVol, 1.0);
		if (oilInvadedVol/TotalVol>1.0e-12)
		{

			vector<int> clusterIndx(elemans.size()+2, 0);

			int numFracWetted(0), clusterIdx(1);
			if(spatialDistrib[0]=='C' || spatialDistrib[0]=='c')    // Spatial correlation approach
			{
				bool clct_more(true);


				targetFrac = oilClustInWat ? targetFrac: 1.0-targetFrac;

//				long long counter(0);
				while(clct_more)
				{
//					if (++counter>nBpPors*100) {out_<<"Warning: could not apply the requested wetting fraction"<<endl; break;}

					set< VoidElem * > frontier;
//					double randNum = comn.rand01();
//					int randPoreIdx = randNum*nBpPors;
					VoidElem* elem = poresToSeed.back(); poresToSeed.pop_back();

					if (clusterIndx[elem->index()]==0)
					{
//						randNum = comn.rand01();
						int clustDiam = comn.weibull(clustDiam1,clustDiam2, deltaClustD, etaClustD) +0.5;

						frontier.insert(elem);
						int frontExpansions(0);
						while(!frontier.empty() && clct_more && frontExpansions < clustDiam)  // Front expansions will go over thraots as well as pores, wheras
						{                                                                           // cluster diam is wrt pores  => Don't divide by 2
							set<VoidElem*> oldFrontier = frontier;
							frontier.clear();
							++frontExpansions;
							for(auto elm:oldFrontier)
							{

								if( clusterIndx[elm->index()] != clusterIdx && clct_more)
								{
									if(clusterIndx[elm->index()] == 0 && elm->exists(OIL))
									{
										clusterIndx[elm->index()] = clusterIdx;                    /// Set flag
										oilwettedVol += elm->eleman()->flowVolumeX();
										if(dynamic_cast<const Pore*>(elm->eleman())) ++numFracWetted;   /// Only count pores

										clct_more = volBased ? oilwettedVol/oilInvadedVol < targetFrac :  double(numFracWetted)/totElem < targetFrac;

										if(!clct_more)			break;
									}

									for(int i = 0; i < elm->eleman()->connectionNum(); ++i)
									{
										if ( clusterIndx[elm->connection(i)->index()] != clusterIdx)
										if(auto el = dynamic_cast<VoidElem*>(elm->ChParent()->connection(i)->ChModel()))
											frontier.insert(el);

									}
								}
							}
						}
						++clusterIdx;
					}
				}
				numFracWetted = 0;
				oilwettedVol = 0.0;
				if(oilClustInWat)
				{
					for(size_t i = 0; i < elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->index()] > 0 && elemans[i]->model()->exists(OIL))
						{
							if(int(i) < nBpPors) { if(auto elem = dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								pores2BAltered.push_back(elem);
							} }
							else if(auto elem = dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								trots2BAltered.push_back(elem);
							}
						}
					}
				}
				else
				{
					for(size_t i = 0; i < elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->index()] == 0 && elemans[i]->model()->exists(OIL))
						{
							if(int(i) < nBpPors) { if(auto elem = dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								pores2BAltered.push_back(elem);
							} }
							else if(auto elem = dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								trots2BAltered.push_back(elem);
							}
						}
					}
				}
				out_<< "Number of correlated regions: " << clusterIdx-1                         << endl;

			}
			else ///. spatialDistrib == rand
			{

				int targetNum(fraction*poresToSeed.size());

				while(!poresToSeed.empty() &&
					((volBased && oilwettedVol/oilInvadedVol < targetFrac) || (!volBased && numFracWetted < targetNum)))
				{
					++numFracWetted;
					VoidElem* elem = poresToSeed.back(); poresToSeed.pop_back();
					oilwettedVol += elem->eleman()->flowVolumeX();
					pores2BAltered.push_back(elem);
					clusterIndx[elem->index()]=1;
					for(int t = 0; t < elem->eleman()->connectionNum(); ++t)
					{
						VoidElem* throat = dynamic_cast<VoidElem*>(elem->ChParent()->connection(t)->ChModel());
						if( throat && throat->exists(OIL) && clusterIndx[throat->index()] == 0)
						{
							double randNum = comn.rand01();
							if(randNum < fraction)
							{
								trots2BAltered.push_back(throat);
								clusterIndx[throat->index()]=1;
								oilwettedVol += throat->eleman()->flowVolumeX();
							}
						}
					}
				}

			}
			out_
				<< "Altered volume (of oil invaded volume): " << oilwettedVol/oilInvadedVol   << endl
				<< "Altered volume (of total net pore volume): " << oilwettedVol/TotalVol   << endl
				<< "Number of pores altered: " << numFracWetted                            << endl
				<< "================================================================="  << endl;

	   }
	   else
			out_
				<< "\nWarning:  fractional wetting was NOT applied, "                  << endl
				<< "too low oil invaded fraction: " << oilInvadedVol/TotalVol             << endl
				<< "=================================================================\n"      << endl;

		setContactAngles(pores2BAltered, trots2BAltered, minAng, maxAng, weiDelta, weiEta, CAMdl, CAMdl2SepAng,  CACrl, nBpPors, comn);

	}
	else if (inp.giv("FRAC_CONT_ANG",data) )
		cout<<"\n******************************************"<<endl
			 <<" Warning keyword FRAC_CONT_ANG is ignored"<<endl
			 <<" because keyword FRAC_CONT_OPT is missing"<<endl
			 <<"******************************************\n"<<endl;

	if(inp.giv("READ_ALTR_CA",data,0))
	{
		readSetCAs(data, elemans, nBpPors, out_);
		return;
	}
	else if(!warn)
	{
		out_<<"\n*******************************"<<endl;
		out_<<"EQUIL_CON_ANG  not found"<<endl;
		out_<<"Leaving contact angles unchanged"<<endl;
		out_<<"********************************"<<endl;
	}
}


void setElemProps(const InputFile& inp, const vector<Element*>& elemans, size_t nBpPors, mstream& out_, const CommonData& comn)
{
	if (comn.dispCycle()==1)
		applyInitWettability(inp, elemans, nBpPors, out_, comn);
	else if (comn.dispCycle()==2)
		applyFWettabilityChange(inp, elemans, nBpPors, out_, comn);
}




