

#include "readSetCAs.h"



void setContactAngles(vector<VoidElem*>& pores, vector<VoidElem*>& trots, const Weibul1& wb, int CAMdl, double CAMdl2SepAng, int nBpPors, const GNMData& comn)
{

 if(wb.cor.size()>4 && wb.cor[4] == 'A') //rMaxAll rMinAll
 {
	const size_t nps= pores.size();
	const size_t nes= pores.size()+trots.size();
	vector<double> conAngles(nes);
	for(size_t i= 0; i < nes; ++i)  conAngles[i]= comn.weibul1(wb);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	vector<VoidElem*> elems(nes);
	for(size_t i= 0; i<pores.size() ; ++i) elems[i]=pores[i];
	for(size_t i= 0; i<trots.size(); ++i) elems[i+nps]=trots[i];

	if(wb.cor[2] == 'a' || wb.cor[2] == 'A')      // rMax
	  sort(elems.begin(), elems.end(), ElemRadCmpRed());
	else if(wb.cor[2] == 'i' || wb.cor[2] == 'I') // rMin
	  sort(elems.begin(), elems.end(), ElemRadCmpInc());
	else
	  shuffle(elems.begin(), elems.end(), comn.randomGenerator());            // Random


	for(size_t i= 0; i<nes; ++i)
	{
		if(auto elem= dynamic_cast<VoidElem*>(elems[i])) elem->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
		//else
		//if(auto elem= dynamic_cast<VoidElem*>(elems[i])) elem->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
	}
 }
 else if(wb.cor.size()>4 && wb.cor[4] == 'E')//rMaxEach rMinEach
 {
  {
	vector<double> conAngles(pores.size());
	for(size_t i= 0; i<conAngles.size(); ++i)
	  conAngles[i]= comn.weibul1(wb);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(wb.cor[2] == 'a' || wb.cor[2] == 'A')      // rMax
	  sort(pores.begin(), pores.end(), ElemRadCmpRed());
	else if(wb.cor[2] == 'i' || wb.cor[2] == 'I') // rMin
	  sort(pores.begin(), pores.end(), ElemRadCmpInc());
	else
	  shuffle(pores.begin(), pores.end(), comn.randomGenerator());            // Random


	for(size_t i= 0; i<pores.size(); ++i)
		pores[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
  }

  {
	vector<double> conAngles(trots.size());
	for(size_t i= 0; i<conAngles.size(); ++i)
	  conAngles[i]= comn.weibul1(wb);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(wb.cor[2] == 'a' || wb.cor[2] == 'A')      // rMax
	  sort(trots.begin(), trots.end(), ElemRadCmpRed());
	else if(wb.cor[2] == 'i' || wb.cor[2] == 'I') // rMin
	  sort(trots.begin(), trots.end(), ElemRadCmpInc());
	else
	  shuffle(trots.begin(), trots.end(), comn.randomGenerator());            // Random


	for(size_t i= 0; i<trots.size(); ++i)
		trots[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
  }
 } else {
	vector<double> conAngles(pores.size());
	for(size_t i= 0; i<conAngles.size(); ++i)
	  conAngles[i]= comn.weibul1(wb);
	sort(conAngles.begin(), conAngles.end(), greater<double>());

	if(wb.cor[2] == 'a' || wb.cor[2] == 'A')      // rMax
	  sort(pores.begin(), pores.end(), ElemRadCmpRed());
	else if(wb.cor[2] == 'i' || wb.cor[2] == 'I') // rMin
	  sort(pores.begin(), pores.end(), ElemRadCmpInc());
	else
	  shuffle(pores.begin(), pores.end(), comn.randomGenerator());            // Random


	vector<double> conAnglesMaped(nBpPors,-1.);
	for(size_t i= 0; i<pores.size(); ++i)
	{
		pores[i]->setContactAngle(conAngles[i], CAMdl, CAMdl2SepAng);
		conAnglesMaped[pores[i]->index()]=conAngles[i];
	}

	for(size_t k= 0; k < trots.size(); ++k)
	{
		const Elem* pOne= trots[k]->neib(0);
		const Elem* pTwo= trots[k]->neib(1);

		if( conAnglesMaped[pOne->index()]>=0 && conAnglesMaped[pTwo->index()]>=0 /*clusterIndx[pOne->index()]>0 && clusterIndx[pTwo->index()]>0*/)
		{
			if(comn.rand01()>0.5)
			{
				 trots[k]->setContactAngle(conAnglesMaped[pOne->index()], CAMdl, CAMdl2SepAng);
				 dbgAsrt(conAnglesMaped[pOne->index()]>=0);
			}
			else
			{
				 trots[k]->setContactAngle(conAnglesMaped[pTwo->index()], CAMdl, CAMdl2SepAng);
				 dbgAsrt(conAnglesMaped[pTwo->index()]>=0);
			}
		}
		else if(conAnglesMaped[pTwo->index()]>=0 /*clusterIndx[pTwo->index()] > 0 &&*/)
			trots[k]->setContactAngle(conAnglesMaped[pTwo->index()], CAMdl, CAMdl2SepAng);
		else
		{
			if (conAnglesMaped[pOne->index()]>=0)
				trots[k]->setContactAngle(conAnglesMaped[pOne->index()], CAMdl, CAMdl2SepAng);
			else
				trots[k]->setContactAngle(comn.weibul1(wb), CAMdl, CAMdl2SepAng);
			//dbgAsrt(conAnglesMaped[pOne->index()]>=0);
		}
		dbgAsrt(trots[k]->conAngleAdv()>0.0001);
//        cout<<trots[k]->conAngleAdv()<<" ";
	}
 }
 if(comn.input().informative) cout<<"Contact Angles ["<<int(wb.minV*180/PI)<<"-"<<int((wb.minV+wb.delV)*180/PI)<<"] set for "<<pores.size()<<" pores and "<<trots.size()<<" trots. Crl:"<<wb.cor<<endl;
}


#ifndef READSETCAS_H
void readSetCAs(istringstream& data, const vector<Elem*>& elemans, int nBpPors, mstream& out_, const GNMData& comn)  { out_<<"\nreadSetCAs: Not supported "<<endl; }
#endif


///. set init contact angles
void applyInitWettability(const InputFile& inp, const vector<Elem*>& elemans, int nBpPors, mstream& out_, const GNMData& comn)
{
	out_<<"\nSetting initial (oil-injection) contact angles:"<<endl;

	istringstream data;
	double CAMdl2SepAng(25.2);
	int       CAMdl= 1;
	Weibul1 wbCA;


	if (inp.giv("INIT_CONT_ANG", data,0))
	{
		data >> CAMdl >> wbCA;
		if (data.good()) data >> CAMdl2SepAng;
	}
	else if(inp.giv("INIT_CON_ANG", data,0)) /// Maintain backward compatibility
	{
		data >> wbCA;
	}
	else if( (!inp.giv("READ_INIT_CA",data,0)) && inp.giv("EQUIL_CON_ANG",data,0) ) /// Maintain backward compatibility
	{
		out_<<"\nINIT_CONT_ANG not found,\n -> Setting contact angles from EQUIL_CON_ANG:"<<endl;

		data >> CAMdl >> wbCA;
		if (data.good()) data >> CAMdl2SepAng;

	}
	wbCA.minV *= acos(-1.)/180.;
	wbCA.scale(acos(-1.)/180.);
	CAMdl2SepAng *= acos(-1.) / 180.;

	vector<VoidElem*> pores2Set;  pores2Set.reserve(nBpPors);
	for(int i= 2; i < nBpPors; ++i)
		if (auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(WTR))  pores2Set.push_back(elem);
	vector<VoidElem*> trots2Set;trots2Set.reserve(elemans.size()-nBpPors);
	for(size_t i= nBpPors; i < elemans.size(); ++i)
		if (auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(WTR))  trots2Set.push_back(elem);

	setContactAngles(pores2Set, trots2Set, wbCA, CAMdl, CAMdl2SepAng, nBpPors, comn);
	if (inp.giv("READ_INIT_CA", data,0)) readSetCAs(data, elemans, nBpPors, out_);
}




/** Change the advancing contact angle depending on several criteria set out in
 * input file. Wettability alterations will only occur in elems that have been
 * drained */
void applyFWettabilityChange(const InputFile& inp, const vector<Elem*>& elemans, int nBpPors, mstream& out_, const GNMData& comn)
{

	out_<<"\nWettability change, setting water-injection contact angles:"<<endl;


	int     CAMdl= 1;
	Weibul1 wbCA; double CAMdl2SepAng(25.2 *acos(-1.)/180.);
	istringstream data;

	int warn=inp.giv("READ_ALTR_CA",data,0) ? 0 : 1;

	if( (	inp.giv("INIT_CONT_ANG",data,0) || inp.giv("INIT_CON_ANG",data,0) || inp.giv("READ_INIT_CA",data,0) ) //! reset the EQUIL_CON_ANG only if INIT_CONT_ANG is also set, otherwise do not as it is already set
		&& inp.giv("EQUIL_CON_ANG",data,warn) )
	{		///set cycle 2 contact angles
		warn=0;
		data  >> CAMdl >> wbCA;  	wbCA.minV *= acos(-1.)/180.;    wbCA.scale(acos(-1.)/180.);
		if (data.good()) {
			data >> CAMdl2SepAng;	CAMdl2SepAng *= acos(-1.)/180.;  inp.checkEndOfData(data, "EQUIL_CON_ANG","EQUIL_CON_ANG  CAModel  minAng  maxAng  weiDelta  weiEta CACrl CAMdl2SepAng ;"); }


		vector<VoidElem*> pores2Set;pores2Set.reserve(nBpPors);
		vector<VoidElem*> trots2Set;trots2Set.reserve(elemans.size()-nBpPors);

		for(int i= 2; i < nBpPors; ++i)
			if(auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(OIL))  pores2Set.push_back(elem);

		for(size_t i= nBpPors; i < elemans.size(); ++i)
			if(auto elem=dynamic_cast<VoidElem*>(elemans[i]->ChModel())) if(elem->exists(OIL))  trots2Set.push_back(elem);

		setContactAngles(pores2Set, trots2Set, wbCA, CAMdl, CAMdl2SepAng, nBpPors, comn);

	}


	if( inp.giv("FRAC_CONT_OPT", data) )//|| inp.giv("FRAC_CON_ANG", data) )
	{
		out_ << "================================================================"    << endl
				<< "Applying fractional wetting"                                       << endl;


		bool volBased(false), totalFraction(false), oilClustInWat(true);
		double fraction; std::string spatialDistrib("rand");
		int totElem(0);
		Weibul1 wbClustr; // wbClustr.cor: cluster seed scheme

		if( inp.giv("FRAC_CONT_OPT", data))
		{

			char volBase('x'), totalFrac('x');
			data >> fraction >> volBase >> totalFrac  >> spatialDistrib;
			volBased= (volBase == 'T' || volBase == 't' || volBase == 'V' || volBase == 'v');
			if (!(volBase=='T' || volBase=='t' || volBase == 'V' || volBase == 'v' || volBase=='F' || volBase=='f' || volBase=='N' || volBase=='n'))	{ out_<<" Error Wrong choice for \"volume-/number- based\" fraction, expected T/t/V/v or F/f/N/n"<<endl; exit(-1); }
			totalFraction= (totalFrac == 'T' || totalFrac == 't');
			if (!(totalFrac=='T' || totalFrac=='t' || totalFrac=='F' || totalFrac=='f' || totalFrac=='O' || totalFrac=='o' || totalFrac=='W' || totalFrac=='w'))	{ out_<<" Error Wrong choice for fraction of \"total/oil-invaded\" elements, expected T/t or F/f/O/o"<<endl; exit(-1);}

			char oInW('T');
			if (spatialDistrib[0]=='C' || spatialDistrib[0]=='c')
			{
				data >> oInW >> wbClustr;                  // Maintain backward compatibility
				oilClustInWat= (oInW == 'W' || oInW == 'w' || oInW == 'T' || oInW == 't');
				inp.Assert(wbClustr.minV>0 && wbClustr.delV>=0, "FRAC_CONT_OPT:  fraction  volBase spatialDistrib  oilInWCluster clustDiam1 clustDiam2 delta  eta wbClustr.cor");
			}
			else wbClustr.cor= spatialDistrib;
			inp.checkEndOfData(data, "FRAC_CONT_OPT","",false);


			if ( inp.giv("FRAC_CONT_ANG", data,2))
			{
				data >> CAMdl >> wbCA; 	wbCA.minV *= acos(-1.)/180.;    wbCA.scale(acos(-1.)/180.);
				if (data.good()) { data >> CAMdl2SepAng; CAMdl2SepAng *= acos(-1.)/180.;}

				inp.checkEndOfData(data, "FRAC_CONT_ANG");

			}
		}
		else if( inp.giv("FRAC_CON_ANG", data) )  {cout<<"*** Error:  old keyword FRAC_CON_ANG ignored ***\n\n\n"; }   ///. Maintain backward compatibility
	   //{

			//char volBase('x');
			//data >> fraction >> volBase  >> wbCA;  	wbCA.minV *= acos(-1.)/180.;    wbCA.maxV *= acos(-1.)/180.;

			//volBased= (volBase == 'T' || volBase == 't');
			//if (!(volBase == 'T' || volBase == 't' || volBase == 'F' || volBase == 'f'))
				//out_<<"      Error Wrong value for volBase, expected T/t or F/f"<<endl;

			//char oInW('T');
			//if (spatialDistrib[0]=='C' || spatialDistrib[0]=='c')
			//{
				//data >> wbClustr.minV >> oInW;
				//oilClustInWat= (oInW == 'T' || oInW == 't');
				//inp.Assert(wbClustr.minV>1, "FRAC_CON_ANG");
				//wbClustr.maxV= wbClustr.minV;
				//wbClustr.cor= "rMax";
			//}
			//else wbClustr.cor= spatialDistrib;
			//inp.checkEndOfData(data, "FRAC_CON_ANG","",false);
	   //}



		vector<VoidElem*> poresToSeed;poresToSeed.reserve(nBpPors);
	   {
			for(int i= 2; i < nBpPors; ++i)
				if (auto elem=dynamic_cast<VoidElem*>(elemans[i]))  if(elem->exists(OIL))  poresToSeed.push_back(elem);

			vector<int> clustDiams(poresToSeed.size());
			for(size_t i= 0; i<clustDiams.size(); ++i)
			  clustDiams[i]= comn.weibul1(wbClustr)+0.5;
			sort(clustDiams.begin(), clustDiams.end(), greater<int>());

			if(wbClustr.cor[2] == 'a' || wbClustr.cor[2] == 'A')      // rMax
				  sort(poresToSeed.begin(), poresToSeed.end(), ElemRadCmpInc());
			else if(wbClustr.cor[2] == 'i' || wbClustr.cor[2] == 'I') // rMin
				  sort(poresToSeed.begin(), poresToSeed.end(), ElemRadCmpRed());
			else  shuffle(poresToSeed.begin(), poresToSeed.end(), comn.randomGenerator());            // Random


	   }
		vector<VoidElem*> pores2BAltered; poresToSeed.reserve(poresToSeed.size());
		vector<VoidElem*> trots2BAltered; trots2BAltered.reserve(poresToSeed.size()*3);


		double oilInvadedVol(0.), TotalVol(0.), oilwettedVol(0.);

		for(size_t i= 2; i < elemans.size(); ++i)
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

		double targetFrac= fraction;
		if (totalFraction)  targetFrac= min(fraction*TotalVol/oilInvadedVol, 1.);
		if (oilInvadedVol/TotalVol>1e-12)
		{

			vector<int> clusterIndx(elemans.size()+2, 0);

			int numFracWetted(0), clusterIdx(1);
			if(spatialDistrib[0]=='C' || spatialDistrib[0]=='c')    // Spatial correlation approach
			{
				bool clct_more(true);


				targetFrac= oilClustInWat ? targetFrac: 1.-targetFrac;

//				long long counter(0);
				while(clct_more)
				{
//					if (++counter>nBpPors*100) {out_<<"Warning: could not apply the requested wetting fraction"<<endl; break;}

					set< VoidElem * > frontier;
//					double randNum= comn.rand01();
//					int randPoreIdx= randNum*nBpPors;
					VoidElem* elem= poresToSeed.back(); poresToSeed.pop_back();

					if (clusterIndx[elem->index()]==0)
					{
//						randNum= comn.rand01();
						int clustDiam= comn.weibul1(wbClustr) +0.5;

						frontier.insert(elem);
						int frontExpansions(0);
						while(!frontier.empty() && clct_more && frontExpansions < clustDiam)  // Front expansions will go over thraots as well as pores, wheras
						{                                                                           // cluster diam is wrt pores  => Don't divide by 2
							set<VoidElem*> oldFrontier= frontier;
							frontier.clear();
							++frontExpansions;
							for(auto elm:oldFrontier)
							{

								if( clusterIndx[elm->index()] != clusterIdx && clct_more)
								{
									if(clusterIndx[elm->index()] == 0 && elm->exists(OIL))
									{
										clusterIndx[elm->index()]= clusterIdx;                    /// Set flag
										oilwettedVol += elm->eleman()->flowVolumeX();
										if(dynamic_cast<const Pore*>(elm->eleman())) ++numFracWetted;   /// Only count pores

										clct_more= volBased ? oilwettedVol/oilInvadedVol < targetFrac :  double(numFracWetted)/totElem < targetFrac;

										if(!clct_more)			break;
									}

									for(int i=0; i<elm->eleman()->nCncts(); ++i)
									{
										if ( clusterIndx[elm->neib(i)->index()] != clusterIdx)
										if(auto el= dynamic_cast<VoidElem*>(elm->ChParent()->neib(i)->ChModel()))
											frontier.insert(el);

									}
								}
							}
						}
						++clusterIdx;
					}
				}
				numFracWetted= 0;
				oilwettedVol= 0.;
				if(oilClustInWat)
				{
					for(size_t i= 0; i<elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->index()] > 0 && elemans[i]->model()->exists(OIL))
						{
							if(int(i) < nBpPors) { if(auto elem= dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								pores2BAltered.push_back(elem);
							} }
							else if(auto elem= dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								trots2BAltered.push_back(elem);
							}
						}
					}
				}
				else
				{
					for(size_t i= 0; i<elemans.size(); ++i)
					{
						if( clusterIndx[elemans[i]->index()] == 0 && elemans[i]->model()->exists(OIL))
						{
							if(int(i) < nBpPors) { if(auto elem= dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
							{
								++numFracWetted;  oilwettedVol+=elemans[i]->flowVolumeX();
								pores2BAltered.push_back(elem);
							} }
							else if(auto elem= dynamic_cast<VoidElem*>(elemans[i]->ChModel()))
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
					VoidElem* elem= poresToSeed.back(); poresToSeed.pop_back();
					oilwettedVol += elem->eleman()->flowVolumeX();
					pores2BAltered.push_back(elem);
					clusterIndx[elem->index()]=1;
					for(int t= 0; t < elem->eleman()->nCncts(); ++t)
					{
						VoidElem* throat= dynamic_cast<VoidElem*>(elem->ChParent()->neib(t)->ChModel());
						if( throat && throat->exists(OIL) && clusterIndx[throat->index()] == 0)
						{
							double randNum= comn.rand01();
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

		setContactAngles(pores2BAltered, trots2BAltered, wbCA, CAMdl, CAMdl2SepAng, nBpPors, comn);

	}
	else if (inp.giv("FRAC_CONT_ANG",data))
		cout<<"\n******************************************"<<endl
			 <<" Warning keyword FRAC_CONT_ANG is ignored"<<endl
			 <<" because keyword FRAC_CONT_OPT is missing"<<endl
			 <<"******************************************\n"<<endl;

	if(inp.giv("READ_ALTR_CA",data,0))  readSetCAs(data, elemans, nBpPors, out_);
	else if(warn)
	{
		out_<<"\n*******************************"<<endl;
		out_<<"EQUIL_CON_ANG  not found"<<endl;
		out_<<"Leaving contact angles unchanged"<<endl;
		out_<<"********************************"<<endl;
	}
}


void setElemProps(const InputFile& inp, const vector<Elem*>& elemans, size_t nBpPors, mstream& out_, const GNMData& comn)
{
	if (comn.dispCycle()==1)
		applyInitWettability(inp, elemans, nBpPors, out_, comn);
	else if (comn.dispCycle()==2)
		applyFWettabilityChange(inp, elemans, nBpPors, out_, comn);
}




