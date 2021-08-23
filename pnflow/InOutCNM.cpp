



/// Creates all the pores and throats
void FlowDomain::initNetworkOld(InputData& input)  {
	/// The network is created and initialized by reading in the connection files that also contains
	/// all other relevant information about the pores/throats. Pointers to all elements are contained
	/// in a vector. Elem 0 and 1 are inlet and outlet. Throats follows after the pores.
	/// Since all rock elements contain pointers to connecting elements rather than vector indicies,
	/// initialization has to be in correct order: throats, pores, in/outlet and finally finishing the
	/// throats. 


		const int     DUMMY_IDX = -99;
		double keepFraction = input.getOr("A_CLOSE_SHAVE",1.);
		input.network(nPors_, nTrots_, box_.x, box_.y, box_.z); //elemData=convertNetworkFromOldFormat(input, nPors_, numThroats, box_);

		vector< pair<int,double> > insidePoreHashs(nPors_);
		int neoNPors(0);// WARNING SOURCE OF ERRORR
		for(int index = 1; index <= nPors_; ++index)  {
			pair<int,double> entry(DUMMY_IDX, 0.);
			input.poreLocation(index, entry.second);
			if(entry.second >= box_.x*(1.-keepFraction)*0.5 && entry.second <= box_.x-box_.x*(1.-keepFraction)*0.5)  {
				entry.first = ++neoNPors;
			}
			insidePoreHashs[index-1] = entry;
		}
		elemans_.resize(neoNPors + 2);

		vector< int > insideThroatHashs(nTrots_, DUMMY_IDX);
		vector< pair<int,int> > neoTrotNeis(nTrots_);
		vector<Elem*> throatsToInlet;								 // All throats connected to the in/outlet are
		vector<Elem*> throatsToOutlet;								// recorded while crating the throats

		int neoNTrots = 0;
		   //readAndCreateThroats(input, neoTrotNeis, throatsToInlet, throatsToOutlet, insidePoreHashs, insideThroatHashs, neoNPors);
		const double clayPorosityToAdd = input.getOr("CLAY_EDIT", 0.);


	/// The data for the throats are read from the link files. Since the pores are not yet created their
	/// indicies are stored in a temporary vector and the pointers will be initialized later. Since
	/// in/outlet does not have separete entries in the data files, we keep track of connecting throats.
	/// The strucure of the link files are as follows:
	///
	/// *_link1.dat:
	/// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
	///
	/// *_link2.dat:
	/// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume

	//int FlowDomain::readAndCreateThroats(InputData& input, vector< pair<int,int> >& neoTrotNeis,
									  //vector<Elem*>& throatsToInlet, vector<Elem*>& throatsToOutlet,
									  //const vector< pair< int, double> >& insidePoreHashs, vector< int >& insideThroatHashs, int neoNPors)  
	{
		cout<<"Reading throats"<<endl;

		int numLengthErrors(0);
		for(int iT = 0; iT < nTrots_; ++iT)  {

			int pore1Idx, pore2Idx;
			double radius, shapeFactor, lenTot, lenPore1, lenPore2, lenThroat, vol, clayVolume;

			input.throatData(iT+1, pore1Idx, pore2Idx, vol, clayVolume, radius, shapeFactor,
				lenPore1, lenPore2, lenThroat, lenTot);


			///. allow initializing pores/throats with zero-length/vol
			radius = max(radius,1e-32); vol = max(vol,1e-96); shapeFactor = max(shapeFactor,1e-32);
			lenThroat = max(lenThroat,1e-32); lenPore1 = max(lenPore1,1e-32); lenPore2 = max(lenPore2,1e-32);
			lenTot = max(lenTot,3.e-32);


			if(pore1Idx > 0)  neoTrotNeis[iT].first = insidePoreHashs[pore1Idx-1].first;
			else				neoTrotNeis[iT].first = pore1Idx;

			if(pore2Idx > 0)  neoTrotNeis[iT].second = insidePoreHashs[pore2Idx-1].first;
			else				neoTrotNeis[iT].second = pore2Idx;

			if(neoTrotNeis[iT].first > 0 || neoTrotNeis[iT].second > 0)  {
				++neoNTrots;
				insideThroatHashs[iT] = neoNTrots;



				clayVolume += clayPorosityToAdd*box_.x*box_.y*box_.z/(nTrots_+nPors_);

				if(fabs(lenPore1+lenPore2+lenThroat-lenTot)/(lenTot+lenPore1+lenPore2+lenThroat) > 0.01)	++numLengthErrors;



				Elem *throat = new Throat(comn_, neoNPors+1+neoNTrots, dbl3(0.,0.,0.), radius, vol, clayVolume, shapeFactor, lenThroat, lenPore1, lenPore2,0);


				elemans_.push_back(throat);

				if(neoTrotNeis[iT].first == DUMMY_IDX)  {
					if(insidePoreHashs[pore1Idx-1].second < box_.x*0.5)
						  neoTrotNeis[iT].first = -1;
					else  neoTrotNeis[iT].first = 0;
				}
				else if(neoTrotNeis[iT].second == DUMMY_IDX)  {
					if(insidePoreHashs[pore2Idx-1].second < box_.x*0.5)
						  neoTrotNeis[iT].second = -1;
					else  neoTrotNeis[iT].second = 0;
				}
				if(neoTrotNeis[iT].first == -1 || neoTrotNeis[iT].second == -1)
					throatsToInlet.push_back(throat);
				else if(neoTrotNeis[iT].first == 0 || neoTrotNeis[iT].second == 0)
					throatsToOutlet.push_back(throat);
			}
		}

		ensure(!numLengthErrors,
				"Warning: For "+ _s(numLengthErrors) + " throats the lengths of the\n"
				"pore-throat-pore did not match the total length.\n"
				"This is generally only an artifact of the network\n"
				"reconstruction process, and is not serious.\n");

	}



		{	cout<<"Reading pores\n";cout.flush();
			/// The pore data is read from the node files. At this point the throats are already created and the pointers
			/// can be set. The strucure of the node files are as follows:
			/// *_node1.dat:
			/// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
			/// *_node2.dat:
			/// index, volume, radius, shape factor, clay volume
			double  shaveOff(box_.x*(1.-keepFraction)*0.5);


			int neoIndex(1);
			for(int iP = 1; iP <= nPors_; ++iP) // WARNING SOURCE OF MISTAKE
			{
				int connNumber;
				double xPos, yPos, zPos, vol, radius, shapeFactor, clayVolume;
				vector< int > connThroats, connPores;
				vector<Elem*> adjTrots;

				input.poreData(iP, xPos, yPos, zPos, connNumber, connThroats, connPores, vol, clayVolume, radius, shapeFactor);

				if(xPos >= shaveOff && xPos <= box_.x-shaveOff)  {
					++neoIndex;

					//double adjustingVol(input.clayEdit()*(vol+clayVolume));
					//adjustingVol = min(adjustingVol, vol);
					//adjustingVol = -min(-adjustingVol, clayVolume);
					//vol -= adjustingVol;
					//clayVolume += adjustingVol;
					clayVolume += clayPorosityToAdd*box_.x*box_.y*box_.z/(nTrots_+nPors_);

					adjTrots.resize(connNumber);
					for(int j = 0; j < connNumber; ++j)   adjTrots[j] = elemans_[ neoNPors+1+insideThroatHashs[connThroats[j]-1] ];

					double initSolvPrs = (outletSolverPrs_ + inletSolverPrs_)*0.5;

					dbl3 nod(xPos-shaveOff, yPos, zPos); ///. shifted to left by half of curtailed fraction
					//bool inSlvrBox(nod.isInsideBox(solverBoxStart_, solverBoxEnd_));
					//bool insideSatBox(nod.isInsideBox(solverBoxStart_, solverBoxEnd_));
					bool inSlvrBox(nod.x >= box_.x*solverBoxStart_ &&  nod.x <=  box_.x*solverBoxEnd_ && iP>=nBSs_);

					elemans_[neoIndex] = new Pore(comn_, neoIndex, nod, radius, vol, clayVolume, shapeFactor, inSlvrBox, inSlvrBox, initSolvPrs, adjTrots,0);
				}
			}
			cout<<" "<<endl;

		}

		input.clearNetworkData(); //functions above create a lot of additional working data not needed anymore, lets clean them to free the memory

		box_.x *= keepFraction;
		nPors_ = neoNPors;
		nBpPors_=nPors_+nBSs_;// not properly implemented yet

		double yMid(box_.y*0.5), zMid(box_.z*0.5);
		/// In and outlet only need to know what throats are connected to it. Their indicies are 0 and 1
		elemans_[0   ] = new InOutBoundary(comn_, 0   , dbl3( -1e-15,      yMid, zMid), throatsToInlet);
		elemans_[OutI] = new InOutBoundary(comn_, OutI, dbl3(box_.x+1e-15, yMid, zMid), throatsToOutlet);


		int runIdx(0);
		for(int i = 0; i < nTrots_; ++i)						// Adding the pore pointers to the throats had
		{															// to be delayed until after having created
			int poreIndex1 =  neoTrotNeis[i].first+1;// oren indices start from -1	  // the pores. The network should now be
			int poreIndex2 =  neoTrotNeis[i].second+1;	 // properly initialized.

			if(poreIndex1 >= 0 && poreIndex2 >= 0)  {
				ensure(poreIndex1 < nBpPors_ && poreIndex2 < nBpPors_);
				Elem* pore1 = elemans_[poreIndex1];
				Elem* pore2 = elemans_[poreIndex2];
				ensure(pore1 != NULL && pore2 != NULL);
				elemans_[nBpPors_+ runIdx]->node()=0.5*(pore1->node()+pore2->node());

				static_cast<Throat*>(elemans_[nBpPors_+ runIdx])->addConnections(pore1, pore2, box_.x*solverBoxStart_, box_.x*solverBoxEnd_, !useAvrPrsAsBdr_, true);
				++runIdx;
			}
		}
		nTrots_ = neoNTrots;





}



/**
// The throat data is written to file in following format:
//
// *_link1.dat (outOne):
// index, pore 1 index, pore 2 index, radius, shape factor, total length (pore center to pore center)
//
// *_link2.dat (outTwo):
// index, pore 1 index, pore 2 index, length pore 1, length pore 2, length throat, volume, clay volume
*/
void Throat::writeNetworkData(ostream& outOne, ostream& outTwo) const
{
	outOne.flags(ios::showpoint);
	outOne.flags(ios::scientific);
	outTwo.flags(ios::showpoint);
	outTwo.flags(ios::scientific);
	double lenPoreOne(poreLength_[0]), lenPoreTwo(poreLength_[1]), lenThroat(length_);
	double lenTotal(poreLength_[0]+poreLength_[1]+length_);

	if(cnctions_.size() != 2)  {
		cerr << endl
			<< "============================================" << endl
			<< "For optimized network to be written to file " << endl
			<< "the option to drain singlets must be enabled" << endl
			<< "============================================" << endl;        exit(-1);
	}

	outOne << setw(7)   << indexOren()
		<< setw(7)      << cnctions_[0]->indexOren()
		<< setw(7)      << cnctions_[1]->indexOren()
		<< setw(15)     << model_->RRR()
		<< setw(15)     << model_->shapeFactor()
		<< setw(15)     << lenTotal
		<< endl;

	outTwo << setw(7)   << indexOren()
		<< setw(7)      << cnctions_[0]->indexOren()
		<< setw(7)      << cnctions_[1]->indexOren()
		<< setw(15)     << lenPoreOne
		<< setw(15)     << lenPoreTwo
		<< setw(15)     << lenThroat
		<< setw(15)     << flowVolume_* (rockIndex_>0 ? 1./model_->porosity() : 1.)
		<< setw(15)     << clayVolume_
		<< endl;
}

void Throat::writeNetworkDataBinary(ostream& out) const
{
	double lenPoreOne(poreLength_[0]), lenPoreTwo(poreLength_[1]), lenThroat(length_);
	double lenTotal(poreLength_[0]+poreLength_[1]+length_);

	ThroatStruct pr;
	pr.index = indexOren();
	pr.poreOne = cnctions_[0]->indexOren();
	pr.poreTwo = cnctions_[1]->indexOren();

	pr.radius = model_->RRR();
	pr.shapeFact = model_->shapeFactor();
	pr.lenPoreOne = lenPoreOne;
	pr.lenPoreTwo = lenPoreTwo;
	pr.lenThroat = lenThroat;
	pr.lenTot = lenTotal;
	pr.volume = flowVolume_* (rockIndex_>0 ? 1./model_->porosity() : 1.);
	pr.clayVol = clayVolume_;
	out.write((char *)(&pr), sizeof(pr));

} 
 
 
 

/**
// The pore data is written to file in following format:
//
// *_node1.dat (outOne):
// index, x_pos, y_pos, z_pos, connection num, connecting nodes..., at inlet?, at outlet?, connecting links...
//
// *_node2.dat (outTwo):
// index, volume, radius, shape factor, clay volume
*/
void Pore::writeNetworkData(ostream& outOne, ostream& outTwo) const
{
	size_t i;
	bool connToIn(false), connToOut(false);
	outOne.flags(ios::showpoint);
	outOne.flags(ios::scientific);
	outTwo.flags(ios::showpoint);
	outTwo.flags(ios::scientific);

	outOne << setw(7)  << indexOren()
		<< node_
		<< setw(5) << nCncts_;

	for(i = 0; i < cnctions_.size(); ++i)  {
		const Throat* throat = dynamic_cast< Throat* >(cnctions_[i]);
		assert(throat);
		const Elem* prj = throat->neighbouringPore(this);
		if(prj->isEntryRes()) connToIn = true;
		if(prj->isExitRes()) connToOut = true;

		outOne << setw(7) << prj->indexOren();                                 // Connecting nodes
	}

	outOne << setw(7) << connToIn << setw(7) << connToOut;                          // In and outlet?

	for(i = 0; i < cnctions_.size(); ++i)
		outOne << setw(7) << cnctions_[i]->indexOren();                     // Connecting throats

	outOne << endl;

	outTwo << setw(7) << indexOren()
		<< setw(15) << flowVolume_* (rockIndex_>0 ? 1./model_->porosity() : 1.)
		<< setw(15) << model_->RRR()
		<< setw(15) << model_->shapeFactor()
		<< setw(15) << clayVolume_
		<< endl;
}

void Pore::writeNetworkDataBinary(ostream& out) const
{
	PoreStruct pr;
	pr.index = indexOren();
	pr.x = node_.x;
	pr.y = node_.y;
	pr.z = node_.z;
	pr.connNum = nCncts_;

	pr.radius = model_->RRR();
	pr.shapeFact = model_->shapeFactor();
	pr.volume = flowVolume_* (rockIndex_>0 ? 1./model_->porosity() : 1.);
	pr.clayVol = clayVolume_;
	out.write((char *)(&pr), sizeof(pr));
	for(size_t i = 0; i < cnctions_.size(); ++i)  {
		const Throat* throat = dynamic_cast< Throat* >(cnctions_[i]);
		assert(throat);
		int idx = throat->neighbouringPore(this)->indexOren();
		out.write((char *)(&idx), sizeof(int));
	}
	for(size_t j = 0; j < cnctions_.size(); ++j)  {
		int idx = cnctions_[j]->indexOren();
		out.write((char *)(&idx), sizeof(int));
	}

} 
 

void InputData::writeNetwork(bool& writeNet, bool& writeBin, string& netName) const
{
	istringstream data;
	string keyword("WRITE_NET");
	char inBinary;

	if(giv(keyword, data))  {
		//if (verbose) cout<< "Reading " << keyword << endl;
		data >> inBinary >> netName;
		writeBin = (inBinary == 'T' || inBinary == 't');
		writeNet = true;
		checkEndOfData(data, keyword);
  }
	else
	{
		writeNet = false;
		writeBin = false;
		netName = "nothing";
	}
}





/**
* Writes the optimized network to file
*/
void FlowDomain::writeNetworkToFile(const InputData& input) const
{

	bool writeNetToFile(false), writeInBinary(false);
	string fileNameBase;

	input.writeNetwork(writeNetToFile, writeInBinary, fileNameBase);

	if(!writeNetToFile) return;


	//if(!drainSinglets)
	//{
	cout<< "==================================================" << endl
		<< "CATION: The option not to drain dangling ends	   " << endl
		<< "DRAIN_SINGLETS not be used together with the option" << endl
		<< "to write the optimized network to file." << endl
		<< "==================================================" << endl;
	//exit(-1);
	//}

	if(writeInBinary)  {
		string poreFileName(fileNameBase + "_node.bin"), throatFileName(fileNameBase + "_link.bin");
		ofstream pOut(poreFileName, ios::binary), tOut(throatFileName, ios::binary);
		ensure(pOut && tOut, "Warning: Could not open " + poreFileName + "for writing. Network is not written to file." );

		pOut.write((char *)(&nPors_), sizeof(int));
		pOut.write((char *)(&box_.x), sizeof(double));
		pOut.write((char *)(&box_.y), sizeof(double));
		pOut.write((char *)(&box_.z), sizeof(double));
		tOut.write((char *)(&nTrots_), sizeof(int));

		for (int i = nBSs_; i < nBpPors_; ++i)
			elemans_[i]->writeNetworkDataBinary(pOut);

		for (size_t j = nBpPors_; j < elemans_.size(); ++j)
			elemans_[j]->writeNetworkDataBinary(tOut);

		pOut.close();
		tOut.close();
	}
	else
	{
		string pOut1FileName(fileNameBase + "_node1.dat"), pOut2FileName(fileNameBase + "_node2.dat");

		ofstream pOut1, pOut2;
		pOut1.open(pOut1FileName);
		pOut2.open(pOut2FileName);
		pOut1.flags(ios::showpoint);
		pOut1.flags(ios::scientific);
		if(!pOut1 || !pOut2)  {
			cout<< endl
				<< "=====================================================" << endl
				<< "Warning: Could not open " << pOut1FileName << endl
				<< "for writing. Network is not written to file." << endl
				<< "=====================================================" << endl
				<< endl;
		}

		pOut1 << nPors_ << "   " << box_.x << "   " << box_.y << "   " << box_.z << '\n';

		for (int i = nBSs_; i < nBpPors_; ++i)
			elemans_[i]->writeNetworkData(pOut1, pOut2);

		pOut1.close();
		pOut2.close();

		string tOut1FileName(fileNameBase + "_link1.dat"), tOut2FileName(fileNameBase + "_link2.dat");
		ofstream tOut1, tOut2;
		tOut1.open(tOut1FileName);
		tOut2.open(tOut2FileName);

		tOut1 << nTrots_ << "   " << '\n';

		for (size_t j = nBpPors_; j < elemans_.size(); ++j)
			elemans_[j]->writeNetworkData(tOut1, tOut2);

		tOut1.close();
		tOut2.close();
	}
}


void FlowDomain::createMatlabLocationData() const
{
	ofstream outp("poreLocation.m");
	outp.flags(ios::scientific);
	outp.precision(3);
	outp << "% The physical dimensions (m) of the network is: " << endl;
	outp << "modSize(1) = " << box_.x << ";" << endl;
	outp << "modSize(2) = " << box_.y << ";" << endl;
	outp << "modSize(3) = " << box_.z << ";" << endl;

	outp << endl << "% The location of individual pores (m). The array is 1-based" << endl;
	outp << "poreLoc = [";
	for(int i = nBSs_; i < nBpPors_; ++i)  {
		outp << elemans_[i]->node().x << ", "
			<< elemans_[i]->node().y << ", "
			<< elemans_[i]->node().z << "; ..."
			<< endl;
	}
	outp << "];" << endl;
	outp.close();

	ofstream outt("throatConnection.m");
	outt << "% Pores to which the throats are connected. The indexing is" << endl
		<< "% 1-based, which is the same as that used by the Oren format" << endl;

	outt << "throatConn = [";
	for(size_t j = nBpPors_; j < elemans_.size(); ++j)  {
		outt << elemans_[j]->neib(0)->indexOren() << ", ";
		if(elemans_[j]->nCncts() == 2)
			outt << elemans_[j]->neib(1)->indexOren() << "; " << endl;
		else
			outt << elemans_[j]->neib(0)->indexOren() << "; " << endl;
	}
	outt << "];" << endl;
	outt.close();
}

