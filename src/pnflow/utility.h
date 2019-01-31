

void loadRockTypeData(const string& fileName, vector<std::string>& poreTypes);


void setShapeFactFromFile(vector<Element*>& elems, const string& fileName, double lowCO, double hiCO);

void modifyInscribedRadii(int model, const string& options, vector<Element*>& elems,
							const string& fileName, bool writeToFile, int numPts, bool forPores = false);
void modifyShapeFactor(int model, const string& options, vector<Element*>& elems,
						const string& fileName, bool writeToFile, int numPts);
void setRadiiFromFile(vector<Element*>& elems, const string& fileName, double lowCO, double hiCO, bool forPores);

void setContactAngles(vector<Element*>& pores, vector<Element*>& throats, double min, double max,
	double delta, double eta,
							int wettingClass, double modelTwoSepAng, string  angDistScheme, bool attractToMarked = true );

inline double tableLookUpX(double yVal, const vector< pair<double, double> >& lookupTable)
{
    if(yVal <= lookupTable.front().second)
        return lookupTable.front().first;
    else if(yVal >= lookupTable.back().second)
        return lookupTable.back().first;

    vector< pair<double, double> >::const_iterator itr = lookupTable.begin()+1;
    while(itr->second < yVal && itr != lookupTable.end()) ++itr;
    assert(itr->second >= yVal && (itr-1)->second < yVal);
    double xOne((itr-1)->first), xTwo(itr->first), yOne((itr-1)->second), yTwo(itr->second);
    double logResX = log10(xOne)+log10(xTwo/xOne)*(yVal-yOne)/(yTwo-yOne);
    return pow(10.0, logResX);
}


inline double tableLookUpY(double xVal, const vector< pair<double, double> >& lookupTable)
{
    vector< pair<double, double> >::const_iterator itr = lookupTable.begin()+1;
    while(itr->first > xVal && itr != lookupTable.end()) ++itr;
    assert(itr->first <= xVal && (itr-1)->first > xVal);
    double xOne((itr-1)->first), xTwo(itr->first), yOne((itr-1)->second), yTwo(itr->second);
    return yOne+log10(xVal/xOne)*(yTwo-yOne)/log10(xTwo/xOne);
}

inline double weibull(double minVal, double maxVal, double delta, double eta) 
{
    double randNum = double(rand()) / double(RAND_MAX);
    //double randNum = 0.5; // delete me
    if(delta < 0.0 && eta < 0.0)                    // Uniform Distribution
        return minVal + (maxVal-minVal)*randNum;
    else                                            // Weibull Distribution
        return (maxVal - minVal) * pow(-delta*log(randNum*(1.0-exp(-1.0/delta))+exp(-1.0/delta)), 1.0/eta) + minVal;
}
