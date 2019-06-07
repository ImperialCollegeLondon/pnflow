#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <map>
using namespace std;

#include "inputData.h"
#include "Element.h"
#include "solver.h"




#include "netsim.h"




void Netsim::recordRes(bool wantRelPerm, bool wantResIdx)
{
	m_resultWaterSat.push_back(m_satWater);
	m_resultCappPress.push_back(m_cappPress);
	//if(m_apexPrsReported)
	//{
	//m_resultBoundingPc.push_back(m_boundPress);
	//}
	if(wantRelPerm)
	{
		pair<double, double> watElem(m_watFlowRate, m_deltaPw), oilElem(m_oilFlowRate, m_deltaPo);
		m_resultWaterFlowRate.push_back(watElem);
		m_resultOilFlowRate.push_back(oilElem);
	}
	if(wantResIdx)
	{

		double resistivityIdx = m_singlePhaseCurrent / (m_current + 1.0e-100);
		if(m_current > 0.0)
		{
			if(m_useAvrPrsAsBdr)
			{
				resistivityIdx = (m_singlePhaseCurrent * m_deltaV) / ((m_current + 1.0e-100) * m_singlePhaseDvolt);
			}
		}

		m_resultResistivityIdx.push_back(resistivityIdx);
	}
	if(m_reportMaterialBal)
	{
		m_resultWaterMass.push_back((m_totalFlowVolume + m_totalClayVolume)*m_satWater*m_water.density());
		m_resultOilMass.push_back((m_totalFlowVolume + m_totalClayVolume)*(1 - m_satWater)*m_oil.density());
	}

	m_vtkWriter.vtuWrite(reinterpret_cast<const std::vector<Element const *>*>(&m_rockLattice), m_numPores, m_cappPress, m_oil.interfacialTen());

}

//JIA
void Netsim::formatResults(bool matlabFormat, bool excelFormat, bool mcpFormat)
{
	bool wantRelPerm(m_resultWaterFlowRate.size() > 0);
	bool wantResIdx(m_resultResistivityIdx.size() > 0);

	string whiteSpace;
	if(matlabFormat) whiteSpace = ", ";
	else if(mcpFormat) whiteSpace = "\t";
	else if(excelFormat) whiteSpace = ",";
	else whiteSpace = "	";

	double krwRefFlow(m_singlePhaseWaterQ), kroRefFlow(m_singlePhaseOilQ);
	double deltaPwRef(m_singlePhaseDprs), deltaPoRef(m_singlePhaseDprs);
	string m_relPermDef = "single";
	m_input.relPermDef(m_relPermDef);
	if(wantRelPerm && (m_relPermDef[2] == 's' || m_relPermDef[2] == 'S'))   // Scale wrt residual cond
	{
		double maxQW(0.0), maxQO(0.0), delPw(0.0), delPo(0.0);
		for (size_t el = 0; el < m_resultWaterSat.size(); ++el)
		{
			if(m_resultWaterFlowRate[el].first > maxQW)
			{
				maxQW = m_resultWaterFlowRate[el].first;
				delPw = m_resultWaterFlowRate[el].second;
			}
			if(m_resultOilFlowRate[el].first > maxQO)
			{
				maxQO = m_resultOilFlowRate[el].first;
				delPo = m_resultOilFlowRate[el].second;
			}
		}
		krwRefFlow = maxQW;
		kroRefFlow = maxQO;
		deltaPwRef = delPw;
		deltaPoRef = delPo;
	}
	else if(wantRelPerm && m_relPermDef[2] != 'n' && m_relPermDef[2] == 'N')
	{
		cout<< "========================================================" << endl
			<< "Warning. Uknown reference flag for relative permeability" << endl
			<< "calculations. Use 'single' or 'residual'. Defaulting" << endl
			<< "to single phase reference." << endl
			<< "========================================================" << endl;
	}

	for (size_t entry = 0; entry < m_resultWaterSat.size(); ++entry)
	{
		ostringstream out;
		out.flags(ios::showpoint);
		out.flags(ios::fixed);

		out << m_resultWaterSat[entry] << whiteSpace;
		out.flags(ios::scientific);
		out << setw(15) << m_resultCappPress[entry] << whiteSpace;


		if(wantRelPerm)
		{
			double krw(0.0), kro(0.0);
			if(m_resultWaterFlowRate[entry].first > 0.0)
			{
				if(m_useAvrPrsAsBdr)
					krw = m_resultWaterFlowRate[entry].first*deltaPwRef / (krwRefFlow*m_resultWaterFlowRate[entry].second);
				else
					krw = m_resultWaterFlowRate[entry].first / krwRefFlow;
			}
			if(m_resultOilFlowRate[entry].first > 0.0)
			{
				if(m_useAvrPrsAsBdr)
					kro = m_resultOilFlowRate[entry].first*deltaPoRef / (kroRefFlow*m_resultOilFlowRate[entry].second);
				else
					kro = m_resultOilFlowRate[entry].first / kroRefFlow;
			}
			out << setw(15) << krw << whiteSpace << setw(15) << kro << whiteSpace;
		}
		else if(mcpFormat)
			m_out << "Error microporosity format requires calc_perm to be true (See SAT_CONTROL keyword)";

		if(wantResIdx)
		{
			out << setw(15) << m_resultResistivityIdx[entry] << whiteSpace;
		}
		else if(mcpFormat)
			out << 1 << whiteSpace;

		if(m_reportMaterialBal && !mcpFormat)
		{
			out << setw(15) << m_resultWaterMass[entry] << whiteSpace
				<< setw(15) << m_resultOilMass[entry] << whiteSpace;
		}


		if(matlabFormat) out << "; ...";

		m_results.push_back(out.str());
	}
}

