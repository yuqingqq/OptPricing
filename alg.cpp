#include "stdafx.h"

void::Alg::set_candidateset(std::vector<int> candidates) {

	__iscandidate = std::vector<bool>(__numV, 0);
	std::vector<double>().swap(__powerof2);
	for (int i = 0; i < candidates.size();i++) {
		__candidate_set.push_back(candidates[i]);
		__iscandidate[candidates[i]] = 1;
	}
	__pHyperG->set_candidate_set(__candidate_set);
	double po0 = 1;
	for (int i = 0; i <= __candidate_set.size(); i++) {
		__powerof2.push_back(po0);
		po0 /= 2;
	}
}

void Alg::build_n_RRsets(int64 numRRsets)
{

	__pHyperG->build_n_RRsets(numRRsets);
	__numRRsets = __pHyperG->get_RRsets_size();
	__pRes->set_RRsets_size(__numRRsets);
}

double Alg::getfC() {
	double fC = 0;
	for (auto RRcand : __pHyperG->__candsinRR) {
		if (RRcand.size()) {
			fC++;
		}
	}
	fC = fC / __pHyperG->get_RRsets_size()  * __numV;
	return fC;
}


void Alg::OPTprice_of_seeds(double fC, double b) {
	int m = __candidate_set.size();
	__numRRsets = __pHyperG->get_RRsets_size();
	double bmfC = 0;
	std::vector<double> prices(__numV, 0);
	for (int i = 0; i < __numRRsets; i++) {
		int Nr = __pHyperG->__candsinRR[i].size();
		if (Nr == 0) continue;
	//	bmfC += 2.0 / (1 + m) * (m - (m - Nr) * __powerof2[Nr]) - 1;
	//	double cons = 1.0 / m * (1 - __powerof2[Nr - 1] * Nr);
	//	double cons = 2.0 / (1 + m) * (1 - (1 + Nr) * __powerof2[Nr]);
		double cons = 1.0 / m * (b/fC - __powerof2[Nr - 1] * Nr);
		std::vector<bool> inRR(__numV);
		for (auto candinRR : __pHyperG->__candsinRR[i]) {
			inRR[candinRR] = 1;
			prices[candinRR] += cons;
			prices[candinRR] += __powerof2[Nr - 1];
		}
		for (auto cand : __candidate_set) {
			if (inRR[cand]) continue;
			prices[cand] += cons;
		}
	}
//	bmfC = bmfC / __numRRsets * __numV;
	for (auto&candi : __candidate_set) {
//		__prices_OPT[candi] = prices[candi]  / __numRRsets* __numV+bmfC/m;
		__prices_OPT[candi] = prices[candi] / __numRRsets * __numV;
	}
}

double Alg::Value_bmfC_OPT(){
	int numrr = __pHyperG->get_RRsets_size();
	double objvalue = 0;
	int nc = __candidate_set.size();
	for (int irr = 0; irr < numrr; irr++) {
		int Nr = __pHyperG->__candsinRR[irr].size();
		if (Nr == 0) continue;
		objvalue += 2.0/(1+nc)*(nc-(nc-Nr)*__powerof2[Nr])-1;
	}
	objvalue = objvalue * __numV / numrr;
	return objvalue;
}

double Alg::Value_bmfC() {
	int numrr = __pHyperG->get_RRsets_size();
	double objvalue = 0;
	for (int irr = 0; irr < numrr; irr++) {
		int Nr = __pHyperG->__candsinRR[irr].size();
		if (Nr == 0) continue;
		objvalue += 1 - __powerof2[Nr - 1];
	}
	objvalue = objvalue * __numV / numrr;
	return objvalue;
}

void Alg::obj_evaluate_1stpart(int numsamples, double sigmaC) {
	std::cout << "evaluating 1st part" << std::endl;
	//int numofrrsets = 1000000;
	//__pHyperG->refresh_hypergraph();
	//__pHyperG->build_n_RRsets(numofrrsets);
	int numofrrsets = __pHyperG->get_RRsets_size();
	std::vector<bool> edgeMark(numofrrsets, false);
	_1st_value = std::vector<double>(_numexp);
	for (int iexp = 0; iexp < _numexp; iexp++) {
		std::vector<std::vector<int>> sets(numsamples);
		for (int iset = 0; iset < numsamples; iset++) {
			std::vector<int> randset;
			for (int ibit = 0; ibit < __candidate_set.size(); ibit++) {
				double randDouble = dsfmt_gv_genrand_open_close();
				if (randDouble > 0.5) {
					randset.push_back(__candidate_set[ibit]);
				}
			}
			sets.push_back(randset);
		}
		double sum_sigmas_squ = 0;
		for (auto subset : sets) {
			double real_inf = 0;
			edgeMark = std::vector<bool>(numofrrsets, false);
			for (auto candi : subset) {
				for (auto edgeIdx : __pHyperG->_vecCover[candi]) {
					if (edgeMark[edgeIdx]) continue;
					edgeMark[edgeIdx] = true;
				}
			}
			real_inf = 1.0*std::count(edgeMark.begin(), edgeMark.end(), true) * __numV / numofrrsets;
			sum_sigmas_squ += pow(real_inf / sigmaC, 2);
		}
		_1st_value[iexp] = sum_sigmas_squ;
		std::cout << "exp " << iexp << " with 1st part = " << sum_sigmas_squ << std::endl;
	}
	std::ofstream fout;
	std::string pricing_error = "first_part.txt";
	fout.open(pricing_error, std::ios::app);
	for (int iexp = 0; iexp < _numexp; iexp++) {
		_1st_value[iexp] *= pow(sigmaC, 2) / numsamples;
		//__obj1stpart[iexp] += beta * pow(sigmaC, 2);
		fout <<_1st_value[iexp] << std::endl;
		std::cout << "first part = " << _1st_value[iexp] << std::endl;
	}
	fout << std::endl;
	fout.close();
	std::cout << "1st part calculation done" << std::endl;
	obj_2nd_sigmai();
}

void Alg::obj_2nd_sigmai() {
	int numrr = __pHyperG->get_RRsets_size();
	for (int irr = 0; irr < numrr; irr++) {
		std::vector<bool> isInrr(__numV, 0);
		int Nr = __pHyperG->__candsinRR[irr].size();
		if (Nr == 0) continue;
		for (auto& candiinR : __pHyperG->__candsinRR[irr]) {
			isInrr[candiinR] = 1;
		}
		for (auto& candi : __candidate_set) {
			if (isInrr[candi]) {
				_sigmai[candi] += 1;
			}
			else {
				_sigmai[candi] += 1 - __powerof2[Nr];
			}
		}
	}
	for (auto& candi : __candidate_set) {
		_sigmai[candi] *= 1.0 * __numV / numrr;
		//std::cout << "sigma " << candi << " = " << sigmai[candi] << std::endl;
	}
	std::cout << "calculated sigma for each seed" << std::endl;
}

void Alg::obj_evaluate_2ndpart(std::vector<double> prices, double b) {
//	std::cout << " using sigmai to compute 2nd part" << std::endl;
	double obj = 0;
	for (auto& candi : __candidate_set) {
		obj += 0.25 * prices[candi] * (prices[candi] + b) - prices[candi] * _sigmai[candi];
	}
	std::cout << "2 nd part is " << obj << '\n';
	_objvalue = std::vector<double>(_numexp);
	for (int iexp = 0; iexp < _numexp; iexp++) {
		_objvalue[iexp] = obj + _1st_value[iexp];
		std::cout << "total is " << _objvalue[iexp] << '\n';
	}
	std::ofstream fout;
	std::string pricing_error = "pricing_error.txt";
	fout.open(pricing_error, std::ios::app);
	for (int iexp = 0; iexp < _numexp; iexp++) {
		std::cout << _objvalue[iexp] << std::endl;
		fout << _objvalue[iexp] << std::endl;
	}
	fout << std::endl;
	fout.close();
}


void Alg::stoppingrule(double epsilon, double delta) {
	Timer t("");
	double ups = 1 + epsilon + (1 + epsilon) * (2 + 2.0 / 3 * epsilon) * log(2.0 / delta) / pow(epsilon, 2);
	std::cout << "ups = " << ups << std::endl;
	int nc = __candidate_set.size();
	double minps = 0;
	int theta = 0;
	std::vector<double> prices(__numV);
	double cons = 2.0 / (1 + nc);
	while (minps < ups) {
		theta = __pHyperG->get_RRsets_size();
	//	if (theta == 1000000) break;
		__pHyperG->build_n_RRsets1by1();
		int nr = __pHyperG->__candsinRR[theta].size();
		if (nr == 0) continue;
		std::vector<bool> inRR(__numV, 0);
		for (auto candinRR : __pHyperG->__candsinRR[theta]) {
			inRR[candinRR] = 1;
		}
		for (auto candi : __candidate_set) {
		//	prices[candi] += 1.0 / nc * (1 - nr * __powerof2[nr - 1]);
			prices[candi] += cons;
			if (inRR[candi]) {
				prices[candi] += __powerof2[nr - 1]*(nc-nr)/(1+nc);
			}
			else {
				prices[candi] -= __powerof2[nr - 1] * (nr + 1) / (nc + 1);
			}
		}
		minps = ups + 1;
		for (auto candi : __candidate_set) {
			minps = min(minps, prices[candi]);
		}
		//std::cout << " minps = " << minps << std::endl;
	}
	std::cout << "build " << theta << " samples \n";
	double costtime = t.get_total_time();
	std::cout << "cost time " << costtime << '\n';
	for (auto candi : __candidate_set) {
		__prices_OPT[candi] = 1.0 * __numV / theta * prices[candi];
	}
	std::ofstream fout;
	std::string tandr = "time.txt";
	fout.open(tandr, std::ios::app);
	fout <<theta<< std::endl;
	fout << costtime << '\n';
	fout.close();
}

void Alg::record_div_opt(double fC) {
	std::ofstream fout;
	std::string valueb = "b.txt";
	fout.open(valueb, std::ios::app);
	int m = __candidate_set.size();
	std::vector<double> newoptprice(__numV);

	double bmfC_OPT = Value_bmfC_OPT();

	for (int i = 0; i < 6; i++) {
		double b = fC + _step * i * fC;
		std::cout << "b is " << b << '\n';
		for (auto cand : __candidate_set) {
			newoptprice[cand] = __prices_OPT[cand]+(b - fC - bmfC_OPT) / m;
		}
		obj_evaluate_2ndpart(newoptprice, b);
		fout << b << '\n';
	}
	double bmfC = Value_bmfC();
	std::cout << " set b is " << bmfC + fC << '\n';
	for (auto cand : __candidate_set) {
		newoptprice[cand] = __prices_OPT[cand] + (bmfC-bmfC_OPT) / m;
	}
	obj_evaluate_2ndpart(newoptprice, fC + bmfC);
	fout << bmfC + fC << '\n';
	std::cout << "opt b is " << bmfC_OPT + fC << '\n';
	for (auto cand : __candidate_set) {
//		newoptprice[cand] = __prices_OPT[cand] + bmfC_OPT / m;
		newoptprice[cand] = __prices_OPT[cand];
	}
	obj_evaluate_2ndpart(newoptprice, fC + bmfC_OPT);
	fout << bmfC_OPT + fC << '\n';
	fout << std::endl;
	fout.close();
}

void Alg::record_div_deg(double fC, std::vector<double> degprice) {
	double sumprice = 0;
	for (auto cand:__candidate_set) {
		sumprice += degprice[cand];
	}
	std::ofstream fout;
	std::string valueb = "b.txt";
	fout.open(valueb, std::ios::app);
	std::vector<double> newdegprice(__numV);
	for (int i = 0; i < 6; i++) {
		double b = fC + _step * i * fC;
		std::cout << "b is " << b << '\n';
		for (auto cand : __candidate_set) {
			newdegprice[cand] = degprice[cand] / sumprice * b;
		}
		obj_evaluate_2ndpart(newdegprice, b);
		fout << b << '\n';
	}
}

void Alg::record_div_uni(double fC) {
	std::vector<double> uniprice(__numV);
	std::ofstream fout;
	std::string valueb = "b.txt";
	fout.open(valueb, std::ios::app);
	for (int i = 0; i < 6; i++) {
		double b = fC + _step * i * fC;
		std::cout << "b is " << b << '\n';
		for (auto cand : __candidate_set) {
			uniprice[cand] = b/__candidate_set.size();
		}
		obj_evaluate_2ndpart(uniprice, b);
		fout << b << '\n';
	}
}

void Alg::record_div_inf(double fC) {
	std::vector<double> infprice(__numV);
	int numRR = __pHyperG->get_RRsets_size();
	double sumprice = 0;
	for (auto cand : __candidate_set) {
		infprice[cand] = 1.0 * __pHyperG->_vecCover[cand].size() / numRR * __numV;
		sumprice += infprice[cand];
	}
	std::ofstream fout;
	std::string valueb = "b.txt";
	fout.open(valueb, std::ios::app);
	std::vector<double> newinfprice(__numV);
	for (int i = 0; i < 6; i++) {
		double b = fC + _step * i * fC;
		std::cout << "b is " << b << '\n';
		for (auto cand : __candidate_set) {
			newinfprice[cand] = infprice[cand] / sumprice * b;
		}
		obj_evaluate_2ndpart(newinfprice, b);
		fout << b << '\n';
	}
}

void Alg::record_div_imrank(double fC) {
	std::vector<double> rankprice(__numV);
	int maxDeg = 0;
	std::vector<int> coverage(__numV);
	for (auto cand : __candidate_set) {
		auto deg = __pHyperG->_vecCover[cand].size();
		coverage[cand] = deg;
		maxDeg = max(deg, maxDeg);
	}
	std::vector<std::vector<int>> degMap(maxDeg + 1);// degMap: map degree to the nodes with this degree
	for (auto cand : __candidate_set) {
		degMap[coverage[cand]].push_back(cand);
	}
	int numRR = __pHyperG->get_RRsets_size();
	std::vector<bool> edgeMark(numRR);
	for(auto deg= maxDeg;deg>0;deg--){
		auto listnodes = degMap[deg];
		for (int i = 0; i < listnodes.size(); i++) {
			int curNode = listnodes[i];
			int curDeg = coverage[curNode];
			if (deg > curDeg) {
				degMap[curDeg].push_back(curNode);
				continue;
			}
			rankprice[curNode] = 1.0 * curDeg /numRR * __numV;
			for (auto edgeId : __pHyperG->_vecCover[curNode]) {
				if (edgeMark[edgeId]) continue;
				edgeMark[edgeId] = true;
				for (auto nodeId : __pHyperG->_vecCoverRev[edgeId]) {
					if (coverage[nodeId] == 0) continue;
					coverage[nodeId]--;
				}
			}
		}
		degMap.pop_back();
	}
	double sumprice = 0;
	for (auto cand : __candidate_set) {
		sumprice += rankprice[cand];
	}
	std::ofstream fout;
	std::string valueb = "b.txt";
	fout.open(valueb, std::ios::app);
	std::vector<double> newrankprice(__numV);
	for (int i = 0; i < 6; i++) {
		double b = fC + _step * i * fC;
		std::cout << "b is " << b << '\n';
		for (auto cand : __candidate_set) {
			newrankprice[cand] = rankprice[cand] / sumprice * b;
		}
		obj_evaluate_2ndpart(newrankprice, b);
		fout << b << '\n';
	}
}

