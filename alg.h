#pragma once

class Alg
{
private:
	int __numV, __numRRsets;
	std::vector<int> __vecSeed;
	const double* __pCost;
	PHyperGraphRef __pHyperG;
	PResult __pRes;
	std::vector<int> __candidate_set;
	std::vector<bool> __iscandidate;
	std::vector<double> __powerof2;
	int __prevsize;
	int _numexp;
	std::vector<double> _sigmai;

	double _step;

public:
	Alg(const double* pCost, const PHyperGraphRef& pHyperG, const PResult& pRes) : __pCost(pCost), __pHyperG(pHyperG), __pRes(pRes)
	{
		__numV = __pHyperG->get_nodes();
		__numRRsets = __pHyperG->get_RRsets_size();
		__iscandidate = std::vector<bool>(__numV, 0);
		__prevsize = 0;
		__prices_OPT = std::vector<double>(__numV, 0);
		_numexp = 1;
		_sigmai = std::vector<double>(__numV, 0);
		_step = 0.2;
	}
	std::vector<double> __prices_OPT;
	
	std::vector<double> _1st_value;
	std::vector<double> _objvalue;

	/// Build a set of n RR sets
	void build_n_RRsets(int64 numRRsets);
	// price for candidate seeds
	void OPTprice_of_seeds(double fC, double b);
	//set candidate set
	void set_candidateset(std::vector<int> candidates);

	double getfC();

	double Value_bmfC_OPT();

	double Value_bmfC();

	void obj_evaluate_1stpart(int numsamples, double sigmaC);

	void obj_2nd_sigmai();

	void obj_evaluate_2ndpart(std::vector<double> prices, double b);

	void Alg::stoppingrule(double epsilon, double delta);

	void record_div_opt(double fC);

	void record_div_deg(double fC,std::vector<double> degprice);

	void record_div_inf(double fC);

	void record_div_uni(double fC);

	void record_div_imrank(double fC);

};

using TAlg = Alg;
using PAlg = std::shared_ptr<TAlg>;

