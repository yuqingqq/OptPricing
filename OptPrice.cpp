
#include "stdafx.h"
#include<chrono>
#define NATURAL_E 2.71828
int numcandidates = 200;
double epsilon = 0.1;
double delta = 0.2;

int main(int argc, char* argv[])
{

	// Randomize the seed for generating random numbers
	dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
	TArgument Arg(argc, argv);
	std::string infilename = Arg._dir + "/" + Arg._graphname;
	if (Arg._func == 0 || Arg._func == 2)
	{// Format the graph
		GraphBase::format_graph(infilename, Arg._mode);
		if (Arg._func == 0) return 1;
	}

	std::cout << "---The Begin of " << Arg._outFileName << "---" << std::endl;
	Timer mainTimer("main");
	// Initialize a result object to record the results
	PResult pResult(new TResult());

	// Load the graph
	//infilename += std::to_string(Arg._probDist);
	Graph graph = GraphBase::load_graph(infilename, true, Arg._probDist);
	Graph dirgraph = GraphBase::load_graph(infilename, false, Arg._probDist);
	std::vector<std::pair<int, int>> nbridx;
	for (int i = 0; i < dirgraph.size(); i++) {
		nbridx.push_back({ i, dirgraph[i].size()});
	}
	std::sort(nbridx.begin(), nbridx.end(), larger());
	std::vector<int> candidates;
	//for (int i = 0; i < numcandidates; i++) {
	//	candidates.push_back(nbridx[i].first);
	//}
	for (int i = 0; i < numcandidates; i+=3) {
		candidates.push_back(nbridx[i].first);
	}
	int numV = (int)graph.size();

	std::vector<double> degprice(numV);
	for (auto cand : candidates) {
		degprice[cand] = dirgraph[cand].size();
	}

	if (Arg._model == TArgument::CascadeModel::LT)
	{// Normalize the propagation probabilities in accumulation format for LT cascade model for quickly generating RR sets
		to_normal_accum_prob(graph);
	}
	double* pCost = (double*)malloc(numV * sizeof(double));
	double* pAccumWeight = nullptr; // Used for non-uniform benefit distribution for activating nodes

	// Create a hyper-graph object to generate/compute RR sets
	PHyperGraphRef pHyperG(new THyperGraphRef(graph, pAccumWeight));
	pHyperG->set_cascade_model(static_cast<THyperGraphRef::CascadeModel>(Arg._model));
	//pHyperG->set_hyper_graph_mode(true);
	TAlg tAlg(pCost, pHyperG, pResult);
	std::cout << "  ==>Graph loaded for RIS! total time used (sec): " + std::to_string(mainTimer.get_total_time()) << std::endl;

	//set candidate set
	tAlg.set_candidateset(candidates);

	//initialize the parameters
	
	delta = 1.0 / numV;
	pHyperG->set_hyper_graph_mode(true);
	
	//tAlg.build_n_RRsets(Arg._numR);
	tAlg.stoppingrule(epsilon, delta);
	
	double fC = tAlg.getfC();

	std::cout << "fC is " << fC << '\n';

	std::ofstream fout;
	std::string fc2 = "beta_sigmac.txt";
	fout.open(fc2, std::ios::app);
	fout << 0.002*pow(fC,2) << std::endl;
	fout << std::endl;
	fout.close();

	int numseedsets = 14144;
	tAlg.obj_evaluate_1stpart(numseedsets, fC);

	tAlg.record_div_opt(fC);
	tAlg.record_div_inf(fC);
	tAlg.record_div_deg(fC, degprice);
	tAlg.record_div_imrank(fC);
	tAlg.record_div_uni(fC);



	std::cout << "---The End of " << Arg._outFileName << "---" << std::endl;
	
	return 0;
}
