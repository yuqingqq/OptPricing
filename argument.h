#pragma once


class Argument
{
public:
	int _func = 1; // Function parameter. 0: format graph, 1: maximize profit, 2: format graph and then maximize profit
	std::string _graphname = "big_twitter"; // Graph name. Default is "facebook".
	/*
	 * three letters: 
	 * 1 -> e: edge file (default), a: adjacent vector,
	 * 2 -> g: graph only (default, WC cascade model), w: with edge property,
	 * 3 -> r: reverse (default), f: forward, b: bilateral
	 */
	std::string _mode = "egr";
//	std::string _dir = "D://InfMaxDistr//graphInfo"; // Directory
//	std::string _dir = "D://shapley value price//graphInfo"; // Directory
	std::string _dir = "F://Yuqing//graph"; // Directory
	std::string _outFileName; // File name of the result
	std::string _resultFolder = "result"; // Result folder. Default is "test".
	std::string _algName = "greedy"; // Algorithm. Default is greedy.
	std::string _probDist = "WC";// Probability distribution for diffusion model. Option: WC, TR. Default is WC.
	std::string _benefitDist = "uniform";// Probability distribution for diffusion model. Option: WC, TR. Default is WC.
	std::string _costDist = "degree"; // Cost type. Default is degree-proportional.
	double _scale = 10.0; // Scale factor: the ratio between total cost and total benefit. Default is 10.
	double _para = 0.0; // Cost parameter: the fraction of base cost. Default is 0 (pure degree-proportional).
	double _eps = 0.1; // Error allowed. Default is 0.5.
	int _sizek = 1000; // Number of seeds. Default is 50.
	int _numR = 10000000; // Number of RR sets to be generated. Default is 1000000.
	enum CascadeModel
	{
		IC,
		LT
	};

	CascadeModel _model = IC; // Diffusion models: IC, LT. Default is IC.

	Argument(int argc, char* argv[])
	{
		std::string param, value;
		for (int ind = 1; ind < argc; ind++)
		{
			if (argv[ind][0] != '-') break;
			std::stringstream sstr(argv[ind]);
			getline(sstr, param, '=');
			getline(sstr, value, '=');
			if (!param.compare("-func")) _func = stoi(value);
			else if (!param.compare("-gname")) _graphname = value;
			else if (!param.compare("-mode")) _mode = value;
			else if (!param.compare("-dir")) _dir = value;
			else if (!param.compare("-outpath")) _resultFolder = value;
			else if (!param.compare("-alg")) _algName = value;
			else if (!param.compare("-pdist")) _probDist = value;
			else if (!param.compare("-bdist")) _benefitDist = value;
			else if (!param.compare("-cdist")) _costDist = value;
			else if (!param.compare("-scale")) _scale = stod(value);
			else if (!param.compare("-para")) _para = stod(value);
			else if (!param.compare("-eps")) _eps = stod(value);
			else if (!param.compare("-size")) _sizek = stoi(value);
			else if (!param.compare("-numR")) _numR = stoi(value);
			else if (!param.compare("-model")) _model = value == "LT" ? LT : IC;
		}
		_outFileName = TIO::get_out_file_name(_graphname);
		if (_model == LT) _outFileName = "LT_" + _outFileName;
	}

	std::string get_outfilename_with_alg(const std::string& algName) const
	{
		return TIO::get_out_file_name(_graphname);
	}
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
