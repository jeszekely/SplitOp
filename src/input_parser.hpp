#ifndef DMRG_INPUT_PARSER
#define DMRG_INPUT_PARSER

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

class programInputs
{
public:
	std::string filename;
	bool DP_Verbose;
	bool TC_BasicMatrix_Run;
	bool TC_Davidson_Run;
	int TC_Davidson_MatrixSize;
	int TC_Davidson_nVecs;
	double TC_Davidson_tolerance;
	bool TC_DMRGSpin_Run;
	int TC_DMRGSpin_TestTrunc1;
	int TC_DMRGSpin_TestTrunc2;
	int TC_DMRGSpin_TestTrunc3;
	int TC_DMRGSpin_ChainLen;

	programInputs(const std::string fn) : filename(fn)
	{
		boost::property_tree::ptree mainIP;
		boost::property_tree::json_parser::read_json(filename,mainIP);

//		[DavidsonParams]
		DP_Verbose = mainIP.get<bool>("DavidsonParams.Verbose",false);

//		[TestCode]
		TC_BasicMatrix_Run = mainIP.get<bool>("TestCode.BasicMatrix.Run",false);

		TC_Davidson_Run = mainIP.get<bool>("TestCode.Davidson.Run",false);
		TC_Davidson_MatrixSize = mainIP.get<int>("TestCode.Davidson.MatrixSize",500);
		TC_Davidson_nVecs = mainIP.get<int>("TestCode.Davidson.nVecs",2);
		TC_Davidson_tolerance = mainIP.get<double>("TestCode.Davidson.tolerance",1.0e-6);

		TC_DMRGSpin_Run = mainIP.get<bool>("TestCode.DMRGSpin.Run",false);
		TC_DMRGSpin_ChainLen = mainIP.get<int>("TestCode.DMRGSpin.ChainLen",10);
		TC_DMRGSpin_TestTrunc1 = mainIP.get<int>("TestCode.DMRGSpin.TestTrunc1",10);
		TC_DMRGSpin_TestTrunc2 = mainIP.get<int>("TestCode.DMRGSpin.TestTrunc2",10);
		TC_DMRGSpin_TestTrunc3 = mainIP.get<int>("TestCode.DMRGSpin.TestTrunc3",10);
	};
};
#endif