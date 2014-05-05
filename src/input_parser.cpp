#include "input_parser.hpp"

using namespace std;

programInputs::programInputs(std::string fn) : filename(fn)
{
	stripComments();
	boost::property_tree::ptree IP;
	boost::property_tree::json_parser::read_json(filename,IP);

	procs   = IP.get<int>("CompCell.Procs",1);
	xmin    = IP.get<double>("CompCell.Xmin",0.0);
	xmax    = IP.get<double>("CompCell.Xmax",100.0);
	ymin    = IP.get<double>("CompCell.Ymin",0.0);
	ymax    = IP.get<double>("CompCell.Ymax",100.0);
	nx      = IP.get<int>("CompCell.Nx",1000);
	ny      = IP.get<int>("CompCell.Ny",100);
	zl      = IP.get<double>("CompCell",0.0);
	zr      = IP.get<double>("CompCell",24.57);
	nphonon = IP.get<int>("CompCell",5);
	nelec   = IP.get<int>("CompCell",3);
	eigen   = IP.get<int>("CompCell",1);

	dt      = IP.get<double>("TimeParameters.Dt",1.0);
	runtime = IP.get<double>("TimeParameters.Runtime",1000.0);
	output  = IP.get<double>("TimeParameters.Output",5.0);

	hbar       = IP.get<double>("PhysicalParameters.hbar",1.0);
	echarge    = IP.get<double>("PhysicalParameters.echarge",-1.0);
	m_electron = IP.get<double>("PhysicalParameters.M_electron",1.0);
	m_C60      = IP.get<double>("PhysicalParameters.M_C60",2.59e6);
	elecvel    = IP.get<double>("PhysicalParameters.ElecVel",0.488756);
	eaffin     = IP.get<double>("PhysicalParameters.EAffin",-2.56);
	lbarrier   = IP.get<double>("PhysicalParameters.LBarrier",120.0);
	rbarrier   = IP.get<double>("PhysicalParameters.RBarrier",120.0);
	xi         = IP.get<double>("PhysicalParameters.Xi",40.0);
	workfxn    = IP.get<double>("PhysicalParameters.WorkFxn",-5.40);
	bias       = IP.get<double>("PhysicalParameters.Bias",0.0);
	elecpos    = IP.get<double>("PhysicalParameters.ElecPos",-150.0);
	epsilon    = IP.get<double>("PhysicalParameters.Epsilon",1.09);
	requil     = IP.get<double>("PhysicalParameters.Requil",12.0);
	coupfxn    = IP.get<int>("PhysicalParameters.CoupFxn",1);
	rad        = IP.get<double>("PhysicalParameters.Rad",10.0);
	sigma      = IP.get<double>("PhysicalParameters.Sigma",1.0);
	del        = IP.get<double>("PhysicalParameters.Del",6.0);
	cconst     = IP.get<double>("PhysicalParameters.CConst",1.0);
	alpha      = 4.0;

	potential    = IP.get<int>("Outputs.Potential",0);
	wvfxninit    = IP.get<int>("Outputs.WvfxnInit",0);
	fullwvfxn    = IP.get<int>("Outputs.FullWvfxn",0);
	elecstates   = IP.get<int>("Outputs.ElecStates",0);
	phononstates = IP.get<int>("Outputs.PhonStates",0);
	animation    = IP.get<int>("Outputs.Animation",0);
	slice        = IP.get<int>("Outputs.Slice",0);
	sliceindex   = IP.get<int>("Outputs.SliceIndex",0);
	eleccorr     = IP.get<int>("Outputs.ElecCorr",0);
	moleccorr    = IP.get<int>("Outputs.MolecCorr",0);
	dcoeffr      = IP.get<int>("Outputs.DCoeffR",0);
	dcoeffl      = IP.get<int>("Outputs.DCoeffL",0);
	phases       = IP.get<int>("Outputs.Phases",0);
	energy       = IP.get<int>("Outputs.Energy",0);
	lifetime     = IP.get<int>("Outputs.Lifetime",0);
	zexp         = IP.get<int>("Outputs.ZExp",0);
	fluxr        = IP.get<int>("Outputs.FluxR",0);
	fluxl        = IP.get<int>("Outputs.FluxL",0);
	fluxint      = IP.get<int>("Outputs.FluxInt",0);
	correlation  = IP.get<int>("Outputs.Correlation",0);

	toAU();

  ofstream ss;
  ss.open("debug.txt");
  boost::property_tree::write_json(ss,IP);
  ss.close();
}

void programInputs::stripComments()
{
	//Open input and output file
	string newfn = "output_data/" + filename;
	fstream inputfile;
	inputfile.open(filename);
	ofstream inputcopy;
	inputcopy.open(newfn);

	//search for '//', delete everything following, print remainder to new file
	string line;
	int found, found2;
	while (getline(inputfile,line))
	{
		found  = line.find('/');
		found2 = line.find('/', found+1);
		if (found != line.npos && found2 == found+1)
			inputcopy << line.erase(found, line.length()) << endl;
		else
			inputcopy << line << endl;
	}
	inputcopy.close();

	//update filename;
	filename = newfn;
}

void programInputs::toAU()
{
	epsilon     /= 27.21;
	eaffin      /= 27.21;
	workfxn     /= 27.21;

	omega       = sqrt(72.0 * epsilon/pow(requil,2)*(1.0/m_C60));
	lbarrieritn = 2*M_PI*M_PI/pow(lbarrier,2);
	rbarrieritn = 2*M_PI*M_PI/pow(rbarrier,2);
}
