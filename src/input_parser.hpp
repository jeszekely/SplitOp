#ifndef SPLITOP_INPUT_PARSER
#define SPLITOP_INPUT_PARSER

#include <string>
#include <iostream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

class programInputs
{
public:
	std::string filename;

	//Computational cell parameters
	int procs;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	int nx;
	int ny;
	double zl;
	double zr;
	int nphonon;
	int nelec;
	int eigen;

	//Time parameters
	double dt;
	double runtime;
	double output;

	//Physical parameters
	double hbar;
	double echarge;
	double m_electron;
	double m_C60;
	double elecvel;
	double eaffin;
	double lbarrier;
	double rbarrier;
	double xi;
	double workfxn;
	double bias;
	double elecpos;
	double epsilon;
	double requil;
	int coupfxn;
	double rad;
	double sigma;
	double del;
	double cconst;
	double alpha;

	//Ouputs
	int potential;
	int wvfxninit;
	int fullwvfxn;
	int elecstates;
	int phononstates;
	int animation;
	int slice;
	int sliceindex;
	int eleccorr;
	int moleccorr;
	int dcoeffr;
	int dcoeffl;
	int phases;
	int energy;
	int lifetime;
	int zexp;
	int fluxr;
	int fluxl;
	int fluxint;
	int correlation;

	//Derived Parameters
	double omega;
	double lbarrieritn;
	double rbarrieritn;

	//Imports the parameters from json file, converts to atomic units
	programInputs(std::string infile);

	//Removes comment lines from the json file
	void stripComments();

	//Convert parameters to atomic units
	void toAU();
};
#endif