/********************************************************************
*
* file: mmAligner.cpp
*
* Copyright (c) 2007, Sittichai Jiampojamarn
* All rights reserved.
* 
* See the file COPYING in the top directory of this distribution
* for more information.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*
* Author: Sittichai Jiampojamarn
* Email: sj@cs.ualberta.ca
*
* Many-to-many alignment : 
* History : 
*
* Sept 9, 2009 : implemented n-best generation model to output
* n-best alignments. It's based on a variation of n-best viterbi algorithm.
*   : added an initialization feature, so we can change 
*   : EM starting point and also inject lingustic knowledge.
*
* Mar 2, 2009 : re-formatting input/output for 
*  NEWS: Shared Task on Transliteration
*
* Dec 10, 2008 : re-implemented the algorithm based on 
*	Sittichai Jiampojamarn, Grzegorz Kondrak and Tarek Sherif. 
*  "Applying Many-to-Many Alignments and Hidden Markov Models 
*	to Letter-to-Phoneme Conversion". Proceedings of the Annual 
*	Conference of the North American Chapter of the Association 
*	for Computational Linguistics (NAACL-HLT 2007), Rochester, 
*	NY, April 2007, pp.372-379
*
* Credits : Tarek Sherif originally proposed this algorithm 
*  based on the Ristad and Yianilos (1997) stochastic transducer
*  as a part of his Master thesis graduated in 2007, U. of Alberta.
*
* Citation : If you use the code for research or commercial purposes,
*  please acknowledge its use with a citation:
*
*  @InProceedings{jiampojamarn2007:main,
*  author    = {Jiampojamarn, Sittichai  and  Kondrak, Grzegorz  and  Sherif, Tarek},
*  title     = {Applying Many-to-Many Alignments and Hidden Markov Models to Letter-to-Phoneme Conversion},
*  booktitle = {Human Language Technologies 2007: The Conference of the North American Chapter of the Association for Computational Linguistics; Proceedings of the Main Conference},
*  month     = {April},
*  year      = {2007},
*  address   = {Rochester, New York},
*  publisher = {Association for Computational Linguistics},
*  pages     = {372--379},
*  url       = {http://www.aclweb.org/anthology/N/N07/N07-1047}
*  }
*
**********************************************************************/

#include <iostream>
#include <string>
#include <list>
#include <ctime>
#include <cstdio>

#include <param.h>
#include <mmEM.h>
#include <util.h>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>

inline void printTimeStamp(time_t startT)
{
	time_t stopT, diffT;

	time(&stopT);
	diffT = difftime(stopT, startT);
	cout << endl << "Time stamp: " << ctime(&stopT);
	cout << diffT << " seconds after started\n\n";
}

namespace boost_po = boost::program_options;

int main(int argc, char** argv)
{

    boost_po::options_description generic("Generic"); // name of help function
        generic.add_options() //detailed specification of command line interface
                ("help,h", "this help") // this option return boolean variable and is used to print command line help
            ;

    boost_po::options_description mandatory("Mandatory"); // name of help function
        mandatory.add_options() //detailed specification of command line interface
                ("inputFilename,i", boost_po::value<std::string>()->required(), "input file")
                ;

    boost_po::options_description optional("Optional"); // name of help function
        optional.add_options() //detailed specification of command line interface
                ("outputFilename,o", boost_po::value<std::string>(), "output file name")
                ("alignerOut", boost_po::value<std::string>()->default_value(""), "aligner model output filename")
                ("alignerIn", boost_po::value<std::string>()->default_value(""), "aligner model input filename")
                ("maxX", boost_po::value<int>()->default_value(2), "Maximum length of substring x")
                ("maxY", boost_po::value<int>()->default_value(2), "Maximum length of substring y")
                ("delX", boost_po::value<bool>()->default_value(false), "Allow deletion of substring x")
                ("delY", boost_po::value<bool>()->default_value(false), "Allow deletion of substring y")
                ("eqMap", boost_po::value<bool>()->default_value(false), "Allow mapping of |x| == |y| > 1")
                ("maxFn", boost_po::value<std::string>()->default_value("conYX"), "Maximization function [conXY, conYX, joint]")
                ("cutOff", boost_po::value<double>()->default_value(0.01), "Training threshold")
                ("printScore", boost_po::value<bool>()->default_value(false), "Report score of each alignment")
                ("prefixProcess", boost_po::value<std::string>()->default_value(""), "Specify prefix output files")
                ("nullChar", boost_po::value<std::string>()->default_value("_"), "Null character used")
                ("sepChar", boost_po::value<std::string>()->default_value("|"), "Separated character used")
                ("sepInChar", boost_po::value<std::string>()->default_value(":"), "Separated character used")
            	("nBest", boost_po::value<int>()->default_value(1), "Generate n-best alignments")
            	("inFormat", boost_po::value<std::string>()->default_value("news"), "Input file format [l2p, news]")
            	("initFile", boost_po::value<std::string>()->default_value(""), "Initial mapping (model) filename")
            	("initProbCut", boost_po::value<long double>()->default_value(0.5), "Cut-off sum prior probability")
            	("errorInFile", boost_po::value<bool>()->default_value(false), "Keep unaligned item in the output file")
            	("limitPair", boost_po::value<bool>()->default_value(false), "Limit the alignment pair to used only from the initFile")
                ;

	boost_po::options_description cmdline_options;
	cmdline_options.add(generic).add(mandatory).add(optional);

	boost_po::variables_map vm;
	boost_po::store(boost_po::parse_command_line(argc, argv, cmdline_options), vm);

    if(vm.count("help"))
    {
            std::cout << cmdline_options << std::endl;
            return 1;
    }

    try {
	    boost_po::notify(vm);
    } catch (boost_po::error& e) {
    	std::cout << "Error: " << e.what() << std::endl ;
    	return 1 ;
    }

    std::string inputFilename=vm["inputFilename"].as<std::string>();
 
    std::string outputFilename="";
    if (vm.count("outputFilename"))
    {
		outputFilename=vm["outputFilename"].as<std::string>();
    }

    std::string alignerOut=vm["alignerOut"].as<std::string>();
 
    std::string alignerIn=vm["alignerIn"].as<std::string>();
 
    int maxX=vm["maxX"].as<int>();
    int maxY=vm["maxY"].as<int>();
    
    bool delX=vm["delX"].as<bool>();
    bool delY=vm["delY"].as<bool>();

    bool eqMap=vm["eqMap"].as<bool>();

	std::array<std::string, 3> allowFnValues = {"conXY", "conYX", "joint"};
	std::string maxFn=vm["maxFn"].as<std::string>();
	if (not boost::algorithm::any_of_equal( allowFnValues.begin(), allowFnValues.end(), maxFn)) {
		std::cerr << "Error: invalid value of maxFn option" << std::endl ;
		return 2 ;
	}

	double cutOff=vm["cutOff"].as<double>();
	bool printScore=vm["printScore"].as<bool>();
    std::string prefixProcess=vm["prefixProcess"].as<std::string>();
    std::string nullChar=vm["nullChar"].as<std::string>();
	std::string sepChar=vm["sepChar"].as<std::string>();
	std::string sepInChar=vm["sepInChar"].as<std::string>();

	std::array<std::string, 2> inFormatValidValues = {"news", "l2p"};
	std::string inFormat=vm["inFormat"].as<std::string>();
	if (not boost::algorithm::any_of_equal( inFormatValidValues.begin(), inFormatValidValues.end(), inFormat)) {
		std::cerr << "Error: invalid value of inFormat option" << std::endl ;
		return 3 ;
	}

	std::string initFile=vm["initFile"].as<std::string>();
	long double initProbCut=vm["initProbCut"].as<long double>();
	bool errorInFile=vm["errorInFile"].as<bool>();
	bool limitPair=vm["limitPair"].as<bool>();
	bool nBest=vm["nBest"].as<int>();
	
	// define program parameters //
	param myParam;
	myParam.inputFilename = inputFilename ;
	myParam.outputFilename = outputFilename ;

	myParam.alignerIn = alignerIn ;
	myParam.alignerOut = alignerOut ;

	myParam.maxX = maxX ;
	myParam.maxY = maxY ;
	myParam.delX = delX ;
	myParam.delY = delY ;
	myParam.eqMap = eqMap ;

	myParam.maxFn = maxFn ;

	myParam.cutOff = cutOff ;

	myParam.printScore = printScore ;
	myParam.prefixProcess = prefixProcess ;

	myParam.nullChar = nullChar ;
	myParam.sepChar = sepChar ;
	myParam.sepInChar = sepInChar ;
	myParam.inFormat = inFormat ;
	myParam.nBest = nBest ;
	
	myParam.initFile = initFile ;
	myParam.initProbCut = initProbCut ;
	myParam.errorInFile = errorInFile ;
	myParam.limitPair = limitPair ;

	cout << "inputFilename" << " : " << inputFilename << endl;
	cout << "outputFilename" << " : " << outputFilename << endl;
	cout << "alignerIn" << " : " << alignerIn << endl;
	cout << "alignerOut" << " : " << alignerOut << endl;
	cout << "maxX" << " : " << maxX << endl;
	cout << "maxY" << " : " << maxY << endl;
	cout << "delX" << " : " << delX << endl;
	cout << "delY" << " : " << delY << endl;
	cout << "eqMap" << " : " << eqMap << endl;
	cout << "maxFn" << " : " << maxFn << endl;
	cout << "cutOff" << " : " << cutOff << endl;
	cout << "printScore" << " : " << printScore << endl;
	cout << "prefixProcess" << " : " << prefixProcess << endl;
	cout << "nullChar" << " : " << nullChar << endl;
	cout << "sepChar" << " : " << sepChar << endl;
	cout << "sepInChar" << " : " << sepInChar << endl;
	cout << "inFormat" << " : " << inFormat << endl;
	cout << "nBest" << " : " << nBest << endl;
	cout << "initFile" << " : " << initFile << endl;
	cout << "initProbCut" << " : " << initProbCut << endl;
	cout << "errorInFile" << " : " << errorInFile << endl;
	cout << "limitPair" << " : " << limitPair << endl;

	if (myParam.prefixProcess == "")
	{
		myParam.prefixProcess = "m-mAlign." + stringify(myParam.maxX) + "-" + stringify(myParam.maxY);
		if (myParam.delX)
		{
			myParam.prefixProcess += ".delX";
		}
		if (myParam.delY)
		{
			myParam.prefixProcess += ".delY";
		}

		myParam.prefixProcess += "." + stringify(myParam.nBest) + "-best";
		myParam.prefixProcess += "." + myParam.maxFn;

		if (myParam.printScore)
		{
			myParam.prefixProcess += ".pScore";
		}

		if (myParam.initFile != "")
		{
			myParam.prefixProcess += ".withInit";
		}

		if (myParam.errorInFile)
		{
			myParam.prefixProcess += ".errorInFile";
		}

		if (myParam.limitPair)
		{
			myParam.prefixProcess += ".limitPair";
		}
	}

	if (outputFilename == "") {
		outputFilename = inputFilename + "." + myParam.prefixProcess + ".align";
    }

	cout << endl;
	time(&myParam.startT);
	cout << "Started at: " << ctime(&myParam.startT) << endl << endl;
	
	mmEM myAligner;

	if (myParam.alignerIn != "")
	{
		// read aligner model from file
		myAligner.readAlignerFromFile(myParam);
	}
	else
	{
		// training aligner //
		myAligner.training(myParam);
		
		// writing aligner model to file //
		if (myParam.alignerOut == "")
		{
			myParam.alignerOut = myParam.inputFilename + "." + myParam.prefixProcess + ".align.model";
		}
		myAligner.writeAlingerToFile(myParam);

		// read model from file //
		myParam.alignerIn = myParam.alignerOut;
		myAligner.readAlignerFromFile(myParam);
	}
	printTimeStamp(myParam.startT);

	// create alignments //
	if (myParam.outputFilename == "")
	{
		myParam.outputFilename = myParam.inputFilename + "." + myParam.prefixProcess + ".align";
	}
	myAligner.createAlignments(myParam);
	printTimeStamp(myParam.startT);		

	return 0;
}
