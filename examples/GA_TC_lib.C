#include <stdio.h>
#include <ga/ga.h>
#include <ga/std_stream.h>

#define cout STD_COUT
#define endl STD_ENDL

#define INSTANTIATE_REAL_GENOME
#include <ga/GARealGenome.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator> 

using namespace std;

float Objective(GAGenome &);

static map<int,string> constraints;
static vector<string> header;
static vector<string> footer;
static map<int,int> group;
static map<int,string> gNames;
static string lib_name;
static string out_msf;
static string gen_name;

GABoolean GATerminateUponMaxScore(GAGeneticAlgorithm &ga) {
  stringstream buffer;
  streambuf * oldCoutStreamBuf = cout.rdbuf(buffer.rdbuf());
  
  cout << ga.statistics().bestIndividual();
  string bestIndividual = buffer.str();
  cout.rdbuf( oldCoutStreamBuf );
  
  cout << bestIndividual << endl;
  
  stringstream ss(bestIndividual);
  istream_iterator<int> begin(ss);
  istream_iterator<int> end;
  vector<int> vstrings(begin, end);
  sort(vstrings.begin(), vstrings.end());
  
  ofstream fout;
  fout.open (gen_name.c_str());
  fout << ga.statistics().bestIndividual() << endl;
  copy(vstrings.begin(), vstrings.end(), ostream_iterator<int>(fout, "\n"));
  fout.close();

  return(ga.generation() < ga.nGenerations() ? gaFalse : gaTrue);
}

int main(int argc, char** argv) {

// See if we've been given a seed to use (for testing purposes).  When you
// specify a random seed, the evolution will be exactly the same each time
// you use that seed number.

  unsigned int seed = 0;
  
  if(argc!=3) {
  	cout << "Execution: ./GA_TC_lib name length" << endl;
  	exit(0);
  }
  string sc_name(argv[1]);
  lib_name = sc_name;
  out_msf = sc_name;
  gen_name = sc_name;
  int length=atoi(argv[2]);

  sc_name += "_bog.dat";
  lib_name += ".lib";
  out_msf += ".msf";
  gen_name += ".gen";
  
  ifstream f;
  string line;
  int nCons = 1;
  int nGroup = 0;
  bool c_header = true;
 
  //consistency library file
  f.open ("TC_DATASET.lib");

  //PARSE FILE IN MEM
  if (f.is_open()) {

        while (! f.eof() ) {
        	getline (f,line);

		if(line[0]=='#') { //new pair
			//cout << line << endl;
			nGroup++;
			gNames[nGroup] = line;
			c_header = false;
		} else if(c_header) { //header
			header.push_back(line);
		} else if(line[0]=='!') { //footer
			footer.push_back(line);
		} else { //constraint
			constraints[nCons] = line;
			group[nCons] = nGroup;
			nCons++;
		}
        }

        f.close();
    }
    else cout << "File error" << endl;

// This genome uses an enumerated list of alleles.  We explictly add each 
// allele to the allele set.  Any element of the genome may assume the value
// of any member of the allele set.

  GARealAlleleSet alleles(1,nCons,1,GAAllele::INCLUSIVE,GAAllele::INCLUSIVE);
  GARealGenome genome(length, alleles, Objective);

// Now that we have the genomes, create a parameter list.

  GAParameterList params;
  GASteadyStateGA::registerDefaultParameters(params);
  params.set(gaNnGenerations, 5000);
  params.set(gaNpopulationSize, 500);
  params.set(gaNscoreFrequency, 10);
  params.set(gaNflushFrequency, 10);
  params.set(gaNselectScores, (int)GAStatistics::AllScores);
  params.parse(argc, argv, gaFalse);

// Now do a genetic algorithm.

  GASteadyStateGA ga(genome);
  ga.terminator(GATerminateUponMaxScore);
  ga.parameters(params);
  ga.set(gaNscoreFilename, sc_name.c_str());
  cout << "\nrunning ga ..."<<endl;
  ga.evolve(seed);
  
  stringstream buffer;
  streambuf * oldCoutStreamBuf = cout.rdbuf(buffer.rdbuf());
  
  cout << ga.statistics().bestIndividual();
  string bestIndividual = buffer.str();
  cout.rdbuf( oldCoutStreamBuf );
  
  cout << bestIndividual << endl;
  
  stringstream ss(bestIndividual);
  istream_iterator<int> begin(ss);
  istream_iterator<int> end;
  vector<int> vstrings(begin, end);
  sort(vstrings.begin(), vstrings.end());
  
  ofstream fout;
  fout.open (gen_name.c_str());
  fout << ga.statistics().bestIndividual() << endl;
  copy(vstrings.begin(), vstrings.end(), ostream_iterator<int>(fout, "\n"));
  fout.close();

  return 0;
}

void generate_lib(vector<int> indexed_lib) {
  ofstream fout;
  fout.open (lib_name.c_str()); 

  for (vector<string>::iterator it = header.begin(); it != header.end(); ++it)
    fout << *it << endl;

  int pgroup = 0;
  for (vector<int>::iterator it = indexed_lib.begin(); it != indexed_lib.end(); ++it) {
    if(pgroup<group[*it]) {
	fout << gNames[group[*it]] << endl;
	pgroup = group[*it];
    }
    fout << constraints[*it] << endl;
    
  }

  for (vector<string>::iterator it = footer.begin(); it != footer.end(); ++it)
    fout << *it << endl;

  fout.close();
}

float run_tc() {
  string cmd("PATH/t_coffee DATASET.tfa -output=msf_aln -dp_mode=myers_miller_pair_wise -outfile=");
  cmd += out_msf;
  cmd += " -lib=";
  cmd += lib_name;
  cmd += " > /dev/null 2>&1";
  system(cmd.c_str());

  string pcmd("PATH/bali_score DATASET.msf ");
  pcmd += out_msf;
  pcmd += " | grep 'SP score=' | cut -d' ' -f3";
  
  FILE* pipe = popen(pcmd.c_str(), "r");
  if (!pipe) 
  	cout << "ERROR\n";
  char buffer[128];
  string result = "";
  while (!feof(pipe)) {
  	if (fgets(buffer, 128, pipe) != NULL)
    	result += buffer;
    }
  pclose(pipe);
  
  return strtod(result.c_str(), 0);
}

float Objective(GAGenome& g) {
  GARealGenome& genome = (GARealGenome&)g;
  vector<int> indexed_lib;
  float value=0.0;

  for(int i=0; i<genome.length(); i++){
    indexed_lib.push_back((int)genome.gene(i));
  }

  sort(indexed_lib.begin(), indexed_lib.end());

  generate_lib(indexed_lib);
  value = run_tc();
  cout << value << endl;

  return value;
}
