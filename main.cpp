

// C++ standard headers
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cstdlib>

#include "boost/regex.hpp"

// C headers
#include <popt.h>        // For command line options

// Unix libraries
#include <unistd.h>
#include <sys/stat.h>

/* THIS PROGRAM REQUIRES THE "GNU Scientific Ligrary"  */
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>


#include "FastaFile.hpp"

using namespace std;

#define log2(A) (log(A)/log(2.0))

static const int alphabet_size = 4;

/* INPUT PARAMETERS */
const char *sequences_file_name = static_cast <const char *>(0); // file containing foreground sequences
const char *motif_consensus = static_cast <const char *>(0);

const char *output_file_name = static_cast<const char *>(0);        // file in which to print output

int use_both_strands = 0;
int allow_overlap_sites = 0;
int ignore_repeat = 0;
int verbose = 0;

char complement(char c) {
  static char btoc[26] = {
  //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
   'T','V','G','H','N','N','C','D','N','N','M','N','K','N','N','N','N','Y','S','A','N','B','W','N','R','N' 
  };
  static char btoc_small[26] = {
  //A,  b,  C,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,  p,  q,  r,  s,  T,  u,  v,  w,  x,  y,  z
   't','v','g','h','n','n','c','d','n','n','m','n','k','n','n','n','n','y','s','a','n','b','w','n','r','n' 
  };
  if (c <= 'Z' && c >= 'A')
  	return btoc[c - 'A'];
  else
	return btoc_small[c - 'a'];
}
string reverse_complement(string &s) {
  string r;
  transform(s.begin(), s.end(), back_inserter(r), complement);
  reverse(r.begin(), r.end());
  return r;
}

int get_sequence_length(vector<string> &s) {
  int count = 0;
  for (vector<string>::iterator i = s.begin(); i != s.end(); ++i)
    for (string::iterator j = i->begin(); j != i->end(); ++j)
      count += (*j != 'N' && *j != 'n');
  return count;
}

/****************** COMMAND LINE OPTIONS ********************/
poptContext optCon = NULL;
static struct poptOption optionsTable[] = {
  { "motif-consensus", 'c', POPT_ARG_STRING, &motif_consensus, 0, "the consensus sequence of motif" },
  { "both-strands", 'b', POPT_ARG_NONE, &use_both_strands, 0, "Search both strands" },
  { "ignore-repeat", 'i', POPT_ARG_NONE, &ignore_repeat, 0, "Ignore repeat (in small letters)"},
  { "allow-overlap-sites", 'r', POPT_ARG_NONE, &allow_overlap_sites, 0, "Allow overlapping sites)"},
  { "verbose", 'v', POPT_ARG_NONE, &verbose, 0, "Verbose mode" },
  { "output-file", 'o', POPT_ARG_STRING, &output_file_name, 0, "output file (default: stdout)" },
   POPT_AUTOHELP { NULL, 0, 0, NULL, 0, }
};

bool  BaseMatch (char q, char c)
{
	//see if the definition of q is included in c
	//the parameter should in in upper case
	
	if (ignore_repeat && q >= 'a' && q<= 'z') //repeat region should be ignored when the switch is on
		return false;
	
	if (c == q || c == 'N')
		return true;
	if ((c == 'R' && (q == 'A' || q == 'G'))
	  ||(c == 'Y' && (q == 'C' || q == 'T')) 
	  ||(c == 'M' && (q == 'A' || q == 'C'))
	  ||(c == 'K' && (q == 'G' || q == 'T'))
	  ||(c == 'S' && (q == 'C' || q == 'G'))
	  ||(c == 'W' && (q == 'A' || q == 'T'))
	  ||(c == 'H' && (q == 'A' || q == 'C' || q == 'T'))
	  ||(c == 'B' && (q == 'C' || q == 'G' || q == 'T'))
	  ||(c == 'V' && (q == 'A' || q == 'C' || q == 'G'))
	  ||(c == 'D' && (q == 'A' || q == 'G' || q == 'T'))
	  )
		return true;
	return false;
}


//////////////////////////////////////////
//////////// THE MAIN METHOD /////////////
//////////////////////////////////////////
int main(int argc, char **argv) {
  char c;
  
  /***************** GET COMMAND LINE ARGUMENTS *******************/
  optCon = poptGetContext("RegExpMatch", argc, (const char **)argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "fasta-file(.gz)");
  if (argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    exit(EXIT_SUCCESS);
  }
  if ((c = poptGetNextOpt(optCon)) >= 0) {
    fprintf(stderr, "Error: bad option ?");
    return EXIT_FAILURE;
  }
  if ((c = poptGetNextOpt(optCon)) < -1) {
    fprintf(stderr, "PatternMatch: bad argument %s: %s\n",
            poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(c));
    return EXIT_FAILURE;
  }
  if (poptPeekArg(optCon) == NULL) {
    poptPrintUsage(optCon, stderr, 0);
    return EXIT_FAILURE;
  }
  else sequences_file_name = poptGetArg(optCon);
  poptFreeContext(optCon);
  /**********************************************************************/

  // READ IN THE SEQUENCES AND THEIR NAMES
  FastaFile faa = FastaFile(const_cast<char *>(sequences_file_name));
  
  //faa.ToUpper(); 
  faa.Clean();
  vector<string> sequence_names = faa.GetNames();
  vector<string> sequences = faa.GetSequences();
  size_t num_seq = sequence_names.size();

  //boost::regex_constants::syntax_option_type flags = ignore_repeat ? match_default : boost::regex_constants::icase;
  //compile the regular expression
  boost::regex motif;
 
  ignore_repeat ? motif.assign (motif_consensus) :  motif.assign (motif_consensus, boost::regex_constants::icase);

  //this is necessary because the subroutine reverse_complement can not handle small letters
  //string motif_rc = reverse_complement (motif);
  
  //ostringstream buffer;
  ofstream out;
  if (output_file_name != static_cast<char *> (0))
  {
	//ofstream out;
	out.open (const_cast<char *>(output_file_name));
	//out << buffer.str ();
	//out.close();	
  }
 
  for (size_t i = 0; i < num_seq; ++i)
  {
	string name = sequence_names[i];
	string seq = sequences[i];
	
	if (verbose & (i % 500 == 0))
		cerr << "processing # " << i << ":" << name << ", size=" << seq.size () << endl; 
	
	size_t len = seq.size ();

	std::string::const_iterator begin = seq.begin();
	std::string::const_iterator end = seq.end();
	std::string::const_iterator iter = begin;
	boost::smatch hit;

	while (boost::regex_search(iter, end, hit, motif))
	{
		size_t hitStart = size_t (hit[0].first-begin);
		size_t hitEnd = size_t (hit[0].second-begin) - 1;
		string hitSeq = string (hit[0].first, hit[0].second);

		if (allow_overlap_sites)
		{
			iter=hit[0].first;
			iter++;
		}
		else
		{
			iter = hit[0].second; //set the next starting point
		}

		//the case is not handled yet
		if (output_file_name != static_cast<char *> (0))
		{
			out << name << "\t" 
				<< hitStart << "\t"
				<< hitEnd + 1 << "\t"
				<< hitSeq << "\t"
				<< 0 << "\t"
				<< "+" << endl;
		}
		else
		{
			cout << name << "\t" 
				<< hitStart << "\t"
				<< hitEnd + 1 << "\t"
				<< hitSeq << "\t"
				<< 0 << "\t"
				<< "+" << endl;
		}
	}


	if (use_both_strands == 0)
		continue;
	
	//the negative strand
	//reverse complement the motif will be difficult in the case of regular expressions
	//so reverse complement the sequence
	string seq_rc = reverse_complement (seq);
	begin = seq_rc.begin();
	end = seq_rc.end();
	iter = begin;
	
	while (boost::regex_search(iter, end, hit, motif))
	{
		size_t hitEnd = len - 1 - size_t (hit[0].first-begin);
		size_t hitStart = len - size_t (hit[0].second-begin);

		string hitSeq = string (hit[0].first, hit[0].second);

		if (allow_overlap_sites)
		{
			iter=hit[0].first;
			iter++;
		}
		else
		{
			iter = hit[0].second; //set the next starting point
		}

		//transform (target.begin (), target.end (), target.begin(), ::toupper);
		if (output_file_name != static_cast<char *> (0))
		{
			out << name << "\t" 
				<< hitStart << "\t"
				<< hitEnd + 1 << "\t"
				<< hitSeq << "\t"
				<< 0 << "\t"
				<< "-" << endl;
		}
		else
		{
			cout << name << "\t" 
				<< hitStart << "\t"
				<< hitEnd + 1 << "\t"
				<< hitSeq << "\t"
				<< 0 << "\t"
				<< "-" << endl;
		}
	}
  }

  if (output_file_name != static_cast<char *> (0))
  {
	out.close();	
  }
  //else
  //{
	//cout << buffer.str ();
  //}
  return EXIT_SUCCESS;
}

