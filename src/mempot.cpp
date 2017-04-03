// Generate transmembrane propensities
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <dirent.h>
#include <sys/stat.h>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

/*
int frequencies[20][32];

int get_slice_index(double& z){

	// 48  Angstrom thick membrane
	// 1.5 Angstrom thick slices
	// 32 slices

	double slice_start = -24.0;
	for(int i = 0;i < 32;i++,slice_start+=1.5){
		//cout << "i= " << i << ", z= " << z << "\t" << slice_start << " - " << slice_start+1.5;
		if(z >= slice_start && z < slice_start+1.5){
			//cout << "\t<---" << endl;
			return(i);
		}else{
			//cout << endl;
		}
	}
	return(-1);
}
*/

int frequencies[20][34];

int get_slice_index(double& z){

	// 48  Angstrom thick membrane
	// 1.5 Angstrom thick slices
	// 34 slices

	double slice_start = -24.0;
	double slice_stop = 24.0;

	if(z < slice_start){
		return(0);
	}else if(z >= slice_stop){
		return(33);
	}
	for(int i = 1;i < 33;i++,slice_start+=1.5){
		//cout << "i= " << i << ", z= " << z << "\t" << slice_start << " - " << slice_start+1.5;
		if(z >= slice_start && z < slice_start+1.5){
			//cout << "\t<---" << endl;
			return(i);
		}
	}
	return(0);
}

static void parse_pdb(std::string& pdbfile, std::vector<string>& chains, map <string, int>& codes){

	
	string chainlist;
	//cout << "Parsing PDB file " << pdbfile << " chains ";
	BOOST_FOREACH(string k, chains){
		chainlist.append(k);
		chainlist.append(",");
	}
	chainlist = chainlist.substr(0,chainlist.length()-1);
	//cout << chainlist << endl;
	

	char buf[1024];
	FILE *fin = fopen (pdbfile.c_str(), "r");
	if (fin != NULL){
		while (fgets (buf, 1024, fin)){
			string line = buf;
			if (strncmp(line.substr(0,4).c_str(),"ATOM",4) == 0){
				int tag = 0;
				BOOST_FOREACH(string k, chains){
					if (strcmp (line.substr(21,1).c_str(),k.c_str()) == 0){
						tag = 1;
					}
				}
				if(tag){
					string atom = line.substr(13,4);
					erase_all(atom, " ");
					string res = line.substr(17,3);
					double z = atof(line.substr(46,8).c_str());
										

					if(strcmp (atom.c_str(),"CA") == 0 && strcmp (res.c_str(),"GLY") == 0){
						//cout << line.substr(21,1) << " -" << atom << "- " << "-" << res << "- " << z << endl;
						int slice = get_slice_index(z);					
						//cout << codes[res] << "\t" << slice << endl;
						//cout << "aa index:\t" << codes[res] << endl;	
						if(slice > -1)frequencies[codes[res]][get_slice_index(z)]++;						
					}else if(strcmp (atom.c_str(),"CB") == 0 && strcmp (res.c_str(),"GLY") != 0){
						//cout << line.substr(21,1) << " -" << atom << "- " << "-" << res << "- " << z << endl;
						int slice = get_slice_index(z);					
						//cout << codes[res] << "\t" << slice << endl;
						//cout << "aa index:\t" << codes[res] << endl;
						if(slice > -1)frequencies[codes[res]][get_slice_index(z)]++;					
					}
				}			
			}			
		}
	}else{
		cout << "Couldn't open file " << pdbfile << endl;
		exit(1);
	}	
}

int main(int argc, const char* argv[]){

	int i = 1, c = 0;
	int precision = 5;
	string target;
	string homologues;
	string pdbchain_list = "TARGETS";
	string pdb_path = "pdbs/";
	string residues[20] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};	
	memset(frequencies, 0, sizeof(frequencies));		
	map <string, int> aacodes_num; 
	map <int, string> aacodes_name;	
	BOOST_FOREACH(string r, residues){
		aacodes_num[r] = c;	
		aacodes_name[c] = r;
		c++;
	}

	if (argc < 3){
		printf("Usage : %s [-p <pdb path>] [-l <pdb chain list>] -t <target chain> -h <homologues>\n", argv[0]);
		exit(1);
	}
	while(i < argc){
  		if( argv[i][0] == '-'){
    			i++;
    			switch(argv[i-1][1]){
	      			case 't' : {target=argv[i]; break;}
					case 'h' : {homologues=argv[i]; break;}
					case 'p' : {pdb_path=argv[i]; break;}
					case 'l' : {pdbchain_list=argv[i]; break;}
   			}   
   		}
    		i++;
  	}

	/*
	cout << "Target:\t\t" << target << endl;
	cout << "Chain list:\t" << pdbchain_list << endl;
	cout << "Homologues:\t" << homologues << endl;
	cout << "PDB path:\t" << pdb_path << endl;
	*/

	// Parse PDB chain list
	map <string, vector<string> > chain_list; 
	char buf[1024];
	FILE *fin = fopen (pdbchain_list.c_str(), "r");
	if (fin != NULL){
		while (fgets (buf, 1024, fin)){
			string line = buf;
			string pdb = line.substr(0,4);
			string chain = line.substr(4,1);
			boost::to_lower(pdb);
			chain_list[pdb].push_back(chain);
		}
	}else{
		cout << "Couldn't open file " << pdbchain_list << endl;
		exit(1);
	}

	// Tokenize list of homologues
	char_separator<char> sep(",");
	tokenizer< char_separator<char> > tokens(homologues,sep);

	// Parse PDB directory
	vector<string> pdb_list;
	DIR *dir;
	struct dirent *ent;
	struct stat st;
	dir = opendir (pdb_path.c_str());
	if (dir != NULL) {
		// Print all file names
	  	while ((ent = readdir (dir)) != NULL) {

			string path = pdb_path;
			path.append(ent->d_name);
			lstat(path.c_str(), &st);

			// PDB File			
			if(S_ISREG(st.st_mode)){
				int tag = 0;							
				BOOST_FOREACH(string t, tokens){
					boost::to_lower(t);
					if (strncmp (t.c_str(),ent->d_name,4) == 0){
						//cout << "Skipping homolog " << path << endl;
						tag++;
					}
				}
				if(!tag) parse_pdb(path,chain_list[path.substr(pdb_path.length(),4)],aacodes_num);
			}
	  	}
	  	closedir (dir);
	}else{
		cout << "Couldn't open PDB path " << pdb_path << endl;
		return (0);
	}

	
	//int z_total[32];
	int z_total[34];
	int res_total[20];
	int total = 0;
	//double mempot[20][32];
	memset(z_total, 0, sizeof(z_total));
	//memset(mempot, 0, sizeof(mempot));
	//cout << endl << "Frequencies" << endl;	
	for(int aa = 0; aa < 20; aa++){
		res_total[aa] = 0;
	}	
	for(int s = 0; s < 34; s++){
		z_total[s] = 0;
	}	


	for(int aa = 0; aa < 20; aa++){
		//cout << aacodes_name[aa] << ":\t";
		for(int s = 0; s < 34; s++){
			// Add pseudocount
			if(!frequencies[aa][s]) frequencies[aa][s]++; 
			//cout << frequencies[aa][s] << "\t";
			z_total[s] += frequencies[aa][s];
			res_total[aa] += frequencies[aa][s];
		}
		//cout << endl;
		//cout << "res_total[" << aa << "] = " << res_total[aa] << endl;
		total += res_total[aa];
		//cout << endl;
	}
	//cout << "total = " << total << endl;

	//cout << endl << "Potential" << endl;		
	for(int aa = 0; aa < 20; aa++){
		//cout << aacodes_name[aa] << ",";
		for(int s = 0; s < 34; s++){
			// Original
			//double pot = -0.582*log((double)frequencies[aa][s]/z_total[s]);
			// Original - no Boltzmann
			double pot = -1*log((double)frequencies[aa][s]/z_total[s]);

			// Fixed
			//double pot = -0.582*log(((double)frequencies[aa][s]/z_total[s])/((double)res_total[aa]/total));

			// Ez3D
			//double pot = -0.582*log((double)(frequencies[aa][s]*total)/(res_total[aa]*z_total[s]));

			/*
			cout << frequencies[aa][s] << "/" << z_total[s] << endl;
			cout << (double)frequencies[aa][s]/z_total[s] << endl;
			cout << log((double)frequencies[aa][s]/z_total[s]) << endl;
			cout << -0.582*log((double)frequencies[aa][s]/z_total[s]) << endl << endl;
			*/

			cout << setprecision(precision) << pot << "\t";
			//double pot = -0.582*log(frequencies[aa][s] * ((double)total/(res_total[aa]*z_total[s])));
			//cout << setprecision(precision) << pot << ",";
			

		}
		cout << endl;
	}

	return(1);
}
