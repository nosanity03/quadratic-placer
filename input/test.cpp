#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<map>
#include<vector>

using namespace std;

typedef vector<vector<int> > vvi;
typedef vector<vector<double> > vvd;
typedef vector<int> vi;
typedef vector<double> vd;

class mothercore{
	int numG,numP,numN;
	map <int, vi > gate;
	map <int, vd > pad;
	map <int, vvi> nets;
	map <int, double> gateX;
	map <int, double> gateY;
	
	public:
	mothercore(){
		numG=0;
		numP=0;
		numN=0;
	}
  
  // Helper function to return the number of gates.
	int get_numG() {
		return this->numG;
  }
  
  // Helper function to return the number of pads
  int get_numP() {
    return this->numP;
  }  

  // Helper function to return the number of nets
  int get_numN() {
    return this->numN;
  }
  
  // Helper function which makes a new gate and adds list of connections
  void add_gate(int gateNum, vi listofconnections) {
    this->gate[gateNum] = listofconnections;
    this->numG++;
    return;
  } 
	
  // Helper function which makes a new pad and adds its connections and location
	void add_pad(int padNum, vd netandlocation) {
		this->pad[padNum] = netandlocation;
		this->numP++;
		return;
  }

  // Helper function which makes a new net, if needed, and appends a connection to the net 'netnum'
	void add_net(int netNum, int connection, int gateorpad) {
		// 0 for gate, 1 for pad

		// if netnum doesn't already exist in dictionary
    if (this->nets.find(netNum) == this->nets.end()) {
      vi gates;
      vi pads;
      vvi empty_netconn;
      empty_netconn.push_back(gates);
      empty_netconn.push_back(pads);
			this->nets[netNum] = empty_netconn; // first vector for gates, second for pads
			this->numN += 1;
    }

		if (gateorpad == 1) this->nets[netNum][1].push_back(connection);
    else this->nets[netNum][0].push_back(connection);
    return;
  }
};

// Creates the "mothercore" class from file
mothercore create(char *filename) {
	cout << "Creating data structure from file ..." << endl;

  /////////////////////////////////////////////////////////////////////
  // Open file
  ifstream in(filename);
  streambuf *cinbuf = std::cin.rdbuf(); //save old buf
  cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

  /////////////////////////////////////////////////////////////////////
  // New object
  mothercore returncore;

  /////////////////////////////////////////////////////////////////////
  // Add gates, nets and pads to object from file

  string line;
  int numG, numN, numP, numnets, gatenum, padnum, netnum, xloc, yloc, i, j;
  vi gatevec;
  vd padvec;

  cin >> numG >> numN;
  //cout << numG << " " << numN << endl;
  for (i = 0; i < numG; i++) {
    cin >> gatenum >> numnets;
    //cout << gatenum << " " << numnets << " ";
    for (j = 0; j < numnets; j++) {
      cin >> netnum;
      gatevec.push_back(netnum);
	    returncore.add_net(netnum, gatenum, 0);
      //cout << netnum << " ";
    }
    returncore.add_gate(gatenum, gatevec);
    gatevec.clear();
    //cout << endl;
  }
  cin >> numP;
  //cout << numP << endl;
  for (i = 0; i < numP; i++) {
    cin >> padnum >> netnum >> xloc >> yloc;
	  returncore.add_net(netnum, padnum, 1);
    padvec.push_back(netnum);
    padvec.push_back(xloc);
    padvec.push_back(yloc);
    returncore.add_pad(padnum, padvec);
    padvec.clear();
    //cout << padnum << " " << netnum << " " << xloc << " " << yloc << endl;
  }

  /////////////////////////////////////////////////////////////////////
  // Close the file
  cin.rdbuf(cinbuf);   //reset to standard input again

  /////////////////////////////////////////////////////////////////////
  // Return
  cout << "Done. Added " << returncore.get_numG() << " gates, " << returncore.get_numP() << " pads, " << returncore.get_numN() << " nets." << endl;
  return returncore;
}

int main(){

  mothercore core = create("toy1");
  return 0;
}
