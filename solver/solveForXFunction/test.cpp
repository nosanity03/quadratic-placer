#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<map>
#include<vector>

#include <cmath>
#include <valarray>
#include <cstdio>
#include "solver.h"

using namespace std;

typedef vector<vector<int> > vvi;
typedef vector<int> vi;
template<class T>
bool IsNaN(T t)
{
    return t != t;
}


class mothercore{
	int numG,numP,numN,side[4];
	map <int, vector <int> > gate;
	map <int, vector <int> > pad;
	map <int, vvi> nets;
	map <int, int> gateX;
	map <int, int> gateY;
	vector <int> sortedgate;
	
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
	
  // Helper function to return gate keys
  vector<int> get_gateKeys() {
    vector<int> v;
    for(map<int,vi>::iterator it = this->gate.begin(); it != this->gate.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }
	
  // Helper function to return the number of pads
	int get_numP() {
		return this->numP;
  }
  
  // Helper function to return pad keys
  vector<int> get_padKeys() {
    vector<int> v;
    for(map<int,vi>::iterator it = this->pad.begin(); it != this->pad.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  // Helper function to return pad coordinates
  vector<int> get_padCoords(int padNum) {
    vector<int> v;
    v.push_back(this->pad[padNum][1]); // x coordinate
    v.push_back(this->pad[padNum][2]); // y coordinate
    return v;
  }
	
  // Helper function to return the number of nets
	int get_numN() {
		return this->numN;
  }

  // Helper function to return net keys
  vector<int> get_netKeys() {
    vector<int> v;
    for(map<int,vvi>::iterator it = this->nets.begin(); it != this->nets.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  // Helper function to return number of net connections
  int get_numNetConns(int netNum) {
    int total = nets[netNum][0].size() + nets[netNum][1].size();
    return total;
  }

  // Helper function to return gate connections to a net
  vector<int> get_netGateConns(int netNum) {
    return nets[netNum][0];
  }

  // Helper function to return pad connections to a net
  vector<int> get_netPadConns(int netNum) {
    return nets[netNum][1];
  }

  // Helper function to modify class using file
  void create (char *filename){
    ifstream file;
    file.open(filename);
    if(!file.is_open()) cout<<"Error in opening file."<<endl;
    else {
      cout<<"File opened."<<endl;
      string line;
      int line_no=0,word_no,key,no_of_nets,num;
      while(file.good()){
        getline(file,line);
        word_no=0;
        for(int count=0;count<=line.size();count++){
          string word="";
          if((line[count]==' ')||(count==0)){
            if (count==0) {
              word.append(&line[count]);
              word_no++;
            }	
            else {
              word.append(&line[count+1]);
              word_no++;
            }	
            num=atoi(word.c_str());
            if (line_no==0){
              if (word_no==1) {
                numG=num;
                //cout<<"numG = "<<numG<<endl;
              }
              else if (word_no==2){
                numN=num;
                //cout<<"numN = "<<numN<<endl;
                vector <int> gatevec;
                vector <int> padvec;
                vector <vector <int> > vec;
                vec.push_back(gatevec);
                vec.push_back(padvec);
                for (int i=0;i<numN;i++) nets[i+1]=vec;
              }	
            }
            else if (line_no<=numG) {
              if (word_no==1) {
                key = num;
                //cout << "Gate no. " << key << " is connected to net no. ";
                vector <int> vec;
                gate[key]=vec;
              }
              else if (word_no==2) {
                no_of_nets=num;
              }
              else if ((word_no>2)&&(word_no<(no_of_nets+3))) {
                gate[key].push_back(num);
                nets[num][0].push_back(key);
                //cout << gate[key][word_no-3] <<' ';
              }
            }
            else if (line_no==(numG+1)) {
              if (word_no==1){
                numP=num;
                //cout << "numP = " << numP <<endl;
              }
            }
            else if ((line_no>(numG+1))&&(line_no<(numG+2+numP))){
              if (word_no==1) {
                key = num;
                //cout << "Pad no. " << key << " is connected to net no. ";
                vector <int> vec;
                pad[key]=vec;
              }
              else if (word_no==2) {
                pad[key].push_back(num);
                //cout << pad[key][0] <<" at coordinates (";
                nets[num][1].push_back(key);
              }
              else {
                pad[key].push_back(num);
                //if (word_no==3) cout << pad[key][1]<<',';
                //else if (word_no==4) cout << pad[key][2]<<')';
              }
            }
          }
        }
        //cout<<endl;
        line_no++;
      }
    }
    file.close();

    /*
    for (int i=0;i<numN;i++){
      cout<<"net "<<i+1<<" is connected to gates ";
      for (int j=0;j<nets[i+1][0].size();j++){
        cout<<nets[i+1][0][j]<<' ';
      }
      if (nets[i+1][1].size()!=0){
        cout<<"and to pads ";
        for (int j=0;j<nets[i+1][1].size();j++){
          cout<<nets[i+1][1][j]<<' ';
        }
      }
      cout<<endl;
    }
    */
  }
};


vector<double> solve(vector<int> R, vector<int> C, vector<double> V, vector<double> ba) {
  coo_matrix A;
  vector<double> aout;

  A.n = ba.size();
  A.nnz = R.size();
  A.row.resize(A.nnz);
  A.col.resize(A.nnz);
  A.dat.resize(A.nnz);

  A.row = valarray<int>(R.data(), R.size());
  A.col = valarray<int>(C.data(), C.size());
  A.dat = valarray<double>(V.data(), V.size());

  valarray<double> x(A.n);
  valarray<double> b(ba.data(), A.n);

  A.solve(b, x);

  // cout << "b = " << endl;
  // print_valarray(b);
  // cout << "x = " << endl;
  // print_valarray(x);

  aout.assign(begin(x), end(x));  // needs -std=c++11 flag 
  return aout;
}

void solveforx(mothercore core) {
  cout << "Solving for locations ..." << endl;
  /////////////////////////////////////////////////////////////////////
  // Initializations

  //gate = core.get_gate();
  //pad = core.get_pad()
  vector<int> keyG = core.get_gateKeys();
  vector<int> keyP = core.get_padKeys();
  vector<int> keyN = core.get_netKeys();
  int G = keyG.size();
  int P = keyP.size();
  int N = keyN.size();
  map <int, double> weights;
  int i, j, k, netval;

  if (G == 0) {
    return;
  }

  //double C[G][G], A[G][G];
  //double bx[G], by[G]; 

  double** C = new double*[G];
  for(i = 0; i < G; ++i)
    C[i] = new double[G];
  
  double** A = new double*[G];
  for(i = 0; i < G; ++i)
    A[i] = new double[G];

  double* bx = new double[G];
  double* by = new double[G];
  
  for (i = 0; i < G; i++) {
    for (j = 0; j < G; j++) {
      C[i][j] = 0;
      A[i][j] = 0;
    }
    bx[i] = 0;
    by[i] = 0;
  }

  /////////////////////////////////////////////////////////////////////
  // Calculating weights using number of connections of nets
  cout << "Calculating weights ..." << endl;
  for(i = 0; i < N; ++i) {
    netval = keyN[i];
    k = core.get_numNetConns(netval);
    //cout << netval << " " << k << endl;
    weights[netval] = 1.0/(k-1);
  }

  /////////////////////////////////////////////////////////////////////
  // Gate numbers may vary, this dictionary keeps them in order
  map <int, int> gateorder;
  for(i = 0; i < G; ++i) {
    gateorder[keyG[i]] = i;
  }

  /////////////////////////////////////////////////////////////////////
  // Calculating C and A matrices, bx and by vectors
  cout << "Calculating valid contributions to the cost function ..." << endl;

  int numgates, numpads;
  double weight;
  vector<int> gates, pads, padCoordinate;
  for(k = 0; k < N; ++k) {
    netval = keyN[k];
    weight = weights[netval];
    gates = core.get_netGateConns(netval);
    numgates = gates.size();
    pads = core.get_netPadConns(netval);
    numpads = pads.size();
    if (numgates > 1) {
      i = 0;
      while (i < numgates-1) {
        j = i+1;
        while (j < numgates) {
          C[gateorder[gates[i]]][gateorder[gates[j]]] = weight;
          A[gateorder[gates[i]]][gateorder[gates[j]]] = -weight;
          C[gateorder[gates[j]]][gateorder[gates[i]]] = weight;
          A[gateorder[gates[j]]][gateorder[gates[i]]] = -weight;
          A[gateorder[gates[i]]][gateorder[gates[i]]] += weight;
          A[gateorder[gates[j]]][gateorder[gates[j]]] += weight;
          j += 1;
        }
        i += 1;
      }
    }

    if (numpads > 0) {
      i = 0;
      while (i < numgates) {
        j = 0;
        while (j < numpads) {
          A[gateorder[gates[i]]][gateorder[gates[i]]] += weight;
          padCoordinate = core.get_padCoords(pads[j]);
          bx[gateorder[gates[i]]] += weight*padCoordinate[0];
          by[gateorder[gates[i]]] += weight*padCoordinate[1];
          j += 1;
        }
        i += 1;
      }
    }
  }

   /* 
      cout << "C:" << endl;
      for (i = 0; i < G; i++) {
      for (j = 0; j < G; j++) {
      cout << C[i][j] << " ";
      }
      cout << endl;
      }

      cout << "A:" << endl;
      for (i = 0; i < G; i++) {
      for (j = 0; j < G; j++) {
      cout << A[i][j] << " ";
      }
      cout << endl;
      }

      cout << "bx:" << endl;
      for (i = 0; i < G; i++) {
      cout << bx[i] << endl;
      }

      cout << "by:" << endl;
      for (i = 0; i < G; i++) {
      cout << by[i] << endl;
      }
      
*/

  /////////////////////////////////////////////////////////////////////
  // Derive R, C, V matrices for sparse matrix generation
  cout << "Forming R, C, V matrices ..." << endl;

  vector<int> Rsm, Csm;
  vector<double> Vsm;
  vector<double> locx;
  vector<double> locy;
  vector<double> bxsm;
  vector<double> bysm;

  for (i = 0; i < G; i++) {
    for (j = 0; j < G; j++) {
      if ((A[i][j] != 0) && isfinite(A[i][j])) {
        Rsm.push_back(i);
        Csm.push_back(j);
        Vsm.push_back(A[i][j]);
      }
    }
    bxsm.push_back(bx[i]);
    bysm.push_back(by[i]);
  }

  /////////////////////////////////////////////////////////////////////
  // Solve for x and y vectors
  cout << "Solving for x and y ..." << endl;
  locx = solve(Rsm, Csm, Vsm, bxsm);
  locy = solve(Rsm, Csm, Vsm, bysm);


  cout << "New Locations = ";
  for (i = 0; i < G; ++i) { 
    cout << endl;
    cout << "Gate " << keyG[i] << ": " << locx[i] << ", " << locy[i];
  }
  cout << endl;


  return;
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "Invalid Command: Run as `executable filename`" << endl;
    return -1;
  }
  mothercore core;
  core.create(argv[1]);
  solveforx(core);
  return 0;
}
