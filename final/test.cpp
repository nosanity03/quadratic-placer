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
  vi get_gateKeys() {
    vi v;
    for(map<int,vi>::iterator it = this->gate.begin(); it != this->gate.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }
	
  // Helper function to return gate coordinates
  vd get_gateCoords(int gateNum) {
    vd v;
    v.push_back(this->gateX[gateNum]); // x coordinate
    v.push_back(this->gateY[gateNum]); // y coordinate
    return v;
  }

  // Helper function to return connections of a gate
  vi get_gateconnections(int gateNum) {
    return this->gate[gateNum];
  }

  // Helper function which makes a new gate and adds list of connections
  void add_gate(int gateNum, vi listofconnections) {
    this->gate[gateNum] = listofconnections;
    this->numG++;
    return;
  }

  // Helper function to return the number of pads
  int get_numP() {
    return this->numP;
  }

  // Helper function to return pad keys
  vi get_padKeys() {
    vi v;
    for(map<int,vd>::iterator it = this->pad.begin(); it != this->pad.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  // Helper function to return pad coordinates
  vd get_padCoords(int padNum) {
    vd v;
    v.push_back(this->pad[padNum][1]); // x coordinate
    v.push_back(this->pad[padNum][2]); // y coordinate
    return v;
  }

	// Helper function which makes a new pad and adds its connections and location
	void add_pad(int padNum, vd netandlocation) {
		this->pad[padNum] = netandlocation;
		this->numP++;
		return;
  }

  // Helper function to return the number of nets
  int get_numN() {
    return this->numN;
  }

  // Helper function to return net keys
  vi get_netKeys() {
    vi v;
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
  vi get_netGateConns(int netNum) {
    return nets[netNum][0];
  }

  // Helper function to return pad connections to a net
  vi get_netPadConns(int netNum) {
    return nets[netNum][1];
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

  // Helper function which adds location values for given gate keys
  bool add_location(vd x, vd y, vi gatekeys, int bound[4]) {
    int l;
    double xloc, yloc;
    if ((gatekeys.size() == x.size()) && (x.size() == y.size())) {
      for (l = 0; l < gatekeys.size(); l++) {
        xloc = x[l];
        yloc = y[l];
        if (x[l] < bound[0]) xloc = bound[0]; // xmin
        if (x[l] > bound[1]) xloc = bound[1]; // xmax
        if (y[l] < bound[2]) yloc = bound[2]; // ymin
        if (y[l] > bound[3]) yloc = bound[3]; // ymax
        this->gateX[gatekeys[l]] = xloc;
        this->gateY[gatekeys[l]] = yloc;
      }
      return true;
    }	else {
      return false;
    }
  }

  vvd get_locations(vi gatekeys) {
    vd xloc, yloc;
    vvd returnvec;
    int l;
    for (l = 0; l < gatekeys.size(); l++) {
      xloc.push_back(this->gateX[gatekeys[l]]);
      yloc.push_back(this->gateY[gatekeys[l]]);
    }
    returnvec.push_back(xloc);
    returnvec.push_back(yloc);
    return returnvec;
  }

  // Helper function which prints the locations of all gates in the present
  // core.
  void print_all_locations() {
    int i;
    vi keyG = this->get_gateKeys();
    cout << "Locations:";
    for (i = 0; i < this->numG; ++i) { 
      cout << endl;
      cout << "Gate " << keyG[i] << ": " << this->gateX[keyG[i]] << ", " << this->gateY[keyG[i]];
    }
    cout << endl;
  }

  // Helper function which prints the pads in the present core.
  void print_all_pads() {
    int i, padnum;
    vi keyP = this->get_padKeys();
    cout << "Pads:";
    for (i = 0; i < this->numP; ++i) { 
      padnum = keyP[i];
      cout << endl;
      cout << "Pad " << padnum << ": Net - " << this->pad[padnum][0] << ", Location - ("<< this->pad[padnum][1] << ", " << this->pad[padnum][2] << ")";
    }
    cout << endl;
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
                vd vec;
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

void writeback(mothercore *core, char *filename) {
  ofstream out(filename);
  streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
  
  vi keyG = core->get_gateKeys();
  vd gateloc;
  int i, gatekey;
  for (i = 0; i < keyG.size(); i++) {
    gatekey = keyG[i];
    gateloc = core->get_gateCoords(gatekey);
    cout << gatekey << " " << gateloc[0] << " " << gateloc[1] << endl;
  }

  cout.rdbuf(coutbuf); //reset to standard output again
  return;
}


vd solve(vi R, vi C, vd V, vd ba) {
  coo_matrix A;
  vd aout;

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

bool solveforx(mothercore *core, int bound[4]) {
  cout << "Solving for locations ..." << endl;
  /////////////////////////////////////////////////////////////////////
  // Initializations

  vi keyG = core->get_gateKeys();
  vi keyP = core->get_padKeys();
  vi keyN = core->get_netKeys();
  int G = keyG.size();
  int P = keyP.size();
  int N = keyN.size();
  map <int, double> weights;
  int i, j, k, netval;

  if (G == 0) {
    return false;
  }

  double** A = new double*[G];
  for(i = 0; i < G; ++i)
    A[i] = new double[G];

  double* bx = new double[G];
  double* by = new double[G];

  for (i = 0; i < G; i++) {
    for (j = 0; j < G; j++) {
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
    k = core->get_numNetConns(netval);
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
  vi gates, pads;
  vd padCoordinate;
  for(k = 0; k < N; ++k) {
    netval = keyN[k];
    weight = weights[netval];
    gates = core->get_netGateConns(netval);
    numgates = gates.size();
    pads = core->get_netPadConns(netval);
    numpads = pads.size();
    if (numgates > 1) {
      i = 0;
      while (i < numgates-1) {
        j = i+1;
        while (j < numgates) {
          A[gateorder[gates[i]]][gateorder[gates[j]]] = -weight;
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
          padCoordinate = core->get_padCoords(pads[j]);
          bx[gateorder[gates[i]]] += weight*padCoordinate[0];
          by[gateorder[gates[i]]] += weight*padCoordinate[1];
          j += 1;
        }
        i += 1;
      }
    }
  }

  /* 
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

  vi Rsm, Csm;
  vd Vsm;
  vd locx;
  vd locy;
  vd bxsm;
  vd bysm;

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

  if (!core->add_location(locx, locy, keyG, bound)) {
    cout<<"Error in adding new locations."<<endl;
    return false;
  } else {
    cout<<"Successfully added new locations."<<endl;
    return true;
  }
}

// Function which sorts given gates according to locations horizontally or
// vertically. Horizontally if hORv = 0; Vertically if hORv = 1
vvi assign(mothercore *core, vi gatekeys, int hORv) { 
  cout << "Assignment ..." << endl;

  /////////////////////////////////////////////////////////////////////
  // Initializations
  int i, gatenum = 0;
  double customlocation = 0;
  int numkeys = gatekeys.size();
  int halfOfTotalGates = numkeys/2;

  vi sortedGates, tempvector;
  vvi returnvectors;
  vd location;
  map< double, int > sortingmap; 
  // map sorts its elements by keys. It doesn't care about the values when sorting. 
  // So, our key is custom location value, and value is gatekey.

  /////////////////////////////////////////////////////////////////////
  // getting locations of given gatekeys, and forming custom location value.
  for (i = 0; i < numkeys; i++) {
    gatenum = gatekeys[i];
    location = core->get_gateCoords(gatenum);

    // merge 2 sort keys into 1 using (100000*x+y) (hORv = 0) or 
    // (100000*y+x) (hORv = 1). This will work as long as there are less than 
    // 100k gates.
    customlocation = (location[hORv]*100000)+location[1-hORv];

    // Checking for duplicate key
    while ( sortingmap.find(customlocation) != sortingmap.end() ) {
      // found a duplicate key, change key
      customlocation += 0.000001;
    }
    sortingmap[customlocation] = gatenum;
  }

  /////////////////////////////////////////////////////////////////////
  // get the sorted gates
  for (map<double, int>::iterator it = sortingmap.begin(); it != sortingmap.end(); it++) {
    // cout << it->second << '\n';
    sortedGates.push_back(it->second);
  }

  /////////////////////////////////////////////////////////////////////
  // return
  for (i = 0; i < halfOfTotalGates; i++) {
    tempvector.push_back(sortedGates[i]);
  }
  returnvectors.push_back(tempvector);
  tempvector.clear();
  for (i = halfOfTotalGates; i < numkeys; i++) {
    tempvector.push_back(sortedGates[i]);
  }
  returnvectors.push_back(tempvector);
  return returnvectors;
}

void update_coordinates(vd *padlocation, int bound[4]) {
  if ((*padlocation)[0] < (double)bound[0]) (*padlocation)[0] = (double)bound[0]; // xmin
  if ((*padlocation)[0] > (double)bound[1]) (*padlocation)[0] = (double)bound[1]; // xmax
  if ((*padlocation)[1] < (double)bound[2]) (*padlocation)[1] = (double)bound[2]; // ymin
  if ((*padlocation)[1] > (double)bound[3]) (*padlocation)[1] = (double)bound[3]; // ymax
  return;
}

vvd containNrun(mothercore *core, vi gatekeys, int bound[4], int hORv, int lORr) {
  cout << "Containment ..." << endl;
  mothercore newcore;
  int i, pgatenum, k, pnet, l, newconn, ppadnum;
  vi pgateconns, donepads, donenets, netconn;
  vd padtemp, padlocation;

  /////////////////////////////////////////////////////////////////////
  // add gates to the new core, which are on appropriate side
  for (i = 0; i < gatekeys.size(); i++) {
    pgateconns = core->get_gateconnections(gatekeys[i]);
    newcore.add_gate(gatekeys[i], pgateconns);
  }

  /////////////////////////////////////////////////////////////////////
  // add pads and nets to the new core

  ppadnum = 0;
  // i recurses through the gates in gatekeys
  for (i = 0; i < gatekeys.size(); i++) {
    pgatenum = gatekeys[i];
    pgateconns = newcore.get_gateconnections(pgatenum);

    // k recurses through the nets in the presentgate
    for (k = 0; k < pgateconns.size(); k++) {
      pnet = pgateconns[k];
      if (find(donenets.begin(), donenets.end(), pnet) != donenets.end()) {
        // processed the present net already
        continue;
      } else {
        donenets.push_back(pnet);
      }

      // for each net get all the pad connections
      netconn = core->get_netPadConns(pnet);

      // process connected pads as new pads
      for (l = 0; l < netconn.size(); l++) {
        newconn = netconn[l];	// pad number
        ppadnum++;
        newcore.add_net(pnet, ppadnum, 1);

        padtemp.push_back((double)pnet);
        padlocation = core->get_padCoords(newconn);

        // correct coordinates which are outside the bounding box
        update_coordinates(&padlocation, bound);
        padtemp.push_back(padlocation[0]);
        padtemp.push_back(padlocation[1]);

        // add new pad
        newcore.add_pad(ppadnum, padtemp);
        padtemp.clear();
      }

      // for each net get all the gate connections
      netconn = core->get_netGateConns(pnet);

      // process connected gates as new pads
      for (l = 0; l < netconn.size(); l++) {
        newconn = netconn[l]; // gate number

        // skipping if newconn is a gate inside the bounding box
        if (find(gatekeys.begin(), gatekeys.end(), newconn) != gatekeys.end()) {
          newcore.add_net(pnet, newconn, 0);
          continue;
        }
        ppadnum++;
        newcore.add_net(pnet, ppadnum, 1);

        padtemp.push_back((double)pnet);
        padlocation = core->get_gateCoords(newconn);

        // correct coordinates of gates if they are outside the bounding box
        update_coordinates(&padlocation, bound);

        // correct coordinates of gates if they are inside the bounding box
        if (hORv == 0) {
          if (lORr == 0)
            // left of a horizontal cut
            padlocation[0] = bound[1];
          else
            // right of a horizontal cut
            padlocation[0] = bound[0];
        } else {
          if (lORr == 0)
            // top of a vertical cut
            padlocation[1] = bound[3];
          else
            // bottom of a vertical cut
            padlocation[1] = bound[2];
        }
        padtemp.push_back(padlocation[0]);
        padtemp.push_back(padlocation[1]);

        // add new pad
        newcore.add_pad(ppadnum, padtemp);
        padtemp.clear();
      }
    }
  }

  cout << "Done. Added " << newcore.get_numG() << " gates, " ;
  cout << newcore.get_numP() << " pads, " << newcore.get_numN() << " nets." << endl;
  //newcore.print_all_locations();
  //newcore.print_all_pads();

  /////////////////////////////////////////////////////////////////////
  // solve "newcore" for locations of gates inside bound
  solveforx(&newcore, bound);
  return newcore.get_locations(gatekeys);
}

void place(mothercore *core, vi gatekeys, int bound[4], int n) {
  if (n >= 8) return;
  else {
    cout << endl;
    cout << endl;
    cout << "n = " << n << endl;
    cout << endl;
    int nnext = n*2;

    // sort overall horizontally
    cout << endl;
    cout << "Horizontal Sort ..." << endl;
    vvi leftrightGates = assign(core, gatekeys, 0);
    int midx = (bound[1]-bound[0])/2;

    // call containNrun
    int left_bound[4] = { bound[0], bound[1]-midx, bound[2], bound[3] };
    int right_bound[4] = { bound[0]+midx, bound[1], bound[2], bound[3] };

    cout << endl;
    cout << "Containing Left Gates:" << endl;
    vvd leftlocs = containNrun(core, leftrightGates[0], left_bound, 0, 0);

    cout << endl;
    cout << "Containing Right Gates:" << endl;
    vvd rightlocs = containNrun(core, leftrightGates[1], right_bound, 0, 1);

    core->add_location(leftlocs[0], leftlocs[1], leftrightGates[0], left_bound);
    core->add_location(rightlocs[0], rightlocs[1], leftrightGates[1], right_bound);

    // sort left half vertically
    cout << endl;
    cout << "Left Half Vertical Sort ..." << endl;
    vvi left_topbottomGates = assign(core, leftrightGates[0], 1);
    int midy = (bound[3]-bound[2])/2;

    // call containNrun
    int leftbottom_bound[4] = { bound[0], bound[1]-midx, bound[2], bound[3]-midy };
    int lefttop_bound[4] = { bound[0], bound[1]-midx, bound[2]+midy, bound[3] };

    cout << endl;
    cout << "Containing Left Bottom Gates:" << endl;
    vvd leftbottomlocs = containNrun(core, left_topbottomGates[0], leftbottom_bound, 1, 0);

    cout << endl;
    cout << "Containing Left Top Gates:" << endl;
    vvd lefttoplocs = containNrun(core, left_topbottomGates[1], lefttop_bound, 1, 1);

    // sort right half vertically
    cout << endl;
    cout << "Right Half Vertical Sort ..." << endl;
    vvi right_topbottomGates = assign(core, leftrightGates[1], 1);

    // call containNrun
    int rightbottom_bound[4] = { bound[0]+midx, bound[1], bound[2], bound[3]-midy };
    int righttop_bound[4] = { bound[0]+midx, bound[1], bound[2]+midy, bound[3] };

    cout << endl;
    cout << "Containing Right Bottom Gates:" << endl;
    vvd rightbottomlocs = containNrun(core, right_topbottomGates[0], rightbottom_bound, 1, 0);

    cout << endl;
    cout << "Containing Right Top Gates:" << endl;
    vvd righttoplocs = containNrun(core, right_topbottomGates[1], righttop_bound, 1, 1);

    core->add_location(leftbottomlocs[0], leftbottomlocs[1], left_topbottomGates[0], leftbottom_bound);
    core->add_location(lefttoplocs[0], lefttoplocs[1], left_topbottomGates[1], lefttop_bound);
    core->add_location(rightbottomlocs[0], rightbottomlocs[1], right_topbottomGates[0], rightbottom_bound);
    core->add_location(righttoplocs[0], righttoplocs[1], right_topbottomGates[1], righttop_bound);

    place(core, left_topbottomGates[0], leftbottom_bound, nnext);
    place(core, left_topbottomGates[1], lefttop_bound, nnext);
    place(core, right_topbottomGates[0], rightbottom_bound, nnext);
    place(core, right_topbottomGates[1], righttop_bound, nnext);

    return;
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "Invalid Command: Run as `executable filename`" << endl;
    return -1;
  }
  mothercore core;
  core.create(argv[1]);
  int initial_bound[4] = { 0, 100, 0, 100 };
  solveforx(&core, initial_bound);
  
  vi keyG = core.get_gateKeys();
  place(&core, keyG, initial_bound, 1);
  
  writeback(&core, "out.txt");
  return 0;
}
