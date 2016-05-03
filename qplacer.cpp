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

/**
 * Class which defines an ASIC and its components.
 */
class mothercore{
	int numG,                 /**< Number of gates in the ASIC */
      numP,                 /**< Number of pads in the ASIC */ 
      numN;                 /**< Number of nets in the ASIC */
	map <int, vi > gate;      /**< A Hash Table storing informations about gates in the ASIC */
	map <int, vd > pad;       /**< A Hash Table storing informations about pads in the ASIC */
	map <int, vvi> nets;      /**< A Hash Table storing informations about nets in the ASIC */
  map <int, double> gateX;  /**< A Hash Table storing x-coordinates of gates in the ASIC */
  map <int, double> gateY;  /**< A Hash Table storing y-coordinates of gates in the ASIC */

  public:
  
  /**
   * Constructor of the class 'mothercore'
   */
  mothercore(){
    numG=0;
    numP=0;
    numN=0;
  }

  /**
   * Helper function to return the number of gates.
   * @see get_numP()
   * @see get_numN()
   * @return The number of gates in the ASIC.
   */
  int get_numG() {
    return this->numG;
  }

  /**
   * Helper function to return gate keys.
   * @see get_padKeys()
   * @see get_netKeys()
   * @return A vector with the gate-id of gates in the ASIC.
   */
  vi get_gateKeys() {
    vi v;
    for(map<int,vi>::iterator it = this->gate.begin(); it != this->gate.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  /**
   * Helper function to return gate coordinates.
   * @param gateNum The gate-id for which coordinates are needed.
   * @see get_padCoords(int padNum)
   * @return The x, y coordinates of the gate with gate-id 'gateNum'.
   */
  vd get_gateCoords(int gateNum) {
    vd v;
    v.push_back(this->gateX[gateNum]); // x coordinate
    v.push_back(this->gateY[gateNum]); // y coordinate
    return v;
  }

  /**
   * Helper function to return connections of a gate.
   * @param gateNum The gate-id of the gate for which connections are needed.
   * @return A vector of the net-ids connected to the gate with gate-id 'gateNum'.
   */
  vi get_gateconnections(int gateNum) {
    return this->gate[gateNum];
  }

  /**
   * Helper function which makes a new gate and adds list of connections.
   * @param gateNum The gate-id of the gate to be added.
   * @param listofconnections A vector of net-ids connected to the gate to
   * be added.
   * @see add_pad(int padNum, vd netandlocation)
   * @see add_net(int netNum, int connection, int gateorpad)
   */
  void add_gate(int gateNum, vi listofconnections) {
    this->gate[gateNum] = listofconnections;
    this->numG++;
    return;
  }

  /**
   * Helper function to return the number of pads.
   * @see get_numG()
   * @see get_numN()
   * @return The number of pads in the ASIC.
   */
  int get_numP() {
    return this->numP;
  }

  /**
   * Helper function to return pad keys.
   * @see get_gateKeys()
   * @see get_netKeys()
   * @return A vector with the pad-id of pads in the ASIC.
   */
  vi get_padKeys() {
    vi v;
    for(map<int,vd>::iterator it = this->pad.begin(); it != this->pad.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  /**
   * Helper function to return pad coordinates.
   * @param padNum The pad-id of the pad for which coordinates are needed.
   * @see get_gateCoords(int gateNum)
   * @return The x, y coordinates of the pad with pad-id 'padNum'.
   */
  vd get_padCoords(int padNum) {
    vd v;
    v.push_back(this->pad[padNum][1]); // x coordinate
    v.push_back(this->pad[padNum][2]); // y coordinate
    return v;
  }

  /**
   * Helper function which makes a new pad and adds its connections and location
   * @param padNum The pad-id of the pad to be added.
   * @param netandlocation A vector with net-id connected to the pad to
   * be added, and its x, y coordinate.
   * @see add_gate(int gateNum, vi listofconnections)
   * @see add_net(int netNum, int connection, int gateorpad)
   */
  void add_pad(int padNum, vd netandlocation) {
    this->pad[padNum] = netandlocation;
    this->numP++;
    return;
  }

  /**
   * Helper function to return the number of nets.
   * @see get_numG()
   * @see get_numP()
   * @return The number of nets in the ASIC.
   */
  int get_numN() {
    return this->numN;
  }

  /**
   * Helper function to return net keys.
   * @see get_gateKeys()
   * @see get_padKeys()
   * @return A vector with the net-id of nets in the ASIC.
   */
  vi get_netKeys() {
    vi v;
    for(map<int,vvi>::iterator it = this->nets.begin(); it != this->nets.end(); ++it) {
      v.push_back(it->first);
    }
    return v;
  }

  /**
   * Helper function to return number of net connections
   * @param netNum The net-id of the net for which the information is desired.
   * @return Total number of gates and pads, combined, connected to a net of
   * net-id 'netNum'.
   */
  int get_numNetConns(int netNum) {
    int total = nets[netNum][0].size() + nets[netNum][1].size();
    return total;
  }

  /**
   * Helper function to return gate connections to a net
   * @param netNum The net-id of the net for which the information is desired.
   * @see get_netPadConns(int netNum)
   * @return A vector of gate-ids of the gates connected to the net of
   * net-id 'netNum'.
   */
  vi get_netGateConns(int netNum) {
    return nets[netNum][0];
  }

  /**
   * Helper function to return pad connections to a net
   * @param netNum The net-id of the net for which the information is desired.
   * @see get_netGateConns(int netNum)
   * @return A vector of pad-ids of the pads connected to the net of
   * net-id 'netNum'.
   */
  vi get_netPadConns(int netNum) {
    return nets[netNum][1];
  }
  
  /**
   * Helper function which makes a new net, if needed, and appends a connection to the net 'netnum'
   * @param netNum The net-id of the net to be added.
   * @param connection The gate-id/ pad-id of the gate/pad to which the net is
   * connected to.
   * @param gateorpad It shows if the given connection is a gate or a pad. It
   * is 0 for gate, 1 for pad.
   * @see add_gate(int gateNum, vi listofconnections)
   * @see add_pad(int padNum, vd netandlocation)
   */
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

  /**
   * Helper function which adds location values for given gate keys
   * @param x A vector containing x-coordinates of gates in the same order as
   * the gate-ids in the vector 'gatekeys'
   * @param y A vector containing y-coordinates of gates in the same order as
   * the gate-ids in the vector 'gatekeys'
   * @param gatekeys A vector containing gate-ids of the gates for which
   * location is given and to be updated.
   * @param bound The minimum and maximum values of the x, y coordinates
   * desired. It is an array of 4 numbers, [x_min, x_max, y_min, y_max]. The 
   * argument is not necessarily used.
   * @see get_locations(vi gatekeys)
   * @return True if the operation is successful, false otherwise. It is false
   * if there is a mismatch amongst sizes of 'x', 'y', and 'gatekeys'.
   */
  bool add_location(vd x, vd y, vi gatekeys, int bound[4]) {
    int l;
    double xloc, yloc;
    if ((gatekeys.size() == x.size()) && (x.size() == y.size())) {
      for (l = 0; l < gatekeys.size(); l++) {
        xloc = x[l];
        yloc = y[l];
        //if (x[l] < bound[0]) xloc = bound[0]; // xmin
        //if (x[l] > bound[1]) xloc = bound[1]; // xmax
        //if (y[l] < bound[2]) yloc = bound[2]; // ymin
        //if (y[l] > bound[3]) yloc = bound[3]; // ymax
        this->gateX[gatekeys[l]] = xloc;
        this->gateY[gatekeys[l]] = yloc;
      }
      return true;
    }	else {
      return false;
    }
  }

  /**
   * Helper function which gets the location values for given gate keys
   * @param gatekeys A vector containing gate-ids of the gates for which
   * location data is needed.
   * @see add_location(vd x, vd y, vi gatekeys, int bound[4])
   * @return A vector of x, y coordinates for each gate-id in vector 'gatekeys'
   * , in the same order as the gate-ids in 'gatekeys'. 
   */
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

  /**
   * Helper function which prints the locations of all gates in the present core.
   * @see print_all_pads()
   */
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

  /**
   * Helper function which prints the pads in the present core.
   * @see print_all_locations()
   */
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
};

/**
 * Function which creates an object of class "mothercore" from a given file.
 * @param filename Pointer to the name of a file where the input is located.
 * @return An object of class mothercore.
 */
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

/**
 * Function which writes the location of gates in an object of class "mothercore" to a file.
 * @param core Pointer to an object of class "mothercore".
 * @param filename Pointer to the name of a file where the output is to be
 * written.
 * @see writebackpads(mothercore *core, char *filename)
 */
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

/**
 * Function which writes the location of pads in an object of class "mothercore" to a file.
 * @param core Pointer to an object of class "mothercore".
 * @param filename Pointer to the name of a file where the output is to be
 * written.
 * @see writeback(mothercore *core, char *filename)
 */
void writebackpads(mothercore *core, char *filename) {
  ofstream out(filename);
  streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!

  vi keyP = core->get_padKeys();
  vd padloc;
  int i, padkey;
  for (i = 0; i < keyP.size(); i++) {
    padkey = keyP[i];
    padloc = core->get_padCoords(padkey);
    cout << padkey << " " << padloc[0] << " " << padloc[1] << endl;
  }

  cout.rdbuf(coutbuf); //reset to standard output again
  return;
}

/**
 * Function which solves a sparse matrix, of the form Ax=b, using coo_matrix 
 * class from solver.h.
 * @param R Vector containing non-zero row values of the matrix A, in order.
 * @param C Vector containing non-zero column values of the matrix A, in order.
 * @param V Vector containing non-zero values of the matrix A, in order.
 * @param ba Vector containing the b vector in the matrix form Ax=b.
 * @return A vector containing values of the solved vector x in the form Ax=b.
 */
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

/**
 * Function which generates the A matrix and b vector for each coordinate, from
 * an object of the class mothercore and sends it to solve for solving. The
 * result is written back to the mothercore object.
 * @param core Pointer to an object of class "mothercore".
 * @param bound The minimum and maximum values of the x, y coordinates
 * desired. It is an array of 4 numbers, [x_min, x_max, y_min, y_max].
 * @see solve(vi R, vi C, vd V, vd ba)
 * @return True if there are no errors, false otherwise.
 */
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
          A[gateorder[gates[i]]][gateorder[gates[j]]] += -weight;
          A[gateorder[gates[j]]][gateorder[gates[i]]] += -weight;
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

  for(i = 0; i < G; ++i) {
    delete [] A[i];
  }
  delete [] A;
  delete [] bx;
  delete [] by;

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

/**
 * Function which sorts given gates according to their locations, horizontally 
 * or vertically. Horizontally if hORv = 0; Vertically if hORv = 1
 * @param core Pointer to an object of class "mothercore".
 * @param gatekeys Vector of gate-ids of gates that need to be sorted.
 * @param hORv It is used to decide if the sorting is done based on x-coordinate or
 * y-coordinate. If hORv = 0, x-coordinate is used, and if it is 1,
 * y-coordinate is used.
 * @return A vector of 2 vectors. First vector contains the gate-ids which are
 * on the lower values of the sorting coordinate. The second vector contains 
 * the gate-ids which are on the higher values of the sorting coordinate.
 */
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

/**
 * Function which updates the coordinates of the given x, y coordinates
 * according to the given bound.
 * @param padlocation Pointer to a vector containing x, y coordinates.
 * @param bound The minimum and maximum values of the x, y coordinates
 * desired. It is an array of 4 numbers, [x_min, x_max, y_min, y_max].
 */
void update_coordinates(vd *padlocation, int bound[4]) {
  if ((*padlocation)[0] < (double)bound[0]) (*padlocation)[0] = (double)bound[0]; // xmin
  if ((*padlocation)[0] > (double)bound[1]) (*padlocation)[0] = (double)bound[1]; // xmax
  if ((*padlocation)[1] < (double)bound[2]) (*padlocation)[1] = (double)bound[2]; // ymin
  if ((*padlocation)[1] > (double)bound[3]) (*padlocation)[1] = (double)bound[3]; // ymax
  return;
}

/**
 * Function which contains the given gate-ids within the given bound, creates
 * virtual pads, and runs the resulting mothercore object.
 * @param core Pointer to an object of class "mothercore".
 * @param gatekeys Vector of gate-ids of gates that need to be contained within
 * the given bound.
 * @param bound The minimum and maximum values of the x, y coordinates
 * desired. It is an array of 4 numbers, [x_min, x_max, y_min, y_max].
 * @param hORv It is used to correct the coordinate of some virtual pads. hORv =
 * 0 means that the present bound is the result of a horizontal cut. hORv = 1
 * means that the present bound is the result of a vertical cut. 
 * @param lORr It is used to correct the coordinate of some virtual pads. lORr =
 * 0 means that the present gate-ids are present in the lower part of a cut. 
 * lORr = 1 means that the present gate-ids are present in the higher part of a cut.
 * @return The new locations of the gates whose gate-ids are in "gatekeys".
 */
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

/**
 * Function which recursively calls itself to place the given gates within the
 * bound. The bound keeps shortening as the depth of the recursion increases.
 * Also the number of gates in each bound decreases as the depth of the 
 * recursion increases. It aims to find a uniform distribution of the gates in
 * all the divisions of the ASIC.
 * @param core Pointer to an object of class "mothercore".
 * @param gatekeys Vector of gate-ids of gates that need to be placed within
 * the given bound.
 * @param bound The minimum and maximum values of the x, y coordinates
 * desired. It is an array of 4 numbers, [x_min, x_max, y_min, y_max].
 * @param n The level of the iteration. A nth iteration means that there are
 * 2^n divisions in the ASIC.
 */
void place(mothercore *core, vi gatekeys, int bound[4], int n) {
  if (n >= 2) return;
  else {
    cout << endl;
    cout << endl;
    cout << "n = " << n << endl;
    cout << endl;
    int nnext = n + 2;

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
    //core->print_all_locations();

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
    //core->print_all_locations();

    place(core, left_topbottomGates[0], leftbottom_bound, nnext);
    place(core, left_topbottomGates[1], lefttop_bound, nnext);
    place(core, right_topbottomGates[0], rightbottom_bound, nnext);
    place(core, right_topbottomGates[1], righttop_bound, nnext);

    return;
  }
}

/**
 * The main function: It creates a new mothercore object from a file and runs
 * place for recursive placing. Final placement is written back to an output
 * file using writeback and writebackpads.
 */
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "Invalid Command: Run as `executable filename`" << endl;
    return -1;
  }
  mothercore core = create(argv[1]);
  int initial_bound[4] = { 0, 100, 0, 100 };
  solveforx(&core, initial_bound);

  vi keyG = core.get_gateKeys();
  place(&core, keyG, initial_bound, 1);
  //core.print_all_locations();

  writeback(&core, "out.txt");
  writebackpads(&core, "outpad.txt");
  return 0;
}
