#include <cstdio>
#include <cstdlib>
#include <valarray>
#include <vector>
using namespace std;

// this is the header file of solver class
#include "solver.h"

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

int main(int argc, char *argv[]) {
  double arrbx[] = {1, 1, 2};
  int arrR[]    = {0, 0, 1, 1, 1, 2, 2};
  int arrC[]    = {0, 1, 0, 1, 2, 1, 2};
  double arrV[] = {4.0, -1.0, -1.0,  4.0, -1.0, -1.0, 4.0};
  
  vector<double> b (arrbx, arrbx + sizeof(arrbx)/sizeof(arrbx[0]));
  vector<int> R (arrR, arrR + sizeof(arrR)/sizeof(arrR[0]));
  vector<int> C (arrC, arrC + sizeof(arrC)/sizeof(arrC[0]));
  vector<double> V (arrV, arrV + sizeof(arrV)/sizeof(arrV[0]));
  vector<double> x;
  
  x = solve(R, C, V, b);

  // should get 
  // x = 
  // 0.375 
  // 0.5 
  // 0.625 
  
  cout << "x = ";
  for (vector<double>::const_iterator i = x.begin(); i != x.end(); ++i) { 
    cout << endl;
    cout << *i << ' ';
  }
  cout << endl;
  return 0;
}
