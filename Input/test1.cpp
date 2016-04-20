#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<map>
#include<vector>

using namespace std;

typedef vector<vector<int> > vvi;

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
								cout<<"numG = "<<numG<<endl;
							}
							else if (word_no==2){
								numN=num;
								cout<<"numN = "<<numN<<endl;
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
								cout << "Gate no. " << key << " is connected to net no. ";
                				vector <int> vec;
								gate[key]=vec;
							}
							else if (word_no==2) {
								no_of_nets=num;
							}
							else if ((word_no>2)&&(word_no<(no_of_nets+3))) {
								gate[key].push_back(num);
								nets[num][0].push_back(key);
								cout << gate[key][word_no-3] <<' ';
							}
						}
						else if (line_no==(numG+1)) {
							if (word_no==1){
								numP=num;
								cout << "numP = " << numP <<endl;
							}
						}
						else if ((line_no>(numG+1))&&(line_no<(numG+2+numP))){
							if (word_no==1) {
								key = num;
								cout << "Pad no. " << key << " is connected to net no. ";
								vector <int> vec;
								pad[key]=vec;
							}
							else if (word_no==2) {
								pad[key].push_back(num);
								cout << pad[key][0] <<" at coordinates (";
								nets[num][1].push_back(key);
							}
							else {
								pad[key].push_back(num);
								if (word_no==3) cout << pad[key][1]<<',';
								else if (word_no==4) cout << pad[key][2]<<')';
							}
						}
					}
				}
				cout<<endl;
				line_no++;
			}
		}
		file.close();
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
	}
};

int main(){

	mothercore core;
	core.create("toy1");
	return 0;
}