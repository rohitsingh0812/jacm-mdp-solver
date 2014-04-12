#include <iostream>
#include "inputs3.h"
using namespace std;
void Assert(bool b,string msg){
	if(!b){
		cout<<"Error :" <<msg<<endl;
		exit(1);
	}
}

int BIGC = NF;//(WEIGHTED_REWARDS)? NF*NF:NF;
int gitr1 = 4000;
int gitr2 = 12000;
#define ulli unsigned long long int

ulli tNF = (ulli) pow(2.0,NF);
map<int, vector <bool> > eimap;
map<int, vector <bool> > simap;

map<int, double> p1probmap;
map<int, double> p2probmap;

template <typename T>
T StringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

//double EPS = 1e-10;
double EPS = 1e-8;
double EPSCHECK = 1e-5;
bool PROBINPUT = false;
int NNZ = 10000000;
string toSTR(int n){
	string s;
	stringstream out;
	out << n;
	return out.str();
}
class sparseMat {
	public:
	vector< pair<int,int> > ij;
	vector<double> v;
	int N;
};

#define WEIGHTED_REWARDS true
#define EXTRA 0
#define DFLAG true
#define pl() if(DFLAG) cout<<__LINE__<<endl

