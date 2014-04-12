#define DBG false
#define DBGP1 false
#define DBGIO false
#define ulli unsigned long long int
#define f(i,n) for(unsigned int i=0;i< (unsigned int)n;i++)
#define fe(i,s) for(vector<edge*>::iterator i=s.begin();i!=s.end();i++)
#define fi(i,s) for(set<int>::iterator i=s.begin();i!=s.end();i++)
#define fs(it,R,T) for(map<int,T*>::iterator it = R.begin();it != R.end();it++)
#include "statevars.h"

class edge {
	public:
	int from,to; //state number ids
	double prob; // for MDP
	int wt; //for Mean Payoff MDP
	//int thread;
	edge(int Id1,int Id2, double p,int w){
		from = Id1;
		to = Id2;
		prob = p;
		wt = w;
		label="";
	}
	string label;
};

class action {
	public:
	int sid; //state to which it belongs
	//int a1,a2;
	vector <edge *> edges;
	double cost;
	action(int SID) {
		sid = SID;
	}
	void addEdge(int Id1,int Id2, double p,int w){
		edges.push_back(new edge(Id1,Id2,p,w));
	}
	void setCost(){
		cost = 0;
		Assert(edges.size() > 0, toSTR((int) __LINE__));
		f(i,edges.size()){
			edge * e = edges[i];
			cost += e->prob * e->wt;
		}
	}
};




class state{
	public:
	int id; //own id
	vector <action *> actions; //assume atleast one action in each state
	int chosen;
	bool bad;
	bool good;
	double value;
	bool badReachable;
	bool goodReachable;
	bool startReachable;
	bool mark,DFSmark;
	//State Variables
    statevars1 sv;//Default is always Env State
    state(int Id,bool b,statevars1 &s){
		id = Id;
		chosen = 0;
		bad =b;
		value =0.0;
		badReachable = b;
		if(id==0) startReachable = true;
		else startReachable = false;
		good = false;
		DFSmark = false;

		//State Vars
		sv.setVals(s);
	}
	void printComp(){
		//cout<<" ("<<OID1<<","<<OID2<<")";
		cout<<sv.str()<<endl;
	}
	vector <edge *> incoming; //to be set before MC computation every time!
	void print(){
		cout<<"STATE "<<id<<sv.str()<<"    "<<actions.size()<<" actions. BadProb = "<<value<<endl;
		f(j,actions.size()){
			if(j==chosen) cout<<"C";
			cout<<"\tACTION "<<j<<endl;
			f(k,actions[j]->edges.size()){
				cout<<"\t\tEDGE "<<actions[j]->edges[k]->to<<" "<<actions[j]->edges[k]->prob<<" "<<actions[j]->edges[k]->wt<<endl;
			}
		}
	}
	void printsucc(){
		cout<<"STATE "<<id<<sv.str()<<" -> [ ";
			f(k,actions[chosen]->edges.size()){
				if(k!=0) cout<<",";
				cout<<actions[chosen]->edges[k]->to;
			}
			cout<<" ]"<<endl;
	}
};


