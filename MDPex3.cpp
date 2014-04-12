#include <iostream>
#include <string>
#include<map>
#include <vector>
#include <stack>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <sstream>
#include <utility>
#include<set>
#include<cmath>
#include <cassert>
#include "sparselib.h"

#include "defs.h"
#include "graph_class.h"
#include "utils.h"
#include "mdp_class.h"


void buildMDP(MDP &m, int probabs){
	map<string,state*> seen;
	//Build Start State *ss
	//Use a Stack for DFS, put ss on the stack
	//As long as stack is non-empty
	//Take the top, 
	//	1. If its an Env State: for each possible map entry in eimap, find the next state and check if its already in seen
	//		Establish the action based on probab info from MC: 1 state, Weight = 0, Also implement the constraint that r_i will 
	//		not become 0 unless c_i is reached and if currFloor is c then r_c should be 0.  
	//	2. If its a player state: For each value of "CurrFloor" find next Env state and check if its already in seen establish all actions 
	//		with weights given by MP Spec: number of active requests
	statevars1 sva1(false);//Env Type
	m.states[0] = new state(0,false,sva1);//Env State
	seen[m.states[0]->sv.idstr()] = m.states[0];
	stack<state*> st;//Only Env states to be put onto the stack
	st.push(m.states[0]);
	int cid=1;//current state id, to be incremented everytime we create a state
	while(!st.empty()){
		state* cs = st.top();
		assert(cs->sv.type == false);//Env state
		st.pop();
		assert(cs->actions.size() == 0);
		action* tact = new action(cs->id);
		cs->actions.push_back(tact);	
		//map<int, vector <bool> >::iterator it;
		//cs->sv.genEnvInpMap();
		//Assert((int)pow(2.0,NF-cs->sv.num1s()) == eimap.size(), "Matching sizes 2^num0s");
		//one should go for only minterms of those requests which are 0

        //Remember that previous step had an error? repair variable! 
        double perr = pall;
        if(cs->sv.repair) perr = pall/2;
            
        f(it,2){
			statevars1 svt;
            svt.setVals(cs->sv);
            if(svt.allscheds()) svt.resetscheds(); 
            svt.type = true;//player type state
            svt.setreq(it==1);//0 - no error, 1 - error
			string chkstr= svt.idstr();
			state *ns;
			
			bool cont=false;
			if(seen.find(chkstr) == seen.end()){
				//didn't find the state already
				//create new state
				ns = new state(cid,false,svt);
				m.states[cid] =ns;
				seen[chkstr] = ns;
				cid++;
			}
			else{
				//found the state already
				ns = seen[chkstr];
				cont =true;
			}
			int rctr=0;
			if(it == 1)
				cs->actions[0]->addEdge(cs->id,ns->id,perr,0);//cs->sv.getwt() for the next edge system to environment
			else 
				cs->actions[0]->addEdge(cs->id,ns->id,1-perr,0);//cs->sv.getwt() for the next edge system to environment
			assert(cs->actions[0]==tact);
			ns->incoming.push_back(cs->actions[0]->edges[it]);
			if(cont) continue;
			//now we have the next Player state
			//Try to do all possible prod-step changes fixing the reqs 
			//if no error: if not last step -> try to take next step for this product, reward = 0 or schedule a new poduct (choice 0 to NF-1)
			//if no error: if last step -> finish producing this product and give the reward and schedule a new product (choice 0 to NF-1)
            vector< sysvars > vsv;
            if(ns->sv.error){
                //if error: give repair and reset the manufacture to prod=-1, reward = 0
			    sysvars sysv(-1,0,true,0);
			    vsv.push_back(sysv);
            }
            else{
                int rew=0;
                if(ns->sv.prod >= 0 && ns->sv.step == steps[ns->sv.prod]){
                    //last step
                    rew = ns->sv.getwt();
                }
                else if(ns->sv.prod >= 0 ){
                    //can also increment the step or start anew
                    if(!ns->sv.sched[ns->sv.prod]){//if not scheduled, try to schedule
                        rew = 0;
                        sysvars sysv(ns->sv.prod,ns->sv.step+1,false,0);
                        vsv.push_back(sysv);
                    }
                }
                f(pctr,NF){
                    if(!ns->sv.sched[pctr]){//if not scheduled, try to schedule
                        sysvars sysv(pctr,0,false,rew);
                        vsv.push_back(sysv);
                    }
                }
            }
            for(int ctr =0; ctr < vsv.size(); ctr++){
                statevars1 svt2(false);//env state
				svt2.setreqp(vsv[ctr].prod,vsv[ctr].step,vsv[ctr].repair);
				if(vsv[ctr].prod >= 0){
				    //cout<<svt2.str()<<endl;
				    svt2.schedadd(vsv[ctr].prod,ns->sv.sched);//add this product as scheduled
                    //cout<<vsv[ctr].prod<<"->"<<svt2.str()<<endl;
				}
				else if(vsv[ctr].prod == -1){
                    //error state
                    svt2.resetscheds();
				}
				state * nns;
				string nstr=svt2.idstr();
				if(seen.find(nstr) != seen.end()){
					nns = seen[nstr];
				}
				else{
					nns = new state(cid,false,svt2);
					st.push(nns);
					m.states[cid] = nns;
					seen[nstr] = nns;
					cid++;
				}
				action* nact = new action(ns->id);
				ns->actions.push_back(nact);
				assert(ns->actions[rctr] ==nact);
				ns->actions[rctr]->addEdge(ns->id,nns->id,1.0,vsv[ctr].reward);
				ns->actions[rctr]->setCost();
				nns->incoming.push_back(ns->actions[rctr]->edges[0]);
				rctr++;
            }
		}
		tact->setCost();
	}
	
}

string sNF = toSTR(NF);

int main(int argc,char *argv[]){
    //Give number of products(np), number of stages for each product and environment MC probability - p 
    
	clock_t start = clock();
	MDP myMDP;
	stringstream fs;
	string sDN = "e3_" +sNF;
	fs<<"_"<<sDN;
	string fapp = fs.str();
	buildMDP(myMDP,MDPTYPE);
	if(DFLAG) myMDP.sanityCheck();
	cout<<"Starting Policy Iteration"<<endl;
	srand((unsigned int) time(NULL));
	myMDP.markStartReachable();
	if(DFLAG) myMDP.printDotFile("./" + sDN + "/mdp" + fapp + ".dot");
	int mdpsize = myMDP.states.size();
	double x = myMDP.ValueStartState();
	if(DFLAG) myMDP.printDotFileColor("./" + sDN + "/mdp" + fapp + "color.dot");
	if(DFLAG) myMDP.printDotFileStratOnly("./" + sDN + "/mdp" + fapp + "only.dot");
	//if(DFLAG) myMDP.print();
	cout<< "Value of Starting State = "<<x<<endl;
	int msize = myMDP.keepStrategy("./" + sDN + "/strategy" + fapp + ".dot");
	clock_t end = clock();
	double duration = ( end - start ) / (double) CLOCKS_PER_SEC;
	ofstream finalout(("./table_" + fapp + ".txt").c_str());
	cout<<"Num_processes"<<'\t'<< "max_states" <<'\t'<<"mdp_size"<<'\t'<<"mealy_machine_size"<<'\t'<<"value"<<'\t'<<"duration"<<endl;
	cout<<NF<<'\t'<<mdpsize<<'\t'<<msize<<'\t'<<x<<'\t'<<duration<<endl;
	finalout<<NF<<'\t'<<mdpsize<<'\t'<<msize<<'\t'<<x<<'\t'<<duration<<endl;
	finalout.close();
	return 0;
}
