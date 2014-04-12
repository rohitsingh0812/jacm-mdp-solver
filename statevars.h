class statevars1{
    public:
	bool type;//0 for environment and 1 for player

	bool sched[NF]; //Already scheduled processes

    //step ranges from 0 to NK-1, the assembly can work on only one product at a time	
    int prod;//which product is being assembled, 0 to NF-1
    int step;//Current state of production,prod i gives NF-1 reward and takes NF-i+1 steps (0 to NF-i)

    bool error; //environent has decided to produce an error, can be true only after an environment state
    bool repair; //repair state: can be made true only after a system state

	//int wtcase; //which product was synthesized last
	//wtcase >= 0 and wtcase <= NF-1
	
	statevars1(){
		statevars1(false);
	}
	statevars1(bool t){
		prod = -1;
		error = false;
		repair = false;
		step = 0;
		type = t;
		f(i,NF) sched[i] = false;
		//wtcase = NF-1;
	}
	statevars1(bool t, int sched_prod){
	Assert(sched_prod >=0 && sched_prod < NF, "invalid product scheduled");
		prod = sched_prod;
		error = false;
		repair = false;
		step = 0;
		type = t;
		f(i,NF) sched[i] = false;

		//wtcase = NF-1;
	}
    bool allscheds(){
        f(i,NF){
            if(!sched[i]) return false;
        }
        return true;
    }
    void resetscheds(){
        //Assert(allscheds(),"Can be called only when all sched array values are true");
        f(i,NF) sched[i] = false; 
    }
    void schedadd(int &p, bool sc[NF]){
        Assert(!sc[p], "Can't schedule product that was already scheduled!");
        f(i,NF) sched[i] = sc[i];
        sched[p] = true;
    } 
	int getwt(){//These are rewards! Not costs! The rewards are maximized
	    //if error in system state, reward = 0
	    if(error) return 0;
	    if(prod == -1) return 0;
	    Assert(prod >= 0 && prod <= NF-1, "prod value not valid");
		if(step == steps[prod]){
            //Reward only when the production ends
            if(sched[prod])
                return NF-rewards[prod]+1;
            else
                return rewards[prod];
		}
		return 0;
	}
	void setreq(bool e){
		error=e;	
	}
    void setreqp(int p, int s, bool r){//TODO include sched?
		prod = p;
		step = s;
		repair = r;
		if(r) error = false; 
	}

	void setVals(const statevars1 &s){
        prod = s.prod;//which product is being assembled, 0 to NF-1
        step = s.step;//Current state of production,There are NK steps (0 to NK-1)
        error = s.error; //environent has decided to produce an error, can be true only after an environment state
        repair = s.repair; //repair state: can be made true only after a system state
	    //wtcase = s.wtcase; //which product was synthesized last
    	type = s.type;
    	f(i,NF) sched[i] = s.sched[i];
	}
	string str(){
        stringstream ss;
		if(type) ss<<'p';
		else ss<<'e';
		f(i,NF) ss<<((int)sched[i]);
		ss<<"("<<prod<<"_"<<step<<")";
		if(error) ss<<"Err";
		if(repair) ss<<"Rep";
		return ss.str();
	}
	string idstr(){
        stringstream ss;
		if(type) ss<<'p';
		else ss<<'e';
        f(i,NF) ss<<((int)sched[i]);
		ss<<"("<<prod<<"_"<<step<<")";
		if(error) ss<<1;
		else if(repair) ss<<2;
        else ss<<0;
		return ss.str();
	}
	string reqstr(){
        stringstream ss;
        f(i,NF) ss<<((int)sched[i]);
		ss<<"("<<prod<<"_"<<step<<")";
        if(error) ss<<1;
		else if(repair) ss<<2;
        else ss<<0;
		return ss.str();
	}
/*	void genEnvInpMap(){
		eimap.clear();
		p1probmap.clear();
		p2probmap.clear();
		
		double sumup1= 0.0,sumup2=0.0;
		int n0s = NF - num1s();
		int p2n0s = (int)pow((double)2,(int)n0s);
		f(i,p2n0s){
			vector<bool> v;
			double prob1=1.0;
			double prob2=1.0;
			int p2k =1;
			int k=0;
			f(j,NF){
				 if(req[j]){
					 v.push_back(true);
					 continue;
				 }
				 if ( i & p2k ){ //kth bit of i is 1
					 v.push_back(true);
					 prob1 *= pall;
					 prob2 *= pdiff[j];
				 }
				 else{
					 v.push_back(false);
					 prob1 *= (1.0-pall);
					 prob2 *= (1.0-pdiff[j]);
				 }
				 p2k *=2;
				 k++;
			}
			eimap[i] =v;
			p1probmap[i] = prob1;
			p2probmap[i] = prob2;
			sumup1+=prob1;
			sumup2+=prob2;
		}
		Assert(abs(sumup1 -1.0) < EPSCHECK,"Total probability is not 1.0");
		Assert(abs(sumup2 -1.0) < EPSCHECK,"Total probability is not 1.0");
	}

	void genSysInpMap(){
		simap.clear();
		
		int n1s = num1s();
		int p2n1s = (int)pow((double)2,(int)n1s);
		f(i,p2n1s){
			vector<bool> v;
			int p2k =1;
			int k=0;
			f(j,NF){
				 if(!req[j]){
					 v.push_back(false);
					 continue;
				 }
				 if ( i & p2k ){ //kth bit of i is 1
					 v.push_back(false);
				 }
				 else{
					 v.push_back(true);
				 }
				 p2k *=2;
				 k++;
			}
			simap[i] =v;
		}
	}*/
};

