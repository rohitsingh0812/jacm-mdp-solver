using namespace std;
/*bool nonConsecutive(bool from[NF], vector<bool> &to){
	Assert(to.size() == NF, "to vector size mismatch");
	//check that 0's in from remained 0's
	set<int> chng;
	f(i,NF){
		if(!from[i]){ 
			Assert(!to[i],"to[i] should be 0 if from[i] was 0");
			continue;
		}
		else{
			if(to[i]) continue;
			else{
				if(i==NF-1 && chng.find(0) != chng.end()) return false;
				if(i>0 && chng.find(i-1) != chng.end()) return false;
				chng.insert(i);
			}
		}
	}
	return true;
}*/


string strchange(statevars1 &from, statevars1 &to){
    //what change happened?
	stringstream ss;
	f(i,NF) ss<<((int)to.sched[i]);
	ss<<"|";
    if(from.prod != to.prod)
        ss<<from.prod<<"p"<<to.prod<<"|";
    else{
        ss<<from.prod<<"p";
        if(from.step != to.step){
            ss<<from.step<<"s"<<to.step;
        }
    }
    if(to.error) ss<<"e";
    if(from.repair) ss<<"r";
	return ss.str();
}

