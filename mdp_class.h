using namespace std;
using namespace ::sparselib_load ;
using           ::sparselib::index_type ;
using namespace ::std ;
using ::sparselib_fun::maxval ;
using ::sparselib_fun::minval ;
using ::sparselib_fun::absval ;

map<edge*,bool> doneEdges;

class MDP {
	public:
/*	string getLabel(edge *e){
		int f = e->from;
		int t = e->to;
		f(i,NF){
			if(states[f]->sv.req[i]
		}
	}*/
	
	map <int,state*> states;
	int N;
	bool valid;
	//vector <int> policy; //stored in actions as in chosen field (default 0)
	MDP(){
		return;
	}
	~MDP(){
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			f(j,states[i]->actions.size()){
				f(k,states[i]->actions[j]->edges.size()){
					delete states[i]->actions[j]->edges[k];
				}
				delete states[i]->actions[j];
			}
			delete states[i];
		}
	}
	void print(){
		f(i,states.size()){
			cout<<"STATE "<<i<<" "<<states[i]->sv.str()<<endl;
			f(j,states[i]->actions.size()){
				cout<<"\tACTION "<<j<<endl;
				f(k,states[i]->actions[j]->edges.size()){
					cout<<"\t\tEDGE "<<states[i]->actions[j]->edges[k]->to<<" "<<states[i]->actions[j]->edges[k]->prob<<" "<<states[i]->actions[j]->edges[k]->wt<<endl;
				}
			}
		}
	}
	void badDebug(){
		f(i,states.size()){
			
			f(j,states[i]->actions.size()){
				
				f(k,states[i]->actions[j]->edges.size()){
					if(states[states[i]->actions[j]->edges[k]->to]->bad){
						cout<<"STATE "<<i<<endl;
						cout<<"\tACTION "<<j<<endl;
						cout<<"\t\tEDGE "<<states[i]->actions[j]->edges[k]->to<<" "<<states[i]->actions[j]->edges[k]->prob<<" "<<states[i]->actions[j]->edges[k]->wt<<endl;
					}
				}
			}
		}
	}
	ulli numStrats(){
		ulli temp=1;
		f(i,states.size()){
			if(states[i]->actions.size() > 1){ cout<<states[i]->actions.size()<<","; temp++; }
		}
		return temp;
	}
	int badStates(){
		int temp=0;
		f(i,states.size()){
			if(states[i]->bad){ temp++; }
		}
		return temp;
	}
	void multAb(double **A,double b[],double res[],int &n){
		f(i,n){
			res[i] =0;
			f(j,n){
				if(A[i][j] != 0 && b[j] != 0) res[i]+=A[i][j]*b[j];
			}
		}
	}
	void Add(double a[],double b[],double res[],int &n){
		f(i,n) res[i] =a[i]+b[i];
	}
	double absDiff(double *x,double *y,int &n){
		double max = 0,temp;
		f(i,n){
			temp=x[i]-y[i];
			if(temp <0) temp = -temp;
			if(temp > max) max=temp;
		}
		return max;
	}
	void solveLinEq(double **A,double b[],double x[],int n){
		//solve x=Ax+b
		//x=0 initialized
		f(i,n) x[i] =0.0;
		double* t1 = new double[n];
		double* t2 = new double[n];

		double val;
		int itr =0;
		while(1){
			itr++;
			multAb(A,x,t1,n);
			Add(t1,b,t2,n);
			val = absDiff(t2,x,n);
			f(i,n) x[i] = t2[i];
			//cout<<itr<<" : "<<val<<endl;
			if(val<EPS) break;			
		}
		delete[] t1;
		delete[] t2;
	}
	void sparseMatMethod(int ctr,map <int,int> &M, map <int,int> &N ){
		ofstream fout("tempmatrix.sce");
		sparseMat Matrix;
		vector <double> b,x;
		f(i,3*ctr){ 
			b.push_back(0.0); x.push_back(0.0);
			//f(j,ctr) Mat[i][j] = 0;
		}
		Matrix.N = 3*ctr;
		f(i,ctr){
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				if(N.find(e->to) != N.end()){
					{pair<int,int> p (i+1,N[e->to]+1);
					Matrix.ij.push_back(p);
					Matrix.v.push_back(e->prob);}
					{pair<int,int> p (ctr+i+1,ctr+N[e->to]+1);
					Matrix.ij.push_back(p);
					Matrix.v.push_back(e->prob);}
					{pair<int,int> p (2*ctr+i+1,2*ctr+N[e->to]+1);
					Matrix.ij.push_back(p);
					Matrix.v.push_back(e->prob);}
				}
				
			}
			b[ctr+i] = states[M[i]]->actions[states[M[i]]->chosen]->cost;
		}
		f(i,ctr){
			{pair<int,int> p (i+1,ctr+i+1);
			Matrix.ij.push_back(p);
			Matrix.v.push_back(-1);}
			{pair<int,int> p (ctr+i+1,2*ctr+i+1);
			Matrix.ij.push_back(p);
			Matrix.v.push_back(-1);}
			
		}
		f(i,3*ctr){
			{pair<int,int> p (i+1,i+1);
			Matrix.ij.push_back(p);
			Matrix.v.push_back(-1);}
		}
		fout<<"mn = [ "<<3*ctr<<" "<<3*ctr<<" ];"<<endl<<endl;
		fout<<"v = [ ";
		f(i,Matrix.v.size()){
			if(i!=0) fout<<";";
			fout<<Matrix.v[i];
		}
		fout<<"];"<<endl<<endl;
		
		fout<<"ij = [ ";
		f(i,Matrix.ij.size()){
			if(i!=0) fout<<";";
			fout<<Matrix.ij[i].first<<" "<<Matrix.ij[i].second;
		}
		fout<<"];"<<endl<<endl;
		
		//fout<<"A = sparse(ij,v,mn)\n"<<endl;
		
		/*fout<<"b = [ ";
		f(i,b.size()){
			if(i!=0) fout<<";";
			fout<<b[i];
		}
		fout<<"];"<<endl<<endl;
		fout<<"[x,kerA] = linsolve(A,b);\n disp('OUTPUT');\n fid = mopen('x.txt', \"w\");\nif (fid == -1)\nerror('cannot open file for writing');\nend\nfor i = x\nmfprintf(fid, \"%f \\n\", i);\nend\nexit\n"<<endl;
		//fout<<"[x,kerA] = linsolve(A,b)\n\n disp('OUTPUT')\n\n disp(x)\n\n disp(kerA)\n\nexit"<<endl<<endl;*/
		//fout.close();
		//exit(0);
		/*system("scilab -nogui -f tempmatrix.sce > tempio.txt");
		//system("cat tempio.txt");
		ifstream fin("x.txt");
		//string c;
		//while(fin>>c && c != "OUTPUT"){}
		//cout<<"[";
		f(i,ctr){ 
			fin>>x[i];
			states[M[i]]->value = x[i];
			//should work!
			//cout<<x[i]<<" ";
			//assert(x[i] >= 0 && x[i] <= 1);
		}
		//cout<<"]"<<endl;
		fin.close();*/
	}
	void fullMatMethod(int &ctr,map <int,int> &M, map <int,int> &N ){
		double** Mat = new double*[ctr];
		f(i,ctr) Mat[i] = new double[ctr];
		//cout<<350<<endl;
		double *b = new double[ctr];
		double *x = new double[ctr];
		f(i,ctr){ 
			b[i]=x[i]=0.0;
			f(j,ctr) Mat[i][j] = 0;
		}
		//cout<<356<<endl;
		f(i,ctr){
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				if(N.find(e->to) != N.end()){
					Mat[i][N[e->to]] = e->prob;
				}
				else if(states[e->to]->bad){
					b[i]+= e->prob;
				}
			}
		}
		//cout<<365<<endl;
		solveLinEq(Mat,b,x,ctr);
		//cout<<367<<endl;
		if(DBG) cout<<"[";
		f(i,ctr){ 
			states[M[i]]->value = x[i];
			//should work!
			if(DBG) cout<<x[i]<<" ";
		}
		if(DBG) cout<<"]"<<endl;
		f(i,ctr) delete [] Mat[i];
		delete[] b;
		delete[] x;
		delete Mat;
	}
	void Identity(CCoorMatrix<double> & Id,int &n){
		Id.new_dim(n,n,n+1);
		f(i,n){
			Id.insert(i,i) = 1.0;
		}
	}
	void Diag(CCoorMatrix<double> &I,CCoorMatrix<double> &A,int &n){
		I.new_dim(n, n, n*2);
		double d;
		f(i,n){
			//cout<<i<<endl;
			d = A(i,i);
			if(d==0){
				cout<<i<<" : "<<n<<endl;
				assert(d!=0);
			}
			I.insert(i,i) = d;
		}
		I.internal_order();
	}
	void solveSparseLinEqSpl(CCoorMatrix<double> &B,CCoorMatrix<double> &A,sparselib::Vector < double> &b,sparselib::Vector < double> &x,int &n,double &alpha){
		
		//DPreco<double> I;
		//I.build(A);
		
		//Diag(I,A,n);
		//solve Ax=b
		//x=0 initialized
		//cout<<n<<endl;
		x=1.;
		//int iter;
		
		//bicgstab(A,b,x,I,EPS,100000,iter);
		//cg(A, b, x, I, EPS, 10000, iter) ;
		//gmres(A,b,x,I,EPS,100,10000,iter);
		//cout<<iter<<endl;
		sparselib::Vector < double> y(n);
		sparselib::Vector < double> z(n);
		double val=1;
		pl();
		while(val > EPS){
			y = A*x;
			y += b;
			z = B*y;
			//y = ((A * x) + b);
			val = sparselib::normi(x-z);
			cout<<setprecision(10)<<"VAL : "<<val<<endl;
			x=z;
		}

	}
	void solveSparseLinEq(CCoorMatrix<double> &A,sparselib::Vector < double> &b,sparselib::Vector < double> &x,int &n){
		
		//solve x=Ax+b
		//x=0 initialized
		//cout<<n<<endl;
		x=0.;
		sparselib::Vector < double> y(n);
		double val=1;
		pl();
		while(val > EPS){
			y = A*x;
			y += b;
			//y = ((A * x) + b);
			val = sparselib::normi(x-y);
			cout<<setprecision(10)<<"VAL : "<<val<<endl;
			x=y;
		}

	}
	void markStartReachable(){
		fs(itr,states,state){
			int i = itr->first;
			states[i]->startReachable = false;
		}
		states[0]->startReachable = true;
		int ctr=1;
		set <int> prevlevel;
		prevlevel.insert(0);
		while(1){
			set <int> thislevel;
			for(set<int>::iterator it = prevlevel.begin();it != prevlevel.end();it++){
				f(k,states[*it]->actions.size()){
					f(j,states[*it]->actions[k]->edges.size()){
						edge * e =states[*it]->actions[k]->edges[j];
						if(!states[e->to]->startReachable){
							states[e->to]->startReachable = true;
							ctr++;
							thislevel.insert(e->to);
						}
					}
				}
			}
			if(thislevel.empty()) break;
			prevlevel = thislevel;
		}
		cout<<"Reachable States = "<<ctr<<endl;
	}
	void markSpecialStartReachable(){
		f(i,states.size()) states[i]->startReachable = false;
		states[0]->startReachable = true;
		set <int> prevlevel;
		prevlevel.insert(0);
		while(1){
			set <int> thislevel;
			for(set<int>::iterator it = prevlevel.begin();it != prevlevel.end();it++){
				//f(k,states[*it]->actions.size()){
					int k = states[*it]->chosen;
					f(j,states[*it]->actions[k]->edges.size()){
						edge * e =states[*it]->actions[k]->edges[j];
						if(!states[e->to]->startReachable){
							states[e->to]->startReachable = true;
							thislevel.insert(e->to);
						}
					}
				//}
			}
			if(thislevel.empty()) break;
			prevlevel = thislevel;
		}
	}
	void sparselibMethod(int ctr,map <int,int> &M, map <int,int> &N ){
		CCoorMatrix<double> A;
		A.new_dim(ctr, ctr,100000);
		//cout<<350<<endl;
		sparselib::Vector < double> b(ctr);
		sparselib::Vector < double> x(ctr);	
		b=0.;
		x=0.;
		f(i,ctr){
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				if(N.find(e->to) != N.end()){
					A.insert(i,N[e->to]) = e->prob;
				}
				else if(states[e->to]->bad){
					b[i] += e->prob;
				}
			}
		}
		A.internal_order();
		solveSparseLinEq(A,b,x,ctr);
		if(DBG) cout<<"[";
		f(i,ctr){ 
			states[M[i]]->value = x[i];
			//should work!
			if(DBG) cout<<x[i]<<" ";
		}
		if(DBG) cout<<"]"<<endl;
	}
	void patRAND(CCoorMatrix<double> & Id,int ctr){
		srand((unsigned int)time(NULL));
		Id.new_dim(ctr,ctr,NNZ);
		f(i,ctr){
			Id.insert(i,i) = (rand()%2 + 1)/2.0;
		}
	}
	void solveSparseEqn(CCoorMatrix<double> &A,sparselib::Vector < double> &c,sparselib::Vector < double> &z,int &n,int ITR1,int ITR2){
		IdPreco<double> I;
		z=0.;
		int iter;
		gmres(A,c,z,I,EPS,ITR1,ITR2,iter);
		//cg(A,c,z,I,EPS,ITR2,iter);
		cout<< "Gmres ITERATIONS : "<<iter<<endl;
	}
	bool testSol(CCoorMatrix<double> &ImP,sparselib::Vector < double> &z,sparselib::Vector < double> &c,sparselib::Vector < double> &x,sparselib::Vector < double> &y,sparselib::Vector < double> &r){
		sparselib::Vector < double> t1(c.size());
		sparselib::Vector < double> t2(c.size());
		t1 = ImP * z;
		t2 = ImP * t1;
		t1 = ImP * t2;
		t2 = t1 -c;
		double val = sparselib::normi(t2);
		if(DBGIO) cout<<"VAL1 : "<<val<<endl;
		if(val > EPSCHECK) return false;
		t1 = ImP*x;
		val = sparselib::normi(t1);
		if(DBGIO) cout<<"VAL2 : "<<val<<endl;
		if(val > EPSCHECK) return false;
		t1 = ImP *y;
		t2 = x + t1 -r;
		val = sparselib::normi(t2);
		if(val > EPSCHECK) return false;
		if(DBGIO) cout<<"VAL3 : "<<val<<endl;
		t1 = y + ImP*z;
		val = sparselib::normi(t1);
		if(val > EPSCHECK) return false;
		if(DBGIO) cout<<"VAL4 : "<<val<<endl;
		return true;
	}
	bool testSol(CCoorMatrix<double> &ImP,sparselib::Vector < double> &x,sparselib::Vector < double> &y,sparselib::Vector < double> &r){
		sparselib::Vector < double> t1(x.size());
		sparselib::Vector < double> t2(x.size());
		t1 = ImP*x;
		double val = sparselib::normi(t1);
		if(DBGIO) cout<<"VAL(a) : "<<val<<endl;
		if(val > EPSCHECK) return false;
		t1 = ImP *y;
		t2 = x + t1 -r;
		val = sparselib::normi(t2);
		if(DBGIO) cout<<"VAL(b) : "<<val<<endl;
		if(val > EPSCHECK) return false;
		
		return true;
	}
	bool findXY(map<int,int> & M,map<int,int> &N,int n,sparselib::Vector < double> &x,sparselib::Vector < double> &y,int MAXITER){
		CCoorMatrix<double> P;
		P.new_dim(n, n,NNZ);
		sparselib::Vector < double> r(n);
		r=0.;
		f(i,n){
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				assert(i != N[e->to]);
				P.insert(i,N[e->to]) = e->prob;
			}
			r[i] = states[M[i]]->actions[states[M[i]]->chosen]->cost;
		}
		P.internal_order();
		//prob matrix and reward vector ready
		sparselib::Vector < double> z;//[MAXITER];
		sparselib::Vector < double> t;//[MAXITER];
		x = 0.;
		y = 0.;
		//f(i,MAXITER) z[i].new_dim(n);
		z.new_dim(n);
		t.new_dim(n);
		f(i,MAXITER){
			if(i%1000 == 0) {cout<<i<<",";
			cout.flush();}
			if(i==0) z =r;
			else z = P * z;
			x += z;
			
		}
		cout<<endl;
		x /= MAXITER;
		f(i,MAXITER){
			if(i%1000 == 0) {cout<<i<<",";
			cout.flush();}
			if(i==0) z =r;
			else { t = P * z; z += t; }
			y += z;
			y -= (i+1) * x; 
		}
		cout<<endl;
		
		//y=0.;
		//f(i,MAXITER){
		//	y += (z[i] - x);
		//}
		//have the value vector and deviation vector
		//test solution
		CCoorMatrix<double> ImP;
		ImP.new_dim(n, n,NNZ);
		f(i,n){
			ImP.insert(i,i) = 1;
			double sum = 1;
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				assert(i != N[e->to]);
				ImP.insert(i,N[e->to]) = -e->prob;
				sum -= e->prob;
			}
			assert(abs(sum) < 0.01);
		}
		ImP.internal_order();
		return testSol(ImP,x,y,r);
		
		if (testSol(ImP,x,y,r)) {
			int allOK = 0;
			f(i,n){
				states[M[i]]->value = x[i];
				if((EPSCHECK+x[i]) <= 0) {
					//cout<<"x["<<i<<"] = "<<x[i]<<endl;
					allOK++;
				}
			
			}
			if(allOK!=0){
				cout<<"Arbitrary Values !"<<endl;
				return false;
			}
			return true;
		} 
		else return false;
		
	}
	bool findValues(map<int,int> & M,map<int,int> &N,int n,sparselib::Vector < double> &x,sparselib::Vector < double> &y){
		CCoorMatrix<double> ImP;
		ImP.new_dim(n, n,NNZ);
		sparselib::Vector < double> r(n);
		sparselib::Vector < double> c(n);
		sparselib::Vector < double> z(n);	
		sparselib::Vector < double> t1(n);
		sparselib::Vector < double> t2(n);
		r=0.;
		z=0.;
		c=0.;
		f(i,n){
			ImP.insert(i,i) = 1;
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				assert(i != N[e->to]);
				ImP.insert(i,N[e->to]) = -e->prob;
			}
			r[i] = states[M[i]]->actions[states[M[i]]->chosen]->cost;
		}
		c = ImP * r;
		c= -c;
		ImP.internal_order();
		if(DBG) cout<<"Started Linear Equation Solving..."<<endl;
		//pl();
		solveSparseEqn(ImP,c,z,n,gitr1,gitr2);
		//pl();
		y = ImP * z;
		f(i,y.size()) y[i] = -y[i];
		t1 = ImP *z;
		t2 = ImP * t1;
		x = t2 + r;
		if(!testSol(ImP,z,c,x,y,r)) 
			return false;
		int allOK = 0;
		f(i,n){
			states[M[i]]->value = x[i];
			if((EPSCHECK+x[i]) <= 0) {
				//cout<<"x["<<i<<"] = "<<x[i]<<endl;
				allOK++;
			}
			
		}
		if(allOK!=0) return false;
		if(DBG) cout<<" Value of Start State = "<<setprecision(10)<<x[0]<<endl;
		return true;
	}
	bool findFirstValues(map<int,int> & M,map<int,int> &N,int &n,sparselib::Vector < double> &x,sparselib::Vector < double> &y){
		CCoorMatrix<double> ImP;
		ImP.new_dim(n, n,NNZ);
		cout<<"States in G: "<< n<<endl;
		sparselib::Vector < double> r(n);
		sparselib::Vector < double> c(n);
		sparselib::Vector < double> z(n);	
		sparselib::Vector < double> t1(n);
		sparselib::Vector < double> t2(n);
		r=0.;
		z=0.;
		c=0.;
		f(i,n){
			ImP.insert(i,i) = 1;
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				assert(i != N[e->to]);
				ImP.insert(i,N[e->to]) = -e->prob;
			}
			r[i] = states[M[i]]->actions[states[M[i]]->chosen]->cost;
		}
		c = ImP * r;
		c= -c;
		ImP.internal_order();
		if(!PROBINPUT) solveSparseEqn(ImP,c,z,n,gitr1,gitr2);
		else solveSparseEqn(ImP,c,z,n,3000,3000);
		y = ImP * z;
		f(i,y.size()) y[i] = -y[i];
		t1 = ImP *z;
		t2 = ImP * t1;
		x = t2 + r;
		
		
		//testSol(ImP,z,c,x,y,r);
		t1 = ImP * z;
		t2 = ImP * t1;
		t1 = ImP * t2;
		t2 = t1 -c;
		bool ret = true;
		double val = sparselib::normi(t2);
		//cout<<"VAL1 : "<<setprecision(10)<<val<<endl;
		if(val > EPSCHECK) ret= false;
		t1 = ImP*x;
		val = sparselib::normi(t1);
		//cout<<"VAL2 : "<<setprecision(10)<<val<<endl;
		if(val > EPSCHECK) ret= false;
		t1 = ImP *y;
		t2 = x + t1 -r;
		val = sparselib::normi(t2);
		//cout<<"VAL3 : "<<setprecision(10)<<val<<endl;
		if(val > EPSCHECK) ret= false;
		t1 = y + ImP*z;
		val = sparselib::normi(t1);
		//cout<<"VAL4 : "<<setprecision(10)<<val<<endl;
		if(val > EPSCHECK) ret= false;
		f(i,n){
			states[M[i]]->value = x[i];
			if((EPSCHECK+x[i]) <= 0) {
				//cout<<"x["<<i<<"] = "<<x[i]<<endl;
				ret= false;
			}
		}
		return ret;
		
	}   
	void setvalue(map<int,int> & M,map<int,int> &N,int n,double &alpha){
		int ctr = 3*n;
		CCoorMatrix<double> A;
		CCoorMatrix<double> B;
		A.new_dim(ctr, ctr,NNZ);
		B.new_dim(ctr, ctr,NNZ);
		//cout<<350<<endl;
		sparselib::Vector < double> b(ctr);
		sparselib::Vector < double> x(ctr);	
		b=0.;
		x=0.;
		f(i,n){
			B.insert(n+i,i) = -1;
			B.insert(2*n+i,n+i) = -1;
			B.insert(2*n+i,i) = 1;
		}
		f(i,n){
			f(j,states[M[i]]->actions[states[M[i]]->chosen]->edges.size()){
				edge * e =states[M[i]]->actions[states[M[i]]->chosen]->edges[j];
				//if( i == N[e->to]) { if(e->prob != 1.0) A.insert(i,N[e->to]) = 1- e->prob; }
				//else { if(e->prob != 0) A.insert(i,N[e->to]) = e->prob; }
				//if( n+i == N[e->to]){ if(e->prob != 1.0) A.insert(n+i,N[e->to]) = 1- e->prob;}
				//else { if(e->prob != 0) A.insert(n+i,n+N[e->to]) = e->prob;}
				//if( 2*n+i == N[e->to]) { if(e->prob != 1.0) A.insert(2*n+i,N[e->to]) = 1- e->prob;}
				//else { if(e->prob != 0) A.insert(2*n+i,2*n+N[e->to]) = e->prob; }
				A.insert(i,N[e->to]) = e->prob;
				A.insert(n+i,n+N[e->to]) = e->prob;
				A.insert(2*n+i,2*n+N[e->to]) = e->prob;
			}
			b[n+i] = states[M[i]]->actions[states[M[i]]->chosen]->cost;
			//if(b[n+i] != ((int)b[n+i])) cout<<b[n+i]<<" ";
		}
		//cout<<endl;
		A.internal_order();
		f(i,ctr){
			//assert(A.position(i,i) == NNZ);
			B.insert(i,i) = 1;
		}
		B.internal_order();
		//CCoorMatrix<double> pat;
		//patRAND(pat,ctr);
		//pat.new_dim(ctr,ctr,NNZ);
		//f(i,ctr) pat.insert(i,i) =1.0;
		//pat.internal_order();
		solveSparseLinEqSpl(B,A,b,x,ctr,alpha);
		if(DBG) cout<<"[";
		f(i,ctr){ 
			states[M[i]]->value = x[i];
			//assert(x[i]>0);
			//should work!
			if(DBG) { if(x[i] != 0) cout<<x[i]<<" "; }
		}
		if(DBG) cout<<"]"<<endl;
	}
	void setIncomingEdges(){
		f(i,states.size()){
			states[i]->incoming.clear();
		}
		f(i,states.size()){
			f(j, states[i]->actions[states[i]->chosen]->edges.size()){
				edge * e = states[i]->actions[states[i]->chosen]->edges[j];
				states[e->to]->incoming.push_back(e);
				//cout<<"("<<e->from<<","<<e->to<<") ";
				//if(e->to == *bads.begin()) cout<<"("<<e->from<<","<<e->to<<") ";
			}
		}
	
	}
	void checkRandomOutOfPreStarB(){
		bool canfind = false;
		f(i,states.size()){
			if(states[i]->badReachable){
				f(j, states[i]->actions[states[i]->chosen]->edges.size()){
					edge * e = states[i]->actions[states[i]->chosen]->edges[j];
					if(states[e->to]->badReachable == false){
					 cout<<"Can Find an outgoing edge :-/ "<<i<<" -> "<<e->to<<endl;
					 canfind = true;
					 break;
					}
				}
			}
			if(canfind) break;
		}
	}
	void checkGoodBadConflict(){
		int ct=0;
		f(i,states.size()){
			if(states[i]->goodReachable && states[i]->badReachable) ct++;
		}
		cout<<"Good Bad Conflicts = "<<ct<<endl;
	}
	bool checkGoodBadConflict(int i){
		if(states[i]->goodReachable && states[i]->badReachable) return true;//cout<<"Can Reach both Good and Bad States From State : "<<i<<endl;
		return false;
	}
	void checkStartReachableGoodBad(){
		f(i,states.size()){
			if(states[i]->startReachable && states[i]->goodReachable) {
				cout<<"Found a Start-Reachable Node leading to the good area! "<<i << " : "<<setprecision(20)<<states[i]->value<<endl;
				break;
			}
		}
	}
	void BFSlocks(vector<int> &dist,vector<int> &cst){
		set<int> prevlevel;
		prevlevel.insert(0);
		states[0]->DFSmark = true;
		while(1){
			set<int> thislevel;
			fi(i,prevlevel){
				states[(*i)]->DFSmark = true;
				fe(it,states[(*i)]->actions[states[(*i)]->chosen]->edges){
					int id = (*it)->to;
					if(dist[id] == -1){
						dist[id] = dist[(*i)] +1;
					}
					else if(dist[id] > dist[(*i)] +1){
						dist[id] = dist[(*i)] +1;
					}
					if(cst[id] == -1){
						cst[id] = cst[(*i)] +(*it)->wt;
					}
					else if(cst[id] < cst[(*i)] +(*it)->wt){
						cst[id] = cst[(*i)] +(*it)->wt;
					}
					
					if(!states[id]->DFSmark){
						thislevel.insert(id);
					}
				}
			}
			if(thislevel.empty()) break;
			//cout<<thislevel.size()<<endl;
			prevlevel = thislevel;
		}
		return;
	
	}
	void testCycles(){
		markSpecialStartReachable();
		vector<int> dist;
		vector<int> cst;
		f(i,states.size()){ 
			states[i]->DFSmark = false;
			dist.push_back(-1);
			cst.push_back(-1);
		}
		dist[0] =0;
		cst[0]=0;
		BFSlocks(dist,cst);
		cout<<"CYCLE SIZES: ";
		double sum = 0;
		int ctr = 0;
		f(id,states.size()){
			if(states[id]->mark) continue;
			bool lastone = false;
			if(dist[id] >= 0){
			fe(it,states[id]->actions[states[id]->chosen]->edges) if((*it)->to == 0) lastone = true;
			}
			if(lastone){
				//if(dist[id] == x) idx = id;
				assert(dist[id]>=0);
				cout<<id<<":("<<cst[id]<<"/"<<(dist[id]+1)<<"),";
				sum += (cst[id]*1.0)/(dist[id]+1);
				ctr++;
			}
		}
		
		cout<<".. AVG = "<<sum/ctr<<endl;
		markStartReachable();
	}
	
	void printStrategy(string fname){
		ofstream fout(fname.c_str());
		f(i,states.size()){
			if(!states[i]->startReachable) continue;
			fout<<states[i]->id<<" "<<states[i]->actions[states[i]->chosen]->sid<<" "<<endl;
		}
	}

	void printDotFile(string fname){//Assumes startReachable is already marked
		//setIncomingEdges();
		//setBadReachable();
		//markSpecialStartReachable();
		ofstream fout(fname.c_str());
		fout<<"digraph {"<<endl;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			if(states[i]->sv.type) fout<<i<<"[shape=box]"<<endl;
			else fout<<i<<"[shape=circle]"<<endl;

		}	
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			f(j, states[i]->actions.size()){
				f(k, states[i]->actions[j]->edges.size()){
					edge * e = states[i]->actions[j]->edges[k];
					string l = states[e->to]->sv.idstr();
					l = l.substr(1,l.length());
					fout<<e->from<<" -> "<<e->to;
					if(states[e->from]->sv.type)
						fout<<"[label=\""<<strchange(states[e->from]->sv,states[e->to]->sv)<<"("<<e->wt<<")"<<"\"];"<<endl;
					else
						fout<<"[label=\""<<states[e->to]->sv.reqstr()<<"("<<e->prob<<")"<<"\"];"<<endl;
				}
			}
		}
		fout<<"}"<<endl;
	}
	void printDotFileColor(string fname){//Assumes startReachable is already marked
		//setIncomingEdges();
		//setBadReachable();
		//markSpecialStartReachable();
		ofstream fout(fname.c_str());
		fout<<"digraph {"<<endl;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			if(states[i]->sv.type) fout<<i<<"[shape=box,label=\""<< states[i]->sv.idstr()<<"\"]"<<endl;
			else fout<<i<<"[shape=circle,label=\""<< states[i]->sv.idstr()<<"\"]"<<endl;

		}	
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			f(j, states[i]->actions.size()){
				f(k, states[i]->actions[j]->edges.size()){
					edge * e = states[i]->actions[j]->edges[k];
					string l = states[e->to]->sv.idstr();
					l = l.substr(1,l.length());
					fout<<e->from<<" -> "<<e->to;
					if(states[e->from]->sv.type){
						if(j == states[i]->chosen)
							fout<<"[color=red,label=\""<<strchange(states[e->from]->sv,states[e->to]->sv)<<"("<<e->wt<<")"<<"\"];"<<endl;
						else
							fout<<"[label=\""<<strchange(states[e->from]->sv,states[e->to]->sv)<<"("<<e->wt<<")"<<"\"];"<<endl;
					}
					else
						fout<<"[label=\""<<states[e->to]->sv.reqstr()<<"("<<e->prob<<")"<<"\"];"<<endl;
				}
			}
		}
		fout<<"}"<<endl;
	}
	int printDotFileStrat(string fname){//Assumes startReachable is already marked
		//setIncomingEdges();
		//setBadReachable();
		//markSpecialStartReachable();
		ofstream fout(fname.c_str());
		fout<<"digraph {"<<endl;
		int stcnt=0;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			if(states[i]->sv.type) fout<<i<<"[shape=box]"<<endl;
			else fout<<i<<"[shape=circle]"<<endl;
			stcnt++;
		}	
		cout<<"States in M: "<<stcnt<<endl;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			f(j, states[i]->actions.size()){
				f(k, states[i]->actions[j]->edges.size()){
					edge * e = states[i]->actions[j]->edges[k];
					fout<<e->from<<" -> "<<e->to<<"[label=\""<<e->label<<"\"];"<<endl;
				
				}
			}
		}
		fout<<"}"<<endl;
		return stcnt;
	}
	void printDotFileStratOnly(string fname){//Assumes startReachable is already marked
		//setIncomingEdges();
		//setBadReachable();
		markSpecialStartReachable();
		ofstream fout(fname.c_str());
		fout<<"digraph {"<<endl;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;
			if(states[i]->sv.type) fout<<i<<"[shape=box]"<<endl;
			else fout<<i<<"[shape=circle]"<<endl;

		}	
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int i =it->first;
			if(!states[i]->startReachable) continue;

			int j = states[i]->chosen;	
			f(k, states[i]->actions[j]->edges.size()){
					edge * e = states[i]->actions[j]->edges[k];
					if(states[e->from]->sv.type)
						fout<<e->from<<" -> "<<e->to<<"[label=\""<<strchange(states[e->from]->sv,states[e->to]->sv)<<"("<<e->wt<<")"<<"\"];"<<endl;
					else
						fout<<e->from<<" -> "<<e->to<<"[label=\""<<states[e->to]->sv.reqstr()<<"("<<e->prob<<")"<<"\"];"<<endl;
			}
		}
		fout<<"}"<<endl;
	}

	void setRandomPolicy(map<int,int> & M,map<int,int> &N,int &n,sparselib::Vector < double> &x,sparselib::Vector < double> &y,int mode){
		int improved = 0;
		if(mode == 0){
		
		f(i,states.size()){
			if(!states[i]->startReachable) continue;
			int bestA1=states[i]->chosen,bestA2=states[i]->chosen;
			double bestS1 = x[N[i]],bestS2=x[N[i]] + y[N[i]];
			f(j,states[i]->actions.size()){
				if(j == states[i]->chosen) continue;
				double s1 = 0,s2=states[i]->actions[j]->cost;
				f(k,states[i]->actions[j]->edges.size()){
					edge * e = states[i]->actions[j]->edges[k];
					s1 += e->prob * x[N[e->to]];
					s2 += e->prob * y[N[e->to]];
				}
				if(bestS1 < s1-EPSCHECK){
					bestA1 = j;
					bestS1 = s1;
				}
				else if ( abs(bestS1 - s1) <= EPSCHECK){
					if(s2-EPSCHECK > bestS2){
						bestA2 = j;
						bestS2 = s2;
					}
				}
			}
			if(bestA1 != states[i]->chosen){
				improved++;
				//if(DBG) cout<<i<<" : (Act1) "<< states[i]->chosen << " -> " <<bestA1 << " -- "<< x[N[i]] << " -> "<< bestS1 << endl;
				states[i]->chosen = bestA1;
			}
			/*else if (bestA2 != states[i]->chosen){
				improved++;
				 cout<<i<<" : (Act2) "<< states[i]->chosen << " -> " <<bestA2 << " -- "<< x[N[i]]+y[N[i]] << " -> "<< bestS2 << endl;
				states[i]->chosen = bestA2;
			}*/
		}
		}
		
		if(improved == 0){
			double tots = 0;
			f(i,states.size()){
				//if(states[i]->actions.size() ==1) continue;
				//int x = rand()%(states[i]->actions.size() -1);
				//if(x == states[i]->chosen) states[i]->chosen = states[i]->actions.size() -1;
				//else 
				states[i]->chosen = rand()%(states[i]->actions.size());

				Assert((states[i]->actions.size()!= 0),"Actions cannot be empty");
				
				tots = tots + log2(states[i]->actions.size());
				Assert(states[i]->actions.size() == 1 || states[i]->actions.size() <= tNF, "Action size constraint violated!");	
			}
			//cout<<"No fake improvement possible..."<<endl;
			cout<<"Total Number of Strategies Possible: 2^" <<tots<<endl;
		}
		else cout<<"fake improvement done..."<<endl;
		//testCycles();
		return;
	}
	
	void setMaps(map<int,int> & M,map<int,int> &N,int &n){
		int ctr=0,rctr=0;
		/*for(set<int>::iterator it = badReach.begin();it != badReach.end();it++){
			if(states[*it]->bad) continue;
			M[ctr] = *it; //sorted
			N[*it] = ctr;
			ctr++;
		}*/
		f(i,states.size()){
			if (states[i]->startReachable){
				M[ctr] =i;
				N[i] = ctr;
				ctr++;
			}
			else rctr++;
		}
		n= ctr;
		//cout<<"RCTR : "<<rctr<<endl;
		//M[0] -> M[ctr-1] are the start reachable states
	}
	void init(map<int,int> & M,map<int,int> &N,int &n,sparselib::Vector < double> &x,sparselib::Vector < double> &y){
		//int ctr = 0;
		int mode = 0;
		while(1){
			setRandomPolicy(M,N,n,x,y,mode);
			//if(findFirstValues(M,N,n,x,y)){
			if(findFirstValues(M,N,n,x,y)){
				break;
			}
			mode  = (mode +1)%2;
			cout<<"RESTARTING..."<<endl;
			//else cout<<ctr<<",";
			//ctr++;
		}
		//cout<<ctr<<endl;
	}
	double ValueStartState(){ //to the states marked bad
		
		map<int,int> M,N;
		int ctr;
		setMaps(M,N,ctr);
		//cout<<"CTR : "<<ctr<<endl;
		sparselib::Vector < double> x(ctr);
		sparselib::Vector < double> y(ctr);
		//init(M,N,ctr,x,y);	
		setRandomPolicy(M,N,ctr,x,y,1);
		//x and y have the desired value and deviation vectors
		//setvalue(M,N,ctr,alpha);
		//sparseMatMethod(ctr,M,N);
		//return states[0]->value;
		//init();	
		int itrs = 0;
		findValues(M,N,ctr,x,y);
		while(1){
			itrs++;
			//testCycles();
			int improved = 0;
			f(i,states.size()){
				
				if(!states[i]->startReachable) continue;
				int bestA1=states[i]->chosen,bestA2=states[i]->chosen;
				double bestS1 = x[N[i]],bestS2=x[N[i]] + y[N[i]];
				f(j,states[i]->actions.size()){
					if(j == states[i]->chosen) continue;
					double s1 = 0,s2=states[i]->actions[j]->cost;
					f(k,states[i]->actions[j]->edges.size()){
						edge * e = states[i]->actions[j]->edges[k];
						s1 += e->prob * x[N[e->to]];
						s2 += e->prob * y[N[e->to]];
					}
					if(bestS1 < s1-EPSCHECK){
						bestA1 = j;
						bestS1 = s1;
					}
					else if ( abs(bestS1 - s1) <= EPSCHECK){
						if(s2-EPSCHECK > bestS2){
							bestA2 = j;
							bestS2 = s2;
						}
					}
				}
				if(bestA1 != states[i]->chosen){
					improved++;
					states[i]->chosen = bestA1;
				}
				else if (bestA2 != states[i]->chosen){
					improved++;
					states[i]->chosen = bestA2;
				}
			}
			cout<<"Value for Strategy at Iteration #"<< itrs <<": "<<setprecision(10)<<states[0]->value<<endl;
			
			if(improved==0){
				cout<< "Strategy Improvement Iterations: "<<itrs<<endl;
			 	return states[0]->value;
			}
			//cout<<" IMPROVED = "<<improved<<endl;
			//pl();
			if (!findValues(M,N,ctr,x,y)){
				pl();
				cout<<"RESTARTING: "<<__LINE__<<endl;
				init(M,N,ctr,x,y);
			}
			//pl();
		}
		
	}
	void sanityCheck(){
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			int i = it->first;
			if(states[i]->sv.type){//player state
				//Assert(states[i]->actions.size() <= tNF, "Enough actions for Player States");
			}
			else{
				Assert(states[i]->actions.size() == 1, "One action for Env States");
				//int p2n0s = (int)pow(2.0,NF-states[i]->sv.num1s());
				//Assert(states[i]->actions[0]->edges.size() == p2n0s, "2^NF Edges for Env State Action");
			}
			f(j,states[i]->actions.size()){
				double tprob=0.0;
				f(k,states[i]->actions[j]->edges.size()){
					int to = states[i]->actions[j]->edges[k]->to;
					Assert( find(states[to]->incoming.begin(),states[to]->incoming.end(),states[i]->actions[j]->edges[k]) != states[to]->incoming.end(), "Edge not found in incoming field of the to state");
					tprob += states[i]->actions[j]->edges[k]->prob;
				}
				Assert(abs(tprob - 1.0) < EPSCHECK, toSTR(__LINE__) + " Total prob not 1.0");
			}
		}
	}
	void sanityCheckStrat(){
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			int i = it->first;
			if(states[i]->sv.type){//player state
				Assert(states[i]->actions.size() == 1, "One action for Player States After Strat Fix");
				Assert(states[i]->actions[0]->edges.size() == 1, "One edge per action for Player States After Strat Fix");
			}
			else{
				Assert(states[i]->actions.size() == 1, "One action for Env States");
				//int p2n0s = (int)pow(2.0,NF-states[i]->sv.num1s());
				//Assert(states[i]->actions[0]->edges.size() == p2n0s, "2^p2n0s Edges for Env State Action");
			}
			f(j,states[i]->actions.size()){
				double tprob=0.0;
				f(k,states[i]->actions[j]->edges.size()){
					int to = states[i]->actions[j]->edges[k]->to;
					Assert( find(states[to]->incoming.begin(),states[to]->incoming.end(),states[i]->actions[j]->edges[k]) != states[to]->incoming.end(), "Edge not found in incoming field of the to state");
					tprob += states[i]->actions[j]->edges[k]->prob;
				}
				Assert(abs(tprob - 1.0) < EPSCHECK, toSTR(__LINE__) + " Total prob not 1.0");
			}
		}
	}
	action* findAction(edge *e){
		for(vector<action*>::iterator it = states[e->from]->actions.begin();it != states[e->from]->actions.end(); it++){
			if(find((*it)->edges.begin(),(*it)->edges.end(),e) != (*it)->edges.end())
				return (*it);
		}
		Assert(false,"No action for this edge!");
		return NULL;
	}
	//removeEdge
	void removeEdge(edge *e){
		if(doneEdges.find(e) != doneEdges.end()) return;
		doneEdges[e] = true;
		action* a = findAction(e);
		//assume edge belongs to the action
		Assert(find(a->edges.begin(),a->edges.end(),e) != a->edges.end(), "Assumption Failure for Edge in removeEdge.1");
		a->edges.erase(find(a->edges.begin(),a->edges.end(),e));

		Assert(find(states[e->to]->incoming.begin(),states[e->to]->incoming.end(),e) != states[e->to]->incoming.end(), "Assumption Failure for Edge in removeEdge.2");
		states[e->to]->incoming.erase(find(states[e->to]->incoming.begin(),states[e->to]->incoming.end(),e));
		
		Assert(find(states[e->from]->actions.begin(),states[e->from]->actions.end(),a) != states[e->from]->actions.end(), "Assumption Failure for Action in removeEdge");
		//Not removing action if there are more edges left in the action
		if(a->edges.size() == 0){
			states[e->from]->actions.erase(find(states[e->from]->actions.begin(),states[e->from]->actions.end(),a));
			delete a;
		}
		int frm = e->from;
		if(states[frm]->actions.size() == 0){//remove state 
			Assert(states[frm]->incoming.size() == 0, "Only such cases have no actions!");
			/*f(k, states[frm]->incoming.size()){
				edge* e2 = states[frm]->incoming[k];
				removeEdge(e2);
			}*/
			//removing all incoming edges should have removed the state
			if(states.find(frm) != states.end()){
				//state still available
				Assert(states[frm]->incoming.size() == 0, "Incoming Empty");
				Assert(states[frm]->actions.size() == 0, "Actions Empty");
				state * s = states[frm];
				states.erase(states.find(frm));
				delete s;
			}
			//Assert(states.find(e->from) == states.end(),"e->from state Id still available!");
			//states.erase(states.find(e->from));
		}
		delete e;
		
	}
	//removeAction -- Taken care of by removing edges, automatically removed when removing last edge for an action
	//removeState -- Taken care of by removing edges, automatically removed when removing last edge from a state
	//cleanup 
	void cleanup(){
		markStartReachable();
		//states are marked startReachable
		vector<edge*> tbre;
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			int aid = cs->chosen;
			if(!cs->startReachable) aid = -1;
			f(j,cs->actions.size()){
				if(j==aid) continue;
				f(k,cs->actions[j]->edges.size()){
					tbre.push_back(cs->actions[j]->edges[k]);
				}
			}
			
			cs->chosen = 0;

		}
		f(i,tbre.size()){
			removeEdge(tbre[i]);
		}
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			Assert(it->second->startReachable == true, "Unreachable state after cleanup!");
		}
		cleanup2();
		//another cleanup needed for minimization
	}
	void cleanup2(){
		sanityCheckStrat();
		vector<edge*> tbre;
					
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			state * cs = it->second;
			Assert(cs->id == it->first,"Id Mismatch");
			if(cs->sv.type){//Player State
					//there is a few incoming edges and one outgoing edge from this state
					//we will remove both but add one edge from the previous state to next state
					//for each incoming edge
					tbre.push_back(cs->actions[cs->chosen]->edges[0]);
					int to = cs->actions[cs->chosen]->edges[0]->to;
					Assert( find(states[to]->incoming.begin(),states[to]->incoming.end(),cs->actions[cs->chosen]->edges[0]) != states[to]->incoming.end(), "Edge not found in incoming field of the to state");
					Assert(states[to]->sv.type == false, "to state not an Env State!");
					f(j,cs->incoming.size()){
						int frm = cs->incoming[j]->from;
						//tbre.push_back(cs->incoming[j]);
						//tbra.push_back(findAction(cs->incoming[j]));
						edge* e = cs->incoming[j];
						e->label = states[e->to]->sv.reqstr() + "/" + strchange(cs->sv,states[to]->sv);
						e->to = to;
						
						Assert(states[frm]->sv.type == false, "frm state not an Env State!");
						

						states[to]->incoming.push_back(e);//remember to delete the P->E edge once out of the loop
					}
					cs->incoming.clear();
					//states[to]->incoming.erase(find(states[to]->incoming.begin(),states[to]->incoming.end(),cs->actions[cs->chosen]->edges[0]));
				
				}
		}
		f(i,tbre.size()){
			removeEdge(tbre[i]);
		}
		for(map<int,state*>::iterator it = states.begin();it != states.end();it++){
			Assert(it->second->startReachable == true, "Unreachable state after cleanup!");
		}
	}
	//keepStrategy
	int keepStrategy(string fname){
		//first remove all edges except the strategy ones
		
		//markStartReachable
		cleanup();
		//cout<<"States in M: " <<states.size()<<endl;
		markStartReachable();
		return printDotFileStrat(fname);
	}
};

