#define NF 2
double pall = 0.1; //MDPTYPE 1
int rewards[NF] = {2,1};
int steps[NF] = {1,2};

class sysvars{
public:
  int prod;
  int step;
  bool repair;  
  int reward;
  sysvars(int p,int s, bool r, int rew){
    prod = p; step= s; repair = r; reward = rew;
  }
};

#define MDPTYPE 2
double pdiff[2] = {0.05,0.12}; //MDPTYPE 2
//#define NK 3

