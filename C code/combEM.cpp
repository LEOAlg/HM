#include "runEM.h"
using namespace std;
using namespace arma;

// g++ -std=c++0x -w -fopenmp combEM.cpp -O1 -larmadillo

// time ./a.out 1 >result.txt -R

#ifndef typedef_emParam_t
#define typedef_emParam_t
typedef struct {
    vec mu;
    mat C;
    mat T;
    double sigma;
} emParam_t;
#endif //typedef_emParam_t

#ifndef typedef_emValue_t
#define typedef_emValue_t
typedef struct {
    vec w;
    emParam_t em_t;
} emReturn_t;
#endif //typedef_emValue_t

#ifndef typedef_sample_t
#define typedef_sample_t
typedef struct {
    unsigned int id; // configuration id
    double performance;
    double power;
} sample_t;
#endif //typedef_sample_t

#ifndef typedef_application_t
#define typedef_application_t
typedef struct {    
    list<sample_t> s;
    int n; // number of samples known
    string name;
} application_t;
#endif //typedef_application_t

void init_EMParam(emParam_t *Old, int n);
void EM(emParam_t *Old, mat *W, mat *y_em, int i, int ep, int iteration_limit, emReturn_t *Appl);
double residualError(vec A, vec B);
long Time_Difference(struct timeval end, struct timeval start);
void loadData(int appId, mat *Dpower, mat *Dperformance, mat *truePower, mat *truePerformance,  mat *W);

void init_EMParam(emParam_t *Old, int n){
	Old->mu = zeros<vec>(n);
	Old->C = eye<mat>(n,n);
	Old->T = zeros<mat>(n,n);
	Old->sigma = 1;
}
void EM(emParam_t *Old, mat *W, mat *y_em, int i, int ep, int iteration_limit, emReturn_t *Appl){
	
	int n = y_em->n_rows;
	int m = y_em->n_cols;	
	int numSamples = sum(sum(*W));	
	double pi = 1;
	double tau = 1;
	double error = INFINITY;
	mat I = diagmat(*W);	
	double sigma = Old->sigma;
	vec mu = Old->mu;
	mat C = Old->C;
	mat wl = zeros<mat>(n,m);	
	
	for(int iterator = 1; iterator <= iteration_limit ; iterator++){		
		//struct timeval t1, t2,t3; long seconds; gettimeofday(&t1, NULL); ///////////////////////
		
		mat Cinv = inv(C);	
		mat Cl0 = inv(I/sigma + Cinv);
		mat Cll= inv(eye<mat>(n,n)/sigma + Cinv);			
		
		wl = (Cll/sigma)*(*y_em) + repmat(Cinv*mu,1,m);		
		wl.col(0) = Cl0*((y_em->col(0))/sigma +Cinv*mu);
		
		//gettimeofday(&t2, NULL); seconds = Time_Difference(t2,t1);
		//cout<<"Part1 milliseconds: "<<seconds<<endl;//////////////////
		
		double normSum = pow(norm(I*(y_em->col(0) -wl.col(0)),2),2) + trace(I* Cl0); 
		vec wlSum = sum(wl, 1);
		mat ClSum = eye<mat>(n,n) + Cl0 +(m-1)*Cll;		
		mat wlSumCov = eye<mat>(n,n) + (wl.col(0)-mu)*trans(wl.col(0)-mu);
		for(int l = 1; l < m; l++)		{
			wlSumCov = wlSumCov + (wl.col(l)-mu)*trans(wl.col(l)-mu);    
			normSum = normSum + pow(norm(y_em->col(l) - wl.col(l),2),2)+ trace( Cll);            
        }
		mu = ((1.0/(double)(pi + m))*wlSum);        
        C = ((1.0/(double)(tau + m))*(pi*(mu*trans(mu)) + tau*eye<mat>(n,n)+ ClSum + wlSumCov));        
        sigma = (1.0/(double)((m-1)*n + numSamples))*normSum;
        error = norm(Old->mu - mu,"fro")+ norm(Old->C - C,"fro") + abs(Old->sigma - sigma);        
        Old->mu = mu;    Old->C = C;    Old->sigma = sigma;  
        
        //gettimeofday(&t3, NULL);  seconds = Time_Difference(t3,t2);////////////////////
		//cout<<"Part2 milliseconds: "<<seconds<<error<<endl;//////////////////
        
        //cout<<"Itr---------------------------------"<<iterator<<endl;
        //cout<<"C ----------------------------------\n"<<C<<endl;            
		//cout<<"mu ----------------------------------\n"<<mu<<endl;       
		//cout<<"sigma: "<<sigma<<endl; 
		//cout<<"wl ----------------------------------\n"<<wl<<endl;   
		   
	}	
	/*
	cout<<"y_em ------------------------------\n"<<*y_em<<endl;
	cout<<"W ------------------------------\n"<<*W<<endl;
	cout<<"C ----------------------------------\n"<<C<<endl;            
    cout<<"mu ----------------------------------\n"<<mu<<endl;       
    cout<<"sigma: "<<sigma<<endl; 
    cout<<"wl ----------------------------------\n"<<wl<<endl;     
           */
    
	Appl->w = wl.col(0);
	Appl->em_t = *Old;	
}		
long Time_Difference(struct timeval end, struct timeval start){
	long mtime, seconds, useconds; 
	seconds  = end.tv_sec  - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;

    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;

    return mtime;
}
double residualError(vec A, vec B){
	double n = A.n_elem;
	double rss = pow(norm(A-B, 2),2);     
    double tss = pow(norm(A-mean(A),2),2);
    double residualsquare = 1-rss/tss;
    return 1-(1-residualsquare)*(n-1)/(n-2);   
}

void loadData(int appId, mat *Dpower, mat *Dperf, mat *truePower, mat *truePerformance, mat *W){
	int n, numSamples;
	mat SupPower, SupPerf;
	
	numSamples = 20;
	Dpower->load("data/AllPower.txt");
	Dperf->load("data/AllPerformance.txt");
	n = Dpower->n_rows;	// # CONFIGURATIONS	
	*W = zeros<vec>(n);	
	for(int i = 0; i < n; i+=(n/numSamples + 1))
		(*W)(i) = 1;
		
	*truePower = Dpower->col(appId);
	*truePerformance = Dperf->col(appId);
	Dpower->col(appId) = (*W)%Dpower->col(appId);
	Dperf->col(appId) = (*W)%Dperf->col(appId);  
	
	Dpower->swap_cols(0,appId);    
	Dperf->swap_cols(0,appId);    
}

int main(int argc, char *argv[]){		
	int appId, n, m;
	double accuracy_power, accuracy_perf;
	string str(argv[1]);	
	string strPow, strPerf;
	mat Dpower, Dperformance, y_em, W, truePower, truePerformance;
	vec pow, perf;
	application_t targetApp;	
	emParam_t Old_power, Old_perf;
    emReturn_t Appl_power, Appl_perf;   
    
	// LOAD DATA AND MISSING VALUES, BS is the name of target application
	appId = stoi(str)-1;
	loadData(appId, &Dpower, &Dperformance,  &truePower, &truePerformance, &W);	
	
    n = Dpower.n_rows;	// # CONFIGURATIONS
    m = Dpower.n_cols;	// # APPLICATIONS	
	
	ifstream infname, inepsilon, initrlimit;
	
    double epsilon[m][2], iteration_limit[m][2];
    string filename[m];

	infname.open("data/filename.txt");
	inepsilon.open("data/epsilon.txt");
	initrlimit.open("data/iteration_limit.txt");
    for (int i = 0; i < m; i++) {
		string s1, s2, s3, s4, s5;
		infname >> filename[i];
		inepsilon >> s1>> s2>> s3;
		epsilon[i][0] = stoi(s1);
		epsilon[i][1] = stoi(s2);	
		initrlimit >> s4 >> s5;
		iteration_limit[i][0] = stoi(s4);
		iteration_limit[i][1] = stoi(s5);
	}

	// START EM ALGORITHM FOR POWER	AND PERFORMANCE  
	#pragma omp parallel sections		
	{
		#pragma omp section			
		{
			 init_EMParam(&Old_perf, n);      
			 EM (&Old_perf, &W , &Dperformance, 0, epsilon[appId][0], iteration_limit[appId][0] , &Appl_perf);     
		} 
		#pragma omp section		
		{    
			 init_EMParam(&Old_power, n);      
			 EM (&Old_power, &W , &Dpower, 0, epsilon[appId][1], iteration_limit[appId][1], &Appl_power);  
		}
	}
    
    // CALCULATE ERROR AND GET OUTPUT TO THE FILES     
    accuracy_perf = residualError( Appl_perf.w, truePerformance.col(0));        
    accuracy_power = residualError( Appl_power.w, truePower.col(0));  
	cout<< filename[appId] <<"\t"<< accuracy_perf <<"\t" << accuracy_power<<endl;	
    return 0;   
}
