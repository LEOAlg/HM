#include "runEM.h"
using namespace std;
using namespace arma;

// g++ -w -fopenmp runEM.cpp -O1 -larmadillo

// time ./runEM>result.txt -R

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
void EM(emParam_t *Old, mat *W, mat *y_em, int i, emReturn_t *Appl);
void mat2app(application_t *targetApp, mat *Dpower, mat *Dperformance, mat *W);
void app2mat(application_t *targetApp, mat *Dpower, mat *Dperformance, mat *W);
void exampleBS(application_t *targetApp );
void example2BS(application_t *targetApp);
double residualError(vec A, vec B);
void exampleMagic(mat *Dpower, mat *W);
long Time_Difference(struct timeval end, struct timeval start);
void init_EMParam(emParam_t *Old, int n){
	Old->mu = zeros<vec>(n);
	Old->C = eye<mat>(n,n);
	Old->T = zeros<mat>(n,n);
	Old->sigma = 1;
}
void EM(emParam_t *Old, mat *W, mat *y_em, int i, emReturn_t *Appl){
	
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
	for(int iterator = 1; iterator <= ITERATION_LIMIT && error > EPSILON; iterator++){		
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
        /*  
        cout<<"Itr---------------------------------"<<iterator<<endl;
        cout<<"C ----------------------------------\n"<<C<<endl;            
		cout<<"mu ----------------------------------\n"<<mu<<endl;       
		cout<<"sigma: "<<sigma<<endl; 
		cout<<"wl ----------------------------------\n"<<wl<<endl;   
		*/   
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
void mat2app(application_t *targetApp, mat *Dpower, mat *Dperformance){
	int n = Dpower->n_rows;	// # CONFIGURATIONS    
	for(int i = 0; i < n; i++){
		sample_t temp;
		temp.id = i;
		temp.power = (*Dpower)(i,0);
		temp.performance = (*Dperformance)(i,0);
		(targetApp->s).push_back(temp);
	}	
	targetApp->n = n;
}
void app2mat(application_t *targetApp, mat *Dpower, mat *Dperformance, mat *W){
	// Load the supporting data for other applications
	mat SupPower, SupPerf;
	
	SupPower.load("data/SupportPower.tsv");	
	SupPerf.load("data/SupportPerf.tsv");		
	int n = SupPower.n_rows;	// # CONFIGURATIONS	
	mat z = zeros<mat>(n,1);		
	*Dpower = join_rows(z, SupPower);
	*Dperformance = join_rows(z, SupPerf);
	
	// Load the target data into the matrix
	// The 0th column of the matrix our target application
	*W = zeros<vec>(n);	
	list<sample_t> s = targetApp->s;
	for(list<sample_t>::iterator it = s.begin(); it != s.end(); ++it){
		sample_t temp = *it;
		(*Dpower)(temp.id, 0) = temp.power;
		(*Dperformance)(temp.id, 0) = temp.performance;
		(*W)(temp.id) = 1;
	}
}
void exampleBS(application_t *targetApp){
	// Loading example BS with 10 samples. This function is only for testing.
	// We write into targetApp.
	int id[] = {62,66,339,639,779,783,803,943,945,1009};	
	int n = sizeof(id)/sizeof(id[0]);
	vec a = zeros<vec>(n);
	vec b = zeros<vec>(n);

	a[0] = 120.472;  a[5] = 174.085; 
	a[1] = 99.605;   a[6] = 135.753; 
	a[2] = 162.0650;  a[7] = 181.127; 
	a[3] = 311.311;  a[8] = 195.00; 
	a[4] = 159.979;  a[9] = 199.756; 
	
	b[0] = 7.48; 	b[5]  = 51.91; 
	b[1] = 4.6;  	b[6]  = 34.35; 
	b[2] = 32.48;	b[7]  = 58.25; 
	b[3] = 74.98;	b[8]  = 63.31; 
	b[4] = 46.72;	b[9]  = 66.73; 
	
	list<sample_t> s; 	
	for(int i = 0; i < n; i++){
		sample_t temp;
		temp.id = id[i]-1;
		temp.power = a[i];
		temp.performance = b[i];
		s.push_back(temp);
	}
	targetApp->s = s;
	targetApp->n = n;
	targetApp->name = "BS";		
}
void example2BS(application_t *targetApp){
	// Loading example BS with 10 samples. This function is only for testing.
	// We write into targetApp.
	int id[] = {62,66,339,639,779,783,803,943,945,1009,10};	
	int n = sizeof(id)/sizeof(id[0]);
	vec a = zeros<vec>(n);
	vec b = zeros<vec>(n);

	a[0] = 120.472;  a[5] = 174.085; a[10] = 99.378; 
	a[1] = 99.605;   a[6] = 135.753; 
	a[2] = 162.0650;  a[7] = 181.127; 
	a[3] = 311.311;  a[8] = 195.00; 
	a[4] = 159.979;  a[9] = 199.756; 
	
	b[0] = 7.48; 	b[5]  = 51.91; b[10] = 2.153; 
	b[1] = 4.6;  	b[6]  = 34.35; 
	b[2] = 32.48;	b[7]  = 58.25; 
	b[3] = 74.98;	b[8]  = 63.31; 
	b[4] = 46.72;	b[9]  = 66.73; 
	
	list<sample_t> s; 	
	for(int i = 0; i < n; i++){
		sample_t temp;
		temp.id = id[i]-1;
		temp.power = a[i];
		temp.performance = b[i];
		s.push_back(temp);
	}
	targetApp->s = s;
	targetApp->n = n;
	targetApp->name = "BS";		
}
void exampleMagic(mat *Dpower, mat *W){
	Dpower->load("data/M.txt");		
	W->load("data/Wmagic.txt");
	Dpower->col(0) = (W->col(0))%(Dpower->col(0));
	if(Dpower->n_rows != W->n_rows){
		cout<<"Error in data"<<endl;
		exit(EXIT_FAILURE);
	}
}

int main(){	
	double accuracy_power, accuracy_perf;
	string strPow, strPerf;
	mat Dpower, Dperformance, y_em, W, truePower, truePerformance;
	vec pow, perf;
	application_t targetApp;	
	emParam_t Old_power, Old_perf;
    emReturn_t Appl_power, Appl_perf;   
    
	// LOAD DATA AND MISSING VALUES, BS is the name of target application
	exampleBS( &targetApp ); // load BS into targetapp	
	app2mat( &targetApp, &Dpower, &Dperformance, &W );// load BS and supporting data into matrices.	
	
	//exampleMagic( &Dpower, &W );	
	
    int n = Dpower.n_rows;	// # CONFIGURATIONS
    int m = Dpower.n_cols;	// # APPLICATIONS	
	
	// START EM ALGORITHM FOR POWER	AND PERFORMANCE  
	#pragma omp parallel sections		
	{
		#pragma omp section			
		{
			cout<<"      ***************     POWER    ****************      \n"<<endl;  init_EMParam(&Old_power, n);      
			EM (&Old_power, &W , &Dpower, 0, &Appl_power);   
		} 
		#pragma omp section		
		{    
			cout<<"      ***************  PERFORMANCE ****************      \n"<<endl; init_EMParam(&Old_perf, n);      
			EM (&Old_perf, &W , &Dperformance, 0, &Appl_perf);       
		}
	}
    
    // INCREMENT A POINT **************************************************************************
    /*
    example2BS( &targetApp ); // load BS into targetapp	
	app2mat( &targetApp, &Dpower, &Dperformance, &W );// load BS and supporting data into matrices.	
	cout<<"      ***************     POWER    ****************      \n"<<endl;   	
    init_EMParam(&Old_power, n);      
    EM (&Old_power, &W , &Dpower, 0, &Appl_power);       
    cout<<"      ***************  PERFORMANCE ****************      \n"<<endl;       
    init_EMParam(&Old_perf, n);      
    EM (&Old_perf, &W , &Dperformance, 0, &Appl_perf); 
    */
    //*********************************************************************************************
	
	
	// LOAD TRUE DATA TO CALCULATE TRUE ERROR
    truePower.load("data/blackscholes_Pow.txt");	
	truePerformance.load("data/blackscholes_Perf.txt");		
	// The estimated values for configurations are saved into targetApp and targetApp can also be taken as output from this program.
    mat2app(&targetApp, &Dpower, &Dperformance);    
    
    // CALCULATE ERROR AND GET OUTPUT TO THE FILES
    accuracy_power = residualError( Appl_power.w, truePower.col(0));   
    accuracy_perf = residualError( Appl_perf.w, truePerformance.col(0));    
    pow = Appl_power.w;
    perf = Appl_perf.w;
    strPow= targetApp.name+ "_"+"Pow.txt";
    strPerf= targetApp.name+ "_"+"Perf.txt";
    pow.save(strPow, raw_ascii);
    perf.save(strPerf, raw_ascii);
    cout<<"------------------RESULT EM ALGORTIHM-----------------------------"<<endl;
    cout<<"Number of Config="<<n<<": Number of App="<<m<<endl;
    cout<<"Number of Iterations="<<ITERATION_LIMIT<<endl;  
	cout<<"Accuracy Power: "<<accuracy_power<<" Perf: "<<accuracy_perf<<endl;
	
    return 0;   
}
