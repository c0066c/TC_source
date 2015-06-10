/*
 * DynamicAlgo.cpp
 *
 *  Created on: Sep 12, 2012
 *      Author: anas
 *      Co-Author: khchen
 */
#define DYNA 1
#include "DynamicAlgo.h"



// Constructor
DynamicAlgo::DynamicAlgo(float alpha,double alphaScale,double RScale,int D,int sigma,int delta) {
	this->alpha=alpha;
	this->D=D;
	this->sigma=sigma;
	this->delta=delta;
	this->alphaScale=alphaScale;
	this->RScale=RScale;
	tableScale=(RScale/50000);


	// read them from data file
	n=5; // <<--
	K=24; // <<--
	DEBUG = 1;

	Dmax=D;
	//Dmax = ceil(D/delta)*delta;
}

// Destructor
DynamicAlgo::~DynamicAlgo() {
	// TODO Auto-generated destructor stub
}


struct order_th{
        double et, R, RT;
        int ind;
}*th_data;
//khchen convex_hull prototype
struct Node{
        double x, y;
        int ind;
}*Point, *CH;

struct Snode{
        int func, ver, skip;
        double R, AVG;
        double WI, MU;
        double addet;
        double density;
}*greedy, *cur;

//for sort et
bool cmp_th_d(order_th a, order_th b) {
    bool result = false;
#if 0
    if(a.R/a.et < b.R/b.et){
        result =true;
    }
#else
    if(a.RT/a.et < b.RT/b.et){
        result =true;
    }
#endif
    return result;
}

int Monotone_Chain(int n);
double cross(Node o, Node a, Node b)   {return (a.x-o.x)*(b.y-o.y)-(b.x-o.x)*(a.y-o.y);}
bool cmp(Node a, Node b) {return a.x < b.x || (a.x == b.x && a.y < b.y);} //first cmp x, next cmp y

int convex(int fuc, int l, double **WI, double **MU, int **sub, int hp) // construct the subitems
{

		int a=0, k=l-1;
        Point = new Node[l];
        CH = new Node[l*2];
//        cout << "Function = " << fuc << " version = "<< n << endl;
        for(int z = 0; z < l; z++){
        	if(z == hp) continue; //dont care the hp version
            //input x = mu execution time, y = effected RT profit
        	if(WI[fuc][z] <= 0){
        		k--;
        		continue;
        	}
        	//if(MU[fuc][z] <= 0) continue;
            Point[a].x = (MU[fuc][z]);
            Point[a].y = WI[fuc][z];
            Point[a].ind = z;
//            cout <<"MU "<< Point[a].x << " WI " << Point[a].y << endl;
//            cout << Point[a].ind << endl;
            a++;
	    }
        int m = Monotone_Chain(k); //kick out hp version
        /*Monotone Chain over, get the whole convex*/
#if 0
        for(int i = 0; i<m+1; i++){
            cout << "Convex hull node ind =" <<CH[i].ind<<endl;
            cout << CH[i].x << " " << CH[i].y << endl;
        }
#endif        

        int r = 1; //r = 0 used to record function versions number
        if(m < 0) return 0;
        if(CH[0].x != 0 && CH[0].y > 0) //prevent the version mu = 0
        sub[fuc][r++]=CH[0].ind;
//        cout << "The first node index = "<< CH[0].ind<< endl;
//        cout<< CH[0].x<< " " << CH[0].y <<endl;
        if(l > 2){
            // after find the convex hull, trace the node
            double cur_x = CH[0].x; //mu
            double cur_y = CH[0].y; //wi
            for(int i = 1; i<m; i++){ //ignore twice start point
//              cout<< CH[i].x<< " " << CH[i].y <<endl;
              if(CH[i].x == 0||CH[i].y < 0) continue;
              if(CH[i].x-cur_x == 0 ||((CH[i].y-cur_y)/(CH[i].x-cur_x)>=0)){ //slope positive

            	  sub[fuc][r++]=CH[i].ind;
//            	  cout << " the positive index " << CH[i].ind << endl;
              }else{ // The slope become negative
//            	  cout << "The slope become negative" <<endl;
            	  break;
            	  //cout << "ignore node index first = "<< CH[i].ind<< endl;
              }
              cur_x = CH[i].x;
              cur_y = CH[i].y;
            }
        }else{	   
           for(int z = 1; z<=m; z++){ //if version less then 2
        	   if(CH[z].x == 0) continue;
        	   if(CH[z].y < 0) continue;
               sub[fuc][r++]=CH[z].ind; //assign the subitems of functions to reduce the burden of construction           
           }
        }
	    sub[fuc][0]=r-1; //function versions number in sub[fuc][0]
        
        delete [] Point;
        delete [] CH;
    
        return 0;
}

int Monotone_Chain(int n){
    sort(Point, Point+n, cmp);
    int m = 0;
    if(n < 3){
        for(int i = 0; i < n; i++){
            CH[m++] = Point[i];
        }
    }else{

	    for(int i = 0; i < n; i++){
	        while(m >= 2 && cross(CH[m-2], CH[m-1], Point[i]) >= 0)   m--;
	        CH[m++] = Point[i];		
	    }

	    for(int i = n-2, t = m+1; i >= 0; i--){     //t = m+1, already cover start point
        	    while(m >= t && cross(CH[m-2], CH[m-1], Point[i]) >= 0)   m--;
            	CH[m++] = Point[i];
	    }
        
    }
    
    return m-1;
}

// Calculate Rmin(i) and Rmax(i) for all i values and store them into minR and maxR arrays
// Note: Rmax is stored in maxR[0]
void DynamicAlgo::minMaxR(){
	double min,max;
	minR = new double[n];
	maxR = new double[n];

	minR[0]=maxR[0]=0;
	for(int l=0;l<n;l++){
		max=R[l][0];
		min=(int)(R[l][0]*tableScale);
		//for(j=1;j<K;j++){
		for(j=1;j<size[l];j++){
			if((int)(R[l][j]*tableScale)<min) min=(int)(R[l][j]*tableScale);
			if(R[l][j]>max) max=R[l][j];
		}
		if(l!=(n-1)){
			minR[(l+1)]=min+minR[l];
			maxR[(l+1)]=max+maxR[l];
		}
		else{
			minR[0]=min+minR[l];       //Rmin is stored in minR[0] because minR[n-1] is used for Rmin(n)
			maxR[0]=max+maxR[l];       //Rmax is stored in maxR[0] because maxR[n-1] is used for Rmin(n)
		}
	}

	for(int l=0;l<n;l++)
			minR[l]=minR[l]/tableScale;

}

// Generates (or reads) the values of R(i,j), C(i,j,e), P(i,j,e)
void DynamicAlgo::getData(int mmode) {
	// Initialize dynamic arrays
	int w;
	special_index = 0;
	H = new double[K];
	theta = new int[n];  // predefined from the user, read it from data file
	garma = new int[n];
	th_data = new order_th[n];
	if(mmode == 0)
	O_theta = new int[n];  // read it from mmode 0;

//	theta[0]=1;theta[1]=1;theta[2]=0;theta[3]=1;theta[4]=1; // <<--
	if(mmode == 0){
//		theta[0]=3;theta[1]=14;theta[2]=0;theta[3]=0;theta[4]=2; // <<--
/*for all data*/
#if DYNA // don't forget to use the data of all input
//		theta[0]=1;theta[1]=17;theta[2]=0;theta[3]=1;theta[4]=4; // prefered version from reliability
		theta[0]=1;theta[1]=11;theta[2]=0;theta[3]=1;theta[4]=4; // prefered version from reliability
		garma[0]=1;garma[1]=11;garma[2]=0;garma[3]=1;garma[4]=4; // for such a case there is no version dominate each other
#else
/*for static input data*/
		theta[0]=0;theta[1]=0;theta[2]=0;theta[3]=0;theta[4]=0;
		garma[0]=0;garma[1]=0;garma[2]=0;garma[3]=0;garma[4]=0; // for such a case there is no version dominate each other
#endif
	}
	else{
		theta[0]=O_theta[0];theta[1]=O_theta[1];theta[2]=O_theta[2];theta[3]=O_theta[3];theta[4]=O_theta[4];
	}
	//theta[0]=1;theta[1]=1;theta[2]=0; // <<--
	//khchen
	subitems = new int*[n];
	for(i=0;i<n;i++)
		subitems[i]=new int[K];

	size = new int[n];
	for(i=0;i<n;i++)
		size[i]=0;



	R = new double*[n];
	for(i=0;i<n;i++)
		R[i]=new double[K];
	O_R = new double*[n];
	for(i=0;i<n;i++)
		O_R[i]=new double[K];

	WI = new double*[n];
		for(i=0;i<n;i++)
			WI[i]=new double[K];

	double **AVG = new double*[n];
	for(i=0;i<n;i++)
		AVG[i]=new double[K];

	O_AVG = new double*[n];
	for(i=0;i<n;i++)
		O_AVG[i]=new double[K];

	MU = new double*[n];
		for(i=0;i<n;i++)
			MU[i]=new double[K];
	// Reads R(i,j) from Rij.txt file
	string rfile;
	if(mmode == 0) special_fun = -1;
#if DYNA
	rfile = "GetSure_CDF_Origin/Rij.txt"; //for all input
#else
	rfile = "GetSure_CDF_Origin_STATIC_AVG/Rij.txt"; //for static
#endif
	if(mmode == 1) //second time test the static ordering given by mmode 0
	rfile = "esweek_input/Rij_test.txt";
	ifstream in(rfile.c_str());
	if(!in){
		cout << "Cannot open Rij.txt file."<<endl;
		//return 0;
  	}
	int iInd,jInd;
	double rValue,avgVAlue;
	while(1){
	//for(i=0;i<36;i++){
		in>>iInd>>jInd>>rValue>>avgVAlue;
		if(in.eof()) break;
		//R[iInd][jInd]=(int)(rValue/RScale);

		R[iInd][jInd]=(rValue/RScale);
		AVG[iInd][jInd]=avgVAlue;
		size[iInd]+=1;
		if(mmode == 0){
			O_R[iInd][jInd]=rValue;
			O_AVG[iInd][jInd]=avgVAlue;
	
		}

	}
	in.close();

	C = new double**[n];
	for(i=0;i<n;i++){
		C[i]=new double*[K];
		for(j=0;j<K;j++){
//			if(Dmax < 100)
//				C[i][j]=new double[100+1];
//			else
//				C[i][j]=new double[Dmax+1];
			C[i][j]=new double[300+1];
		}
	}

	for(i=0;i<n;i++){
		for(j=0;j<K;j++){
			for(int m=0;m<=300;m++)
				C[i][j][m]=0;
		}
	}

	P = new double**[n];
	for(i=0;i<n;i++){
		P[i]=new double*[K];
		for(j=0;j<K;j++)
			P[i][j]=new double[300+1];
	}

	for(i=0;i<n;i++){
		for(j=0;j<K;j++){
			for(int m=0;m<=300;m++)
					P[i][j][m]=0;
		}
	}


	// Reads C(i,j,e) and P(i,j,e) Old
	string fileName;
	int eSample;
	double cValue,cValueOld=0;
	for(i=0;i<n;i++){
		for(w=0;w<size[i];w++){
			//cout << "number of version" <<size[i] <<endl;
			fileName="GetSure_CDF_Origin/";
//			fileName="GetSure_CDF_Results/";
			if(mmode == 1) //second time test the static ordering given by mmode 0
			fileName="GetSure_CDF_hand/";
			fileName+=('0'+i);
			fileName+="/";
			if((w/10)!=0) fileName+=((48+w/10));
			fileName+=((48+w%10));
			fileName+=".txt";
			if(mmode==0){
				// Reads the file
				if(i!=1){ // <<--
					in.open(fileName.c_str());
					if(!in){
						cout << "Cannot open "<<fileName<<" file."<<endl;
					}
					t=0;
					cValueOld=0;
					//for(int m=0;m<20;m++){
					for(;;){
						in>>eSample>>cValue;
						if(in.eof()) break;
						/*KHCHEN TESTING*/
						for(;t<(int)(eSample/100000.0) && t <= 300;t++){
								C[i][w][t]=cValueOld;
								//cout << "t=" << t << " cVo=" << cValueOld <<endl;
						}
						cValueOld=cValue;

					}
					/*KHCHEN TESTING*/
					for(;t<=300;t++){
						C[i][w][t]=1;
					}
					in.close();
				} // <<--
				else{ //C for greater than the AVG is equal to 1 // <<--
					/*KHCHEN TESTING*/
					for(t=0;t<=300;t++){ // <<--
						if(t>=(int)(AVG[i][w]/100000.0)) C[i][w][t]=1; // <<--
					} // <<--
				} // <<--
			}else{ //mmode == 1
				if(i!=special_fun){ // <<--
						in.open(fileName.c_str());
						if(!in){
							cout << "AAACannot open "<<fileName<<" file."<<endl;
						}
						t=0;
						cValueOld=0;
						//for(int m=0;m<20;m++){
						for(;;){
							in>>eSample>>cValue;
							if(in.eof()) break;
							/*KHCHEN TESTING*/
							for(;t<(int)(eSample/100000.0) && t <= 300;t++){
									C[i][w][t]=cValueOld;
									//cout << "t=" << t << " cVo=" << cValueOld <<endl;
							}
							cValueOld=cValue;

						}
						/*KHCHEN TESTING*/
						for(;t<=300;t++){
							C[i][w][t]=1;
						}
						in.close();
					} // <<--
					else{ //C for greater than the AVG is equal to 1 // <<--
						/*KHCHEN TESTING*/
						for(t=0;t<=300;t++){ // <<--
							//if(t>=round(AVG[i][w]/100000.0)) C[i][w][t]=1; // <<--
							if(t>=(int)(AVG[i][w]/100000.0)) C[i][w][t]=1; // <<--
						} // <<--
					} // <<--
			}
		}
	}

	// Assign P(i,j,e) values
	for(i=0;i<n;i++){
		for(j=0;j<size[i];j++){
			P[i][j][0] = 0;
			for(int m=1;m<=300;m++)
				P[i][j][m]=C[i][j][m]-C[i][j][(m-1)];
		}
	}

	/*khchen: after get the P() and C(), now we can calculate the CDFz(), PDFz()*/
	/*calculate the summation of function execution time*/
	if(mmode == 1){
		for(i=0;i<n;i++){
			theta[i] = O_theta[i]; //no need to calculate again.
			//cout << theta[i] <<endl;
		}
	}

	else if(mmode == 0){

	int hp = -1, hp_ver[n], check = 0;
	double sum_fun = 0.0;
	double exp[n][K]; //temp expected gap

	//計算完所有點的w_i跟mu_i

//seg 11
		//把所有的convex hull的點合並進greedy list to sort.

	}


//	/* Display the gap of hp and lp
//	 * khchen */
//	if(mmode == 0){
//		for(i=0;i<n;i++){
////				for(w=0;w<size[i];w++){
//					//cout<<endl<<"Function: "<<i<<", Version: "<<w<<endl;
//					for(t=0;t<=Dmax;t++){
//
//						if(i > 0){
//						cout<<t<<": "<<C[i][1][t]-C[i][0][t]<<endl;
//						cout<<t<<"hp: "<<C[i][1][t]<<endl;
//						cout<<t<<"lp: "<<C[i][0][t]<<endl;
//						}
//					}
//
////				}
//		}
//	}

/*
	for(i=0;i<n;i++){
		for(j=0;j<size[i];j++){
			double prob_t = 0;
			for(int m=1;m<=Dmax;m++){
				P[i][j][m]= 0;
			}
			for(int m=;m<=Dmax;m+=){
				P[i][j][m]=C[i][j][m]-prob_t;
				prob_t = C[i][j][m];
			}
		}
	}   
*/

	// Prints the gap of the different versions array P or C
#if 0
	if(mmode==0){
		int u = 0, i=4;//i = function
		for(w=0;w<size[i];w++){ //assume 0 is hp version for each function
			cout<<"Gap of Function: "<<i<<", Version hp & lp "<<w<<endl;
#if 1
			for(t=-150;t<=200;t++){
				if(w == 2) continue;
				cout<<integ_G(i, w, 2, t)<<endl; //for z
			}
#else
			for(u=0; u<size[i]; u++){
				if(w == u) continue;
//				cout<<"version "<<u <<" "<<integ_G(i, w, u, 0)<<endl; //for z
				for(t=-0;t<=200;t++){
//					if(t == 0) cout << "here"<<endl;
					cout<<"t = "<<t<<" "<<"version "<<u <<" "<<integ_G(i, w, u, t)<<endl;
				}
			}
#endif
		}
	}
#endif

#if 0
	// Prints the contents of the array P or C
	if(mmode==0){
	i=4;
//		for(i=0; i<n; i++){
			for(w=0;w<size[i];w++){
				cout<<endl<<"Function: "<<i<<", Version: "<<w<<endl;
#if 1
				for(t=0;t<=100;t++){
//					cout<<t<<": "<<C[i][w][t]<<endl;
					cout<<C[i][w][t]<<endl;
				}
#else
				for(t=0;t<=Dmax;t++){
					cout<<P[i][w][t]<<endl;
				}
#endif

			}
		}
//	}
#endif
	double RTvalue = 100;
	int e;
	/*
	for(i=0; i<n; i++){
				for(w=0;w<size[i];w++){
					if((alpha/alphaScale)*(R[i][w])+(1-alpha)*(1-C[i][w][t])<RTvalue){
						e = w;
						RTvalue = (alpha/alphaScale)*(R[i][w])+(1-alpha)*(1-C[i][w][t]);
					}
				}
				cout<<"func : "<<i<<" version: "<< e <<endl;
	}
	*/

	// Prints the contents of the array R
//	 for(i=0;i<n;i++)
//		for(w=0;w<size[i];w++)
//		//	cout<<i<<" "<<w<<": "<<R[i][w]<<endl;
//			cout<<i<<" "<<w<<": "<<AVG[i][w]<<endl;

	////////////////////////////////////////////////////until here, we have convex hull subitems and P, C, Rmax, Rmin arrays.
	minMaxR();
	Rmax = ceil((maxR[0])/sigma)*sigma;
	Rmax = ceil( (maxR[0]*tableScale) );



	//j* and G initial here
	G = new double**[n];
	for(i=0;i<n;i++){
		G[i]=new double*[Rmax+1];
		//for(j=0;j<(maxR[0]);j++)
		for(j=0;j<=Rmax;j++)
			G[i][j]=new double[(Dmax+1)];
	}

	jStar = new int**[n];
	for(i=0;i<n;i++){
		jStar[i]=new int*[Rmax+1];
		//for(j=0;j<(maxR[0]);j++)
		for(j=0;j<=Rmax;j++)
			jStar[i][j]=new int[(Dmax+1)];
	}

	//j* and G initial here
}
//khchen

double DynamicAlgo::integ_Z(int f, int v, int y, int e){ //r = function, v = version
	double sum=0;
	for(int x=0;x<=y+e;x++){
//		cout<<"X "<<x<<endl;
//		if(x > 100){
//			sum+=0;
//		}else{
			sum+=(P[f][v][x]);
//		}
	}
//	cout<<sum<<endl;
	return sum;
}

// Returns the gap of x and y
double DynamicAlgo::integ_G(int f, int x, int y, int e){ //Here r==function y=Y q==X version  e == t, calculate (Q-P)
	double sum=0;
	for(int p=0;p<=200;p++){ //p = dy
			sum = sum + (P[f][y][p] * integ_Z(f, x, p, e));
	}
//	cout<<sum<<endl;
	return sum;
}

// Returns the integration
double DynamicAlgo::integ_P(int r, int s){
	double sum=0;
    t=Dmax; //t global
	for(int x=0;x<=t;x++){
		sum+=((double)x*(P[r][s][x])); //expected value, multiply x
	}

	return sum;
}


// Returns the integration
double DynamicAlgo::integ(int j){
	double sum=0;

	for(int x=0;x<=t;x++){
		//if(x>0 && x<t)
			//sum+=(P[i][j][x]*G[(i+1)][(int)(floor((r+R[i][j])/sigma)*sigma)][(int)(floor((t-x)/delta)*delta)]);
		sum+=(P[i][j][x]*G[(i+1)][(int)( floor(((r+R[i][j])*tableScale)) ) ][(int)(floor((t-x)/delta)*delta)]);
		//else
			//sum+=((P[i][j][x]*G[(i+1)][(int)(floor((r+R[i][j])/sigma)*sigma)][(int)(floor((t-x)/delta)*delta)])/2);
	}

	return sum;
}

// Returns the sum of reliability levels for the functions from u to n with versions predetermined by the user
double DynamicAlgo::rho(int u){
	double sum=0;

	for(int l=u;l<n;l++){
		sum+=R[l][(theta[l])];
	}
	//cout<<u<<" :"<<sum<<endl;
	return sum;
}

// Calculates Hj (for j=1,2,...,Ki), returns the maximum value
// and stores the index of the minimum value into minInd variable
double DynamicAlgo::minH(){ //i is global variable, to record the function item. j is version, t is deadline.
	double min;

	//for(j=0;j<K;j++){
	for(j=0;j<size[i];j++){
		if(i!=n-1)
			H[j]=integ(j)+(1-C[i][j][t])*( (alpha/alphaScale)*(r+R[i][j]+rho(i+1)) +(1-alpha));
		else
			H[j]=(alpha/alphaScale)*(r+R[i][j])+(1-alpha)*(1-C[i][j][t]); // for i=n (here n-1) and may be i=1 (here 0)

		if(j==0){
			min=H[j];
			minInd=0;
		}
		else
			if(H[j]<min){
				min=H[j];
				minInd=j;
			}
	}


	return min;
}

/*
// Stores the index of the version that achieves the minimum expected RT penalty into minInd
// for the current values of i(n and 1),t and r, and returns the maximum value
double DynamicAlgo::minJ(){
	double min;



	return min;
}
*/


bool cmp_density(Snode a, Snode b) {
    bool result = false;
    if(a.density < b.density){ 
        result =true;
    }
    return result;
} 


// Fills the order table and returns the RT penalty
void DynamicAlgo::fillOrder(){
	/* 目前已經有一個greedy的list
	 * 接下來要做的是sort跟決定好ordering
	 * 丟給RTAS algorithm跑
	 * */
    // construct subitems properties
	int fi=0, current_fuc=-1; //assume first function's first version is 0
	double current_et=0, current_val= 0;
	reverse(greedy, greedy+subindex);
	for(i=0; i < subindex; i++){

	//計算每個greedy list裡面的item對於自己function所增加的add_et, 用以計算slope
		if(current_fuc != greedy[i].func){
			current_fuc = greedy[i].func; //found new function
			fi=i; //now is the function first version
	        greedy[i].addet=greedy[i].MU;
	        current_et = greedy[i].MU;
	        current_val = greedy[i].WI;
	        greedy[i].density = greedy[i].WI/greedy[i].MU;
		}
		else{		

	        greedy[i].addet= greedy[i].MU - current_et;
	        greedy[i].density = (greedy[i].WI - current_val)/greedy[i].addet; //Addet已經存好增加的et
	        if(greedy[i].addet == 0){
	        	greedy[i].addet=greedy[i].MU;
				greedy[i].density = greedy[i].WI/greedy[i].MU;
	        }
	        current_et = greedy[i].MU;
	        current_val = greedy[i].WI;
//	        cout << "Density " << greedy[i].density << endl;
		}
		greedy[i].R = O_R[greedy[i].func][greedy[i].ver];
		greedy[i].AVG = O_AVG[greedy[i].func][greedy[i].ver];
//		cout<<"greedy sorted list "<<i<<" density =" <<greedy[i].density<<endl;
//		cout << "==========================" <<endl;
	}

	////////到這之前，已經把所有的relative slope算好存進去了

    //order 1~s with cmp_density()
	sort(greedy, greedy+subindex, cmp_density); //in Ascending
	reverse(greedy, greedy+subindex);


//    cout << "---------Greedy ordering: " <<endl;
//    for(i=0; i < subindex; i++){
//    	cout<<"greedy sorted list "<<i<<" density =" <<greedy[i].density<<endl;
//    	cout<<"greedy sorted list "<<i<<" function =" <<greedy[i].func<<endl;
//    	cout << "==========================" <<endl;
//    }

    /* 接下來要把同樣的function出現不同version的時候的處理，
     * 選最後出現的version去算slope，決定先執行的是哪一個
     * */

	double r_sum = 0, best_r = 0;
	double best_value=0, best_sum=0, sum_val=0;
	double curfunc_WI[n];
	int curfunc_ind[n], curfunc_func[n];
	int best_skip[subindex];
	int r = 1; //index
	int kver = 0, best_first=0, indd=0;




	/*實驗區塊 去掉iteration且用density相加，越高者越優先選擇 for N1 ordering*/
	cur = new Snode[n];
	for(i=0;i<n; i++){ //init
		cur[i].density = 0;
		cur[i].func=i;
	}
	r=0;
//	cout << "subindex " << subindex << " deadline " << Dper <<endl;
	for(i=0; i<subindex; i++){ //這邊計算r可以到多少順便照greedy array裡面的順序放到cur array準備排序
		r_sum += greedy[i].addet;
		sum_val += greedy[i].WI;
//		cout << "r_sum "<< r_sum <<endl;
		if(r_sum >= Dper){ //r_sum已經超過可用的slack 超過的話把最後一個加上去的資訊去掉 break for loop. (雖然這數值沒有用到
			sum_val = sum_val - greedy[i].WI;
			r_sum = r_sum - greedy[i].addet;
			break;
		}

//		cur[greedy[i].func].density=greedy[i].density;
		cur[greedy[i].func].density=greedy[i].WI/greedy[i].MU;
		r++;

	}
//	cout<<"R index ="<<r<<endl;
	//把已經到達總和限制的function version排列從小到大
	sort(cur, cur+n, cmp_density);
	reverse(cur, cur+n);


//    cout << "---------curfunc_density ordering: " <<endl;
//    for(i=0; i < n; i++){
//    	cout<<"curfunc_density sorted list "<< i <<" "<<cur[i].func<<" density =" <<cur[i].density<<endl;
//    	cout << "==========================" <<endl;
//    }

    /*Selection ordering*/
	ofstream fileout;
	string outputfile;
	outputfile="esweek_input/Rij_test";
//	outputfile+=('0'+i);
	outputfile+=".txt";
	fileout.open(outputfile.c_str());
	int cur_index, exist_num=0, nk = 0;
	int record_func[n];
	string ofs;
	string ifs;
	int eSample;
	double cValue;
	for(i=0;i<n;i++) {
		record_func[i]=0;
	}

	/*
	 *  如果D'小於0，全部都直接用High performance.
	 *  N1->N2, N1用Cur[i].func找哪個function 絕對w/mu最高排前面
	 *  這邊透過重新建立餵給DP的input找到最佳version
	 * */

	//沒有density的function 表示是paper中的Theta = -1的function continue等一下在N2時直接用Theta.
	for(i=0;i<n; i++){ //index the function
		if(cur[i].density == 0) {
//				cout << cur[i].func <<endl;
			continue;
		}

		if(cur[i].func == 1){
				O_theta[nk]=14;
//			O_theta[nk]=23;
			special_fun = nk;
			for(int w=0; w<size[1]; w++){
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[1][w]<<" ";
				fileout<<O_AVG[1][w];
				fileout<<"\n";
			}
		}else{
			O_theta[nk]=theta[cur[i].func];
			for(int w=0; w<size[cur[i].func]; w++){ //put all the version in function nk folders;
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[cur[i].func][w]<<" ";
				fileout<<O_AVG[cur[i].func][w];
				fileout<<"\n";

				ifs="GetSure_CDF_Origin/";
				ifs+=('0'+cur[i].func);
				ifs+="/";
				if((w/10)!=0) ifs+=((48+w/10));
				ifs+=((48+w%10));
				ifs+=".txt";
//				cout<< ifs <<endl;
				ofs="GetSure_CDF_hand/";  //0 function
				ofs+=('0'+nk);
				ofs+="/";
				if((w/10)!=0) ofs+=((48+w/10));
				ofs+=((48+w%10));
				ofs+=".txt";
//				cout<< ofs <<endl;
				ifstream file1;
				ofstream file2;
				file1.open(ifs.c_str(), ios::in |ios::binary);
				file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
				for(;;){
					file1>>eSample>>cValue;
	//				cout << eSample << " "<< cValue <<endl;
					if(file1.eof()) break;
					file2<<eSample<<" "<<cValue;
					file2<<"\n";
				}
				file1.close();
				file2.close();
			}
		}
		record_func[cur[i].func] = 1;
		nk++;
	}

	//這是second type的部分 N2
	int k=0;
	for(i=nk;i<n;i++){//remaining function
		for(k=0; k<n; k++){
			//cout <<th_et[k].ind<< endl;
			if(record_func[th_data[k].ind] !=1) break;
		}
//			cout << "static function "<<nk<<" = origin " <<th_data[k].ind <<endl;
		if(th_data[k].ind == 1){
			special_fun = nk;
				O_theta[nk]=14;
//			O_theta[nk]=23;
			for(int w=0; w<size[1]; w++){
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[1][w]<<" ";
				fileout<<O_AVG[1][w];
				fileout<<"\n";
				record_func[th_data[k].ind] = 1;
			}
		}
		else{
			O_theta[nk]=theta[th_data[k].ind ];
			for(int w=0; w<size[th_data[k].ind ]; w++){
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[th_data[k].ind ][w]<<" ";
				fileout<<O_AVG[th_data[k].ind][w];
				fileout<<"\n";
				ifs="GetSure_CDF_Origin/";
				ifs+=('0'+th_data[k].ind);
				ifs+="/";
				if((w/10)!=0) ifs+=((48+w/10));
				ifs+=((48+w%10));
				ifs+=".txt";
				ofs="GetSure_CDF_hand/";  //0 function
				ofs+=('0'+nk);
				ofs+="/";
				if((w/10)!=0) ofs+=((48+w/10));
				ofs+=((48+w%10));
				ofs+=".txt";
				ifstream file1;
				ofstream file2;
				file1.open(ifs.c_str(), ios::in |ios::binary);
				file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
				for(;;){
					file1>>eSample>>cValue;
	//				cout << eSample << " "<< cValue <<endl;
					if(file1.eof()) break;
					file2<<eSample<<" "<<cValue;
					file2<<"\n";
				}
				file1.close();
				file2.close();
			}
		}
		record_func[th_data[k].ind] = 1;
		nk++;
	}
	fileout.close();


    /* 原本做iteration 跟 亂做同function加法*/
#if 0
	for(int kver=0; kver<subindex; kver++){ //init
		best_skip[kver] = 0; //use to record the best skip way
	}
	for(kver=0; kver<subindex; kver++){ //handle different k first version (Check k sorted list for different possibility)
		for(int ver=0; ver<subindex; ver++){
			greedy[ver].skip = 0;
		}
		for(i=0;i<n; i++){ //init
			curfunc_WI[i]=0;
			curfunc_ind[i]=0;
			curfunc_density[i]=0;
		}
		sum_val = 0;
		r = 1;

		/*先記錄第一個item的val跟mu, sum_val用來評估優劣，r_sum用來評估r*/
		r_sum = greedy[kver].addet;
//		r_sum = greedy[kver].MU;
		sum_val = greedy[kver].WI;

		curfunc_WI[greedy[kver].func]=WI[greedy[kver].func][greedy[kver].ver];
		curfunc_ind[greedy[kver].func]=kver;
		cout << "==============="<<endl;
		cout << "First K index : " <<kver<<" Profit density "<<greedy[kver].density<<" MU "<<greedy[kver].MU<<endl;
		cout << "function: "<<greedy[kver].func<<" ver: "<<greedy[kver].ver<<" Reliability profit "<<  WI[greedy[kver].func][greedy[kver].ver] <<endl;
		/*******************************************************/

		/*接著要計算後續的部分，要用last appearance的slope去重新排序*/
		for(i=0; i<subindex; i++){
			//find the minimum index r such that sum mu 1~r<D
			if(kver == i) continue; //表示算到重複的K


			if(curfunc_WI[greedy[i].func]>WI[greedy[i].func][greedy[i].ver]){
//				cout << "skip subitem  : " <<i<<" Profit density "<<greedy[i].density<<" MU "<<greedy[i].MU<<endl;
//				cout << "function: "<<greedy[i].func<<" ver: "<<greedy[i].ver<<" Reliability profit "<<  WI[greedy[i].func][greedy[i].ver] <<endl;
				greedy[i].skip = 1;
				r++;
				continue;
			//if current wi > now, no need to add it.
			}else{
				indd = curfunc_ind[greedy[i].func];
				if(greedy[i].addet < 0){
//					cout << "skip subitem  : " <<i<<" Profit density "<<greedy[i].density<<" MU "<<greedy[i].MU<<endl;
//					cout << "function: "<<greedy[i].func<<" ver: "<<greedy[i].ver<<" Reliability profit "<<  WI[greedy[i].func][greedy[i].ver] <<endl;

					greedy[i].skip = 1;
					r++;
					continue;
				}
				else{
					r_sum += greedy[i].addet;

//							cout << "subitem order : " <<i<<" Profit density "<<greedy[i].density<<" MU "<<greedy[i].MU<<endl;
//							cout << "function: "<<greedy[i].func<<" ver: "<<greedy[i].ver<<" Reliability profit "<<  WI[greedy[i].func][greedy[i].ver] <<endl;

					/*special case, when the later one become before one*/

					sum_val += greedy[i].WI;
		//					if(r_sum <= 0){
		//						cout << "CCCCC"<<greedy[i].MU<<endl;
		//						cout << "BBBBBBB"<<greedy[i].addet<<endl;
		//						cout << "AAAAAAAAAAA" <<endl;
		//					}
					if(curfunc_WI[greedy[i].func] != 0){
						sum_val = sum_val - curfunc_WI[greedy[i].func];
						greedy[indd].skip = 1; //just choose the latest version
						if(r_sum >= Dper){
							sum_val = sum_val - greedy[i].WI;
							r_sum = r_sum - greedy[i].addet;
							greedy[indd].skip = 0;
							break;
						}
					}else if(r_sum >= Dper){
						sum_val = sum_val - greedy[i].WI;
						r_sum = r_sum - greedy[i].addet;
						break;
					}
					curfunc_WI[greedy[i].func] = WI[greedy[i].func][greedy[i].ver];
					curfunc_ind[greedy[i].func] = i;
					r++;
				}
			}
		}
//			cout <<"value " << sum_val <<endl;
//			cout <<"mu_sum " << r_sum <<endl;
		if(best_value < sum_val){
			best_sum = r_sum;
			best_value = sum_val;
			best_first = kver;
			best_r = r;
			for(int ver=0; ver<subindex; ver++){
				best_skip[ver] = greedy[ver].skip;
//					cout <<ver << " skip " << best_skip[ver] <<endl;
			}
		}

	}

////////////////以下是建立測試用的文件檔

//    	cout <<"Best_value " << best_value <<" Best k " <<best_first<< " Best sum "<<best_sum<<endl;
//    	cout << "R is "<< best_r <<endl;
    /*Selection ordering*/
	ofstream fileout;
	string outputfile;
	outputfile="esweek_input/Rij_test";
//	outputfile+=('0'+i);
	outputfile+=".txt";
	fileout.open(outputfile.c_str());
	int cur_index, exist_num=0, nk = 0;
	int record_func[n];
	string ofs;
	string ifs;
	int eSample;
	double cValue;
	for(i=0;i<n;i++) {
		record_func[i]=0;
	}

	/*
	 *  如果D'小於0，全部都直接用High performance.
	 * */
	if(Dper < 0){
		int k=0;
		for(i=nk;i<n;i++){//remaining function
			for(k=0; k<n; k++){
				//cout <<th_et[k].ind<< endl;
				if(record_func[th_data[k].ind] !=1) break;
			}
			if(th_data[k].ind == 1){ //如果是special case, 就不去CDF資訊
				special_fun = nk;
				O_theta[nk]=14;
				for(int w=0; w<size[1]; w++){
					fileout<<nk<<" "<<w<<" ";
					fileout<<O_R[1][w]<<" ";
					fileout<<O_AVG[1][w];
					fileout<<"\n";
					record_func[th_data[k].ind] = 1;
				}
			}
			else{
				O_theta[nk]=theta[th_data[k].ind ];
				for(int w=0; w<size[th_data[k].ind ]; w++){  //如果不是special case, 就直接印出對應的version資料
					fileout<<nk<<" "<<w<<" ";
					fileout<<O_R[th_data[k].ind ][w]<<" ";
					fileout<<O_AVG[th_data[k].ind][w];
					fileout<<"\n";
					ifs="GetSure_CDF_Origin/";
					ifs+=('0'+th_data[k].ind);
					ifs+="/";
					if((w/10)!=0) ifs+=((48+w/10));
					ifs+=((48+w%10));
					ifs+=".txt";
					ofs="GetSure_CDF_hand/";  //0 function
					ofs+=('0'+nk);
					ofs+="/";
					if((w/10)!=0) ofs+=((48+w/10));
					ofs+=((48+w%10));
					ofs+=".txt";
					ifstream file1;
					ofstream file2;
					file1.open(ifs.c_str(), ios::in |ios::binary);
					file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
					for(;;){
						file1>>eSample>>cValue;
						if(file1.eof()) break;
						file2<<eSample<<" "<<cValue;
						file2<<"\n";
					}
					file1.close();
					file2.close();
				}
			}
			record_func[th_data[k].ind] = 1;
			nk++;
		}
	}else{//如果D'大於0 表示要看我們的ordering了

	/*Best first should be consider first*/
	/*find one function*/
//	cout << "first function 0 = origin " <<greedy[best_first].func <<endl;
	if(greedy[best_first].func == 1){
//				O_theta[nk]=23;
				O_theta[nk]=14;
				special_fun = nk;
				for(int w=0; w<size[1]; w++){
					fileout<<nk<<" "<<w<<" ";
					fileout<<O_R[1][w]<<" ";
					fileout<<O_AVG[1][w];
//					cout<<nk<<"A ";
//					cout<<O_R[1][w]<<" ";
//					cout<<O_AVG[1][w]<<endl;
					fileout<<"\n";
				}
	}else{
		O_theta[nk]=theta[greedy[best_first].func];
		for(int w=0; w<size[greedy[best_first].func]; w++){ //put all the version in function nk folders;

			fileout<<nk<<" "<<w<<" ";
			fileout<<O_R[greedy[best_first].func][w]<<" ";
			fileout<<O_AVG[greedy[best_first].func][w];
//			cout<<nk<<"A ";
//			cout<<O_R[greedy[best_first].func][w]<<" ";
//	 	 	cout<<O_AVG[greedy[best_first].func][w]<<endl;
			fileout<<"\n";

			ifs="GetSure_CDF_Origin/";
			ifs+=('0'+greedy[best_first].func);
			ifs+="/";
			if((w/10)!=0) ifs+=((48+w/10));
			ifs+=((48+w%10));
			ifs+=".txt";
//			cout<< ifs <<endl;
			ofs="GetSure_CDF_hand/";  //0 function
			ofs+=('0'+nk);
			ofs+="/";
			if((w/10)!=0) ofs+=((48+w/10));
			ofs+=((48+w%10));
			ofs+=".txt";
//			cout<< ofs <<endl;
			ifstream file1;
			ofstream file2;
			file1.open(ifs.c_str(), ios::in |ios::binary);
			file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
			for(;;){
				file1>>eSample>>cValue;
//				cout << eSample << " "<< cValue <<endl;
				if(file1.eof()) break;
				file2<<eSample<<" "<<cValue;
				file2<<"\n";
			}
			file1.close();
			file2.close();
		}
	}
	record_func[greedy[best_first].func] = 1;
	best_skip[best_first] = 1;
	//cout << greedy[best_first].func<<endl;
	nk++;

	for(i=0;i<best_r; i++){ //index the function
		if(best_skip[i] == 1) continue; //this one skip //
		if(greedy[i].func == greedy[best_first].func) continue;
		/*find one function*/
//		cout << "Preorder function "<<nk<<" = origin " <<greedy[i].func <<endl;
		if(greedy[i].func == 1){
//					O_theta[nk]=23;
					O_theta[nk]=14;
					special_fun = nk;
					for(int w=0; w<size[1]; w++){
						fileout<<nk<<" "<<w<<" ";
						fileout<<O_R[1][w]<<" ";
						fileout<<O_AVG[1][w];
//						cout<<nk<<"B ";
//						cout<<O_R[1][w]<<" ";
//						cout<<O_AVG[1][w]<<endl;
						fileout<<"\n";
					}
		}else{
			O_theta[nk]=theta[greedy[i].func];
			for(int w=0; w<size[greedy[i].func]; w++){ //put all the version in function nk folders;
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[greedy[i].func][w]<<" ";
				fileout<<O_AVG[greedy[i].func][w];
//				cout<<nk<<"B ";
//				cout<<O_R[greedy[i].func][w]<<" ";
//		 	 	cout<<O_AVG[greedy[i].func][w]<<endl;
				fileout<<"\n";

				ifs="GetSure_CDF_Origin/";
				ifs+=('0'+greedy[i].func);
				ifs+="/";
				if((w/10)!=0) ifs+=((48+w/10));
				ifs+=((48+w%10));
				ifs+=".txt";
//				cout<< ifs <<endl;
				ofs="GetSure_CDF_hand/";  //0 function
				ofs+=('0'+nk);
				ofs+="/";
				if((w/10)!=0) ofs+=((48+w/10));
				ofs+=((48+w%10));
				ofs+=".txt";
//				cout<< ofs <<endl;
				ifstream file1;
				ofstream file2;
				file1.open(ifs.c_str(), ios::in |ios::binary);
				file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
				for(;;){
					file1>>eSample>>cValue;
	//				cout << eSample << " "<< cValue <<endl;
					if(file1.eof()) break;
					file2<<eSample<<" "<<cValue;
					file2<<"\n";
				}
				file1.close();
				file2.close();
			}
		}
		record_func[greedy[i].func] = 1;
		nk++;
	}

	//這是second type的部分
	int k=0;
	for(i=nk;i<n;i++){//remaining function
		for(k=0; k<n; k++){
			//cout <<th_et[k].ind<< endl;
			if(record_func[th_data[k].ind] !=1) break;
		}
//		cout << "static function "<<nk<<" = origin " <<th_data[k].ind <<endl;
		if(th_data[k].ind == 1){
//			special_fun = 4;
			special_fun = nk;
//			cout <<special_fun<< "HEREERERER"<<endl;
//			O_theta[4]=23;
//			O_theta[nk]=23;
			O_theta[nk]=14;
			for(int w=0; w<size[1]; w++){
//				fileout<<4<<" "<<w<<" ";
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[1][w]<<" ";
				fileout<<O_AVG[1][w];
//				cout<<nk<<"C ";
//				cout<<O_R[1][w]<<" ";
//				cout<<O_AVG[1][w]<<endl;
				fileout<<"\n";
				record_func[th_data[k].ind] = 1;
			}
//			nk--;
		}
		else{
			O_theta[nk]=theta[th_data[k].ind ];
			for(int w=0; w<size[th_data[k].ind ]; w++){
				fileout<<nk<<" "<<w<<" ";
				fileout<<O_R[th_data[k].ind ][w]<<" ";
				fileout<<O_AVG[th_data[k].ind][w];
//				cout<<nk<<"C ";
//				cout<<O_R[k][w]<<" ";
//				cout<<O_AVG[k][w]<<endl;
				fileout<<"\n";
				ifs="GetSure_CDF_Origin/";
				ifs+=('0'+th_data[k].ind);
				ifs+="/";
				if((w/10)!=0) ifs+=((48+w/10));
				ifs+=((48+w%10));
				ifs+=".txt";
				ofs="GetSure_CDF_hand/";  //0 function
				ofs+=('0'+nk);
				ofs+="/";
				if((w/10)!=0) ofs+=((48+w/10));
				ofs+=((48+w%10));
				ofs+=".txt";
				ifstream file1;
				ofstream file2;
				file1.open(ifs.c_str(), ios::in |ios::binary);
				file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
				for(;;){
					file1>>eSample>>cValue;
	//				cout << eSample << " "<< cValue <<endl;
					if(file1.eof()) break;
					file2<<eSample<<" "<<cValue;
					file2<<"\n";
				}
				file1.close();
				file2.close();
			}
		}
		record_func[th_data[k].ind] = 1;
		nk++;
	}
	fileout.close();
//	cout << nk << endl;

}
#endif
}



 	/*Used to make the special data for test*/



// Fills the dynamic programming table and returns the RT penalty
double DynamicAlgo::fillTable(){
int count=0;
	// The base case i=n (here index n-1)
	i=n-1;
	for(r=minR[i];r<=maxR[i];r+=(1/tableScale)){
//	for(r=0;r<=maxR[i];r+=(1/tableScale)){
		//for(r=0;r<=(ceil(maxR[i]/sigma)*sigma);r+=sigma){
		for(t=0;t<=Dmax;t+=delta){
			if(t==0){
				jStar[i][(int)(r*tableScale)][t]=theta[i];
				G[i][(int)(r*tableScale)][t]=(alpha/alphaScale)*(r+R[i][(theta[i])])+(1-alpha)*(1-C[i][(theta[i])][t]);
			}
			else{
				G[i][(int)(r*tableScale)][t]=minH();
				jStar[i][(int)(r*tableScale)][t]=minInd;
			}
/*
			if (DEBUG) {
				cout << "G["<<i<<"]"<<"["<<r<<"]"<<"["<<t<<"]"<<(int)(r*tableScale)<<endl;
			}
*/
			
		}
	}
	// For i=n-1 to 2 (here index from n-2 to 1)
	for(i=n-2;i>0;i--){//cout<<i<<endl;
//		for(r=0;r<=maxR[i];r+=(1/tableScale)){
		for(r=minR[i];r<=maxR[i];r+=(1/tableScale)){
			for(t=0;t<=Dmax;t+=delta){
				if(t==0){
					jStar[i][(int)(r*tableScale)][t]=theta[i];
					G[i][(int)(r*tableScale)][t]=(alpha/alphaScale)*(r+rho(i))+(1-alpha);
				}
				else{
					G[i][(int)(r*tableScale)][t]=minH();
					jStar[i][(int)(r*tableScale)][t]=minInd;
				}
			}
		}
	}

	//cout<<i<<endl;
	// For i=1 (here index i=0), r=0, t=Dmax
    //testing   start from here
	i=r=0;t=Dmax;
	G[i][0][t]=minH();
	jStar[i][0][t]=minInd;


/*
	int t0=88,t1=24,t2=3,t3=24,t4=32;
	t=0;
	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
	t+=t0;
	i++;
cout<<" ,t: "<<t<<endl;
//Function 2
	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
	t+=t1;
	i++;
cout<<" ,t: "<<t<<endl;
//Function 3
	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
	t+=t2;
	i++;
cout<<" ,t: "<<t<<endl;
//Function 4
	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
	t+=t3;
	i++;
cout<<" ,t: "<<t<<endl;
//Function 5
	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]]<<endl;
	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
	t+=t4;

	cout<<"Total R= "<<r<<endl;
	cout<<"Total t= "<<t<<endl;
	*/
	return G[0][0][Dmax];
}

//private:
//	float alpha;
//	int n,i,j,r,t,e,sigma,delta,K,D,minInd,Rmax,Dmax; //K:max no. of versions, minInd:the index of the
//	double ***G,***C,***P,**R,*H,*minR,*maxR;
//	int *theta,***jStar;
