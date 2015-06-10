/*
 * DynamicAlgo.cpp
 *
 *  Created on: Sep 12, 2012
 *      Author: anas
 */

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
	DEBUG = 0;

	Dmax=D;
	//Dmax = ceil(D/delta)*delta;
}

// Destructor
DynamicAlgo::~DynamicAlgo() {
	// TODO Auto-generated destructor stub
}

/*
double DynamicAlgo::minR(int i){
	double sum=0,min;

	for(int l=0;l<i;l++){
		min=R[l][0];
		for(int j=1;j<K;j++)
			if(R[l][j]<min) min=R[l][j];
		sum+=min;
	}
	return sum;
}*/

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
			minR[0]=min+minR[l];       //Rmin is stored in minR[0] because minR[n-1] is used for Rmin(n), I think we don't need it!
			maxR[0]=max+maxR[l];       //Rmax is stored in maxR[0] because maxR[n-1] is used for Rmin(n)
		}
	}

	for(int l=0;l<n;l++)
			minR[l]=minR[l]/tableScale;

	//for(int l=0;l<n;l++)
	//		cout<<maxR[l]<<" "<<minR[l]<<endl;
	//cout<<"Max. R: "<<maxR[0]<<endl;
	//cout<<"Min. R: "<<minR[0]<<endl;

}

// Generates (or reads) the values of R(i,j), C(i,j,e), P(i,j,e)
void DynamicAlgo::getData() {
	// Initialize dynamic arrays
	int w;
	H = new double[K];
	theta = new int[n];  // predefined from the user, read it from data file
	//theta[0]=1;theta[1]=1;theta[2]=0;theta[3]=1;theta[4]=1; // <<--
//	theta[special_fun0]=0;theta[special_fun]=14;theta[special_fun2]=0;theta[special_fun3]=0;theta[special_fun4]=4;
//	theta[special_fun0]=3;theta[special_fun]=14;theta[special_fun2]=0;theta[special_fun3]=0;theta[special_fun4]=2;
	theta[special_fun0]=1;theta[special_fun]=17;theta[special_fun2]=0;theta[special_fun3]=1;theta[special_fun4]=4;
	//theta[0]=1;theta[1]=1;theta[2]=0; // <<--


	size = new int[n];
	for(i=0;i<n;i++)
		size[i]=0;

	R = new double*[n];
	for(i=0;i<n;i++)
		R[i]=new double[K];

	double **AVG = new double*[n];
	for(i=0;i<n;i++)
		AVG[i]=new double[K];

	// Reads R(i,j) from Rij.txt file
	//ifstream in("GetSure_CDF_Results/Rij.txt");
	ifstream in("estimedia_input/Rij_test.txt");
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
	}
	in.close();

	//in.open("GetSure_CDF_Results/Rij.txt");
	/*
	for(i=0;i<n;i++){
		cout<<size[i]<<endl;
		for(int w=0;w<size[i];w++)
			cout<<i<<" "<<w<<" "<<R[i][w]<<endl;
	}*/

	/*
	// Generates random numbers between 0 and 1 for R(i,j), C(i,j,e), P(i,j,e), read them from data file
	for(i=0;i<n;i++){
		for(j=0;j<K;j++){
			R[i][j]=(rand()/(double)RAND_MAX);
			for(e=0;e<Dmax;e++){
				C[i][j][e]=(rand()/(double)RAND_MAX);
				P[i][j][e]=(rand()/(double)RAND_MAX);
			}
		}
	}*/

	C = new double**[n];
	for(i=0;i<n;i++){
		C[i]=new double*[K];
		for(j=0;j<K;j++){
			C[i][j]=new double[Dmax+1];
		}
	}

	for(i=0;i<n;i++){
		for(j=0;j<K;j++){
			for(int m=0;m<=Dmax;m++)
			C[i][j][m]=0;}}

	P = new double**[n];
	for(i=0;i<n;i++){
		P[i]=new double*[K];
		for(j=0;j<K;j++)
			P[i][j]=new double[Dmax+1];
	}

	for(i=0;i<n;i++){
		for(j=0;j<K;j++){
			for(int m=0;m<=Dmax;m++)
			P[i][j][m]=0;}}


	// Reads C(i,j,e) and P(i,j,e) Old
	string fileName;
	int eSample;
	double cValue,cValueOld=0;
	for(i=0;i<n;i++){
		for(w=0;w<size[i];w++){
			fileName="GetSure_CDF_hand/";
			fileName+=('0'+i);
			fileName+="/";
			if((w/10)!=0) fileName+=((48+w/10));
			fileName+=((48+w%10));
			fileName+=".txt";

			// Reads the file
			if(i!=special_fun){ // <<--
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
					//if ((int)(round(eSample/100000.0)) > Dmax) {cout<<"eSample exceeded the deadline"<<endl;break;}
					//P[i][w][(int)(round(eSample/100000.0))]=cValue-cValueOld;
					//for(;t<round(eSample/100000.0) && t <= Dmax;t++){
					for(;t<(int)(eSample/100000.0) && t <= Dmax;t++){
						C[i][w][t]=cValueOld;
					}
					cValueOld=cValue;

				}
				for(;t<=Dmax;t++){
					C[i][w][t]=1;
				}
				in.close();
			} // <<--
			else{ //C for greater than the AVG is equal to 1 // <<--
				for(t=0;t<=Dmax;t++){ // <<--
					//if(t>=round(AVG[i][w]/100000.0)) C[i][w][t]=1; // <<--
					if(t>=(int)(AVG[i][w]/100000.0)) C[i][w][t]=1; // <<--
				} // <<--
			} // <<--
		}
	}

	// Assign P(i,j,e) values
	for(i=0;i<n;i++){
		for(j=0;j<size[i];j++){
			for(int m=1;m<=Dmax;m++)
			P[i][j][m]=C[i][j][m]-C[i][j][(m-1)];}}

/*
	// Prints the contents of the array P or C
	 for(i=0;i<n;i++){
		for(w=0;w<size[i];w++){
			cout<<endl<<"Function: "<<i<<", Version: "<<w<<endl;
			for(t=0;t<=Dmax;t++){
				cout<<t<<": "<<P[i][w][t]<<endl;
			}

		}
	}*/

//	// Prints the contents of the array R
//	 for(i=0;i<n;i++)
//		for(w=0;w<size[i];w++)
//			cout<<i<<" "<<w<<": "<<R[i][w]<<endl;
//cout <<"AAAAAAA" <<endl;

	minMaxR();
	//Rmax = ceil((maxR[0])/sigma)*sigma;
	Rmax = ceil( (maxR[0]*tableScale) );

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

	// Reads data from file
	/* //---
	// Open the data file
	ifstream in("GetSure_CDF_Results/CDF_ADPCM/CDF_ADPCM_V1.csv");
	if(!in){
		cout << "Cannot open file."<<endl;
		//return 0;
  	}

	double temp1,temp2,comma;

	// Read data from data.txt file
	for(i=0;i<20;i++){
		in>>temp1>>comma>>temp2;
		cout<<temp1<<" "<<temp2<<endl;
	}

	in.close();
	//--- */
}

// Retruns the integration
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

// Retruns the sum of reliability levels for the functions from u to n with versions predetermined by the user
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
double DynamicAlgo::minH(){
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
			if (DEBUG) {
				cout << "G["<<i<<"]"<<"["<<r<<"]"<<"["<<t<<"]"<<(int)(r*tableScale)<<endl;
			}
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
				if (DEBUG) {
					cout << "G["<<i<<"]"<<"["<<r<<"]"<<"["<<t<<"]"<<(int)(r*tableScale)<<" "<<G[i][(int)(r*tableScale)][t]<<endl;
				}
			}
		}
	}
	//cout<<i<<endl;
	// For i=1 (here index i=0), r=0, t=Dmax
	i=r=0;t=Dmax;
	G[i][0][t]=minH();
	jStar[i][0][t]=minInd;

//	cout << "output:" << G[0][0][Dmax] << endl;


//	int t0=88,t1=24,t2=3,t3=24,t4=32;
//	t=0;
//	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
//	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
//	t+=t0;
//	i++;
//cout<<" ,t: "<<t<<endl;
//	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
//	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
//	t+=t1;
//	i++;
//cout<<" ,t: "<<t<<endl;
//	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
//	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
//	t+=t2;
//	i++;
//cout<<" ,t: "<<t<<endl;
//	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]];
//	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
//	t+=t3;
//	i++;
//cout<<" ,t: "<<t<<endl;
//	cout<<"F"<<i<<", Version: "<<jStar[i][(int)(r*tableScale)][Dmax-t]<<", R["<<i<<"]["<<jStar[i][(int)(r*tableScale)][Dmax-t]<<"]: "<<R[i][jStar[i][(int)(r*tableScale)][Dmax-t]]<<endl;
//	r+=R[i][(jStar[i][(int)(r*tableScale)][Dmax-t])];
//	t+=t4;
//
//	cout<<"Total R= "<<r<<endl;
//	cout<<"Total t= "<<t<<endl;

	return G[0][0][Dmax];
}

//private:
//	float alpha;
//	int n,i,j,r,t,e,sigma,delta,K,D,minInd,Rmax,Dmax; //K:max no. of versions, minInd:the index of the
//	double ***G,***C,***P,**R,*H,*minR,*maxR;
//	int *theta,***jStar;
