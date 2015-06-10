/*
 * ReliabilitySim.cpp
 *
 *  Created on: Sep 13, 2012
 *      Author: anas
 */
#include "DynamicAlgo.h"
#include <time.h>
#include <sys/time.h>
#include <mach/clock.h>
#include <mach/mach.h>

long dtime(timespec start, timespec end)
{
	long sec, nsec;
	sec = end.tv_sec  - start.tv_sec;
 	nsec = end.tv_nsec - start.tv_nsec;
    	//return ((sec) * 1000 + nsec/1000000); // time in milliseconds
// 	cout<<sec<<endl;
// 	cout<<nsec<<endl;
 	return ((sec) * 1000000 + nsec/1000); // time in microseconds
}

int proc_notime(int start_time, int during){
	double best_result = 2;
	double worst_result = 0;
	ofstream fileout;
	string outputfile;
	string ofs;
	string ifs;
	int RScale = 100000;
	int index = 0;
	int w_size[5];
	int i=0, n=5, K=24;
	int record_func[120][n];
	int record_order[n];
	int wrecord_order[n];
	double **AVG = new double*[n];
		for(i=0;i<n;i++)
			AVG[i]=new double[K];
	double **R = new double*[n];
		for(i=0;i<n;i++)
			R[i]=new double[K];

	for(i=0; i<n; i++) {
		w_size[i] = 0;
		record_order[i] = -1;
	}
	//w_size[0]=4;w_size[1]=24;w_size[2]=1;w_size[3]=2;w_size[4]=5;

	for(i=0; i<120; i++){ //times
		for(int j=0; j<n; j++){ //function order
			record_func[i][j] = 0;
		}
	}

	// Reads R(i,j) from Rij.txt file
	ifstream in("GetSure_CDF_Results/Rij.txt");
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
		R[iInd][jInd]=(rValue);
		AVG[iInd][jInd]=avgVAlue;
		w_size[iInd]+=1;
	}
	in.close();

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			for(int k=0; k<n; k++){
				for(int l=0; l<n; l++){
					for(int m=0; m<n; m++){
						if(i == j || i == k || i == l || i == m) continue;
						if(j == k || j == l || j == m) continue;
						if(k == l || k == m) continue;
						if(l == m ) continue;
						record_func[index][0] = i;
						record_func[index][1] = j;
						record_func[index][2] = k;
						record_func[index][3] = l;
						record_func[index][4] = m;
						//cout << record_func[index][0] << record_func[index][1] << record_func[index][2] << record_func[index][3] << record_func[index][4] <<endl;
						index++;
					}
				}
			}
		}
	}
	for (int i = start_time; i <= 300; i+=during) {
		worst_result = 0;
		best_result = 2;
				DynamicAlgo DAImp(0.5,1000,RScale,i,1,1); //6
				/* For example: 0.1, 10, 10^5 or 0.1, 100, 10^4
				 * it means 0.1 is alpha, alpha scale=10, RScale = 10^5 -> 0.1 alpha, 10^-6
				 * */
#if 1
				for(int z = 0; z<120; z++){
#else
				for(int z = 0; z<1; z++){// origin version, new order
					record_func[z][0]=4;
					record_func[z][1]=0;
					record_func[z][2]=3;
					record_func[z][3]=2;
					record_func[z][4]=1;
#endif
					/*Selection ordering*/
//					cout << "AAA" <<endl;
					outputfile="estimedia_input/Rij_test";
					//outputfile="GetSure_CDF_Results/Rij_aaa";
					outputfile+=".txt";
					fileout.open(outputfile.c_str());
					int cur_index, exist_num=0, nk = 0;
					int eSample;
					double cValue;
					DAImp.special_fun0 = record_func[z][0];
					DAImp.special_fun = record_func[z][1];
					DAImp.special_fun2 = record_func[z][2];
					DAImp.special_fun3 = record_func[z][3];
					DAImp.special_fun4 = record_func[z][4];

					for(int j=0;j<n;j++){// function index
//						cout << "BBB" <<endl;
					if(j == 1){
						for(int w=0; w<w_size[1]; w++){
							fileout<<record_func[z][j]<<" "<<w<<" ";
							fileout<<R[1][w]<<" ";
							fileout<<AVG[1][w];
//							cout<<record_func[z][j]<<"C ";
//							cout<<R[1][w]<<" ";
//							cout<<AVG[1][w]<<endl;
							fileout<<"\n";
						}
					}
					else{
						for(int w=0; w<w_size[j]; w++){
							fileout<<record_func[z][j]<<" "<<w<<" ";
							fileout<<R[j][w]<<" ";
							fileout<<AVG[j][w];
//							cout<<record_func[z][j]<<"C ";
//							cout<<R[j][w]<<" ";
//							cout<<AVG[j][w]<<endl;
							fileout<<"\n";
							ifs="GetSure_CDF_Results/";
							ifs+=('0'+j);
							ifs+="/";
							if((w/10)!=0) ifs+=((48+w/10));
							ifs+=((48+w%10));
							ifs+=".txt";
//							cout<< ifs <<endl;
							ofs="GetSure_CDF_hand/";  //0 function
							ofs+=('0'+record_func[z][j]);
							ofs+="/";
							if((w/10)!=0) ofs+=((48+w/10));
							ofs+=((48+w%10));
							ofs+=".txt";
//							cout<< ofs <<endl;
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

				}//5
				fileout.close();
				DAImp.getData();
				double result=DAImp.fillTable();
//					cout<<result<<endl;
				if(result>worst_result) {
					worst_result = result;
					wrecord_order[0] = record_func[z][0];
					wrecord_order[1] = record_func[z][1];
					wrecord_order[2] = record_func[z][2];
					wrecord_order[3] = record_func[z][3];
					wrecord_order[4] = record_func[z][4];
				}
				if(result<best_result) {
					best_result = result;
					record_order[0] = record_func[z][0];
					record_order[1] = record_func[z][1];
					record_order[2] = record_func[z][2];
					record_order[3] = record_func[z][3];
					record_order[4] = record_func[z][4];
				}
			}//120

//			cout << i<<endl;
//			cout << best_result <<endl;
//			cout << record_order[0] << record_order[1] << record_order[2] <<record_order[3] << record_order[4] <<endl;
//
			cout << worst_result <<endl;
//			cout << wrecord_order[0] << wrecord_order[1] << wrecord_order[2] <<wrecord_order[3] << wrecord_order[4] <<endl;
				}
			return 0;
}


int proc_time(int start_time){
/*
	DynamicAlgo DAImp(0.005,120,1,1);
	DAImp.getData();
	cout << DAImp.fillTable()<<endl;
*/
	double best_result = 2;
	double worst_result = 0;
	ofstream fileout;
	string outputfile;
	string ofs;
	string ifs;
	int RScale = 100000;
	int index = 0;
	int w_size[5];
	int i=0, n=5, K=24;
	int record_func[120][n];
	int record_order[n];
	int wrecord_order[n];
	double **AVG = new double*[n];
		for(i=0;i<n;i++)
			AVG[i]=new double[K];
	double **R = new double*[n];
		for(i=0;i<n;i++)
			R[i]=new double[K];

	for(i=0; i<n; i++) {
		w_size[i] = 0;
		record_order[i] = -1;
	}
	//w_size[0]=4;w_size[1]=24;w_size[2]=1;w_size[3]=2;w_size[4]=5;

	for(i=0; i<120; i++){ //times
		for(int j=0; j<n; j++){ //function order
			record_func[i][j] = 0;
		}
	}

	// Reads R(i,j) from Rij.txt file
	ifstream in("GetSure_CDF_Results/Rij.txt");
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
		R[iInd][jInd]=(rValue);
		AVG[iInd][jInd]=avgVAlue;
		w_size[iInd]+=1;
	}
	in.close();

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			for(int k=0; k<n; k++){
				for(int l=0; l<n; l++){
					for(int m=0; m<n; m++){
						if(i == j || i == k || i == l || i == m) continue;
						if(j == k || j == l || j == m) continue;
						if(k == l || k == m) continue;
						if(l == m ) continue;
						record_func[index][0] = i;
						record_func[index][1] = j;
						record_func[index][2] = k;
						record_func[index][3] = l;
						record_func[index][4] = m;
						//cout << record_func[index][0] << record_func[index][1] << record_func[index][2] << record_func[index][3] << record_func[index][4] <<endl;
						index++;
					}
				}
			}
		}
	}

//	for (int i = 160; i <= 300; i+=10) {
//		worst_result = 0;
//		best_result = 2;
//				DynamicAlgo DAImp(0.1,10,RScale,i,1,1); //6
//				/* For example: 0.1, 10, 10^5 or 0.1, 100, 10^4
//				 * it means 0.1 is alpha, alpha scale=10, RScale = 10^5 -> 0.1 alpha, 10^-6
//				 * */
//#if 1
//				for(int z = 0; z<120; z++){
//#else
//				for(int z = 0; z<1; z++){// origin version, new order
//					record_func[z][0]=0;
//					record_func[z][1]=1;
//					record_func[z][2]=2;
//					record_func[z][3]=3;
//					record_func[z][4]=4;
//#endif
//				    /*Selection ordering*/
////					cout << "AAA" <<endl;
//					outputfile="estimedia_input/Rij_test";
//					//outputfile="GetSure_CDF_Results/Rij_aaa";
//					outputfile+=".txt";
//					fileout.open(outputfile.c_str());
//					int cur_index, exist_num=0, nk = 0;
//					int eSample;
//					double cValue;
//					DAImp.special_fun0 = record_func[z][0];
//					DAImp.special_fun = record_func[z][1];
//					DAImp.special_fun2 = record_func[z][2];
//					DAImp.special_fun3 = record_func[z][3];
//					DAImp.special_fun4 = record_func[z][4];
//
//					for(int j=0;j<n;j++){// function index
////						cout << "BBB" <<endl;
//					if(j == 1){
//						for(int w=0; w<w_size[1]; w++){
//							fileout<<record_func[z][j]<<" "<<w<<" ";
//							fileout<<R[1][w]<<" ";
//							fileout<<AVG[1][w];
////							cout<<record_func[z][j]<<"C ";
////							cout<<R[1][w]<<" ";
////							cout<<AVG[1][w]<<endl;
//							fileout<<"\n";
//						}
//					}
//					else{
//						for(int w=0; w<w_size[j]; w++){
//							fileout<<record_func[z][j]<<" "<<w<<" ";
//							fileout<<R[j][w]<<" ";
//							fileout<<AVG[j][w];
////							cout<<record_func[z][j]<<"C ";
////							cout<<R[j][w]<<" ";
////							cout<<AVG[j][w]<<endl;
//							fileout<<"\n";
//							ifs="GetSure_CDF_Results/";
//							ifs+=('0'+j);
//							ifs+="/";
//							if((w/10)!=0) ifs+=((48+w/10));
//							ifs+=((48+w%10));
//							ifs+=".txt";
////							cout<< ifs <<endl;
//							ofs="GetSure_CDF_hand/";  //0 function
//							ofs+=('0'+record_func[z][j]);
//							ofs+="/";
//							if((w/10)!=0) ofs+=((48+w/10));
//							ofs+=((48+w%10));
//							ofs+=".txt";
////							cout<< ofs <<endl;
//							ifstream file1;
//							ofstream file2;
//							file1.open(ifs.c_str(), ios::in |ios::binary);
//							file2.open(ofs.c_str(), ios::out | ios::trunc |ios::binary);
//							for(;;){
//								file1>>eSample>>cValue;
//				//				cout << eSample << " "<< cValue <<endl;
//								if(file1.eof()) break;
//								file2<<eSample<<" "<<cValue;
//								file2<<"\n";
//							}
//							file1.close();
//							file2.close();
//						}
//					}
//
//				}//5
//				fileout.close();
//				DAImp.getData();
//				double result=DAImp.fillTable();
////					cout<<result<<endl;
//				if(result>worst_result) {
//					worst_result = result;
//					wrecord_order[0] = record_func[z][0];
//					wrecord_order[1] = record_func[z][1];
//					wrecord_order[2] = record_func[z][2];
//					wrecord_order[3] = record_func[z][3];
//					wrecord_order[4] = record_func[z][4];
//				}
//				if(result<best_result) {
//					best_result = result;
//					record_order[0] = record_func[z][0];
//					record_order[1] = record_func[z][1];
//					record_order[2] = record_func[z][2];
//					record_order[3] = record_func[z][3];
//					record_order[4] = record_func[z][4];
//				}
//			}//120
//
////			cout << i<<endl;
//			cout << best_result <<endl;
////			cout << record_order[0] << record_order[1] << record_order[2] <<record_order[3] << record_order[4] <<endl;
////
////			cout << worst_result <<endl;
////			cout << wrecord_order[0] << wrecord_order[1] << wrecord_order[2] <<wrecord_order[3] << wrecord_order[4] <<endl;
//
//}


	//---for mac
	// Timer code
	timespec start, end;
	double sum=0;
	clock_serv_t cclock;
	mach_timespec_t mts;

	for (int i = start_time; i <= 300; i+=10) {
		for(int j=0;j<10;j++){

			//---
			host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
			clock_get_time(cclock, &mts);
			mach_port_deallocate(mach_task_self(), cclock);
			start.tv_sec = mts.tv_sec;
			start.tv_nsec = mts.tv_nsec;

			worst_result = 0;
				best_result = 2;
				DynamicAlgo DAImp(0.01,1000,RScale,i,1,1); //6
				/* For example: 0.1, 10, 10^5 or 0.1, 100, 10^4
				 * it means 0.1 is alpha, alpha scale=10, RScale = 10^5 -> 0.1 alpha, 10^-6
				 * */
#if 1
				for(int z = 0; z<120; z++){
#else
				for(int z = 0; z<1; z++){// origin version, new order
					record_func[z][0]=0;
					record_func[z][1]=1;
					record_func[z][2]=2;
					record_func[z][3]=3;
					record_func[z][4]=4;
#endif
					/*Selection ordering*/
//					cout << "AAA" <<endl;
					outputfile="estimedia_input/Rij_test";
					//outputfile="GetSure_CDF_Results/Rij_aaa";
					outputfile+=".txt";
					fileout.open(outputfile.c_str());
					int cur_index, exist_num=0, nk = 0;
					int eSample;
					double cValue;
					DAImp.special_fun0 = record_func[z][0];
					DAImp.special_fun = record_func[z][1];
					DAImp.special_fun2 = record_func[z][2];
					DAImp.special_fun3 = record_func[z][3];
					DAImp.special_fun4 = record_func[z][4];

					for(int j=0;j<n;j++){// function index
//						cout << "BBB" <<endl;
					if(j == 1){
						for(int w=0; w<w_size[1]; w++){
							fileout<<record_func[z][j]<<" "<<w<<" ";
							fileout<<R[1][w]<<" ";
							fileout<<AVG[1][w];
//							cout<<record_func[z][j]<<"C ";
//							cout<<R[1][w]<<" ";
//							cout<<AVG[1][w]<<endl;
							fileout<<"\n";
						}
					}
					else{
						for(int w=0; w<w_size[j]; w++){
							fileout<<record_func[z][j]<<" "<<w<<" ";
							fileout<<R[j][w]<<" ";
							fileout<<AVG[j][w];
//							cout<<record_func[z][j]<<"C ";
//							cout<<R[j][w]<<" ";
//							cout<<AVG[j][w]<<endl;
							fileout<<"\n";
							ifs="GetSure_CDF_Results/";
							ifs+=('0'+j);
							ifs+="/";
							if((w/10)!=0) ifs+=((48+w/10));
							ifs+=((48+w%10));
							ifs+=".txt";
//							cout<< ifs <<endl;
							ofs="GetSure_CDF_hand/";  //0 function
							ofs+=('0'+record_func[z][j]);
							ofs+="/";
							if((w/10)!=0) ofs+=((48+w/10));
							ofs+=((48+w%10));
							ofs+=".txt";
//							cout<< ofs <<endl;
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

				}//5
				fileout.close();
				DAImp.getData();
				double result=DAImp.fillTable();
//					cout<<result<<endl;
				if(result>worst_result) {
					worst_result = result;
					wrecord_order[0] = record_func[z][0];
					wrecord_order[1] = record_func[z][1];
					wrecord_order[2] = record_func[z][2];
					wrecord_order[3] = record_func[z][3];
					wrecord_order[4] = record_func[z][4];
				}
				if(result<best_result) {
					best_result = result;
					record_order[0] = record_func[z][0];
					record_order[1] = record_func[z][1];
					record_order[2] = record_func[z][2];
					record_order[3] = record_func[z][3];
					record_order[4] = record_func[z][4];
				}
			}//120

//			cout << i<<endl;
//			cout << best_result <<endl;
//			cout << record_order[0] << record_order[1] << record_order[2] <<record_order[3] << record_order[4] <<endl;
//
			cout << worst_result <<endl;
//			cout << wrecord_order[0] << wrecord_order[1] << wrecord_order[2] <<wrecord_order[3] << wrecord_order[4] <<endl;

			//---
			host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
			clock_get_time(cclock, &mts);
			mach_port_deallocate(mach_task_self(), cclock);
//			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
			end.tv_sec = mts.tv_sec;

			end.tv_nsec = mts.tv_nsec;
			sum+=dtime(start,end);
//			cout<<sum<<endl;
		}
		cout<<i<<" "<<(sum/10)<<endl;
		sum=0;
		}
		return 0;
}
int main(void){
		proc_notime(200,5);
//		proc_time(230);
		return 0;
}



