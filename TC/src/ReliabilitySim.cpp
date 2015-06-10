/*
 * ReliabilitySim.cpp
 *
 *  Created on: Sep 13, 2012
 *      Author: anas
 *      Co-Author: Kuan-Hsun
 */
#include "DynamicAlgo.h"
#include <time.h>
#include <sys/time.h>
#include <mach/clock.h>
#include <mach/mach.h>
#include <unistd.h>
#include <sys/resource.h>

long dtime(timespec start, timespec end)
{
	long sec, nsec;
	sec = end.tv_sec  - start.tv_sec;
 	nsec = end.tv_nsec - start.tv_nsec;
    	//return ((sec) * 1000 + nsec/1000000); // time in milliseconds
 	return ((sec) * 1000000 + nsec/1000); // time in microseconds
}
size_t getPeakRSS()
{
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
	return (size_t) rusage.ru_maxrss;
}
size_t getCurrentRSS()
{
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if(task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) return (size_t)0L;
	return (size_t)info.resident_size;
}

int main(void){
/*
	DynamicAlgo DAImp(0.005,120,1,1);
	DAImp.getData();
	cout << DAImp.fillTable()<<endl;
*/
#if 1
// i = relative deadline

	for (int i =0; i <= 300; i+=20) {
//		for (int i =108; i > 107; i-=1) {

				DynamicAlgo DAImp(0.5,1000,1000,i,1,1); //6~8
				DAImp.getData(0);
				//cout<<"Origin "<<i<<" "<<DAImp.fillTable()<<endl;
				cout<<DAImp.fillTable()<<endl;
//<---- DATE ----->//
//				DAImp.fillOrder(); //put in subitem set and print out in new file
//				DAImp.getData(1); //check the new file with one ordering
//				cout<<DAImp.fillTable()<<endl;

				cout<<"PEAKRSS "<<8*getPeakRSS()/1000/1000<<" MB"<<endl;
				cout<<sizeof(size_t)<<" Bytes"<<endl;
				cout<<sizeof(int)<<" Bytes"<<endl;
//				cout<<"CURRENT RSS "<<getCurrentRSS()<<endl;
	}
#else


	//--- for mac
		// Timer code
		timespec start, end;
		double sum=0;
		clock_serv_t cclock;
		mach_timespec_t mts;

		for (int i = 0; i <= 300; i+=20) {
			for(int j=0;j<10;j++){

				//---
	//			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);

				host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
				clock_get_time(cclock, &mts);
				mach_port_deallocate(mach_task_self(), cclock);
				start.tv_sec = mts.tv_sec;
				start.tv_nsec = mts.tv_nsec;

				DynamicAlgo DAImp(0.01,1000,10000,i,1,1); //6
				DAImp.getData(0);

				DAImp.fillOrder(); //put in subitem set and print out in new file
				DAImp.getData(1); //check the new file with one ordering
				DAImp.fillTable();
//				cout<<DAImp.fillTable()<<endl;
				//---
				host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
				clock_get_time(cclock, &mts);
				mach_port_deallocate(mach_task_self(), cclock);
	//			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
				end.tv_sec = mts.tv_sec;
				end.tv_nsec = mts.tv_nsec;
				sum+=dtime(start,end);

			}
			cout<<i<<" "<<(sum/10)<<endl;
			sum=0;
		}

#endif
}



