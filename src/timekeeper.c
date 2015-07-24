#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "init.h"
#include "initialize.h"
#include "timekeeper.h"

extern struct timekeeper_t timekeeper;
extern struct parameters run_parameters;
const int days_in_month[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
const double seconds_in_year = (double) 365 * 24 * 60 * 60;
const double seconds_in_day = (double) 24 * 60 * 60;

void initialize_timekeeper( void ) {

	// Note that here the start time is actually one interval before the read in start time.
	// This is because the layer thicknesses from HIM/GOLD/MOM6 are actually end of month snapshots.
	// For example: the start interval for a monthly timestep run is 0 (January). Layer thickness from
	// December will be read in to hstart and on the first integration step, hstart will be copied
	// to hend, and average January transports will be used to integrate from the start to the end of
	// January


	if (!strcmp(run_parameters.timestep,"5day")) timekeeper.num_intervals_year = 73;
	else if (!strcmp(run_parameters.timestep,"month")) timekeeper.num_intervals_year = 12;
	else { printf("Unknown timestep\n"); exit(-1); };
	printf("Timekeeper.num_intervals: %d\n",timekeeper.num_intervals_year);

	// Calculate the total number of intervals to integrate over
	timekeeper.total_intervals = (run_parameters.eyear - run_parameters.syear)*timekeeper.num_intervals_year +
			(run_parameters.einterval-run_parameters.sinterval);
	printf("Timekeeper.total_intervals: %d\n",timekeeper.total_intervals);

	timekeeper.write_intervals = ceil(timekeeper.total_intervals/run_parameters.wrint);
	if (run_parameters.sinterval == 0) {
		run_parameters.syear = run_parameters.syear-1;
	}
	// Check to see if we should begin by reading a hindcast year
	if (run_parameters.use_hindcast) {

		timekeeper.read_hind_flag = run_parameters.syear >= BEGHIND && run_parameters.syear <= ENDHIND;
		if (timekeeper.read_hind_flag)
			timekeeper.current_interval = mod(run_parameters.sinterval - 1,timekeeper.num_intervals_year);

	}
	// Always keep track of the climatology index
	timekeeper.climatology_index = mod(run_parameters.sinterval - 1,run_parameters.ntime_climatology);


	// Calculate starting year fraction
	// Note that 1 is added to the interval since the timestamp more rightfully refers to the end of the timestep
	timekeeper.current_time = run_parameters.syear +
			(double) (run_parameters.sinterval+1)/timekeeper.num_intervals_year;

	timekeeper.current_interval = mod(run_parameters.sinterval-1,timekeeper.num_intervals_year);
	timekeeper.current_year = run_parameters.syear;
	

		
}

void update_timekeeper( void ) {

	struct timespec wallclock;
	timekeeper.current_interval++;
	timekeeper.averaging_counter++;
	timekeeper.current_time = timekeeper.current_year + (double) (timekeeper.current_interval+1)/timekeeper.num_intervals_year;

	if ( (timekeeper.current_interval == timekeeper.num_intervals_year) )	{
		timekeeper.current_year++;
		timekeeper.current_interval = 0;
	}

	if (run_parameters.use_hindcast) {
		timekeeper.read_hind_flag = run_parameters.syear >= BEGHIND && run_parameters.syear <= ENDHIND;
	}

	// The climatology index should always increment so that it is always on the correct time interval
	timekeeper.climatology_index++;
	timekeeper.climatology_index = timekeeper.climatology_index % run_parameters.ntime_climatology;


	if (timekeeper.num_intervals_year==73) { timekeeper.dt = 5.0 * seconds_in_day; }
	else if (timekeeper.num_intervals_year==12) { timekeeper.dt = (double) days_in_month[timekeeper.current_interval] * seconds_in_day; }
	else {
		printf("Invalid number of time intervals in year, choose either 5-day or monthly\n");
		exit(-1);
	}
	timekeeper.accumulated_time_since_writing+=timekeeper.dt;


}

