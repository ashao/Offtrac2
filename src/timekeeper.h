struct timekeeper_t {

	int iteration_counter; 	// Begins at 0, increases each iteration
	int averaging_counter;// Keeps track of how many intervals between averaging
	int climatology_index; 	//Index to read from climatology
	int hindcast_index; 	// Index to read from hindcast

	int total_intervals; 	// How many intervals to integrate over
	int write_intervals; 	// How many intervals to write
	int num_intervals_year; // Number of intervals in a year (12, monthly, 73 for 5-day)

	int read_hind_flag;	// 1 if the current month should read from hindcast or 0 if from climatology

	double current_time;	// Timestamp in year fraction format
	int current_year;
	int current_interval; 	// How far along in a year we are (e.g. 5 is June for monthly intervals, or January 30 for 5-day)

	double accumulated_time_since_writing; // How much time has elapsed since the writing, used for averaging
	double dt;				// Timestep length

	int num_records;	// Number of output records

};

void initialize_timekeeper( void );
void update_timekeeper( void );
