

#ifndef POWER_H
#define POWER_H

#include "XML_Parse.h"


namespace Core {
    struct BaseCore;
};



namespace Memory {
    struct Controller;
    struct Interconnect;
    struct MemoryHierarchy;
};

// platform specific config
static int  _num_cores = 1;
static int _num_l2_caches = 1; 
static bool _uncore_LLC = false; 
static bool _OoO_cores = true; 
static int _cpu_clock_rate =  3200 ; /* in MHz*/



/* core configurations */
//l1 cache config
static int _l1_line_size = 64;
static int _l1_num_sets = 256; 
static int _l1_assoc = 8 ; 
static int _l1_latency =2; 
static int _l1_bank_width = 64; 
static short  _l1_write_policy = 0; 

//l2 cache config
static int _l2_num_sets = 4096;
static int _l2_line_size = 64;
static int _l2_assoc = 8; 
static int _l2_latency  = 5; 
static int _l2_bank_width = 64; 
static short _l2_write_policy = 1; /* 0 -> write-through  1 -> write_back */



/* LLC configurations */
//static short _llc_write_policy = 1; /* 0 -> write-through  1 -> write_back */
// default LLC is 2 MB (2048 sets, 64 b line size, 16 associativity)
static int _llc_num_sets = 4096;
static int _llc_line_size = 64;
static int _llc_assoc = 8; 
static int _llc_latency  = 5; 
static int _llc_bank_width = 64; 
static short _llc_write_policy = 1; /* 0 -> write-through  1 -> write_back */

static int _uncore_freq = 800 ; /* MHZ */


//static const char * CacheType[4] = {"L1_I_","L1_D_","L2_","L3_"};// hardwire since we cant get the type from the controller class.
							  // Alternative is to add the type to the controller. not necessary


 
void init_power(const char *filename);
void kill_power();
int calc_core_power(Core::BaseCore *core, bool print_power);
void translate_params(system_core *core_params, system_L2 * L2_params); 
void get_OoO_params(system_core *core_params, system_L2* L2_params);// add out of order core parameters
void get_IO_params(system_core *core_params, system_L2* L2_params);// add in order core parameters

//get stats
void translate_core_stats(Core::BaseCore *core, system_core * core_stats, unsigned long sim_cycles);
void translate_L1Cache_stats(Memory::Controller *IL1, Memory::Controller *DL1, system_core *core_stats);
void translate_L2Cache_stats(Memory::Controller *L2,system_L2 *L2_stats);
void translate_UncoreCache_stats(Memory::Controller *LLC, unsigned long sim_cycles, root_system * stats); 

void getcore_stats(int coreid, Core::BaseCore *core, Memory::Controller * IL1, Memory::Controller * DL1, Memory::Controller *L2,unsigned long sim_cycles); 
int test_main();

struct core_power_t
{
	double rt_power; 
	double lk_power;
	double default_vdd; 
};

static struct core_power_t *cores_power; 

#endif // POWER_H
