
/*
 * MARSSx86 : A Full System Computer-Architecture Simulator
 *
 * This code is released under GPL.
 *
 * Copyright 2011 Avadh Patel <apatel@cs.binghamton.edu>
 *
 */

#include <basecore.h>
#include <globals.h>
#include <decode.h>

using namespace Core;

BaseCore::BaseCore(BaseMachine& machine, const char* name)
    : Statable(name, &machine)
      , machine(machine)
{

	 pow_total_uops = 0 ; 
	 pow_committed_uops = 0 ; 
 	 pow_committed_insns = 0 ; 
	 pow_total_branches = 0;
	 pow_branch_mispredicts = 0 ; 
	 pow_total_loads = 0 ; 
	 pow_total_stores  = 0 ; 
	 pow_cycles = 0;
	 pow_idle_cycles = 0 ;
	 pow_int_regreads= 0; 
	 pow_fp_regreads = 0 ; 
	 pow_int_regwrites = 0 ; 
	 pow_fp_regwrites =0 ; 
	 pow_alu_accesses = 0;  
	 pow_fpu_accesses  =0; 
	 pow_mul_accesses = 0 ; 
	 pow_itlb_misses = 0; 
	 pow_itlb_accesses  = 0; 
	 pow_dtlb_misses = 0; 
	 pow_dtlb_accesses = 0 ; 
	 pow_btb_lookups = 0;
	 pow_btb_updates =0;
	 // additional stats for OoO core
 	pow_fetch_uops = 0 ; 
	pow_rob_reads = 0 ; 
 	pow_rob_writes = 0 ; 
 	pow_dispatch_uops = 0 ; //allocated uops
 	pow_ctx_switches = 0; 

}
void BaseCore::reset_pow_stats() {

	 pow_total_uops = 0 ; 
	 pow_committed_uops = 0 ; 
 	 pow_committed_insns = 0 ; 
	 pow_total_branches = 0;
	 pow_branch_mispredicts = 0 ; 
	 pow_total_loads = 0 ; 
	 pow_total_stores  = 0 ; 
	 pow_cycles = 0;
	 pow_idle_cycles = 0 ;
	 pow_int_regreads= 0; 
	 pow_fp_regreads = 0 ; 
	 pow_int_regwrites = 0 ; 
	 pow_fp_regwrites =0 ; 
	 pow_alu_accesses = 0;  
	 pow_fpu_accesses  =0; 
	 pow_mul_accesses = 0 ; 
	 pow_itlb_misses = 0; 
	 pow_itlb_accesses  = 0; 
	 pow_dtlb_misses = 0; 
	 pow_dtlb_accesses = 0 ; 
	 pow_btb_lookups = 0;
	 pow_btb_updates =0;
	 // additional stats for OoO core
 	pow_fetch_uops = 0 ; 
	pow_rob_reads = 0 ; 
 	pow_rob_writes = 0 ; 
 	pow_dispatch_uops = 0 ; //allocated uops
 	pow_ctx_switches = 0; 

}

void BaseCore::update_memory_hierarchy_ptr() {
    memoryHierarchy = machine.memoryHierarchyPtr;
}


extern "C" void ptl_flush_bbcache(int8_t context_id) {
    if(in_simulation) {
      foreach(i, NUM_SIM_CORES) {
        bbcache[i].flush(context_id);
        // Get the current ptlsim machine and call its flush tlb
        PTLsimMachine* machine = PTLsimMachine::getcurrent();

        if(machine) {
            Context& ctx = machine->contextof(context_id);
            machine->flush_tlb(ctx);
        }
      }
    }
}
