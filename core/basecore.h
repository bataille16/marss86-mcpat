
/*
 * MARSSx86 : A Full System Computer-Architecture Simulator
 *
 * This code is released under GPL.
 *
 * Copyright 2011 Avadh Patel <apatel@cs.binghamton.edu>
 *
 */

#ifndef BASE_CORE_H
#define BASE_CORE_H

#include <ptlsim.h>
#include <ptl-qemu.h>
#include <machine.h>
#include <statsBuilder.h>
#include <memoryHierarchy.h>

namespace Core {

    struct BaseCore : public Statable {
        BaseCore(BaseMachine& machine, const char* name);
        virtual ~BaseCore() {}

        virtual void reset() = 0;
        virtual void check_ctx_changes() = 0;
        virtual void flush_tlb(Context& ctx) = 0;
        virtual void flush_tlb_virt(Context& ctx, Waddr virtaddr) = 0;
        virtual void dump_state(ostream& os) = 0;
        virtual void update_stats() = 0;
        virtual void flush_pipeline() = 0;
        virtual W8 get_coreid() = 0;
		virtual void dump_configuration(YAML::Emitter &out) const = 0;

        void update_memory_hierarchy_ptr();

        BaseMachine& machine;
        Memory::MemoryHierarchy* memoryHierarchy;


	void reset_pow_stats(); 
	// stats collected for dynamic power calcualtion
	unsigned long pow_total_uops; // commited + misspred (issued)
	unsigned long pow_committed_uops;
	unsigned long pow_committed_insns; 
	unsigned long pow_total_branches;
	unsigned long pow_branch_mispredicts;
	unsigned long pow_total_loads; 
	unsigned long pow_total_stores ;
	unsigned long pow_cycles;
	unsigned long pow_idle_cycles; 
	unsigned long pow_int_regreads;
	unsigned long pow_fp_regreads;
	unsigned long pow_int_regwrites; 
	unsigned long pow_fp_regwrites; 
	unsigned int pow_alu_accesses; 
	unsigned int pow_fpu_accesses ;
	unsigned int pow_mul_accesses; 
	unsigned long pow_itlb_misses; 
	unsigned long pow_itlb_accesses;
	unsigned long pow_dtlb_misses; 
	unsigned long pow_dtlb_accesses;
	unsigned long pow_btb_lookups;
	unsigned long pow_btb_updates;

	// additional stats for OoO core
	unsigned long pow_fetch_uops; 
	unsigned long pow_rob_reads;
	unsigned long pow_rob_writes; 
	unsigned long pow_dispatch_uops; //allocated uops
	unsigned long pow_ctx_switches; 


    };


};

#endif // BASE_CORE_H
