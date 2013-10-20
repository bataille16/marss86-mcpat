#include "power.h"
#include <stdio.h>
#include <basecore.h>

#include <ooo-const.h>
//#include "XML_Parse.h"
#include "mcpat.h"

using namespace Core;

class ParseXML *XML = NULL;

bool private_l2 = false; 

double *cores_leakage;
double uncore_leakage; 

double *cores_rtp;
double uncore_rtp;

FILE *pow_out_trace;


//int testing = OOO_CORE_MODEL::FETCH_WIDTH;



void init_power(const char * filename)
{


  XML = new ParseXML();
  XML->initialize();

  /* Translate uncore params */
  XML->sys.number_of_cores = _num_cores;
  XML->sys.number_of_L1Directories = 0;
  XML->sys.number_of_L2Directories = 0;
  XML->sys.number_of_NoCs = 0;
  XML->sys.number_of_dir_levels = 0;
  XML->sys.homogeneous_cores = 0;
  XML->sys.homogeneous_L1Directories = 0;
  XML->sys.homogeneous_L2Directories = 0;
  XML->sys.core_tech_node = 45;
  XML->sys.target_core_clockrate = _cpu_clock_rate;
  XML->sys.temperature = 380; // K
  XML->sys.interconnect_projection_type = 0; // aggressive
  XML->sys.longer_channel_device = 1; // use when appropirate
  XML->sys.machine_bits = 64;
  XML->sys.virtual_address_width = 64;
  XML->sys.physical_address_width = 52;
  XML->sys.virtual_memory_page_size = 4096;

  int num_l2 = 1;
  bool has_hp_core = false;
    // If any core has a private L2 and we have more than one core, we assume L3 is LLC
  if(_num_l2_caches == _num_cores && _num_cores > 1)
  {
      private_l2 = true;
      num_l2 = _num_cores;
      _uncore_LLC = true;    
  }

  if (_OoO_cores)
     has_hp_core = true;
   
  XML->sys.device_type = has_hp_core ? 0 /*HP*/: 2 /*LOP*/;

  XML->sys.Private_L2 = private_l2;
  XML->sys.number_of_L2s = private_l2 ? num_l2 : 1;
  XML->sys.number_cache_levels = private_l2 ? 3 : 2;
  XML->sys.number_of_L3s = private_l2 ? 1 : 0;
  

  if (_uncore_LLC && ! private_l2)
 { 


    XML->sys.L2[0].L2_config[0] = _llc_num_sets * _llc_assoc * _llc_line_size;
    XML->sys.L2[0].L2_config[1] = _llc_line_size;
    XML->sys.L2[0].L2_config[2] = _llc_assoc;
    XML->sys.L2[0].L2_config[3] = 1; // single bank  
    XML->sys.L2[0].L2_config[4] = 1;
    XML->sys.L2[0].L2_config[5] = _llc_latency;
    XML->sys.L2[0].L2_config[6] = _llc_bank_width;
    XML->sys.L2[0].L2_config[7] = _llc_write_policy; 
    XML->sys.L2[0].ports[0] = 1;
    XML->sys.L2[0].ports[1] = 0;
    XML->sys.L2[0].ports[2] = 0;

    XML->sys.L2[0].buffer_sizes[0] = 1;
    XML->sys.L2[0].buffer_sizes[1] = 2;
    XML->sys.L2[0].buffer_sizes[2] = 2;
    XML->sys.L2[0].buffer_sizes[3] = 2;

    XML->sys.L2[0].clockrate = _uncore_freq; ;
    XML->sys.L2[0].device_type = 2;
  } else if (_uncore_LLC) // LLC is L3
  {
    XML->sys.L3[0].L3_config[0] = _llc_num_sets * _llc_assoc * _llc_line_size;
    XML->sys.L3[0].L3_config[1] = _llc_line_size;
    XML->sys.L3[0].L3_config[2] = _llc_assoc;
    XML->sys.L3[0].L3_config[3] = 1; // single bank  
    XML->sys.L3[0].L3_config[4] = 1;
    XML->sys.L3[0].L3_config[5] = _llc_latency;
    XML->sys.L3[0].L3_config[6] = _llc_bank_width;
    XML->sys.L3[0].L3_config[7] = _llc_write_policy; 

    XML->sys.L3[0].ports[0] = 1;
    XML->sys.L3[0].ports[1] = 0;
    XML->sys.L3[0].ports[2] = 0;

    XML->sys.L3[0].buffer_sizes[0] = 1;
    XML->sys.L3[0].buffer_sizes[1] = 2;
    XML->sys.L3[0].buffer_sizes[2] = 2;
    XML->sys.L3[0].buffer_sizes[3] = 2;

    XML->sys.L3[0].clockrate = (int)_uncore_freq;
    XML->sys.L3[0].device_type = 0;
  }

  XML->sys.mc.number_mcs = 0;
  XML->sys.flashc.number_mcs = 0;
  XML->sys.niu.number_units = 0;
  XML->sys.pcie.number_units = 0;


 // initiate the core power structure for all cores
  cores_power = (core_power_t*)calloc(_num_cores, sizeof(core_power_t));
  if (cores_power == NULL)
	exit(-1); 

  for (int i=0; i<_num_cores; i++)
  	translate_params(&XML->sys.core[i], &XML->sys.L2[i]);

  for (int i=0; i<_num_cores; i++)
  	get_OoO_params(&XML->sys.core[i], &XML->sys.L2[i]);


  cores_leakage = (double*)calloc(_num_cores, sizeof(*cores_leakage));
  if (cores_leakage == NULL)
    exit(-1);

  mcpat_initialize(XML, &std::cerr, cores_leakage, &uncore_leakage, 5);


  cores_rtp = (double*)calloc(_num_cores, sizeof(*cores_rtp));
  if (cores_rtp == NULL)
    exit(-1); ;


  if (filename !=NULL)
  {  
      pow_out_trace = fopen(filename, "w"); 
      if (pow_out_trace == NULL)
        exit(-1); 
  }

/*
 // Initialize Vdd from mcpat defaults
  core_power_t::default_vdd = g_tp.peri_global.Vdd;
  for (int i=0; i<num_cores; i++)
    if (cores[i]->vf_controller) {
      cores[i]->vf_controller->vdd = g_tp.peri_global.Vdd;
      cores[i]->vf_controller->vf_controller_t::change_vf();
    }

*/
	
}

void translate_params(system_core *core_params, system_L2 * L2_params)
{

  (void) L2_params;
  core_params->clock_rate = _cpu_clock_rate; 
  core_params->opt_local = false;
  core_params->x86 = true;
  core_params->machine_bits = 64;
  core_params->virtual_address_width = 64;
  core_params->physical_address_width = 52; //XXX
  core_params->opcode_width = 16;
  core_params->micro_opcode_width = 8;
  core_params->instruction_length = 32;

  core_params->fetch_width =OOO_CORE_MODEL::FETCH_WIDTH;

  core_params->decode_width = OOO_CORE_MODEL::FRONTEND_WIDTH;
  core_params->issue_width = OOO_CORE_MODEL::MAX_ISSUE_WIDTH;
  core_params->peak_issue_width = OOO_CORE_MODEL::MAX_ISSUE_WIDTH;
  core_params->commit_width = OOO_CORE_MODEL::COMMIT_WIDTH; 

  core_params->ALU_per_core = OOO_CORE_MODEL::ALU_FU_COUNT;
  core_params->MUL_per_core = OOO_CORE_MODEL::ALU_FU_COUNT;
  core_params->FPU_per_core = OOO_CORE_MODEL::FPU_FU_COUNT;
  core_params->instruction_buffer_size = OOO_CORE_MODEL::FETCH_QUEUE_SIZE;  
  core_params->decoded_stream_buffer_size = OOO_CORE_MODEL::ISSUE_QUEUE_SIZE;  

  core_params->ROB_size = OOO_CORE_MODEL::ROB_SIZE;  
  core_params->load_buffer_size = OOO_CORE_MODEL::LDQ_SIZE; 
  core_params->store_buffer_size = OOO_CORE_MODEL::STQ_SIZE; 
  core_params->RAS_size =  1024; //default RAAS size obtained from ptlsim code review

//if (dl1)
   core_params->dcache.dcache_config[0] =  _l1_num_sets * _l1_assoc * _l1_line_size;
    core_params->dcache.dcache_config[1] = _l1_line_size;
    core_params->dcache.dcache_config[2] = _l1_assoc;
    // Hardcode banks to 1, McPAT adds big overhead for multibanked caches
    core_params->dcache.dcache_config[3] = 1;
    core_params->dcache.dcache_config[4] = 1;
    core_params->dcache.dcache_config[5] = _l1_latency;
    core_params->dcache.dcache_config[6] = _l1_bank_width;
    core_params->dcache.dcache_config[7] = _l1_write_policy;  



    core_params->dcache.buffer_sizes[0] = 4;//core->memory.DL1->MSHR_size;
    core_params->dcache.buffer_sizes[1] = 4;//core->memory.DL1->fill_num[0]; //XXX
    core_params->dcache.buffer_sizes[2] = 4;//core->memory.DL1->PFF_size;
    core_params->dcache.buffer_sizes[3] = 4;//core->memory.DL1->WBB_size;
 

//  if (il1)
  
    core_params->icache.icache_config[0] = _l1_num_sets * _l1_assoc * _l1_line_size;
    core_params->icache.icache_config[1] = _l1_line_size;
    core_params->icache.icache_config[2] = _l1_assoc;           
    // Hardcode banks to 1, McPAT adds big overhead for multiba 
    core_params->icache.icache_config[3] = 1;                   
    core_params->icache.icache_config[4] = 1;                   
    core_params->icache.icache_config[5] = _l1_latency;         
    core_params->icache.icache_config[6] = _l1_bank_width;      
    core_params->icache.icache_config[7] = _l1_write_policy;  

    core_params->icache.buffer_sizes[0] = 2;//core->memory.IL1->MSHR_size;
    core_params->icache.buffer_sizes[1] = 2;//core->memory.IL1->fill_num[0]; //XXX
    core_params->icache.buffer_sizes[2] = 2;//core->memory.IL1->PFF_size;
    core_params->icache.buffer_sizes[3] = 2;//core->memory.IL1->WBB_size;
  

 // if (itlb)
  
    core_params->itlb.number_entries = OOO_CORE_MODEL::ITLB_SIZE;
    core_params->itlb.cache_policy = 1;
  

  //if (dtlb)
  
    core_params->dtlb.number_entries = OOO_CORE_MODEL::DTLB_SIZE;;
    core_params->dtlb.cache_policy = 1;
  

//if (bpred && btb)
    core_params->BTB.BTB_config[0] = 4096; // BTB entries obtained from ptlsim code review 
    core_params->BTB.BTB_config[1] = 8 / 8;// harcoded from ptlsim code review
    core_params->BTB.BTB_config[2] = 4;// hardcoded from pltsim code review
    core_params->BTB.BTB_config[3] = 1; //# banks
    core_params->BTB.BTB_config[4] = 1; //troughput
    core_params->BTB.BTB_config[5] = 1; //latency
   

    // standard 2-level predictor 
    core_params->predictor.local_predictor_entries  = 4096;
    core_params->predictor.local_predictor_size[0]  = 20;
    core_params->predictor.local_predictor_size[1]  = 0;
    core_params->predictor.global_predictor_entries  = 1;
    core_params->predictor.global_predictor_bits  = 12;
    core_params->predictor.chooser_predictor_entries = 0; 
    core_params->predictor.chooser_predictor_bits  =  0 ;
  


  // AF for max power computation
  core_params->IFU_duty_cycle = 1.0;
  core_params->LSU_duty_cycle = 0.5;
  core_params->MemManU_I_duty_cycle = 1.0;
  core_params->MemManU_D_duty_cycle = 0.5;
  core_params->ALU_duty_cycle = 1.0;
  core_params->MUL_duty_cycle = 0.3;
  core_params->FPU_duty_cycle = 0.3;
  core_params->ALU_cdb_duty_cycle = 1.0;
  core_params->MUL_cdb_duty_cycle = 0.3;
  core_params->FPU_cdb_duty_cycle = 0.3;


}


void get_OoO_params(system_core *core_params, system_L2* L2_params)
{

  core_params->machine_type = 0; // OoO
  core_params->number_hardware_threads = 2;
  core_params->number_instruction_fetch_ports = 2;
  core_params->fp_issue_width = 2;
  core_params->prediction_width = 1;
  core_params->pipelines_per_core[0] = 1;
  core_params->pipelines_per_core[1] = 1;
  core_params->pipeline_depth[0] = 15;
  core_params->pipeline_depth[1] = 15;

  core_params->instruction_window_scheme = 1; //RSBASED 0; // PHYREG
  core_params->instruction_window_size = 20;
  core_params->archi_Regs_IRF_size = 16;
  core_params->archi_Regs_FRF_size = 32;
  core_params->phy_Regs_IRF_size = 256;
  core_params->phy_Regs_FRF_size = 256;
  core_params->rename_scheme = 0; //RAM-based
  core_params->register_windows_size = 0;
  strcpy(core_params->LSU_order, "inorder");
  core_params->memory_ports = 2;

  L2_params->L2_config[0] = _l2_num_sets * _l2_assoc * _l2_line_size;
  L2_params->L2_config[1] =_l2_line_size;
  L2_params->L2_config[2] = _l2_assoc;
  L2_params->L2_config[3]  = 1; 
  L2_params->L2_config[4] = 1;
  L2_params->L2_config[5] = _l2_latency;
  L2_params->L2_config[6] = _l2_bank_width;
  L2_params->L2_config[7] = _l2_write_policy;
  L2_params->device_type = XML->sys.device_type;

    L2_params->ports[0] = 1;
    L2_params->ports[1] = 0;
    L2_params->ports[2] = 0;

    L2_params->buffer_sizes[0] = 1;
    L2_params->buffer_sizes[1] = 2;
    L2_params->buffer_sizes[2] = 2;
    L2_params->buffer_sizes[3] = 2;

}


//get stats

void translate_core_stats(Core::BaseCore *core, system_core * core_stats, unsigned long sim_cycles)
{
  core_stats->total_instructions = core->pow_total_uops; 
  core_stats->branch_instructions = core->pow_total_branches;
  core_stats->branch_mispredictions = core->pow_branch_mispredicts;  
  core_stats->load_instructions = core->pow_total_loads;  
  core_stats->store_instructions = core->pow_total_stores; 
  core_stats->committed_instructions = core->pow_committed_uops;
 
  // core cycles at potentially variable frequency
  core_stats->total_cycles = core->pow_cycles; 
  // get average frequency for this period
  core_stats->clock_rate = (int) ceil(sim_cycles* _cpu_clock_rate/ (double) core_stats->total_cycles);

  core_stats->idle_cycles = core->pow_idle_cycles;
  core_stats->busy_cycles = core_stats->total_cycles - core_stats->idle_cycles;

  core_stats->int_regfile_reads = core->pow_int_regreads;
  core_stats->float_regfile_reads = core->pow_fp_regreads;
  core_stats->int_regfile_writes = core->pow_int_regwrites;
  core_stats->float_regfile_writes = core->pow_fp_regwrites;


  core_stats->function_calls = core->pow_ctx_switches; 
  core_stats->cdb_alu_accesses = core->pow_alu_accesses;
  core_stats->cdb_fpu_accesses =  core->pow_fpu_accesses;
  core_stats->cdb_mul_accesses =  core->pow_mul_accesses;
  core_stats->ialu_accesses = core_stats->cdb_alu_accesses;
  core_stats->fpu_accesses = core_stats->cdb_fpu_accesses;
  core_stats->mul_accesses = core_stats->cdb_mul_accesses;

  core_stats->pipeline_duty_cycle = (double)core-> pow_committed_insns;
  core_stats->pipeline_duty_cycle /= (double) core->pow_cycles;
  core_stats->pipeline_duty_cycle /= (double)OOO_CORE_MODEL::COMMIT_WIDTH; 

  core_stats->itlb.total_accesses = core->pow_itlb_accesses;
  core_stats->itlb.total_misses = core->pow_itlb_misses;

  core_stats->dtlb.total_accesses = core->pow_dtlb_accesses; 
  core_stats->dtlb.total_misses = core->pow_dtlb_misses;  
  
  
   core_stats->BTB.read_accesses = core->pow_btb_lookups;
   core_stats->BTB.write_accesses = core->pow_btb_updates;
   /// OoO based stats
  core_stats->int_instructions = core->pow_fetch_uops; 
  core_stats->fp_instructions = 0;

  core_stats->committed_int_instructions = core->pow_committed_uops; 
  core_stats->committed_fp_instructions = 0;

  core_stats->ROB_reads = core->pow_rob_reads;
  core_stats->ROB_writes = core->pow_rob_writes;  
  core_stats->rename_reads = core->pow_int_regreads; 
  core_stats->rename_writes = core->pow_int_regwrites;
  core_stats->fp_rename_reads = core->pow_fp_regreads; 
  core_stats->fp_rename_writes = core->pow_fp_regwrites;

  core_stats->inst_window_reads = core->pow_dispatch_uops; 
  core_stats->inst_window_writes = core->pow_dispatch_uops;  
  core_stats->inst_window_wakeup_accesses = 0;
  core_stats->fp_inst_window_reads = 0;
  core_stats->fp_inst_window_writes = 0;
  core_stats->fp_inst_window_wakeup_accesses = 0;
  core_stats->context_switches = core->pow_ctx_switches;
  
}
void translate_L1Cache_stats(Memory::Controller *IL1, Memory::Controller *DL1, system_core *core_stats)
{
    core_stats->icache.read_accesses = IL1->pow_il1_hits+ IL1->pow_il1_misses; 
    core_stats->icache.read_misses = IL1->pow_il1_misses; 
   
    core_stats->dcache.read_accesses = DL1->pow_dl1_load_hits + DL1->pow_dl1_load_misses; 
    core_stats->dcache.read_misses = DL1->pow_dl1_load_misses; 
    core_stats->dcache.write_accesses = DL1->pow_dl1_store_hits + DL1->pow_dl1_store_misses; 
    core_stats->dcache.write_misses = DL1->pow_dl1_store_misses; 
  
}

void translate_L2Cache_stats(Memory::Controller *L2,system_L2 *L2_stats)
{

    L2_stats->read_accesses = L2->pow_l2_load_hits + L2->pow_l2_load_misses; 
    L2_stats->read_misses = L2->pow_l2_load_misses; 
    L2_stats->write_accesses = L2->pow_l2_store_hits + L2->pow_l2_store_misses; 
    L2_stats->write_misses = L2->pow_l2_store_misses; 
}

void translate_UncoreCache_stats(Memory::Controller *LLC, unsigned long sim_cycles, root_system * stats)
{

     // hack 
     stats->total_cycles  = double(double(_uncore_freq)/double(_cpu_clock_rate)) * (double)sim_cycles; 

   if (!_uncore_LLC)
  	return; 
   if (!private_l2)//FIXME: for all simulated platforms L3 is uncore and L2 is private
   {  
 
      stats->L2[0].read_accesses = LLC->pow_l2_load_hits + LLC->pow_l2_load_misses;
      stats->L2[0].read_misses =  LLC->pow_l2_load_misses;  
      stats->L2[0].write_accesses = LLC->pow_l2_store_hits + LLC->pow_l2_store_misses;  
      stats->L2[0].write_misses = LLC->pow_l2_store_misses;
   }  
   else{


      stats->L3[0].read_accesses = LLC->pow_l3_load_hits + LLC->pow_l3_load_misses;
      stats->L3[0].read_misses =  LLC->pow_l3_load_misses;  
      stats->L3[0].write_accesses = LLC->pow_l3_store_hits + LLC->pow_l3_store_misses;  
      stats->L3[0].write_misses = LLC->pow_l3_store_misses;

  }

}



void getcore_stats(int coreid, Core::BaseCore *core, Memory::Controller * IL1, Memory::Controller * DL1, Memory::Controller *L2,unsigned long sim_cycles)
{   

//	translate_core_stats(Core::BaseCore *core, system_core * core_stats, unsigned long sim_cycles);
	if (_num_cores > 1 && _uncore_LLC)//in multicore platforms all L2 are shared
	{
		assert(private_l2 & (_num_cores == _num_l2_caches));
		translate_core_stats(core, &XML->sys.core[coreid], sim_cycles);
		//void translate_L1Cache_stats(Memory::Controller *IL1, Memory::Controller *DL1, system_core *core_stats);
		translate_L1Cache_stats(IL1, DL1, &XML->sys.core[coreid]);
		//void translate_L2Cache_stats(Memory::Controller *L2,system_L2 *L2_stats);
		translate_L2Cache_stats(L2,&XML->sys.L2[coreid]); 
	} 

	else
	{
		//FIXME: why  does assertion fails
		//assert(coreid = 0);
		translate_core_stats(core, &XML->sys.core[coreid],sim_cycles);
		translate_L1Cache_stats(IL1, DL1, &XML->sys.core[coreid]);
		// FIXME mechanism to make sure that single core has L2 is not uncore
		translate_L2Cache_stats(L2, &XML->sys.L2[coreid]);
	}

}

void get_uncore_stats(Memory::Controller *LLC, unsigned long sim_cycles)
{
	if (LLC == NULL)
		assert(!_uncore_LLC); 
	translate_UncoreCache_stats(LLC, sim_cycles, &XML->sys); 

}



void calc_power(bool print_power)
{

	mcpat_compute_energy(print_power, cores_rtp, &uncore_rtp); 
	if(pow_out_trace)
	{

		for (int i =0 ; i <_num_cores; i++)
			fprintf(pow_out_trace, "%.4f, %.4f\n", cores_rtp[i], cores_leakage[i]);
		fprintf(pow_out_trace, "%.4f, %.4f\n", uncore_rtp, uncore_leakage);
			
	
	}
}

void dump_pow_stats()
{
	FILE *dump_stats = fopen("/home/prism/dump","w");
	
	fprintf(dump_stats,"total uops: %f\n", XML->sys.core->total_instructions); 	
	fprintf(dump_stats,"total branches: %f\n", XML->sys.core->branch_instructions); 	
	fprintf(dump_stats,"mispredictions: %f\n", XML->sys.core->branch_mispredictions); 	
	fprintf(dump_stats,"commited uops: %f\n", XML->sys.core->committed_instructions); 	
	fprintf(dump_stats,"load insns.: %f\n", XML->sys.core->load_instructions); 	
	fprintf(dump_stats,"store insns.: %f\n", XML->sys.core->store_instructions); 	
	fprintf(dump_stats,"clock cycles.: %f\n", XML->sys.core->total_cycles); 	
	fprintf(dump_stats,"clock rate: %d\n", XML->sys.core->clock_rate); 	
	fprintf(dump_stats,"idle cycles: %f\n", XML->sys.core->idle_cycles); 	
	fprintf(dump_stats,"busy cycles: %f\n", XML->sys.core->busy_cycles); 	
	fprintf(dump_stats,"int reg reads: %f\n", XML->sys.core->int_regfile_reads); 	
	fprintf(dump_stats,"int reg writes: %f\n", XML->sys.core->int_regfile_writes); 	
	fprintf(dump_stats,"fp reg reads: %f\n", XML->sys.core->float_regfile_reads); 	
	fprintf(dump_stats,"fp reg writes: %f\n", XML->sys.core->float_regfile_writes); 	
	fprintf(dump_stats,"fn calls: %f\n", XML->sys.core->function_calls); 	
	fprintf(dump_stats,"alu accesses: %f\n", XML->sys.core->ialu_accesses); 	
	fprintf(dump_stats,"fpu accesses: %f\n", XML->sys.core->fpu_accesses); 	
	fprintf(dump_stats,"mul accesses: %f\n", XML->sys.core->mul_accesses); 	
	fprintf(dump_stats,"duty cycle: %f\n", XML->sys.core->pipeline_duty_cycle); 	


	fprintf(dump_stats,"int fetch uops: %f\n", XML->sys.core->int_instructions); 	
	fprintf(dump_stats,"fp fetch uops: %f\n", XML->sys.core->fp_instructions); 	
	fprintf(dump_stats,"ROB reads: %f\n", XML->sys.core->ROB_reads); 	
	fprintf(dump_stats,"ROB writes: %f\n", XML->sys.core->ROB_writes); 	
	fprintf(dump_stats,"rename reads: %f\n", XML->sys.core->rename_reads); 	
	fprintf(dump_stats,"rename writes: %f\n", XML->sys.core->rename_writes); 	
	fprintf(dump_stats,"fp rename reads: %f\n", XML->sys.core->fp_rename_reads); 	
	fprintf(dump_stats,"fp rename writes: %f\n", XML->sys.core->fp_rename_writes); 	
	fprintf(dump_stats,"alloc uops: %f\n", XML->sys.core->inst_window_reads); 	
	fprintf(dump_stats,"context switches: %f\n", XML->sys.core->context_switches); 	






	fprintf(dump_stats,"itlb accesses %f\n", XML->sys.core->itlb.total_accesses); 	
	fprintf(dump_stats,"itlb misses: %f\n", XML->sys.core->itlb.total_misses); 	
	fprintf(dump_stats,"dtlb accesses: %f\n", XML->sys.core->dtlb.total_accesses); 	
	fprintf(dump_stats,"dtlb misses: %f\n", XML->sys.core->dtlb.total_misses); 	
	fprintf(dump_stats,"il1 lookups.: %f\n", XML->sys.core->icache.read_accesses); 	
	fprintf(dump_stats,"il1 misses: %f\n", XML->sys.core->icache.read_misses); 	
	fprintf(dump_stats,"dl1 load lookups: %f\n", XML->sys.core->dcache.read_accesses); 	
	fprintf(dump_stats,"dl1 load misses: %f\n", XML->sys.core->dcache.read_misses); 	
	fprintf(dump_stats,"dl1 store lookups: %f\n", XML->sys.core->dcache.write_accesses); 	
	fprintf(dump_stats,"dl1 store misses: %f\n", XML->sys.core->dcache.write_misses); 	
	fprintf(dump_stats,"dl2 load lookups: %f\n",  XML->sys.L2->read_accesses); 	
	fprintf(dump_stats,"dl2 load misses: %f\n",  XML->sys.L2->read_misses); 	
	fprintf(dump_stats,"dl2 store lookups : %f\n", XML->sys.L2->write_accesses); 	
	fprintf(dump_stats,"dl2 store misses: %f\n",  XML->sys.L2->write_misses); 	
	fprintf(dump_stats,"btb lookups: %f\n", XML->sys.core->BTB.read_accesses); 	
	fprintf(dump_stats,"btb updates: %f\n", XML->sys.core->BTB.write_accesses); 	
	

	fprintf(dump_stats,"uncore cycles: %f\n", XML->sys.total_cycles); 

	fclose(dump_stats);
//	return XML->sys.L2->read_accesses ; 
	
}
void kill_power()
{
	free(cores_rtp); 
	if(pow_out_trace)
		fclose(pow_out_trace); 
	mcpat_finalize(); 
}
