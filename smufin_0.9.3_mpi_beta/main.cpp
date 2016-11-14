#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <queue>
#include <mpi.h>
#include <string>
#include <sstream>
#include <algorithm>
#include "samtools-0.1.19/sam_header.h"
#include "samtools-0.1.19/sam.h"
#include "samtools-0.1.19/faidx.h"
#include "samtools-0.1.19/kstring.h"
#include "samtools-0.1.19/khash.h"
#include "util_structs.h" 

#include <sys/time.h>
#include <sys/resource.h>
#include "benchmark.h"
#include <chrono>


#define PACKAGE_VERSION "0.9.2 Beta"

#define BUFFER_SIZE 1<<22

char* buffer;

// MPI define functions

#define MASTER_ID 0
#define MPI_EXIT_CMD 1
#define MPI_INIT_BWA_INDX_CMD 2
#define MPI_DO_TREE_BLOCK 3
#define MPI_GET_READS_FROM_REGION 4
#define MPI_GET_BAM_CHRS 5
#define MPI_GET_MUTS_FROM_BLOCK 6
#define MPI_GET_MUTS_FROM_BLOCK_DONE 7
#define MPI_GET_SEQUENCES_FROM_BLOCK 8
#define MPI_GET_SEQS_FROM_PARTITION 9
#define MPI_GET_MUTS_FROM_BLOCK_FASTQ 10
#define MPI_SET_SUPP_READS_CMD 11

int NUM_NODES;
int CPUS_PER_NODE = 1;
int min_sup_reads = 4;
int tumor_cont_perc = 0;
int mut_count = 1;

char* patient_id;
int chr_id = 0, pos = 0;

extern "C" {
int my_bwa_index(char* ref_file);
void init_bwa_idx(char* ref_file);
int gsuffix_load_fastq_sequences(char* mem_fastq_1, size_t size_mem_fastq_1, char* mem_fastq_2, size_t size_mem_fastq_2, FILE* fout);

int my_main_mem(char* ref_file, char* mem_seqs, size_t size_mem_seqs, FILE* fout);

int gsuffix_co(char* mem, size_t size_mem, FILE* out);
int gsuffix_so(char* mem, size_t size_mem, FILE* stream_out);
int gsuffix_me(FILE* stream_out, int n, char** mem_files, size_t* size_mem_files, int n_threads);
int gsuffix_fi(char* mem, size_t size_mem, FILE* stream_out);
int gsuffix_fget_size(FILE* mem);
void* gsuffix_routine(void* param);

int get_partitions_from_header(char* mem, size_t size_mem, char*** parts_id, int* num_parts, int** parts_sizes);

int get_seqs_from_region(char* mem_block, size_t size_mem_block, char* mem_iblock, size_t size_mem_iblock, int part_id, int offset_ini, int offset_end, FILE* fout);

int gsuffix_somatic_bp(char* ref_file, FILE* normal_mem, FILE* cancer_mem, FILE* out_file);

}

std::pair<char*, size_t> gsuffix_get_reads_from_blockfile(char* block_file, char* c, int pini, int pend);
void gsuffixfilteraux(std::map<std::string, int>& hash_cs, char* mem, size_t size, char* mem_normal, size_t size_normal, char* mem_inormal, size_t size_inormal, std::vector<SNV>* vsnvs);

static int use_normal_bam = 0;
static int use_cancer_bam = 0;
static int num_threads = 1;
static int print_brkp_ext = 0;
static int usage(FILE* stream_out)
{
        fprintf(stream_out, "\n");
        fprintf(stream_out, "Program: SMuFin (Somatic Mutation Finder)\n");
        fprintf(stream_out, "Version: %s\n", PACKAGE_VERSION);
        fprintf(stream_out, "Contact: http://cg.bsc.es/smufin\n");
        fprintf(stream_out, "Usage:   SMuFin <command> [options]\n\n");
        fprintf(stream_out, "Command: --ref                 <FILE>	Reference genome\n");
        fprintf(stream_out, "         --normal_fastq_1      <FILE>	Normal FASTQ 1st Paired-End file\n");
	fprintf(stream_out, "         --normal_fastq_2      <FILE>	Normal FASTQ 2nd Paired-End file\n");
        fprintf(stream_out, "         --tumor_fastq_1       <FILE>	Tumor FASTQ 1st Paired-End file\n");
        fprintf(stream_out, "         --tumor_fastq_2       <FILE>	Tumor FASTQ 2nd Paired-End file\n");
        fprintf(stream_out, "         --min_sup_reads       <NUM>	Minimum number of supporting tumor reads to accept breakpoint <4>\n");
        fprintf(stream_out, "         --tumor_cont_perc     <NUM>	Expected percentage of tumor contamination in normal dataset 0-100 <0>\n");
        fprintf(stream_out, "         --cpus_per_node       <NUM>	Use NUM threads(cores) per node <1>\n");
        fprintf(stream_out, "\n\n");
        return 1;
}

int is_int(char* p)
{
	while (*p && isdigit(*p))
		p++;
	return !*p;
}

int get_n_lines_from_gzip(gzFile file, char** mem_out, size_t* size_mem_out, int num_lines)
{
	char* buffer = new char[BUFFER_SIZE];
	FILE* fph_out = open_memstream(mem_out, size_mem_out);
	int i = 0;
	while (i < num_lines)
	{	
		char* line = gzgets(file, buffer, BUFFER_SIZE);
		if (line)		
			fwrite(line, strlen(line), 1, fph_out);
		else 
			break;
		i++;
	}
	fclose(fph_out);
	delete[] buffer;
	return i;
}

void read_mem_from_file(char* filename, char** mem, size_t* size)
{
	FILE* fh = fopen(filename, "rb");
	FILE* out_fh = open_memstream(mem, size);
	size_t bytes_read;
	void* buffer_aux = malloc(BUFFER_SIZE);
	while ((bytes_read = fread(buffer_aux, 1, BUFFER_SIZE, fh)) > 0)
	{
		fwrite(buffer_aux, 1, bytes_read, out_fh);
	}
	free(buffer_aux);
	fclose(fh);
	fclose(out_fh);
}

void print_point_mutation(int id, const char* chr, int pos, const char org_bp, const char mut_bp, int num_normal_in_normal, int num_tumor_in_normal, int num_normal_in_tumor, int num_tumor_in_tumor, FILE* fph_out)
{
		fprintf(fph_out, "%d\tSNV\t%s\t%d\t%c\t%c\n", id, chr, pos, org_bp, mut_bp);
}

void print_deletion(int id, const char* chr, int pos, int size, BRK_EXT* ext, int num_normal_in_normal, int num_tumor_in_normal, int num_normal_in_tumor, int num_tumor_in_tumor, FILE* fph_out)
{
	if (!ext) fprintf(fph_out, "%d\tDEL\t%s\t%d\t%d\n", id, chr, pos, size);
}

void print_short_insertion(int id, const char* chr, int pos, int size, const char* seq, int num_normal_in_normal, int num_tumor_in_normal, int num_normal_in_tumor, int num_tumor_in_tumor, FILE* fph_out)
{
	fprintf(fph_out, "%d\tINS\t%s\t%d\t%d\t%s\n", id, chr, pos, size, seq);
}

void print_long_insertion(char* chr, char* pos, BRK_EXT* ext_1, BRK_EXT* ext_2, FILE* fph_out)
{
	//fprintf(fph_out, "%d\tINS\t%s\t%s\tND\t#BreakPoint_1\t%d\t#BreakPoint_2\t%d\n", mut_count++, chr, pos, ext_1->_id, ext_2->_id);
}

void print_inversion(int id, const char* chr, int pos, int size, BRK_EXT* ext_1, BRK_EXT* ext_2, int num_normal_in_normal, int num_tumor_in_normal, int num_normal_in_tumor, int num_tumor_in_tumor, FILE* fph_out)
{
	if (!ext_1 && !ext_2) fprintf(fph_out, "%d\tINV\t%s\t%d\t%d\n", id, chr, pos, size);
	//else if (ext_1 && ext_2) fprintf(fph_out, "INV\t%s\t%s\t%d\t%d\t#BreakPoint_1\t%d\t#BreakPoint_2\t%d\n", chr, pos, size, supp_reads, ext_1->_id, ext_2->_id);
}

void init_det_muts_block(det_muts_block* block)
{
	block->_task_assigned = 0;
	block->_nblock_file = 0;
	block->_tblock_file = 0;	
	block->_mem_control_part = 0;
	block->_mem_control_index = 0;
	block->_mem_tumor_part = 0;
	block->_mem_tumor_index = 0;
	block->_vsnvs = new std::vector<SNV>;
	block->_vdels = new std::vector<DEL>;
	block->_vins = new std::vector<INS>;
	block->_vinvs = new std::vector<INV>;
	block->_vbps = new std::vector<BKP>;
	block->_vbpse = new std::vector<BKPE>;
}

void free_det_muts_block(det_muts_block* block)
{
	delete block->_vsnvs;
	delete block->_vdels;
	delete block->_vins;
	delete block->_vinvs;
	delete block->_vbps;
	delete block->_vbpse;
	if (block->_mem_control_part) delete[] block->_mem_control_part;
	if (block->_mem_control_index) delete[] block->_mem_control_index;
	if (block->_mem_tumor_part) delete[] block->_mem_tumor_part;
	if (block->_mem_tumor_index) delete[] block->_mem_tumor_index;
	init_det_muts_block(block);
}

void* process_det_mut_block(void* param)
{
	det_muts_block* b = (det_muts_block*)param;

	char* mem_inormal;
	char* mem_itumor;
	size_t size_inormal;
	size_t size_itumor;
	std::pair<char*, size_t> mem_normal_reads;
	std::pair<char*, size_t> mem_tumor_reads;

	if (b->_nblock_file && b->_tblock_file)
	{
		mem_normal_reads = gsuffix_get_reads_from_blockfile(b->_nblock_file, b->_c, b->_pini, b->_pend);
		mem_tumor_reads = gsuffix_get_reads_from_blockfile(b->_tblock_file, b->_c, b->_pini, b->_pend);
	
		FILE* stream_inormal = open_memstream(&mem_inormal, &size_inormal);
		FILE* stream_itumor = open_memstream(&mem_itumor, &size_itumor);

		gsuffix_fi(mem_normal_reads.first, mem_normal_reads.second, stream_inormal);
		gsuffix_fi(mem_tumor_reads.first, mem_tumor_reads.second, stream_itumor);
	}
	else
	{
		mem_normal_reads.first = b->_mem_control_part;
		mem_normal_reads.second = b->_size_control_part;
		mem_inormal = b->_mem_control_index;
		size_inormal = b->_size_control_index;
		mem_tumor_reads.first = b->_mem_tumor_part;
        	mem_tumor_reads.second = b->_size_tumor_part;
        	mem_itumor = b->_mem_tumor_index;
        	size_itumor = b->_size_tumor_index;

	}

        sfile_t* block_file = blockfopen(fmemopen(mem_normal_reads.first, mem_normal_reads.second, "rb"), 0);
	std::map<std::string, int> hash_cs;
	for (int i=0; i<block_file->header->n_targets; i++)
	{
		hash_cs[std::string(block_file->header->target_name[i])] = i;
	}

	det_params* pdet = new det_params;
	pdet->_block_id = b->_block_id;
	pdet->_ref_filename = b->_ref_genome_filename;
	pdet->_mem_normal_reads = mem_normal_reads.first;
	pdet->_size_normal_reads = mem_normal_reads.second;
	pdet->_mem_tumor_reads = mem_tumor_reads.first;
	pdet->_size_tumor_reads = mem_tumor_reads.second;
	pdet->_mem_inormal = mem_inormal;
	pdet->_size_inormal = size_inormal;
	pdet->_mem_itumor = mem_itumor;
	pdet->_size_itumor = size_itumor;
	if (b->_nblock_file && b->_tblock_file)
	{
		pdet->_c = b->_c;
		pdet->_pini = b->_pini;
		pdet->_pend = b->_pend;
	}
	else
	{
		pdet->_c = block_file->header->target_name[b->_partition._part_id];
		pdet->_pini = b->_partition._offset_ini;
		pdet->_pend = b->_partition._offset_end;
	}
	pdet->_hash_cs = &hash_cs;
	pdet->_vdels = b->_vdels;
	pdet->_vins = b->_vins;
	pdet->_vinvs = b->_vinvs;
	pdet->_vbps = b->_vbps;
	pdet->_vbpse = b->_vbpse;
	pdet->_nsize = b->_nsize;
	pdet->_tsize = b->_tsize;
	gsuffix_routine(pdet);

	FILE* normal_block = fmemopen(mem_normal_reads.first, mem_normal_reads.second, "rb");
	FILE* tumor_block = fmemopen(mem_tumor_reads.first, mem_tumor_reads.second, "rb");
	
	char* mem_snvs;
	size_t size_snvs;
	FILE* stream_snvs = open_memstream(&mem_snvs, &size_snvs);

	gsuffix_somatic_bp(b->_ref_genome_filename, tumor_block, normal_block, stream_snvs);
	

	fclose(normal_block);
	fclose(tumor_block);
	fclose(stream_snvs);
	
        gsuffixfilteraux(hash_cs, mem_snvs, size_snvs, mem_normal_reads.first, mem_normal_reads.second, mem_inormal, size_inormal, b->_vsnvs);

	blockclose(block_file);
	delete pdet;

	if (b->_nblock_file && b->_tblock_file)
	{
		free(mem_normal_reads.first);
		free(mem_tumor_reads.first);
		free(mem_inormal);
		free(mem_itumor);

		delete[] b->_ref_genome_filename;
		delete[] b->_c;
		delete[] b->_nblock_file;
		delete[] b->_tblock_file;
	}
	else
	{
		delete[] b->_ref_genome_filename;
	}
	b->_task_finished = 1;
	pthread_exit(0);
}

void* process_block_func(void* param)
{
	gsufix_tree_block* b = (gsufix_tree_block*)param;

	char* mem;
	size_t size_mem;
	FILE* fph = open_memstream(&mem, &size_mem);
	gsuffix_load_fastq_sequences(b->_mem_fastq_1, b->_size_fastq_1, b->_mem_fastq_2, b->_size_fastq_2, fph);
	fclose(fph);

	char* mem_block;
	size_t size_block;

	FILE* block = open_memstream(&mem_block, &size_block);	
	gsuffix_co(mem, size_mem, block);
	FILE* blocks = open_memstream(&b->_mem_block_out, &b->_size_block_out);
	gsuffix_so(mem_block, size_block, blocks);
	
	free(b->_mem_fastq_1);
	free(b->_mem_fastq_2);
	free(mem);
	free(mem_block);
	
	FILE* iblock = open_memstream(&b->_mem_iblock_out, &b->_size_iblock_out);
	gsuffix_fi(b->_mem_block_out, b->_size_block_out, iblock);
	b->_task_finished = 1;
	pthread_exit(0);
}

void recv_tree_block(std::queue<std::pair<int, int> >& free_cpus, std::vector<TreeBlock*>& blocks)
{
	MPI_Status status;
        int msg_recv;

	MPI_Iprobe(MPI_ANY_SOURCE, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &msg_recv, &status);
	if (msg_recv)
	{
		TreeBlock* btree = new TreeBlock;
		int node_id = status.MPI_SOURCE;
		int cpu_id;
		btree->_num_node = node_id;
		MPI_Recv(&cpu_id, 1, MPI_INT, node_id, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
		btree->_num_cpu = cpu_id;
		MPI_Recv(&btree->_block_address, sizeof(char*), MPI_CHAR, node_id, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
		MPI_Recv(&btree->_block_address_size, sizeof(size_t), MPI_CHAR, node_id, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
		MPI_Recv(&btree->_iblock_address, sizeof(char*), MPI_CHAR, node_id, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
                MPI_Recv(&btree->_iblock_address_size, sizeof(size_t), MPI_CHAR, node_id, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
		blocks.push_back(btree);
		free_cpus.push(std::pair<int, int>(node_id, cpu_id));
		//fprintf(stderr, "Received Tree Block from node %d cpu %d\n", node_id, cpu_id);
	}
}

void generate_quadtree(char* fastq_1_file, char* fastq_2_file, std::vector<TreeBlock*>& blocks)
{
	FILE* fqs_cont_1 = fopen(fastq_1_file, "r");
	FILE* fqs_cont_2 = fopen(fastq_2_file, "r");
	char* filename_1 = new char[FILENAME_MAX];
	char* filename_2 = new char[FILENAME_MAX];

	std::queue<std::pair<int, int> > free_cpus;

	for (int i=1; i<=NUM_NODES; i++)
		for (int j=0; j<CPUS_PER_NODE; j++)
			free_cpus.push(std::pair<int, int>(i, j));

	while (fgets(filename_1, FILENAME_MAX, fqs_cont_1))
	{
		char* mem_1, *mem_2;
		size_t size_mem_1, size_mem_2;

		strtok(filename_1, "\n\t");
		gzFile fp_1 = gzopen(filename_1, "r");
		fgets(filename_2, FILENAME_MAX, fqs_cont_2);
		strtok(filename_2, "\n\t");
		gzFile fp_2 = gzopen(filename_2, "r");

		while (1)
		{
			recv_tree_block(free_cpus, blocks);
			
			if (!free_cpus.empty())
			{
				if (get_n_lines_from_gzip(fp_1, &mem_1, &size_mem_1, 40000) > 0 && 
				    get_n_lines_from_gzip(fp_2, &mem_2, &size_mem_2, 40000))
				{
					MPI_Send(&free_cpus.front().second, 1, MPI_INT, free_cpus.front().first, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD);
					MPI_Send(mem_1, size_mem_1, MPI_CHAR, free_cpus.front().first, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD);
                		        MPI_Send(mem_2, size_mem_2, MPI_CHAR, free_cpus.front().first, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD);
					free(mem_1);
					free(mem_2);
					free_cpus.pop();
				}
				else
				{
					gzclose(fp_1);
					gzclose(fp_2);
					break;
				}
			}
		}
	}
	while (free_cpus.size() < NUM_NODES*CPUS_PER_NODE)
		recv_tree_block(free_cpus, blocks);

	fclose(fqs_cont_1);
	fclose(fqs_cont_2);
	delete[] filename_1;
	delete[] filename_2;
}

FILE* stdout_backup;
FILE* stdin_backup;
FILE* stderr_backup;

int check_fastqs_access(char* filename, FILE* fout)
{
	FILE* f = fopen(filename, "r");
	char* buf_filename = 0;
	size_t filen_size = 0;
	while (getline(&buf_filename, &filen_size, f) > 0)
	{
		char* new_line = strchr(buf_filename, '\n');
		if (new_line) *new_line = '\0';
		if (access(buf_filename, R_OK))
		{
			fprintf(fout, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", buf_filename);
			return 0;
		
		}
		free(buf_filename);
		buf_filename=0;
	}
	fclose(f);
	return 1;
}

void slave(int node_id)
{
        struct rusage slave_rs;
        getrusage(RUSAGE_SELF, &slave_rs);
	FILE* stdnull = fopen("/dev/null", "w");
	FILE* stderr_backup = stderr;
	//stderr = stdnull;


	MPI_Status status;
	MPI_Request request;
	char* buffer;
	int num_elems_recv;
	int msg_recv;
	int cpu;
	std::vector<gsufix_tree_block*> gst_b;
	std::vector<det_muts_block*> det_b;

	
	while (1)
	{
		usleep(1000);
		MPI_Iprobe(MASTER_ID, MPI_ANY_TAG, MPI_COMM_WORLD, &msg_recv, &status);
		
		if (!msg_recv)
		{
			for (int i=0; i<gst_b.size(); i++)
			{
				if (gst_b[i]->_task_assigned && gst_b[i]->_task_finished)
				{
					MPI_Isend(&gst_b[i]->_cpu, 1, MPI_INT, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					MPI_Isend(&gst_b[i]->_mem_block_out, sizeof(char*), MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					MPI_Isend(&gst_b[i]->_size_block_out, sizeof(size_t), MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					MPI_Isend(&gst_b[i]->_mem_iblock_out, sizeof(char*), MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
                                        MPI_Isend(&gst_b[i]->_size_iblock_out, sizeof(size_t), MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					gst_b[i]->_task_assigned = 0;
				}
			}
			for (int i=0; i<det_b.size(); i++)
			{
				if (det_b[i]->_task_assigned && det_b[i]->_task_finished)
				{
					MPI_Isend(&det_b[i]->_cpu_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
			        	MPI_Isend(&det_b[i]->_block_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					MPI_Isend(&det_b[i]->_vsnvs->front(), det_b[i]->_vsnvs->size()*sizeof(SNV), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
			                MPI_Isend(&det_b[i]->_vdels->front(), det_b[i]->_vdels->size()*sizeof(DEL), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
                	                MPI_Isend(&det_b[i]->_vins->front(), det_b[i]->_vins->size()*sizeof(INS), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
                        	        MPI_Isend(&det_b[i]->_vinvs->front(), det_b[i]->_vinvs->size()*sizeof(INV), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
                                	MPI_Isend(&det_b[i]->_vbps->front(), det_b[i]->_vbps->size()*sizeof(BKP), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					MPI_Isend(&det_b[i]->_vbpse->front(), det_b[i]->_vbpse->size()*sizeof(BKPE), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
					det_b[i]->_task_assigned = 0;
				}
			}
		}
		
		else if (status.MPI_TAG == MPI_EXIT_CMD)
		{
			MPI_Recv(0, 0, MPI_CHAR, MASTER_ID, MPI_EXIT_CMD, MPI_COMM_WORLD, &status);
			break;
		}
		else if (status.MPI_TAG == MPI_INIT_BWA_INDX_CMD)
		{
			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			char* ref_file = new char[num_elems_recv];
			MPI_Recv(ref_file, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_INIT_BWA_INDX_CMD, MPI_COMM_WORLD, &status);
			init_bwa_idx(ref_file);
			MPI_Send(0, 0, MPI_CHAR, MASTER_ID, MPI_INIT_BWA_INDX_CMD, MPI_COMM_WORLD);
			delete[] ref_file;
		}
		else if (status.MPI_TAG == MPI_SET_SUPP_READS_CMD)
		{
			MPI_Recv(&min_sup_reads, 1, MPI_INT, MASTER_ID, MPI_SET_SUPP_READS_CMD, MPI_COMM_WORLD, &status);
			MPI_Send(0, 0, MPI_CHAR, MASTER_ID, MPI_SET_SUPP_READS_CMD, MPI_COMM_WORLD);
		}
		else if (status.MPI_TAG == MPI_DO_TREE_BLOCK)
		{
			MPI_Recv(&cpu, 1, MPI_INT, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
			if (cpu+1 > gst_b.size()) gst_b.push_back(new gsufix_tree_block);
			gst_b[cpu]->_cpu = cpu;
                        MPI_Probe(MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
                        gst_b[cpu]->_mem_fastq_1 = new char[num_elems_recv];
                        gst_b[cpu]->_size_fastq_1 = num_elems_recv;
                        MPI_Recv(gst_b[cpu]->_mem_fastq_1, gst_b[cpu]->_size_fastq_1, MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Probe(MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
                        gst_b[cpu]->_mem_fastq_2 = new char[num_elems_recv];
                        gst_b[cpu]->_size_fastq_2 = num_elems_recv;
                        MPI_Recv(gst_b[cpu]->_mem_fastq_2, gst_b[cpu]->_size_fastq_2, MPI_CHAR, MASTER_ID, MPI_DO_TREE_BLOCK, MPI_COMM_WORLD, &status);
			gst_b[cpu]->_task_assigned = 1;
			gst_b[cpu]->_task_finished = 0;
                        pthread_create(&gst_b[cpu]->_t, 0, process_block_func, gst_b[cpu]);
		}
		else if (status.MPI_TAG == MPI_GET_BAM_CHRS)
		{
			TreeBlock block;
			MPI_Recv(&block, sizeof(TreeBlock), MPI_CHAR, MASTER_ID, MPI_GET_BAM_CHRS, MPI_COMM_WORLD, &status);
			char** chrs;
			int num_chrs;
			int* chr_sizes;
			get_partitions_from_header(block._block_address, block._block_address_size, &chrs, &num_chrs, &chr_sizes);
			MPI_Ssend(&num_chrs, 1, MPI_INT, MASTER_ID, MPI_GET_BAM_CHRS, MPI_COMM_WORLD);
			MPI_Ssend(chr_sizes, num_chrs, MPI_INT, MASTER_ID, MPI_GET_BAM_CHRS, MPI_COMM_WORLD);
			for (int i=0; i<num_chrs; i++)
			{
				MPI_Ssend(chrs[i], strlen(chrs[i])+1, MPI_CHAR, MASTER_ID, MPI_GET_BAM_CHRS, MPI_COMM_WORLD);
				free(chrs[i]);
			}
			free(chrs);
		}
		else if (status.MPI_TAG == MPI_GET_MUTS_FROM_BLOCK)
		{
			int cpu_id;
			int block_id;
			char* ref_filename;
			char* c;
			int p_ini;
			int p_end;
			char* nblock_file;
			char* tblock_file;
			int nsize;
			int tsize;

                        MPI_Recv(&cpu_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Recv(&block_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
                        ref_filename = new char[num_elems_recv];
			MPI_Recv(ref_filename, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			c = new char[num_elems_recv];
			MPI_Recv(c, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Recv(&p_ini, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Recv(&p_end, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			nblock_file = new char[num_elems_recv];
			MPI_Recv(nblock_file, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
                        tblock_file = new char[num_elems_recv];
                        MPI_Recv(tblock_file, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Recv(&nsize, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);
                        MPI_Recv(&tsize, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK, MPI_COMM_WORLD, &status);

			if (cpu_id+1 > det_b.size()) 
			{
				det_b.push_back(new det_muts_block);
        			init_det_muts_block(det_b[cpu_id]);
			}
			free_det_muts_block(det_b[cpu_id]);
		
			det_b[cpu_id]->_cpu_id = cpu_id;
			det_b[cpu_id]->_block_id = block_id;
			det_b[cpu_id]->_ref_genome_filename = ref_filename;
			det_b[cpu_id]->_c = c;
			det_b[cpu_id]->_pini = p_ini;
			det_b[cpu_id]->_pend = p_end;
			det_b[cpu_id]->_nblock_file = nblock_file;
			det_b[cpu_id]->_tblock_file = tblock_file;
			det_b[cpu_id]->_task_assigned = 1;
			det_b[cpu_id]->_task_finished = 0;
			det_b[cpu_id]->_nsize = nsize;
			det_b[cpu_id]->_tsize = tsize;
			pthread_create(&det_b[cpu_id]->_t, 0, process_det_mut_block, det_b[cpu_id]);
		}
		else if (status.MPI_TAG == MPI_GET_MUTS_FROM_BLOCK_FASTQ)
		{
			int cpu_id;
			int block_id;
			char* ref_filename;
			TreePartition partition;
			char* mem_control_part;
			size_t size_control_part;
			char* mem_tumor_part;
			size_t size_tumor_part;
			char* mem_control_index;
			size_t size_control_index;
			char* mem_tumor_index;
			size_t size_tumor_index;
			int nsize;
			int tsize;

                        MPI_Recv(&cpu_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Recv(&block_id, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
                        ref_filename = new char[num_elems_recv];
			MPI_Recv(ref_filename, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Recv(&partition, sizeof(TreePartition), MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Recv(&nsize, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Recv(&tsize, 1, MPI_INT, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			size_control_part = num_elems_recv;
			mem_control_part = new char[num_elems_recv];
			MPI_Recv(mem_control_part, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			size_control_index = num_elems_recv;
			mem_control_index = new char[num_elems_recv];
			MPI_Recv(mem_control_index, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			size_tumor_part = num_elems_recv;
                        mem_tumor_part = new char[num_elems_recv];
                        MPI_Recv(mem_tumor_part, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
			MPI_Probe(MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);
                        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			size_tumor_index = num_elems_recv;
			mem_tumor_index = new char[num_elems_recv];
			MPI_Recv(mem_tumor_index, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD, &status);

			if (cpu_id+1 > det_b.size()) 
			{
				det_b.push_back(new det_muts_block);
        			init_det_muts_block(det_b[cpu_id]);
			}

                        //fprintf(stderr, "MPI_GET_MUTS_FROM_BLOCK_FASTQ_3\n");
			free_det_muts_block(det_b[cpu_id]);
		
                        //fprintf(stderr, "MPI_GET_MUTS_FROM_BLOCK_FASTQ_4\n");
			det_b[cpu_id]->_cpu_id = cpu_id;
			det_b[cpu_id]->_block_id = block_id;
			det_b[cpu_id]->_ref_genome_filename = ref_filename;
			det_b[cpu_id]->_partition = partition;
			det_b[cpu_id]->_task_assigned = 1;
			det_b[cpu_id]->_task_finished = 0;
			det_b[cpu_id]->_nsize = nsize;
			det_b[cpu_id]->_tsize = tsize;
			det_b[cpu_id]->_mem_control_part = mem_control_part;
			det_b[cpu_id]->_size_control_part = size_control_part;
			det_b[cpu_id]->_mem_tumor_part = mem_tumor_part;
			det_b[cpu_id]->_size_tumor_part = size_tumor_part;
			det_b[cpu_id]->_mem_control_index = mem_control_index;
			det_b[cpu_id]->_size_control_index = size_control_index;
			det_b[cpu_id]->_mem_tumor_index = mem_tumor_index;
			det_b[cpu_id]->_size_tumor_index = size_tumor_index;
			pthread_create(&det_b[cpu_id]->_t, 0, process_det_mut_block, det_b[cpu_id]);
		}
		else if (status.MPI_TAG == MPI_GET_SEQUENCES_FROM_BLOCK)
		{
			TreeBlock block;
                	MPI_Recv(&block, sizeof(TreeBlock), MPI_CHAR, MASTER_ID, MPI_GET_SEQUENCES_FROM_BLOCK, MPI_COMM_WORLD, &status);
			MPI_Isend(block._block_address, block._block_address_size, MPI_CHAR, MASTER_ID, MPI_GET_SEQUENCES_FROM_BLOCK, MPI_COMM_WORLD, &request); MPI_Request_free(&request);
		}
		else if (status.MPI_TAG == MPI_GET_SEQS_FROM_PARTITION)
                {
                        TreeBlock* blocks;
			TreePartition part;
			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
			int num_blocks = num_elems_recv/sizeof(TreeBlock);
			blocks = new TreeBlock[num_blocks];
                        MPI_Recv(blocks, num_elems_recv, MPI_CHAR, MASTER_ID, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
                        MPI_Recv(&part, sizeof(TreePartition), MPI_CHAR, MASTER_ID, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
			char* mems[num_blocks];
			size_t sizes[num_blocks];
			for (int i=0; i<num_blocks; i++)
			{
				FILE* stream_out = open_memstream(&mems[i], &sizes[i]);
                        	get_seqs_from_region(blocks[i]._block_address, blocks[i]._block_address_size, blocks[i]._iblock_address, blocks[i]._iblock_address_size, part._part_id, part._offset_ini, part._offset_end, stream_out);
			}
			char* merged_mem;
			size_t merged_size;
			FILE* merged_partition = open_memstream(&merged_mem, &merged_size);
			gsuffix_me(merged_partition, num_blocks, mems, sizes, 1);
			MPI_Send(merged_mem, merged_size, MPI_CHAR, MASTER_ID, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD);
			for (int i=0; i<num_blocks; i++)
			{
				free(mems[i]);
			}
			delete[] blocks;
			free(merged_mem);
                }
	}
	for (int i=0; i<gst_b.size(); i++)
	{
		delete gst_b[i];
	}

	for (int i=0; i<det_b.size(); i++)
	{
		delete det_b[i];
	}

	stderr = stderr_backup;
	fclose(stdnull);
        getrusage(RUSAGE_SELF, &slave_rs);
        std::cout << "Slave " << node_id << " mem usage: " << slave_rs.ru_maxrss << std::endl;
}

MutationBlock* recv_mutations_detected_block(int* num_node, int* cpu_id)
{
	MPI_Status status;
	MPI_Request request;
	int num_elems_recv;
	int msg_recv;
	int block_id;

	MPI_Iprobe(MPI_ANY_SOURCE, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &msg_recv, &status);
        if (!msg_recv)	return NULL;

	MutationBlock* block = new MutationBlock;
	*num_node = status.MPI_SOURCE;
	MPI_Recv(cpu_id, 1, MPI_INT, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
	MPI_Recv(&block->_block_id, 1, MPI_INT, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	block->_snvs = new SNV[num_elems_recv/sizeof(SNV)];
	MPI_Recv(block->_snvs, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
	block->_num_snvs = num_elems_recv/sizeof(SNV);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	block->_dels = new DEL[num_elems_recv/sizeof(DEL)];
        MPI_Recv(block->_dels, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
	block->_num_dels = num_elems_recv/sizeof(DEL);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	block->_ins = new INS[num_elems_recv/sizeof(INS)];
        MPI_Recv(block->_ins, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        block->_num_ins = num_elems_recv/sizeof(INS);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	block->_invs = new INV[num_elems_recv/sizeof(INV)];
        MPI_Recv(block->_invs, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        block->_num_invs = num_elems_recv/sizeof(INV);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	block->_bps = new BKP[num_elems_recv/sizeof(BKP)];
        MPI_Recv(block->_bps, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        block->_num_bps = num_elems_recv/sizeof(BKP);
	MPI_Probe(*num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
        block->_bpse = new BKPE[num_elems_recv/sizeof(BKPE)];
        MPI_Recv(block->_bpse, num_elems_recv, MPI_CHAR, *num_node, MPI_GET_MUTS_FROM_BLOCK_DONE, MPI_COMM_WORLD, &status);
        block->_num_bpse = num_elems_recv/sizeof(BKPE);
	return block;
}


void get_sequences_from_block(TreeBlock* block, char** mem, size_t* mem_size)
{
	MPI_Status status;
	int num_elems_recv;
	MPI_Send(block, sizeof(TreeBlock), MPI_CHAR, block->_num_node, MPI_GET_SEQUENCES_FROM_BLOCK, MPI_COMM_WORLD);
	MPI_Probe(block->_num_node, MPI_GET_SEQUENCES_FROM_BLOCK, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
	*mem = new char[num_elems_recv];
	*mem_size = num_elems_recv;
	MPI_Recv(*mem, num_elems_recv, MPI_CHAR, block->_num_node, MPI_GET_SEQUENCES_FROM_BLOCK, MPI_COMM_WORLD, &status);
}

std::vector<MutationBlock*> get_breakpoints(char* ref_file, std::vector<TreeBlock*>& control_blocks, std::vector<TreeBlock*>& tumor_blocks)
{
	MPI_Status status;
	int msg_recv;
	int num_elems_recv;
	std::queue<std::pair<int, int> > free_nodes;
	char** parts;
	int* parts_sizes;
	int num_parts;

	char* control_blocks_reg[NUM_NODES];
        char* tumor_blocks_reg[NUM_NODES];
        size_t control_sizes_reg[NUM_NODES];
        size_t tumor_sizes_reg[NUM_NODES];

	char* block_mem_aux[100];
	size_t block_size_aux[100];
	for (int i=0; i<100; i++)
	{
		get_sequences_from_block(control_blocks[i], &block_mem_aux[i], &block_size_aux[i]);
	}
	char* block_mem;
	size_t block_size;
	FILE* fmem = open_memstream(&block_mem, &block_size);
	gsuffix_me(fmem, 100, block_mem_aux, block_size_aux, 1);

	FILE* fb = fmemopen(block_mem, block_size, "rb");

	int nsize = gsuffix_fget_size(fb);

        get_partitions_from_header(block_mem, block_size, &parts, &num_parts, &parts_sizes);

	free(block_mem);
	for (int i=0; i<100; i++) free(block_mem_aux[i]);

	for (int i=0; i<100; i++)
        {
                get_sequences_from_block(tumor_blocks[i], &block_mem_aux[i], &block_size_aux[i]);
        }

        fmem = open_memstream(&block_mem, &block_size);
        gsuffix_me(fmem, 100, block_mem_aux, block_size_aux, 1);
        fb = fmemopen(block_mem, block_size, "rb");
        int tsize = gsuffix_fget_size(fb);
        free(block_mem);
        for (int i=0; i<100; i++) free(block_mem_aux[i]);

	std::vector<TreeBlock> control_blocks_nodes[NUM_NODES];
        std::vector<TreeBlock> tumor_blocks_nodes[NUM_NODES];

        for (int i=0; i<control_blocks.size(); i++)
                control_blocks_nodes[control_blocks[i]->_num_node-1].push_back(*control_blocks[i]);
        for (int i=0; i<tumor_blocks.size(); i++)
                tumor_blocks_nodes[tumor_blocks[i]->_num_node-1].push_back(*tumor_blocks[i]);

	std::vector<MutationBlock*> breakpoints;

	for (int i=1; i<=NUM_NODES; i++)
		for (int j=0; j<1; j++)
			free_nodes.push(std::pair<int, int>(i, j));

	int BLOCK_SIZE = 10000000;
	TreePartition part;
	part._part_id = 0;
	part._offset_ini = 0;
	part._offset_end = BLOCK_SIZE - 1;
	int block_id = 0;
	
	while (part._part_id < num_parts)
	{
		std::pair<int, int> node_cpu;
		MutationBlock* block;
		while ((block = recv_mutations_detected_block(&node_cpu.first, &node_cpu.second)))
		{
			breakpoints.push_back(block);
			free_nodes.push(node_cpu);
		}
		if (!free_nodes.empty())
		{
			std::pair<int, int> next_node = free_nodes.front();
			int nodes_with_blocks = 0;
			for (int i=0; i<NUM_NODES; i++)
			{
				if (control_blocks_nodes[i].size() > 0)
				{
                        		MPI_Send(&control_blocks_nodes[i][0], control_blocks_nodes[i].size()*sizeof(TreeBlock), MPI_CHAR,  i+1, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD);
					MPI_Send(&part, sizeof(TreePartition), MPI_CHAR, i+1, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD);
					nodes_with_blocks++;
				}
			}
			char* partition_mems[NUM_NODES];
			size_t partition_sizes[NUM_NODES];

			for (int i=0; i<nodes_with_blocks; i++)
			{
				MPI_Probe(MPI_ANY_SOURCE, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
        			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
        			partition_mems[i] = new char[num_elems_recv];
				partition_sizes[i] = num_elems_recv;
        			MPI_Recv(partition_mems[i], partition_sizes[i], MPI_CHAR, status.MPI_SOURCE, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
			}

			char* control_merged;
			size_t size_control_merged;
			FILE* fmerged = open_memstream(&control_merged, &size_control_merged);
			gsuffix_me(fmerged, nodes_with_blocks, partition_mems, partition_sizes, 1);

			for (int i=0; i<nodes_with_blocks; i++)
				delete[] partition_mems[i];
			char* control_index;
			size_t size_control_index;
			FILE* findex = open_memstream(&control_index, &size_control_index);
			gsuffix_fi(control_merged, size_control_merged, findex);

			nodes_with_blocks = 0;
			for (int i=0; i<NUM_NODES; i++)
			{
				if (tumor_blocks_nodes[i].size() > 0)
				{
                        		MPI_Send(&tumor_blocks_nodes[i][0], tumor_blocks_nodes[i].size()*sizeof(TreeBlock), MPI_CHAR,  i+1, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD);
					MPI_Send(&part, sizeof(TreePartition), MPI_CHAR, i+1, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD);
					nodes_with_blocks++;
				}
			}
			for (int i=0; i<nodes_with_blocks; i++)
			{
				MPI_Probe(MPI_ANY_SOURCE, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
        			MPI_Get_count(&status, MPI_CHAR, &num_elems_recv);
        			partition_mems[i] = new char[num_elems_recv];
				partition_sizes[i] = num_elems_recv;
        			MPI_Recv(partition_mems[i], partition_sizes[i], MPI_CHAR, status.MPI_SOURCE, MPI_GET_SEQS_FROM_PARTITION, MPI_COMM_WORLD, &status);
			}
			char* tumor_merged;
			size_t size_tumor_merged;
			fmerged = open_memstream(&tumor_merged, &size_tumor_merged);
			gsuffix_me(fmerged, nodes_with_blocks, partition_mems, partition_sizes, 1);
			for (int i=0; i<nodes_with_blocks; i++)
				delete[] partition_mems[i];
			char* tumor_index;
			size_t size_tumor_index;
			findex = open_memstream(&tumor_index, &size_tumor_index);
			gsuffix_fi(tumor_merged, size_tumor_merged, findex);
			MPI_Send(&next_node.second, 1, MPI_INT,  next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
                        MPI_Send(&block_id, 1, MPI_INT,  next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
                        MPI_Send(ref_file, strlen(ref_file)+1, MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
                        MPI_Send(&part, sizeof(TreePartition), MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
                        MPI_Send(&nsize, 1, MPI_INT, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
                        MPI_Send(&tsize, 1, MPI_INT, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
			MPI_Send(control_merged, size_control_merged, MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
			MPI_Send(control_index, size_control_index, MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
			MPI_Send(tumor_merged, size_tumor_merged, MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);
			MPI_Send(tumor_index, size_tumor_index, MPI_CHAR, next_node.first, MPI_GET_MUTS_FROM_BLOCK_FASTQ, MPI_COMM_WORLD);

			block_id++;
			free_nodes.pop();
			part._offset_ini += BLOCK_SIZE;
			if (part._offset_ini >= parts_sizes[part._part_id])
			{
				part._offset_ini = 1;
				part._part_id++;
			}
			part._offset_end = part._offset_ini + BLOCK_SIZE - 1;
		}
	}

	while (free_nodes.size() < NUM_NODES)
	{
		std::pair<int, int> node_cpu;
		MutationBlock* block;
		block = recv_mutations_detected_block(&node_cpu.first, &node_cpu.second);
		if (block)
		{
			breakpoints.push_back(block);
			free_nodes.push(node_cpu);
		}
	}
	
	for (int i=0; i<num_parts; i++)
		free(parts[i]);
	free(parts);
	free(parts_sizes);

	return breakpoints;
}

bool sort_breakpoints(MutationBlock* a, MutationBlock* b)
{
	if (a->_block_id < b->_block_id) return true;
	return false;
}

int master(int argc, char* argv[])
{      
  
        struct rusage master_rs;
        getrusage(RUSAGE_SELF, &master_rs);
	FILE* stdout_backup = stdout;
	FILE* stdin_backup = stdin;
	FILE* stderr_backup = stderr;
	FILE* stdnull = fopen("/dev/null", "w");

        int c;
	char* ref_file = 0;
	char* normal_bam_file = 0;
	char* cancer_bam_file = 0;
	char* normal_fastq_1_file = 0;
	char* normal_fastq_2_file = 0;
	char* cancer_fastq_1_file = 0;
	char* cancer_fastq_2_file = 0;
	char* mem_normal_bam_sorted = 0;
	char* brkp_ext_file = 0;
	size_t size_normal_bam_sorted;
	char* mem_normal_bai;
	size_t size_normal_bai;
	char* mem_cancer_bam_sorted;
	size_t size_cancer_bam_sorted;
	char* mem_cancer_bai;
	size_t size_cancer_bai;
	int cpus_per_node_spec = 0;

	FILE* fph_brkp_ext = 0;

	while (1)
	 {
	   static const struct option long_options[] =
	     {
	       {"ref",     required_argument, 0, 'r'},
	       {"normal_bam",  required_argument,     &use_normal_bam, 1},
	       {"tumor_bam",  required_argument,     &use_cancer_bam, 1},
	       {"normal_fastq_1",  required_argument,     &use_normal_bam, 0},
	       {"tumor_fastq_1",  required_argument,     &use_cancer_bam, 0},
	       {"normal_fastq_2",  required_argument,     &use_normal_bam, 0},
	       {"tumor_fastq_2",  required_argument,     &use_cancer_bam, 0},
	       {"min_sup_reads", required_argument,	0, 0},
	       {"print_brkp_ext",  required_argument,     &print_brkp_ext, 1},
	       {"chr",  required_argument,     0, 0},
	       {"pos",  required_argument,     0, 0},
               {"patient_id",  required_argument,     0, 0},
	       {"cpus_per_node",    required_argument, &cpus_per_node_spec, 1},
	       {0, 0, 0, 0}
	     };
	   /* getopt_long stores the option index here. */
	   int option_index = 0;

	   c = getopt_long (argc, argv, "r:n:c:r:",
		            long_options, &option_index);

	   /* Detect the end of the options. */
	   if (c == -1)
	     break;
	   switch (c)
	     {
	     case 0:
	       /* If this option set a flag, do nothing else now. */
	       //if (long_options[option_index].flag != 0)
	       // break;
	       
	       if (strcmp(long_options[option_index].name, "normal_bam") == 0)
	       {
	       		if (optarg)
	       		{
				normal_bam_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "tumor_bam") == 0)
	       {
	       		if (optarg)
	       		{
				cancer_bam_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "normal_fastq_1") == 0)
	       {
	       		if (optarg)
	       		{
				normal_fastq_1_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "normal_fastq_2") == 0)
	       {
	       		if (optarg)
	       		{
				normal_fastq_2_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "tumor_fastq_1") == 0)
	       {
	       		if (optarg)
	       		{
				cancer_fastq_1_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "tumor_fastq_2") == 0)
	       {
	       		if (optarg)
	       		{
				cancer_fastq_2_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "print_brkp_ext") == 0)
	       {
	       		if (optarg)
	       		{
				brkp_ext_file = strdup(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "min_sup_reads") == 0)
	       {
	       		if (optarg)
	       		{
				min_sup_reads = atoi(optarg);
			}
	       }
	       if (strcmp(long_options[option_index].name, "chr") == 0)
               {
                        if (optarg)
                        {
                                chr_id = atoi(optarg);
                        }
               }
	       if (strcmp(long_options[option_index].name, "pos") == 0)
               {
                        if (optarg)
                        {
                                pos = atoi(optarg);
                        }
               }
	       if (strcmp(long_options[option_index].name, "cpus_per_node") == 0)
               {
                        if (optarg)
                        {
                                CPUS_PER_NODE = atoi(optarg);
                        }
               }
	       if (strcmp(long_options[option_index].name, "patient_id") == 0)
               {
                        if (optarg)
                        {
                                patient_id = strdup(optarg);
                        }
               }
	       break;

	     case 'r':
	       if (!optarg)
			return usage(stderr_backup);
	       ref_file = strdup(optarg);
	       break;
	     
	     case 't':
		if (is_int(optarg))
		{
			num_threads = atoi(optarg);
		}
		else
			return usage(stderr_backup);
	       break;

	     case '?':
	       /* getopt_long already printed an error message. */
	       break;

	     default:
	     fprintf(stderr, "Aborta!!!\n");
	       abort ();
	     }
	 }

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	 {
	   fprintf (stderr_backup, "non-option ARGV-elements: ");
	   while (optind < argc)
	     fprintf (stderr_backup, "%s ", argv[optind++]);
	   fprintf(stderr_backup, "\n");
	 }

	if (!ref_file || !cpus_per_node_spec ||
	    !use_normal_bam && (!normal_fastq_1_file || !normal_fastq_2_file) || 
	    !use_cancer_bam && (!cancer_fastq_1_file || !cancer_fastq_2_file) || 
	    use_normal_bam && !normal_bam_file || 
	    use_cancer_bam && !cancer_bam_file
	    )
		return usage(stderr_backup);

	
	/* Checking if BWT index was previously build for reference genome */
	int gen_bwt_index = 0;

	char* str_buff = (char*)malloc(FILENAME_MAX);
	strcpy(str_buff, ref_file);
	strcat(str_buff, ".bwt");
	if (access(str_buff, R_OK | W_OK))
		gen_bwt_index = 1;
	strcpy(str_buff, ref_file);
	strcat(str_buff, ".amb");
	if (access(str_buff, R_OK | W_OK))
		gen_bwt_index = 1;
	strcpy(str_buff, ref_file);
	strcat(str_buff, ".ann");
	if (access(str_buff, R_OK | W_OK))
		gen_bwt_index = 1;
	strcpy(str_buff, ref_file);
	strcat(str_buff, ".pac");
	if (access(str_buff, R_OK | W_OK))
		gen_bwt_index = 1;
	strcpy(str_buff, ref_file);
	strcat(str_buff, ".sa");
	if (access(str_buff, R_OK | W_OK))
		gen_bwt_index = 1;


	fprintf(stderr_backup, "Checking input data files...\n");
	if (gen_bwt_index)
	{
		fprintf(stderr_backup, "No BWT index found for %s. Generating index...\n", ref_file);
		if (access(ref_file, R_OK))
		{
			fprintf(stderr_backup, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", ref_file);
			return 1;
		}		
		my_bwa_index(ref_file);
	}
	else
	{
		fprintf(stderr_backup, "BWT index found for %s\n", ref_file);
	}

	if (print_brkp_ext)
	{
		/*		
		if (access(brkp_ext_file, W_OK))
		{
			fprintf(stderr_backup, "It was not possible to create %s. Check if you entered the correct path to the file and you have the write permissions.\n", brkp_ext_file);
			exit(1);
		}
		*/
		fph_brkp_ext = fopen(brkp_ext_file, "w");
	}

if (!use_normal_bam)
{
	if (access(normal_fastq_1_file, R_OK))
	{ 
		fprintf(stderr_backup, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", normal_fastq_1_file);
			return 1;
	}
	else
	{
		if (!check_fastqs_access(normal_fastq_1_file, stderr_backup))
			return 1;
		
	}
	if (access(normal_fastq_2_file, R_OK))
	{ 
		fprintf(stderr_backup, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", normal_fastq_2_file);
			return 1;
	}
	else
	{
		if (!check_fastqs_access(normal_fastq_2_file, stderr_backup))
			return 1;
		
	}

}
if (!use_cancer_bam)
{
	if (access(cancer_fastq_1_file, R_OK))
	{ 
		fprintf(stderr_backup, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", cancer_fastq_1_file);
			return 1;
	}
	else
	{
		if (!check_fastqs_access(cancer_fastq_1_file, stderr_backup))
			return 1;
		
	}

	if (access(cancer_fastq_2_file, R_OK))
	{ 
		fprintf(stderr_backup, "It was not possible to read %s. Check if you entered the correct path to the file and you have the read permissions.\n", cancer_fastq_2_file);
			return 1;	
	}
	else
	{
		if (!check_fastqs_access(cancer_fastq_2_file, stderr_backup))
			return 1;
	}

}

	MPI_Status status;

	for (int i=1; i<=NUM_NODES; i++)
		MPI_Send(ref_file, strlen(ref_file)+1, MPI_CHAR, i, MPI_INIT_BWA_INDX_CMD, MPI_COMM_WORLD);

	for (int i=1; i<=NUM_NODES; i++)
		MPI_Recv(0, 0, MPI_CHAR, MPI_ANY_SOURCE, MPI_INIT_BWA_INDX_CMD, MPI_COMM_WORLD, &status);

        for (int i=1; i<=NUM_NODES; i++)
                MPI_Send(&min_sup_reads, 1, MPI_INT, i, MPI_SET_SUPP_READS_CMD, MPI_COMM_WORLD);

        for (int i=1; i<=NUM_NODES; i++)
                MPI_Recv(0, 0, MPI_CHAR, MPI_ANY_SOURCE, MPI_SET_SUPP_READS_CMD, MPI_COMM_WORLD, &status);

	std::vector<TreeBlock*> normal_tree_blocks;
	std::vector<TreeBlock*> tumor_tree_blocks;

	char** chrs;
        int* chr_sizes;
        int num_chrs;
	if (!use_normal_bam)
		generate_quadtree(normal_fastq_1_file, normal_fastq_2_file, normal_tree_blocks);

	if (!use_cancer_bam)
                generate_quadtree(cancer_fastq_1_file, cancer_fastq_2_file, tumor_tree_blocks);

	std::vector<MutationBlock*> breakpoints_detected;


	breakpoints_detected = get_breakpoints(ref_file, normal_tree_blocks, tumor_tree_blocks);
	char* mem_block;
	size_t size_mem_block;
        get_sequences_from_block(normal_tree_blocks[0], &mem_block, &size_mem_block);
	get_partitions_from_header(mem_block, size_mem_block, &chrs, &num_chrs, &chr_sizes);
	delete[] mem_block;

	std::sort(breakpoints_detected.begin(), breakpoints_detected.end(), sort_breakpoints);
	
	char filename_buffer[1000];
	sprintf(filename_buffer, "somatic_SNVs.txt.%s", patient_id);
	FILE* fsnvs = fopen(filename_buffer, "w");
        sprintf(filename_buffer, "somatic_small_SVs.%s", patient_id);
	FILE* findels = fopen(filename_buffer, "w");
        sprintf(filename_buffer, "somatic_large_SVs.%s", patient_id);
	FILE* flarge = fopen(filename_buffer, "w");

        fprintf(fsnvs, "Mut_ID\tType\tChr\tPos\tNormal_NT\tTumor_NT\n");
	fprintf(findels, "Mut_ID\tType\tChr\tPos\tSize\tSequence\n");
        fprintf(flarge, "Mut_ID\tType\tChr_BKP_1\tPos_BKP_1\tChr_BKP_2\tPos_BKP_2\tlocal(left_strand left_ini..left_end,right_strand right_ini..right_end)\tExt_Sequence\n");
	int snv_counter = 1;
	int smallSV_counter = 1;
	int largeSV_counter = 1;
	int extensions_counter = 1;
	
	for (int i=0; i<breakpoints_detected.size(); i++)
	{
		for (int im=0; im<breakpoints_detected[i]->_num_snvs; im++)
                {
                        print_point_mutation(snv_counter, chrs[breakpoints_detected[i]->_snvs[im]._chr_id], breakpoints_detected[i]->_snvs[im]._pos, breakpoints_detected[i]->_snvs[im]._org_base, breakpoints_detected[i]->_snvs[im]._mut_base, breakpoints_detected[i]->_snvs[im]._num_normal_in_normal, breakpoints_detected[i]->_snvs[im]._num_tumor_in_normal, breakpoints_detected[i]->_snvs[im]._num_normal_in_tumor, breakpoints_detected[i]->_snvs[im]._num_tumor_in_tumor, fsnvs);
			snv_counter++;
                }
		for (int im=0; im<breakpoints_detected[i]->_num_dels; im++)
		{
			print_deletion(smallSV_counter, chrs[breakpoints_detected[i]->_dels[im]._chr_id], breakpoints_detected[i]->_dels[im]._pos, breakpoints_detected[i]->_dels[im]._size, 0, breakpoints_detected[i]->_dels[im]._num_normal_in_normal, breakpoints_detected[i]->_dels[im]._num_tumor_in_normal, breakpoints_detected[i]->_dels[im]._num_normal_in_tumor, breakpoints_detected[i]->_dels[im]._num_tumor_in_tumor, findels);
			smallSV_counter++;
		}
		for (int im=0; im<breakpoints_detected[i]->_num_ins; im++)
                {
                        print_short_insertion(smallSV_counter, chrs[breakpoints_detected[i]->_ins[im]._chr_id], breakpoints_detected[i]->_ins[im]._pos, breakpoints_detected[i]->_ins[im]._size, breakpoints_detected[i]->_ins[im]._sequence, breakpoints_detected[i]->_ins[im]._num_normal_in_normal, breakpoints_detected[i]->_ins[im]._num_tumor_in_normal, breakpoints_detected[i]->_ins[im]._num_normal_in_tumor, breakpoints_detected[i]->_ins[im]._num_tumor_in_tumor, findels);
			smallSV_counter++;
                }
		for (int im=0; im<breakpoints_detected[i]->_num_invs; im++)
                {
                        print_inversion(smallSV_counter, chrs[breakpoints_detected[i]->_invs[im]._chr_id], breakpoints_detected[i]->_invs[im]._pos, breakpoints_detected[i]->_invs[im]._size, 0, 0, breakpoints_detected[i]->_invs[im]._num_normal_in_normal, breakpoints_detected[i]->_invs[im]._num_tumor_in_normal, breakpoints_detected[i]->_invs[im]._num_normal_in_tumor, breakpoints_detected[i]->_invs[im]._num_tumor_in_tumor, findels);
			smallSV_counter++;
                }
		for (int im=0; im<breakpoints_detected[i]->_num_bps; im++)
                {
			fprintf(flarge, "%d\tBKP\t%s\t%d\t%s\t%d\tlocal(%c%d..%d,%c%d..%d)\t%s\n", largeSV_counter, chrs[breakpoints_detected[i]->_bps[im]._chr_id_1], breakpoints_detected[i]->_bps[im]._pos_1, chrs[breakpoints_detected[i]->_bps[im]._chr_id_2], breakpoints_detected[i]->_bps[im]._pos_2, breakpoints_detected[i]->_bps[im]._strand_1, breakpoints_detected[i]->_bps[im]._q1_ini, breakpoints_detected[i]->_bps[im]._q1_end, breakpoints_detected[i]->_bps[im]._strand_2, breakpoints_detected[i]->_bps[im]._q2_ini, breakpoints_detected[i]->_bps[im]._q2_end, breakpoints_detected[i]->_bps[im]._extension);
			largeSV_counter++;
                }
	}
	fclose(fsnvs);
	fclose(findels);
	fclose(flarge);
        fprintf(stderr, "Finish!!!\n");
        getrusage(RUSAGE_SELF, &master_rs);
        std::cout <<  " Master res usage : " << master_rs.ru_maxrss << std::endl;
	
	return 0;
}

int main(int argc, char* argv[])
{
  CSTART(mainly);
  int myrank;

  int size_node_name;
  char node_name[1000];
  
  /* Initialize MPI */

  MPI_Init(&argc, &argv);

  buffer = new char[BUFFER_SIZE];

  MPI_Get_processor_name(node_name, &size_node_name);
  MPI_Comm_size(MPI_COMM_WORLD, &NUM_NODES);
  NUM_NODES--;

  /* Find out my identity in the default communicator */

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0) {
    master(argc, argv);
    /* Sending Ending Message to all Slaves */
    for (int i=1; i<=NUM_NODES; i++)
	MPI_Send(0, 0, MPI_INT, i, MPI_EXIT_CMD, MPI_COMM_WORLD);
  } else {
    slave(myrank);
  }

  delete[] buffer;

  /* Shut down MPI */

  MPI_Finalize();
  CEND(mainly);
  CTIME(mainly);
  CPRINT(mainly);
  return 0;
}
