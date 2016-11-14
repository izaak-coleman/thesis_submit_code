#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "kseq.h"
KSEQ_DECLARE(FILE*)
extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

int main_mem(int argc, char *argv[])
{
}

bwaidx_t *idx;

int my_main_mem(char* ref_file, char* seqs_mem, size_t size_seqs_mem, FILE* fout)
{

        mem_opt_t *opt;
        int fd, fd2, i, c, n, copy_comment = 0;
        FILE* fp, *fp2 = 0;
        kseq_t *ks, *ks2 = 0;
        bseq1_t *seqs;
        //bwaidx_t *idx;
        char *rg_line = 0;
        void *ko = 0, *ko2 = 0;
        
        fp = fmemopen(seqs_mem, size_seqs_mem, "rb");
	if (!fp) return 0;
	opt = mem_opt_init();
        opt->flag |= MEM_F_ALL;
        if (opt->n_threads < 1) opt->n_threads = 1;

        bwa_fill_scmat(opt->a, opt->b, opt->mat);
        //if ((idx = bwa_idx_load(ref_file, BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
        bwa_print_sam_hdr(idx->bns, rg_line, fout);

        ks = kseq_init(fp);
        while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
                int64_t size = 0;
                if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
                        if (bwa_verbose >= 2)
                                fprintf(stderr, "");
                        n = n>>1<<1;
                }
                if (!copy_comment)
                        for (i = 0; i < n; ++i) {
                                free(seqs[i].comment); seqs[i].comment = 0;
                        }
                for (i = 0; i < n; ++i) size += seqs[i].l_seq;
                if (bwa_verbose >= 3)
                        fprintf(stderr, "");
                mem_process_seqs_fout(opt, idx->bwt, idx->bns, idx->pac, n, seqs, fout);
                free(seqs);
        }

        free(opt);
        //bwa_idx_destroy(idx);
        kseq_destroy(ks);
        fclose(fp);
        if (ks2) {
                kseq_destroy(ks2);
                fclose(fp2);
        }
        return 0;
}


void init_bwa_idx(char* ref_file)
{
	idx = bwa_idx_load(ref_file, BWA_IDX_ALL);
}

void destroy_bwa_idx()
{
	bwa_idx_destroy(idx);
}

int my_main_mem_paired(char* mem_fastq_1, size_t size_mem_fastq_1, char* mem_fastq_2, size_t size_mem_fastq_2, FILE* fout)
{
	mem_opt_t *opt;
	int fd, fd2, i, c, n, copy_comment = 0;
	FILE* fp, *fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	char *rg_line = 0;
	
	fp = fmemopen(mem_fastq_1, size_mem_fastq_1, "r"); 
	fp2 = fmemopen(mem_fastq_2, size_mem_fastq_2, "r");

	opt = mem_opt_init();
	opt->flag |= MEM_F_ALL;
	if (opt->n_threads < 1) opt->n_threads = 1;
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	bwa_print_sam_hdr(idx->bns, rg_line, fout);
	ks = kseq_init(fp);
        ks2 = kseq_init(fp2);
        opt->flag |= MEM_F_PE;

	while ((seqs = bseq_read(opt->chunk_size * opt->n_threads, &n, ks, ks2)) != 0) {
		int64_t size = 0;
		if ((opt->flag & MEM_F_PE) && (n&1) == 1) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "");
			n = n>>1<<1;
		}
		if (!copy_comment)
			for (i = 0; i < n; ++i) {
				free(seqs[i].comment); seqs[i].comment = 0;
			}
		for (i = 0; i < n; ++i) size += seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "");
		mem_process_seqs_fout(opt, idx->bwt, idx->bns, idx->pac, n, seqs, fout);
		free(seqs);
	}

	free(opt);
	//bwa_idx_destroy(idx);
	kseq_destroy(ks);
	fclose(fp);
	if (ks2) {
		kseq_destroy(ks2);
		fclose(fp2);
	}
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
}
