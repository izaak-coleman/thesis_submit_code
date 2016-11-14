#include <map>

typedef struct _TreeBlock {

	int	_num_node;
	int	_num_cpu;
	char*	_block_address;
	size_t	_block_address_size;
	char*	_iblock_address;
	size_t	_iblock_address_size;

} TreeBlock;

typedef struct _TreePartition {

	int   _part_id;
	int   _offset_ini;
	int   _offset_end;

} TreePartition;

typedef struct _SNV {

	int _chr_id;
	int _pos;
	char _org_base;
	char _mut_base;
	int _num_normal_in_normal;
	int _num_tumor_in_normal;
	int _num_normal_in_tumor;
	int _num_tumor_in_tumor;

} SNV;

typedef struct _DEL {
	
	int	_chr_id;
	int	_pos;
	int	_size;
	int 	_num_normal_in_normal;
        int 	_num_tumor_in_normal;
        int 	_num_normal_in_tumor;
        int 	_num_tumor_in_tumor;
} DEL;

typedef struct _INS {

	int	_chr_id;
	int	_pos;
	int	_size;
	char	_sequence[1000];
	int 	_num_normal_in_normal;
        int 	_num_tumor_in_normal;
        int 	_num_normal_in_tumor;
        int 	_num_tumor_in_tumor;
} INS;

typedef struct _INV {

	int	_chr_id;
	int	_pos;
	int	_size;
	int 	_num_normal_in_normal;
        int 	_num_tumor_in_normal;
        int 	_num_normal_in_tumor;
        int 	_num_tumor_in_tumor;
} INV;

typedef struct _BKP {

	int	_ext_id;
	int	_chr_id_1;
	int	_chr_id_2;
	int	_pos_1;
	int	_pos_2;
	int	_q1_ini;
	int	_q1_end;
	int	_q2_ini;
	int	_q2_end;
	char	_strand_1;
	char	_strand_2;
	char	_extension[10000];
	int	_cov;

} BKP;

typedef struct _BKPE {

        char    _extension[10000];
	int	_loffset;
        int     _cov;

} BKPE;

typedef struct _MutationBlock {

	int	_block_id;
	int	_num_snvs;
	int	_num_dels;
	int	_num_ins;
	int	_num_invs;
	int	_num_bps;
	int	_num_bpse;
	SNV*	_snvs;
	DEL*	_dels;
	INS*	_ins;
	INV*	_invs;
	BKP*	_bps;
	BKPE*	_bpse;
	
} MutationBlock;

typedef struct _break_point {
        char* _ext_id;
        char* _chr_1;
        char* _chr_2;
        int _pos_1;
        int _pos_2;
        int _pos_1_ini;
        int _pos_1_end;
        int _pos_2_ini;
        int _pos_2_end;
        int _q0, _q1, _q2, _q3;
        char _q1_strand;
        char _q2_strand;
        char* _ext;
        int _num_sup_reads;

} BreakPoint;

typedef struct _break_point_extension {

        int _id;
        char* _chr;
        int _pos;
        char* _ext;
        int _loffset;
        int _num_C;
        int _used;

} BRK_EXT;


typedef struct _gsufix_tree_block {

        int     _task_assigned;
        int     _task_finished;
        int     _cpu;
        char*   _mem_fastq_1;
        size_t  _size_fastq_1;
        char*   _mem_fastq_2;
        size_t  _size_fastq_2;
        
        char*   _mem_block_out;
        size_t  _size_block_out;
        char*   _mem_iblock_out;
        size_t  _size_iblock_out;

        struct gsuffix *gsuf;

        pthread_t _t;

} gsufix_tree_block;

typedef struct _det_muts_block {

        int     _task_assigned;
        int     _task_finished;
        int     _block_id;
        int     _cpu_id;

        char*   _ref_genome_filename;
        char*   _c;
        int     _pini;
        int     _pend;
        char*   _nblock_file;
        char*   _tblock_file;
        int     _nsize;
        int     _tsize;
        TreePartition   _partition;
        char*   _mem_control_part;
        size_t  _size_control_part;
        char*   _mem_control_index;
        size_t  _size_control_index;
        char*   _mem_tumor_part;
        size_t  _size_tumor_part;
        char*   _mem_tumor_index;
        size_t  _size_tumor_index;

        std::vector<SNV>* _vsnvs;
        std::vector<DEL>* _vdels;
        std::vector<INS>* _vins;
        std::vector<INV>* _vinvs;
        std::vector<BKP>* _vbps;
        std::vector<BKPE>* _vbpse;

        pthread_t _t;

} det_muts_block;

typedef struct _det_params {

        char* _mem_normal_reads;
        size_t _size_normal_reads;
        char* _mem_tumor_reads;
        size_t _size_tumor_reads;
        char* _mem_inormal;
        size_t _size_inormal;
        char* _mem_itumor;
        size_t _size_itumor;
        char* _ref_filename;
        char* _c;
        int _pini;
        int _pend;
        int _nsize;
        int _tsize;
        std::map<std::string, int>* _hash_cs;

        std::vector<DEL>* _vdels;
        std::vector<INS>* _vins;
        std::vector<INV>* _vinvs;
        std::vector<BKP>* _vbps;
        std::vector<BKPE>* _vbpse;

        int _block_id;

} det_params;

