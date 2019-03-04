//
//  assess_fold.c
//  
//
//  Created by Eric Strobel on 8/9/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define BLOCKSIZE 30000	//block size for memory allocation
#define MAXLINE 2048	//maximum line length
#define MAXSEQS 2048	//maximum number of sequences to evaluate
#define MAXLEN 64		//maximum input sequence length

/* randomization mode values; 0 = no randomization*/
#define RNDM 1
#define CNSNS 2
#define NBAL 3

/* base definitions */
#define A 0
#define G 1
#define U 2
#define C 3
#define TOT 4

/* errors */
#define EMPTY_CT_LINE -1

/* structure declarations */
struct fasta_entry {			//Structure to store fasta entries
    char * nm;					//sequence name
    char * sq;					//sequence with multiple sequence alignment characters preserved
    int len;					//sequence length
};

struct RNA {					//structure to hold information about individual RNA sequences
    struct fasta_entry *fa;		//pointer to fasta entry
    char *sq;					//RNA sequence without alignment characters
    char typ;					//flag to distinguish sequences sent to eval_ct from main
    int s_cnt;					//total number of structures
    struct bp_tbl *bpt;			//pointer to root base pair table
    struct bp_tbl *cr_bpt;		//pointer to current base pair table
};

struct bp_tbl {					//structure containing information about a predicted RNA structure
    int *bp;					//values describing predicted RNA structure
    float dg;					//deltaG of predicted RNA structure
    int numpr;					//number of pairs in predicted structure
    int mxhel;					//longest contiguous helix
    int gcp;					//number of G:C pairs
    int aup;					//number of A:U pairs
    int gup;					//number of G:U pairs
    struct bp_tbl *bpt;			//pointer to next base pair table for sequences with multiple structure predictions
};

struct bp_tbl_bank {			//bank for dynamic allocation of bp_tbl memory
    struct bp_tbl *bpt;			//pointer for bp_tbl allocation
    struct bp_tbl_bank *nxt;	//pointer to next bp_tbl_bank
    int count;					//number of allocated bp_tbls from bank currently in use
};

struct deltaG_frequency {		//structure for tracking deltaG counts
    float val;					//deltaG value
    int cnt;					//number of times deltaG value has been observed.
};

/* function prototypes */
// file management functions
int mk_out_dir();                        														//make output directory
int read_fasta(struct fasta_entry *fa_ipt, char *ipt); 											//open input fasta file
int get_nt_freq(struct fasta_entry * fa_ipt, float (*base_cnt)[5]); 							//get positional nt frequency
int mk_rnd_tbl(float (*base_cnt)[5], int rnd_tbl[MAXLEN][4], int mode);							//make randomization table
int mk_ipt_fa_rnd(struct fasta_entry * fa_ipt, int rnd_tbl[MAXLEN][4], char * rnd_filename);	//make fasta file using randomization table
int mk_tmp_fa(struct fasta_entry * fa_ipt, struct RNA *RNAvar);									//make fasta input for RNAstructure
int isbase(char c);                        														//test for IUPAC base notation

// ct file functions
int eval_ct(struct RNA *RNAvar, struct bp_tbl_bank *crrnt_bpt_bank);							//evaluate ct file
int parse_ct_header(char *line, struct RNA *RNAvar);											//parse ct file header line
int parse_ct_seq(char *line, struct RNA *RNAvar, int i);										//parse ct file sequence lines

// metrics functions
int test_dg_uniq(float dg, struct deltaG_frequency * dg_freq_tbl);								//track deltaG distribution
int test_hel_contig(struct RNA *RNAvar, int hel, int i);										//test for contiguous helix
int identify_bps(struct RNA *RNAvar);															//track base pair distribution
int ispair(char nt1, char nt2);																	//test for base pair

int print_output(struct RNA *RNAvar);
/* end of function prototypes */

/* global variables */
// options variables
int rndm = 0;                   				//randomization mode, default is 0 for no randomization
int verbose = 0;								//flag for verbose output
char consensus[MAXLEN] = {0};                    //user supplied consensus sequence for randomization mode 2

// file name variables
char out_dir[MAXLINE] = {0};    				//output directory
char filename[MAXLINE] = {0};					//input fasta prefix
int run_id = 0;             					//unique run identifier

// totals variables
int max_seq_len = 0;							//maximum observed sequence length
int max_seq_indx = 0;							//maximum observed sequence index
int total_variants = 0;         				//total sequences analyzed
int total_structures = 0;						//total structures predicted

// metrics variables
int len_dist[MAXSEQS] = {0};            		//distribution of input sequence lengths
int bp_cnt[MAXLEN] = {0}; 						//tracks base pair distribution
int pr_cnt_tbl[MAXLEN][MAXLEN] = {0};			//matrix of base pair distributions
int tot_pr_cnt[MAXLEN] = {0};					//distribution of pairs observed at each sequence position
int tot_obs_cnt[MAXLEN] = {0};
int helix_contig[MAXLEN] = {0};					//distribution of contiguous helix lengths
struct deltaG_frequency dg_freq[MAXSEQS] = {0}; //main deltaG frequency table
/* end of global variables */

int main(int argc, char *argv[])
{
    extern int total_variants;
    extern int total_structures;
    extern int rndm;
    extern int verbose;
    extern char consensus[MAXLEN];
    
    int i = 0;
    int j = 0;
    
    struct fasta_entry * fa_ipt  =  NULL;

    if ((fa_ipt = calloc(MAXSEQS, sizeof(*fa_ipt))) == NULL) {
        printf("main: error - fa_ipt memory allocation failed\n");
        return 0;
    }
    
    
    for (i = 0; i < MAXLEN; i++) {
        tot_obs_cnt[i] = 0;
    }
    
    
    /* parse options using getopt_long */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"sequences",       required_argument,  0,  's'},	//input fasta file
            {"randomize", 		required_argument,  0,  'r'},	//randomization option
            {"consensus",       required_argument,  0,  'c'},	//consensus sequence input for randomization mode 2
            {"verbose",			no_argument,		0,	'v'},	//flag for verbose output
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "s:r:c:v", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 's':													//input fasta file
                total_variants = read_fasta(fa_ipt, argv[optind-1]);	//get input fasta
                for (i = 0; argv[optind-1][i] != '.'; i++) { 			//get input filename for output files TODO:grab up to last '.'
                    filename[i] = argv[optind-1][i];
                }
                filename[i] = '\0';
                break;
            case 'r':													//randomization option
                if ((rndm = atoi(argv[optind-1])) > 3 || rndm == 0) { 	//TODO add specific errors for 0 or >2 conditions
                    printf("error: main - invalid value (%d) for option r\n", rndm);
                }
                break; //TODO: add error if atoi fails
            case 'c': strcpy(consensus, argv[optind-1]); break;			//set consensus for randomization mode 2
            case 'v': verbose = 1; break;										//set verbose output on
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    printf("%s\n", filename);					//print input file name to screen
    printf("randomization mode %d\n", rndm);	//print randomization mode to screen
    /* end of options parsing */
    
    srand(time(NULL));		//set randomization starting point
    run_id = rand();   		//generate random run id
    printf("%X\n", run_id);	//print run id to screen
    mk_out_dir(); 			//make output directory
    
    /* if randomization option used, construct randomized sequences */
    char rnd_filename[MAXLINE] = {0};		//file name prefix for randomization table and fasta outputs
    float base_cnt[MAXLEN][5] = {0};    	//number of times each base is observed at each position
    int rnd_tbl[MAXLEN][4] = {0};			//distribution of values used for randomization
    struct fasta_entry * rnd_ipt =  NULL;	//pointer for randomized fasta entries
    if (rndm) {
        get_nt_freq(fa_ipt, base_cnt);         							//get positional nucleotide frequency
        mk_rnd_tbl(base_cnt, rnd_tbl, rndm); 							//make randomization table
        mk_ipt_fa_rnd(fa_ipt, rnd_tbl, rnd_filename);					//use randomization table to make input fasta
        
        printf("returned from making randomization fasta\n");
        if ((rnd_ipt = calloc(MAXSEQS, sizeof(*rnd_ipt))) == NULL) {	//allocate memory for randomized fasta
            printf("main: error - fa_ipt memory allocation failed\n");
            return 0;
        }
        printf("attempting to read randomization fasta\n");
        if	(read_fasta(rnd_ipt, rnd_filename) != total_variants) {		//open randomized input fasta
            printf("main: error - variant count when making randomized fasta does not equal input variant count");
        }
    }
    /* end of randomized sequence construction */
    
    /* initialize structure metrics memory */
    struct bp_tbl_bank bpt_bank; 					//bank for bp_tbl structures TODO: add function that allows for expansion of bank
    struct bp_tbl_bank *crrnt_bpt_bank = &bpt_bank;	//pointer to current bp_tbl bank
    bpt_bank.count = 0;
    if ((bpt_bank.bpt = calloc(BLOCKSIZE, sizeof(*(bpt_bank.bpt)))) == NULL) {
        printf("main: error - bp_tbl_bank memory allocation failed\n");
        return 0;
    }
    
    struct RNA *RNAvar; //RNAvar stores extracted ct file data and useful values for runtime
    if ((RNAvar = calloc(BLOCKSIZE, sizeof(*RNAvar))) == NULL) {
        printf("main: error - RNAvar memory allocation failed\n");
        return 0;
    }
    /* end of structure metrics memory initialization */
    
    struct fasta_entry * fa= NULL;	//pointer to fasta entries used for structure prediction
    if (rndm) {
        fa = rnd_ipt;	//predict structures from randomized fasta
    } else {
        fa = fa_ipt;	//predict structures from user input fasta
    }
    
    int count = 0;		//number of structures predicted for a given fasta entry
    for(i = 0; i < total_variants; i++) {
        mk_tmp_fa(&fa[i], &RNAvar[i]);
        system("/Users/ericstrobel/Documents/RNAstructure/exe/Fold -mfe ./tmp.fasta ./tmp.ct > out.log"); //call RNAstructure 'Fold'
        if	((count = eval_ct(&RNAvar[i], crrnt_bpt_bank)) == EMPTY_CT_LINE) {	//read and evaluate ct file
            printf("main: error - found blank ct file line\n");
        } else {
            total_structures += count;
        }
        system("rm tmp.fasta tmp.ct out.log");	//remove temporary files
    }
    print_output(RNAvar);
}

/* mk_out_dir: make output directory using input filename and run id */
int mk_out_dir() {
    extern char filename[MAXLINE];	//input fasta prefix
    extern char out_dir[MAXLINE];	//output directory
    extern int run_id;            	//unique run identifier
    
    switch (rndm) {                //make string for output directory;
        case 0: sprintf(out_dir, "%s_wt_%X", filename, run_id); break;
        case 1: sprintf(out_dir, "%s_rnd_%X", filename, run_id); break;
        case 2: sprintf(out_dir, "%s_cnsns_%X", filename, run_id); break;
        case 3: sprintf(out_dir, "%s_nbal_%X", filename, run_id); break;
        default:
            break;
    }
    
    struct stat st = {0};
    
    if (stat(out_dir, &st) == -1) {    //check that output directory does not exist
        mkdir(out_dir, 0755);        //make output directory
    } else {
        printf("mk_out_dir: error - output directory %s already exists\n", out_dir);
    }
    chdir(out_dir);                    //change directory to output directory
    
    return 1;
}

/* read_fasta: open file input and assign to file pointer */
int read_fasta(struct fasta_entry * fa_ipt, char *ipt)
{
    printf("reading fasta\n");
    FILE * ifp = NULL;
    
    //open input fasta
    if ((ifp = fopen(ipt, "r")) == NULL) {
        printf("error: could not open %s as read one file. Aborting program...\n", ipt);
        abort();
    }

    int i = 0;
    int j = 0;
    
    char line[MAXLINE] = {0};	//input line
    char *nm_strt = NULL;		//pointer to bypass '>' in fasta name line
    int line_len = 0;			//length of input line
    int seq_len = 0;			//length of input sequence without multiple sequence alignment characters
    
    for(i = 0; (fgets(line, MAXLINE, ifp)) != NULL && i < MAXSEQS; i++) {
        /* get fasta name line */
        if (line[0] == '>') {					//check that name line begins with '>'
            line_len = strlen(line);			//get line length for memory allocation
            nm_strt = &line[1]; 				//bypass leading '>'
            line[line_len-1] = '\0'; 			//remove trailing newline
            for (j = 0; line[j]; j++) {   		//change pipes to underscores
                if (line[j] == '|') {
                    line[j] = '_';
                }
            }
            fa_ipt[i].nm = malloc((line_len) * sizeof(*(fa_ipt[i].nm)));	//allocate memory for name TODO: check success
            strcpy(fa_ipt[i].nm, nm_strt);									//copy name to fasta entry
        } else {
            printf("read_fasta: error - expected leading '>' on name line\n");
            abort();
        }
        /* end of getting fasta name line */
        
        /* get fasta sequence line and positional nt frequency */
        if ((fgets(line, MAXLINE, ifp)) == NULL) {
            return 0;
        } else if (line[0] != '>') {		//check that line does not begin with '>'
            line_len = strlen(line);		//get line length for memory allocation
            line[line_len-1] = '\0';		//remove trailing newline
            fa_ipt[i].sq = malloc((line_len) * sizeof(fa_ipt[i].sq));    	//allocate memory for sequence TODO: check success
            strcpy(fa_ipt[i].sq, line);										//copy sequence (includes non-base chars);
            for (j = 0, seq_len = 0; line[j] && j < MAXLEN; j++) {			//determine sequence length
                if (isbase(line[j])) {
                    seq_len++;
                }
            }
            fa_ipt[i].len = seq_len;		//copy seq len to fasta entry length
            if (seq_len > max_seq_len) { 	//record maximum sequence length
                max_seq_len = seq_len;
            }
            if (j == MAXLEN && line[j]) {	//test that full sequence was is within MAXLEN bounds
                printf("mk_rnd_tbl: reached MAXLEN without reaching line end\n");
                abort();
            }
            if (j > max_seq_indx) {			//record max index that contains a base
                max_seq_indx = j;
            }
            len_dist[i] = (seq_len) ? seq_len : -1; //len 0 is recorded as -1, needed if input is trimmed
        } else {
            printf("read_fasta: error - sequence line contains leading '>'\n");
            abort();
        }
    }
    
    return i;
}

/* get_nt_freq: get positional nt frequency */
int get_nt_freq(struct fasta_entry * fa_ipt, float (*base_cnt)[5])
{
    printf("getting nucleotide frequencies\n");
    int i = 0;
    int j = 0;
    
    for (i = 0; i < total_variants; i++) {
        for (j = 0; fa_ipt[i].sq[j] && j < MAXLEN; j++) {
            if (isbase(fa_ipt[i].sq[j])) {
                base_cnt[j][TOT]++;
                switch (fa_ipt[i].sq[j]) { //length is only incremented if line[j] != '-'
                    case 'A': base_cnt[j][A]++; break;
                    case 'G': base_cnt[j][G]++; break;
                    case 'U': base_cnt[j][U]++; break;
                    case 'C': base_cnt[j][C]++; break;
                    default:
                        printf("mk_rnd_tbl: error - base %c not currently supported\n", fa_ipt[i].sq[j]);
                        abort();
                        break;
                }
            }
        }
    }
    return 1;
}

/* mk_rnd_tbl: construct randomization table based on randomization mode */
int mk_rnd_tbl(float (*base_cnt)[5], int rnd_tbl[MAXLEN][4], int mode)
{
    printf("making randomization table\n");
    extern int max_seq_indx;    //max observed sequence length
    
    int i = 0;
    int j = 0;
    int len = 0;				//length of current sequence
    
    if (mode == RNDM) {	//unbiased randomization: equal probability for each nucleotide
        for (i = 0; i < max_seq_indx; i++) {
            rnd_tbl[i][A] = 249;
            rnd_tbl[i][G] = 499;
            rnd_tbl[i][U] = 749;
            rnd_tbl[i][C] = 999;
        }
    } else if (mode == CNSNS) {	//consensus randomization: probability biased by input consensus sequence
        if (consensus[0]) {
            for (i = 0; i < max_seq_indx; i++) {
                if (consensus[i]) {
                    switch (consensus[i]) {
                        case 'A': rnd_tbl[i][A] = 999; rnd_tbl[i][G] = 999; rnd_tbl[i][U] = 999; rnd_tbl[i][C] = 999; break;
                        case 'G': rnd_tbl[i][A] = -1;  rnd_tbl[i][G] = 999; rnd_tbl[i][U] = 999; rnd_tbl[i][C] = 999; break;
                        case 'U': rnd_tbl[i][A] = -1;  rnd_tbl[i][G] = -1;  rnd_tbl[i][U] = 999; rnd_tbl[i][C] = 999; break;
                        case 'C': rnd_tbl[i][A] = -1;  rnd_tbl[i][G] = -1;  rnd_tbl[i][U] = -1;  rnd_tbl[i][C] = 999; break;
                        case 'N': rnd_tbl[i][A] = 249; rnd_tbl[i][G] = 499; rnd_tbl[i][U] = 749; rnd_tbl[i][C] = 999; break;
                        case 'R': rnd_tbl[i][A] = 499; rnd_tbl[i][G] = 999; rnd_tbl[i][U] = 999; rnd_tbl[i][C] = 999; break;
                        case 'Y': rnd_tbl[i][A] = -1;  rnd_tbl[i][G] = -1;  rnd_tbl[i][U] = 499; rnd_tbl[i][C] = 999; break;
                        case  0 : rnd_tbl[i][A] = 249; rnd_tbl[i][G] = 499; rnd_tbl[i][U] = 749; rnd_tbl[i][C] = 999; break;
                        default: break;
                    }
                }
            }
        } else {
            printf("mk_rnd_tbl: error - attempting to make consensus randomization table without providing consensus sequence\n");
        }
    } else if (mode == NBAL) {	//frequency randomization: probability biased by observed nt frequency
        printf("max seq len = %d\n", max_seq_len);
        for (i = 0; i < max_seq_indx; i++) {
            if (base_cnt[i][TOT]) {
                rnd_tbl[i][A] = (int)(((base_cnt[i][A]/base_cnt[i][TOT])*1000)-1);
                rnd_tbl[i][G] = rnd_tbl[i][A] + (int)((base_cnt[i][G]/base_cnt[i][TOT])*1000);
                rnd_tbl[i][U] = rnd_tbl[i][G] + (int)((base_cnt[i][U]/base_cnt[i][TOT])*1000);
                rnd_tbl[i][C] = 999; //always set C to max value
            } else {
                rnd_tbl[i][A] = rnd_tbl[i][G] = rnd_tbl[i][U] = rnd_tbl[i][C] = -1;
            }
            
        }
    }
    
    /* output record of randomization table */
    FILE *tmp_fp;
    char rnd_tbl_filename[MAXLEN];
    sprintf(rnd_tbl_filename, "%s_%X_rand_tbl.txt", filename, run_id);
    tmp_fp = fopen(rnd_tbl_filename, "w");
    
    fprintf(tmp_fp, "randomizationt table %X\tmode %d\n", run_id, mode);
    fprintf(tmp_fp, "pos\tA\tG\tU\tC\tTot\t|\tpA\tpG\tpU\tpC\n");
    for (i = 0; i < max_seq_indx; i++) {
        fprintf(tmp_fp, "%d\t", i+1);
        for (j = 0; j < 5; j++) {
            fprintf(tmp_fp, "%3.0f\t", base_cnt[i][j]);
        }
        fprintf(tmp_fp, "|\t");
        for (j = 0; j < 4; j++) {
            fprintf(tmp_fp, "%d\t", rnd_tbl[i][j]);
        }
        fprintf(tmp_fp, "\n");
    }
    
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    return 0;
}

/* mk_ipt_fa_rnd: use randomization table to make fasta file containing randomized sequences */
int mk_ipt_fa_rnd(struct fasta_entry * fa_ipt, int rnd_tbl[MAXLEN][4], char * rnd_filename)
{
    printf("making randomized input fasta\n");
    FILE *tmp_fp = NULL;
    int i = 0;
    int j = 0;
    int k = 0;
	int rnd = 0;
    char rnd_seq[MAXLEN] = {0};
    
    sprintf(rnd_filename, "%s_%X_rnd.fasta", filename, run_id);
    tmp_fp = fopen(rnd_filename, "w");
    
    for (i = 0; i < total_variants; i++) {
        for (j = 0, k = 0; fa_ipt[i].sq[j]; j++) { //TODO add error checking to ensure not leaving table bounds
            rnd = rand() % 1000;
            if (fa_ipt[i].sq[j] != '-') {
                if (rnd <= rnd_tbl[j][0]) {
                    rnd_seq[k++] = 'A';
                } else if (rnd <= rnd_tbl[j][1]) {
                    rnd_seq[k++] = 'G';
                } else if (rnd <= rnd_tbl[j][2]) {
                    rnd_seq[k++] = 'U';
                } else if (rnd <= rnd_tbl[j][3]) {
                    rnd_seq[k++] = 'C';
                } else {
                    printf("mk_ipt_fa_rnd: error - modulus (%d) at postion %d exceeds randomization table bounds (%d-%d)\n", rnd, k+1, rnd_tbl[j][0], rnd_tbl[j][3]);
                }
            } else {
                rnd_seq[k++] = '-';
            }
        }
        rnd_seq[k] = '\0';
        fprintf(tmp_fp, ">Rnd%d\n%s\n", i, rnd_seq);
    }
    
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    printf("made randomization input fasta\n");
    return 1;
}

/* mk_tmp_fa: make fasta file containing single sequence */
int mk_tmp_fa(struct fasta_entry * fa_ipt, struct RNA *RNAvar)
{
    int i = 0;
    int j = 0;
    
    RNAvar->fa = fa_ipt;	//set RNAvar->fa to point to corresponding fasta entry
    RNAvar->typ = 'm';		//set typ to indicate call from main TODO: make this a function argument
    
    RNAvar->sq = malloc((fa_ipt->len+1) * sizeof(*RNAvar->sq));	//allocate sequence memory
    for (i = 0, j = 0; fa_ipt->sq[i]; i++) {					//copy base characters to RNAvar
        if (isbase(fa_ipt->sq[i])) {
            RNAvar->sq[j++] = fa_ipt->sq[i];
        }
    }
    RNAvar->sq[j] = '\0';
    
    /* make single entry fasta file using RNAvar */
    FILE *tmp_fp;
    tmp_fp = fopen("tmp.fasta", "w");
    fprintf(tmp_fp, ">%s\n%s\n", RNAvar->fa->nm, RNAvar->sq);
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    return 1;
}


/* eval_ct: parse ct file and evaluate predicted structures */
int eval_ct(struct RNA *RNAvar, struct bp_tbl_bank *crrnt_bpt_bank)
{
    FILE *fp_ct;

    char line[MAXLINE];
    int i, hel;
    int align_indx = 0;
    
    //printf("opening ct\n");
    if ((fp_ct = fopen("tmp.ct", "r")) == NULL) {
        printf("error: could not open tmp.ct. Aborting program...\n");
        abort();
    }
    
    
    int lcl_bp_cnt[MAXLEN] = {0};
    int structures = 0;
    while ((fgets(line, MAXLINE, fp_ct)) != NULL) {
        for (i = 0; i < MAXLEN; i++) {
            lcl_bp_cnt[i] = 0;
        }
        
        //assign bp_table pointer
        if (!structures) {
            RNAvar->s_cnt++;
            RNAvar->bpt = &crrnt_bpt_bank->bpt[crrnt_bpt_bank->count];
            RNAvar->cr_bpt = &crrnt_bpt_bank->bpt[crrnt_bpt_bank->count++];
            
        } else {
            RNAvar->s_cnt++;
            RNAvar->cr_bpt->bpt = &crrnt_bpt_bank->bpt[crrnt_bpt_bank->count];
            RNAvar->cr_bpt = &crrnt_bpt_bank->bpt[crrnt_bpt_bank->count++];
        }
        
        //TODO: check that crrnt_bpt_bank isn't full
        
		//parse ct file header line
        parse_ct_header(line, RNAvar);
		
    	//parse ct file body
        hel = 0;
        RNAvar->cr_bpt->bp = malloc((RNAvar->fa->len+1) * sizeof(*RNAvar->cr_bpt->bp)); //allocate memory for base pair data
        for (i = 0, align_indx = 0; i < RNAvar->fa->len; i++, align_indx++) {
            if (!isbase(RNAvar->fa->sq[align_indx])) {
                //printf("bumping align index for %s\n", RNAvar->fa->nm);
                align_indx++;	//tracks alignment position for base pair totals
            }

            if	((fgets(line, MAXLINE, fp_ct)) == NULL) {
                return EMPTY_CT_LINE; //return -1 to signal error
            }
            RNAvar->cr_bpt->numpr += parse_ct_seq(line, RNAvar, i); //increments if nt is paired
            hel = test_hel_contig(RNAvar, hel, i);					//test if pair adds to contiguous helix
            if (RNAvar->typ == 'm') { 								//only consider sequences from
                //pair count table does not handle alignment spacing - fix this before using
                //pr_cnt_tbl[i+1][RNAvar->cr_bpt->bp[i]]++;			//increment pair table entry
                if (RNAvar->cr_bpt->bp[i]) {
                    tot_pr_cnt[align_indx]++;                //increment total pairs
                }
                tot_obs_cnt[align_indx]++;
                lcl_bp_cnt[align_indx]++;
                
            }
            
        }
        RNAvar->cr_bpt->numpr /= 2;		//divide numpr by 2 because pairs are counted twice
        identify_bps(RNAvar);     		//evaluate base pair types
        if (RNAvar->typ == 'm') { 		//ensures main dg table counts come from main and not peripheral calls to eval_ct
            test_dg_uniq(RNAvar->cr_bpt->dg, &dg_freq[0]);
            helix_contig[RNAvar->cr_bpt->mxhel]++;    //track contiguous helix totals
            bp_cnt[RNAvar->bpt->numpr]++;             //track total bp distribution
        }
        structures++;
        if (verbose) {
            printf("%25s\t%d\t%3d nts\t%3d bps\tmxhlx %3d\tdeltaG %4.1f\t%s\n", RNAvar->fa->nm, RNAvar->s_cnt, RNAvar->fa->len, RNAvar->cr_bpt->numpr, RNAvar->cr_bpt->mxhel, RNAvar->cr_bpt->dg, RNAvar->sq);
        }
    }
    
    if (fclose(fp_ct) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    return structures;
}

/* parse_ct_header: parse header line of ct file for deltaG value */
int parse_ct_header(char * line, struct RNA *RNAvar) {

    int i = 0;
    int j = 0;
    char tmp_str[16] = {0};
    
    //bypass values preceeding deltaG; if no '=', proceed to null character
    while (line[i] != '=' && line[i]) {i++;}
    
    //get deltaG value
    if (line[i]) { //confirms that there are pairs to evaluate, no '=' if no pairs
        i+=2; //increment by 2 to get to number TODO: change to be affirmative identification of dg value
        for (j = 0; !isspace(line[i]); j++) {
            tmp_str[j] = line[i++];
        }
        tmp_str[j] = '\0';
        RNAvar->cr_bpt->dg = atof(tmp_str);
    } else {
        RNAvar->cr_bpt->dg = 0; //code for no pairs TODO: give this a definition
    }

    return 1;
}
/* parse_ct_seq: parse a ct file non-header line */
int parse_ct_seq(char *line, struct RNA *RNAvar, int i)
{
    int j = 0;
    int k = 0;
    char tmp_str[16];
    
    //bypass unneeded values
    while (!isalpha(line[j])) {j++;} //bypass spaces
    j++; //bypass base identity
    while (isspace(line[j])) {j++;} //bypass spaces
    while (isdigit(line[j])) {j++;} //bypass digits
    while (isspace(line[j])) {j++;} //bypass spaces
    while (isdigit(line[j])) {j++;} //bypass digits
    while (isspace(line[j])) {j++;} //bypass spaces
    
    //get base pair value
    for (k = 0; isdigit(line[j]); k++) { //get base pair value
        tmp_str[k] = line[j++];
    }
    tmp_str[k] = '\0';

    return (RNAvar->cr_bpt->bp[i] = atoi(tmp_str)) ? 1 : 0; //1 if base pair, 0 if unpaired
}

/* test_dg_uniq: build table that tracks frequency of predicted structure free energy */
int test_dg_uniq(float dg, struct deltaG_frequency * dg_freq_tbl) //TODO: pass DG table in function call
{
    int dg_indx;
    float tmp_dg_indx;
    
    //produce integer index from float deltaG value
    tmp_dg_indx = (fabs(dg)*10);
    dg_indx = (int)(tmp_dg_indx);
    
    if (dg_freq_tbl[dg_indx].cnt == 0) {
        dg_freq_tbl[dg_indx].val = dg;
        dg_freq_tbl[dg_indx].cnt++;
        return 1; //TODO use return value?
    } else {
        dg_freq_tbl[dg_indx].cnt++;
        return 0;
    }
}

/* test_hel_contig: tests whether a base pair adds to a contiguous helix */
int test_hel_contig(struct RNA *RNAvar, int hel, int i)
{
    if (RNAvar->cr_bpt->bp[i]) {	//non-zero value indicates base pair
        //increment hel if no current helix or if new bp is contiguous (no bulge), else hel = 0
        hel = (!hel || abs(RNAvar->cr_bpt->bp[i] - RNAvar->cr_bpt->bp[i-1]) == 1) ? hel+1 : 0;
    } else {						//no base pair
        hel = 0;
    }

    if (hel > RNAvar->cr_bpt->mxhel) {
        RNAvar->cr_bpt->mxhel = hel;
    }
    
    return hel;
}

/* identify_bps: counts number of each base pair type (G:C, A:U/T, G:U/T) in predicted structure */
int identify_bps(struct RNA *RNAvar)
{
    int i = 0;
    int pair = 0;
    
    for (i = 0; i < RNAvar->fa->len; i++) {
        if (RNAvar->cr_bpt->bp[i] != 0) {
            pair = ispair(RNAvar->sq[i], RNAvar->sq[RNAvar->cr_bpt->bp[i]-1]);
            switch (pair) {
                case 3: RNAvar->cr_bpt->gcp++; break;
                case 2: RNAvar->cr_bpt->aup++; break;
                case 1: RNAvar->cr_bpt->gup++; break;
                case 0:
                    printf("identify_bps: error - called base pair failed to match partner\n");
                    printf("%c:%c\n", RNAvar->sq[i], RNAvar->sq[RNAvar->cr_bpt->bp[i]-1]);
                    abort();
                    break;
                default:
                    printf("identify_bps: error - unexpected return value from ispair(). pair = %d\n", pair);
                    abort();
                    break;
            }
        }
    }
    
    RNAvar->cr_bpt->gcp /= 2; 	//divide pairs in half to get total bp
    RNAvar->cr_bpt->aup /= 2;	//divide pairs in half to get total bp
    RNAvar->cr_bpt->gup /= 2; 	//divide pairs in half to get total bp
    return 0;
}

/* ispair: test whether two nucleotides constitute a G:C, A:U/T, or G:U/T pair.
   returns a distinct value for each type of pair and 0 for no pair. */
int ispair(char nt1, char nt2)
{
    if ((nt1 == 'G' && nt2 == 'C') ||  (nt1 == 'C' && nt2 == 'G')){
        return 3;
    } else if ((nt1 == 'A' && nt2 == 'T') ||  (nt1 == 'T' && nt2 == 'A')){
        return 2;
    } else if ((nt1 == 'G' && nt2 == 'T') ||  (nt1 == 'T' && nt2 == 'G')){
        return 1;
    } else if ((nt1 == 'A' && nt2 == 'U') ||  (nt1 == 'U' && nt2 == 'A')){
        return 2;
    } else if ((nt1 == 'G' && nt2 == 'U') ||  (nt1 == 'U' && nt2 == 'G')){
        return 1;
    } else {
        return 0;
    }
}

int print_output(struct RNA *RNAvar)
{
    printf("in print output\n");
    FILE *tmp_fp;
    
    int i, j, gcp, aup, gup, dg_sum;
    
    i = j = 0;
    
    char cdist_filename[MAXLINE];
    sprintf(cdist_filename, "%s_%X_dg_cdist.txt", filename, run_id);
    tmp_fp = fopen(cdist_filename, "w");
    
    //print deltaG distribution
    fprintf(tmp_fp, "%s_dG\t%sFrac >= dG\n", out_dir, out_dir);
    for (i = MAXSEQS-1, dg_sum = 0; i >= 0; i--) {
        if (dg_freq[i].cnt > 0) {
            dg_sum += dg_freq[i].cnt;
            fprintf(tmp_fp, "%f\t%1.4f\n", dg_freq[i].val, ((float)(dg_sum))/((float)(total_structures)));
        }
    }
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    char pairs_filename[MAXLINE];
    sprintf(pairs_filename, "%s_%X_pair_totals.txt", filename, run_id);
    tmp_fp = fopen(pairs_filename, "w");
    
    fprintf(tmp_fp, "\ntotal variants: %d\n", total_variants);
    //print base pair totals
    gcp = aup = gup = 0;
    fprintf(tmp_fp, "\nbase pair totals\n");
    for (i = 0; i < total_variants; i++) {
        gcp += RNAvar[i].bpt->gcp;
        aup += RNAvar[i].bpt->aup;
        gup += RNAvar[i].bpt->gup;
    }
    fprintf(tmp_fp, "total gc pairs:\t%d\n", gcp);
    fprintf(tmp_fp, "total au pairs:\t%d\n", aup);
    fprintf(tmp_fp, "total gu pairs:\t%d\n", gup);
    
    //print base pair distribution
    fprintf(tmp_fp, "\nbase pair distribution\n");
    for (i = 0; i <= max_seq_len; i++) {
        fprintf(tmp_fp, "%3d bps\t%5d\n", i, bp_cnt[i]);
    }
    
    fprintf(tmp_fp, "\nmax helix distribution\n");
    for (i = 0; i <= max_seq_len; i++) {
        fprintf(tmp_fp, "%3d bps\t%5d\n", i, helix_contig[i]);
    }
    
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    char pos_bp_cnt_filename[MAXLINE];
    sprintf(pos_bp_cnt_filename, "%s_pos_bp_cnt.txt", outdir);
    tmp_fp = fopen(pos_bp_cnt_filename, "w");
    fprintf(tmp_fp, "pos\t%s_bp_cnt\ttotobs\n", outdir);
    for (i = 0; i <= max_seq_indx; i++) {
        fprintf(tmp_fp, "%d\t%d\t%d\n", i+1, tot_pr_cnt[i], tot_obs_cnt[i]);
    }
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
    return 1;
}

/* isbase: test that character matches IUPAC notation */
int isbase(char c) {
    switch (c) {
        case 'A': return 1; break;
        case 'T': return 1; break;
        case 'G': return 1; break;
        case 'C': return 1; break;
        case 'U': return 1; break;
        case 'N': return 1; break;
        case 'R': return 1; break;
        case 'Y': return 1; break;
        case 'K': return 1; break;
        case 'M': return 1; break;
        case 'S': return 1; break;
        case 'W': return 1; break;
        case 'B': return 1; break;
        case 'D': return 1; break;
        case 'H': return 1; break;
        case 'V': return 1; break;
        case 'a': return 1; break;
        case 't': return 1; break;
        case 'g': return 1; break;
        case 'c': return 1; break;
        case 'u': return 1; break;
        case 'n': return 1; break;
        case 'r': return 1; break;
        case 'y': return 1; break;
        case 'k': return 1; break;
        case 'm': return 1; break;
        case 's': return 1; break;
        case 'w': return 1; break;
        case 'b': return 1; break;
        case 'd': return 1; break;
        case 'h': return 1; break;
        case 'v': return 1; break;
        default: break;
    }
    return 0;
}





