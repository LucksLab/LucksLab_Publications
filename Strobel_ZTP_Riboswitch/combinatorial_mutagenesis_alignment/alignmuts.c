//
//  alignmuts_v0.0.10.c
//  
//
//  Created by Eric Strobel on 5/11/18.
//
//

//v0.0.10

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>

#define READ1 0					//index for read 1
#define READ2 1					//index for read 2
#define MAX_LINE 1024			//maximum line length, input files should never have lines this long
#define MAX_3p_LEN 8			//maximum length for string containing end number
#define TABLE_SIZE 16001    	//size of hash table
#define BLOCK_SIZE 16384    	//memory allocation block size
#define DEF_ANCHOR_LEN 31		//default anchor length, probably longer than it ever needs to be
#define MAX_RND_NT 32			//maximum number of variable nucleotides
#define MAX_HASH_IPT 16 		//maximum number of chars to hash
#define READ1_TRIM 4			//number of nts to trim from read1 to remove random prefix, could make option but should be constant
#define MIN_TERM_LEN 130		//minimum length for calling an end as terminated TODO: value should come from 3pEnd targets file
#define MAX_TERM_LEN 134		//maximum length for calling an end as terminated TODO: value should come from 3pEnd targets file

/* structure declarations */
struct target {					//structure containing values for read1 and read 2 targets
    char id[MAX_LINE];			//target identifier
    char sq[MAX_LINE];			//target sequence
    int val;            		//if read1: end value; if read2: total observations
    float term;         		//used for read2: number of reads identified as terminated
    float antiterm;     		//used for read2: number of reads identified as antiterminated
};

struct rd2_anchor {				//contains values for read 2 anchor sequence
    char sq[MAX_LINE];			//anchor sequence
    int pos;					//anchor position in read 2 reference target
    int len;					//length of anchor sequence
    int rnd[MAX_RND_NT];		//position of variable nucleotides relative to .pos
    int rnd_cnt;				//number of variable nucleotides in read 2 reference target
};

struct h_node {     			//hash table node
    struct target *trg;			//points to read 1 target
    struct h_node *nxt;			//points to next node in linked list, used to handle collisions
};

struct h_node_bank {      		//hash table node bank for memory management
    struct h_node *hn;      	//pointer for h_node memory allocation
    struct h_node_bank *nxt;   	//pointer to next bank
    int count;           		//number of h_nodes used in current bank
} ;

/* global variables */
//variables below are global to either avoid passing them constantly
//or to allow functions to be easily generalized to both reads without
//passing read 2 specific variables when performing operations on read 1
struct rd2_anchor anchor;		//anchor for read 2 mapping
int trgt_cnt[2];            	//number of read 1/2 targets - indx0:rd1, indx1:rd2

/* function prototypes */
//initialization
int get_file(FILE **ifp, char *ipt);

//target parsing
int get_line(char *line, FILE *ifp);
void parse_targets_rd1(FILE *ifp, struct target *trgts);
int parse_targets_rd2(FILE *ifp, struct target *trgts);
int find_anchor(char *ref, struct target *trgts_rd2);

//hashing and read mapping
void mk_htbl(struct h_node **htbl, struct h_node_bank *bank, struct target *trgts, int read_indx);
struct h_node** srch_htbl(char *read, struct h_node **htbl, int read_indx);
uint32_t hash_rd1(char *hash_ipt);
uint32_t hash_rd2(char *hash_ipt, int *offset);
uint32_t seq2bin_hash(char *hash_seq);
void extend_h_bank(struct h_node_bank *crrnt_hn_rd1_bank);
int map_reads(FILE *fp_rd1, FILE *fp_rd2, struct h_node **htbl_rd1, struct h_node **htbl_rd2, struct target *trgts);
int lnkrMM_srch(char *read1, struct target *trgts);

//report results
void print_outpt(struct target *trgts_rd2);

int main(int argc, char *argv[])
{
    FILE *fp_targets_rd1 = NULL;   	//pointer for read 1 targets file
    FILE *fp_targets_rd2 = NULL; 	//pointer for read 2 targets file
    FILE *fp_rd1 = NULL;   			//pointer for read 1 fastq
    FILE *fp_rd2 = NULL;   			//pointer for read 2 fastq
    
    int libf1 = 0;		//track read 1 file input
    int libf2 = 0;		//track read 2 file input
    int trgf1 = 0;		//track targets 1 file input
    int trgf2 = 0;		//track targets 2 file input
    
    /* parse options using getopt_long */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"read1",       required_argument,  0,  't'},
            {"read2",       required_argument,  0,  'b'},
            {"m_targets",   required_argument,  0,  'm'},
            {"e_targets",   required_argument,  0,  'e'},
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "t:b:m:e:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
            case 't': get_file(&(fp_rd1), argv[optind-1]); libf1++; break;
            case 'b': get_file(&(fp_rd2), argv[optind-1]); libf2++; break;
            case 'm': get_file(&fp_targets_rd2, argv[optind-1]); trgf1++; break;
            case 'e': get_file(&fp_targets_rd1, argv[optind-1]); trgf2++; break;
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    
    if (optind < argc) {
        printf("\nnon-option ARGV-elements:\n");
        while (optind < argc)
            printf("\n%s \n", argv[optind++]);
        putchar('\n');
        printf("Aborting program.\n\n");
        abort();
    }
    /* end of option parsing */
    
    /* check that essential files are supplied */
    if ((libf1 + libf2) < 2) {
        printf("main: error - too few fastq files supplied\n");
        abort();
    }
    
    if ((trgf1 + trgf2) < 2) {
        printf("main: error - too few targets files supplied\n");
        abort();
    }
    /* end check */
    
    /* targets initialization and parsing */
    struct target *trgts_rd1;   //pointer for array of read 1 targets
    struct target *trgts_rd2;   //pointer for array of read 2 targets
    
    if ((trgts_rd1 = calloc(BLOCK_SIZE, sizeof(*trgts_rd1))) == NULL) {
        printf("main: error - hash table memory allocation failed\n");
        return 1;
    }
    
    if ((trgts_rd2 = calloc(BLOCK_SIZE, sizeof(*trgts_rd2))) == NULL) {
        printf("main: error - hash table memory allocation failed\n");
        return 1;
    }
    
    extern int trgt_cnt[2];		//number of targets for from each targets file
    trgt_cnt[0] = trgt_cnt[1] = 0;
    
    parse_targets_rd1(fp_targets_rd1, trgts_rd1);
    parse_targets_rd2(fp_targets_rd2, trgts_rd2);
    /* end of targets initialization and parsing*/
    
    /* hash table initialization and construction*/
    struct h_node **htbl_rd1;    		//read 1 hash table root
    struct h_node **htbl_rd2;     		//read 2 hash table root
    struct h_node_bank hn_rd1_bank;		//bank for read 1 hash table nodes
    struct h_node_bank hn_rd2_bank;		//bank for read 2 hash table nodes
    
    hn_rd1_bank.count = 0;
    hn_rd2_bank.count = 0;
    
    if ((htbl_rd1 = calloc(TABLE_SIZE, sizeof(*htbl_rd1))) == NULL) {
        printf("main: error - hash table memory allocation failed\n");
        return 1;
    }
    
    if ((htbl_rd2 = calloc(TABLE_SIZE, sizeof(*htbl_rd2))) == NULL) {
        printf("main: error - hash table memory allocation failed\n");
        return 1;
    }
    
    if ((hn_rd1_bank.hn = calloc(BLOCK_SIZE, sizeof(*(hn_rd1_bank.hn)))) == NULL) {
        printf("main: error - hash table node memory allocation failed\n");
        return 1;
    }
    
    if ((hn_rd2_bank.hn = calloc(BLOCK_SIZE, sizeof(*(hn_rd2_bank.hn)))) == NULL) {
        printf("main: error - hash table node memory allocation failed\n");
        return 1;
    }
    
    struct h_node_bank *crrnt_hn_rd1_bank = &hn_rd1_bank;	//pointer for handling hash node read 1 bank
    struct h_node_bank *crrnt_hn_rd2_bank = &hn_rd2_bank;	//pointer for handling hash node read 2 bank
    
    mk_htbl(htbl_rd1, crrnt_hn_rd1_bank, trgts_rd1, READ1);
    mk_htbl(htbl_rd2, crrnt_hn_rd2_bank, trgts_rd2, READ2);
    /* end of hash table initialization and construction*/
    
    /* map reads */
    map_reads(fp_rd1, fp_rd2, htbl_rd1, htbl_rd2, trgts_rd1);	//map reads
    print_outpt(trgts_rd2);											//report results
}

/* get_file: open file input and assign to file pointer */
int get_file(FILE **ifp, char *ipt)
{
    if ((*ifp = fopen(ipt, "r")) == NULL) {
        printf("error: could not open %s as read one file. Aborting program...\n", ipt);
        abort();
    }
    return 1;
}

/* get_line: get line from file, place into array, and return line length if successful */
int get_line(char *line, FILE *ifp) {
    /* 	function gets line and replaces terminal newline with null character.
     	there is no need for buffering in this case because lines that exceed
     	MAX_LINE should not exist and if they do, they are an error.
     	the only acceptable mode of failure is to reach the end of the file
     	without getting any preceeding characters */
    
    int i = 0;
    char c = 0;
    
    for (i = 0; (c = fgetc(ifp)) != '\n' && c != EOF && c &&  i < MAX_LINE; i++) {
        line[i] = c;
    }
    if (c == '\n' && i != MAX_LINE) {
        line[i] = '\0'; 		//remove trailing newline
        return i;				//success
    } else if (c == EOF) {  	//reached end of file
        if (i == 0) {			//EOF is expected at the start of a line
            return 0;
        } else {				//unexpected EOF
            printf("get_line: error -last line in file not terminated by newline\n");
            abort();
        }
	}else if (!c) {				//unexpected null character
        printf("get_line: error - unanticipated null character\n");
        abort();
    } else if (i == MAX_LINE) {	//unexpected long line
        printf("get_line: error - unanticipated long line\n");
        abort();
    } else {
        return 0; 				//TODO: need an other error code? this should not be reachable
    }
}

/* parse_targets_rd1: parse read 1 targets file */
void parse_targets_rd1(FILE *ifp, struct target *trgts)
{
    /* 	each target has an id line and a sequence line.
     	the id line contains the end number preceeded by '3p. */
    extern int trgt_cnt[2];	//global variable tracking read1/2 target counts

    int i = 0;
    char *ptr_3p = NULL; 		//pointer to '3p' in character id
    char end_3p[MAX_3p_LEN];	//array for end string

    while (get_line(trgts[trgt_cnt[READ1]].id, ifp)) {	//get target id
        if ((ptr_3p = strstr(trgts[trgt_cnt[READ1]].id, "3p")) != NULL && ptr_3p[2] && strlen(&ptr_3p[2]) < (MAX_3p_LEN-1)) {
            for (i = 0; isdigit(ptr_3p[i+2]) && ptr_3p[i+2] && i < (MAX_3p_LEN-1); i++) {
                end_3p[i] = ptr_3p[i+2];
            }
            end_3p[i] = '\0';
            trgts[trgt_cnt[READ1]].val = atoi(end_3p); 	//set target end value
        } else {
            printf("parse_targets_rd1: error incorrect read 1 targets file format\n");
            abort();
        }
        if (!get_line(trgts[trgt_cnt[READ1]++].sq, ifp)) {	//get target sequence
            printf("parse_targets_rd1: target missing sequence\n");
            abort();
        }
    }
}

/* parse_targets_rd2*/
int parse_targets_rd2(FILE *ifp, struct target *trgts)
{
    extern int trgt_cnt[2];	//global variable tracking read1/2 target counts
    struct target rd2_ref; 	//reference sequence anchor search
	
    //confirm that targets file starts with /REF target
    //then evaluate /REF target to find anchor
    if (get_line(rd2_ref.id, ifp) && rd2_ref.id[0] == '/') {
        get_line(rd2_ref.sq, ifp);
    } else {
        printf("parse_targets_rd2: error - incorrect read 2 targets file format\n");
        abort();
    }
	
    while (get_line(trgts[trgt_cnt[READ2]].id, ifp)) {
        if (!get_line(trgts[trgt_cnt[READ2]++].sq, ifp)) {
            printf("parse_targets_rd2: target missing sequence\n");
            abort();
        }
    }
    
    if (!find_anchor(rd2_ref.sq, trgts)) {	//identify anchor sequence
        printf("parse_targets_rd2: failed to find anchor sequence\n");
        abort();
    }
    
    return trgt_cnt[READ2];
}

/* find anchor: finds anchor sequence for mapping read 2 */
int find_anchor(char *ref, struct target *trgts_rd2)
{
    extern struct rd2_anchor anchor;	//global variable for anchor values
    
    int i = 0;
    int anch_indx = 0;	//index for anchor sequence
    int trgt = 0;		//index for targets
    int pass = 0;		//tracks whether candidate anchor is valid
    int len_match =0;	//tracks whether initial anchor search yielded long enough string
    int crrnt_len= 0;	//current sequence length of anchor search
    int prev_len = 0; 	//previous sequence length of anchor search
    char *uniq = NULL;	//tracks whether anchor is unique within read 2 reference target
    char *mtch1 = NULL;	//pointer to first occurrence of anchor in target sequence
    char *mtch2 = NULL;	//pointer to second occurrence of anchor in target sequence
    char *test2 = NULL;	//pointer for testing if more than one occurrence of anchor in target sequence
    
    
    crrnt_len = DEF_ANCHOR_LEN;
    while (1) {
        //find candidate anchor
        for (i = 0, anch_indx = 0; anch_indx < crrnt_len && ref[i]; i++) {
            if (ref[i] == 'N' || ref[i] == 'R' || ref[i] == 'Y') {
                anch_indx = 0;
            } else {
                anchor.sq[anch_indx++] = ref[i];
            }
        }
        anchor.sq[anch_indx] = '\0';
        if ((len_match = (anch_indx == crrnt_len) ? 1 : 0) && (uniq = strstr(&ref[i-crrnt_len+1], anchor.sq)) == NULL) {
            for (trgt = 0, pass = 1, mtch1 = mtch2 = NULL; trgt < trgt_cnt[READ2] && pass; trgt++) {
                mtch1 = strstr(trgts_rd2[trgt].sq, anchor.sq);
                test2 = &mtch1[1];
                mtch2 = strstr(test2, anchor.sq);
                if (mtch1 == NULL || mtch2 != NULL) {
                    pass = 0;
                }
            }
            if (pass) {
                anchor.pos = i-crrnt_len;
                anchor.len = crrnt_len;
                for (i = 0, anchor.rnd_cnt = 0; ref[i]; i++) { //identify variable positions relative to anchor.pos
                    if (ref[i] == 'N'|| ref[i] == 'R' || ref[i] == 'Y') {
                        anchor.rnd[anchor.rnd_cnt++] = i - anchor.pos; //subtract anchor.pos from i to get relative location of variable position
                    }
                }
                //report anchor values
                printf("anchor = %s\npos = %d\nlen = %d\n%d variable bases\n", anchor.sq, anchor.pos, anchor.len, anchor.rnd_cnt);
                for (i = 0; i < anchor.rnd_cnt; i++) {
                    printf("rnd%d = %d\n", i, anchor.rnd[i]);
                }
                return 1;
            } else {
                uniq = mtch2;
            }
        }
        
        if (!len_match) { //could not find unique sequence that does not overlap random nt, shorten anchor
            crrnt_len--;
        } else if (uniq != NULL) { //could not find sufficiently unique sequence, lengthen anchor
            crrnt_len++;
        }
        if (crrnt_len == prev_len) { //no sufficient anchor exists
            return 0;
        } else {
            prev_len = crrnt_len;
        }
    }
}

/* mk_htbl: construct hash table from target structure */
void mk_htbl(struct h_node **htbl, struct h_node_bank *bank, struct target *trgts, int read_indx)
{
    /* hash table has linked list buckets for possible collisions */
    extern int trgt_cnt[2];	//global variable tracking read1/2 target counts
	
    struct h_node **p_rdnd = NULL; //pointer for h_node handling
    int i = 0;
    
    for (i = 0; i < trgt_cnt[read_indx]; i++) {
        p_rdnd = srch_htbl(trgts[i].sq, htbl, read_indx);	//search hash table for duplicate entries
        if ((*p_rdnd) == NULL) {							//target is unique
            (*p_rdnd) = &bank->hn[bank->count++];			//assign node from hash node bank
            if (bank->count == BLOCK_SIZE) {				//check that bank was not filled
                extend_h_bank(bank);						//extend bank if neededt
            }
            (*p_rdnd)->trg = &(trgts[i]);					//set node to point to target
        } else {
            printf("mk_htbl: error - non-unique target found in read %d targets file\n", read_indx);
            printf("new: %s\nprev: %s\n", trgts[i].sq, (*p_rdnd)->trg->sq);
            abort();
        }
    }
    printf("made htbl: %d entries\n", i);
}

/* srch_htbl: search hash table for sequence match */
struct h_node** srch_htbl(char *read, struct h_node **htbl, int read_indx)
{
    static struct h_node *dummy = NULL; //TODO: set this to point to an identifiable target to signal anchor failure
    
    struct h_node **p_rdnd = NULL; 		//pointer for h_node handling
    uint32_t hash = 0xFFFFFFFF;			//hash value
    int offset = 0;						//offset for anchor position in read2 string comparison, 0 for read1
    
    //determine appropriate hash based on read
    if (read_indx == READ1) {
        hash = hash_rd1(read);
    } else if (read_indx == READ2) {
        if ((hash = hash_rd2(read, &offset)) == 0xFFFFFFFF) {
            p_rdnd = &dummy;
            return p_rdnd;
        }
    }
    
    //hash table search
    p_rdnd = &htbl[hash]; //initialize rd1 node pointer to hash table entry
    
    int i = 0;
    int j = 0;
    int mtch = 0;	//tracks whether comparison is currently a match
   
    while ((*p_rdnd) != NULL) {
        j = (read_indx == READ2) ? (0 + offset - anchor.pos) : 0; //j is conditionally set depending on read_indx
        for (i = 0, mtch = 1;
             (*p_rdnd)->trg->sq[i] && read[j] && mtch; i++, j++) {
            if ((*p_rdnd)->trg->sq[i] != read[j]) {
                mtch = 0;
            }
        }
        if (mtch && !(*p_rdnd)->trg->sq[i]) { //for loop reached end of target string without a mismatch
            return p_rdnd; //match, return pointer to target match node
        } else {
            p_rdnd = &(*p_rdnd)->nxt; //proceed to next node in linked list
        }
    }
    return p_rdnd; //no match
}

/* hash_rd1 */
uint32_t hash_rd1(char *hash_ipt)
{
    /* the hash is simply encoding the expected non-linker segment of read 1
     into a 2 bit notation and taking the modulus of the (prime number) table size */
    
    int len = 0;			//length of input string
    char *hash_ptr = NULL;	//pointer to sequence that will be hashed
    
    len = strlen(hash_ipt); //need to start at end of string to get non-linker sequence
    hash_ptr = (len > MAX_HASH_IPT) ? &hash_ipt[len-MAX_HASH_IPT] : &hash_ipt[0];

    return (seq2bin_hash(hash_ptr) % TABLE_SIZE);
}

/* hash_rd2 */
uint32_t hash_rd2(char *hash_ipt, int *offset)
{
    /* the hash is simply encoding the variable nts of read 2
     into a 2 bit notation and taking the modulus of the (prime number) table size */
    
    extern struct rd2_anchor anchor;
    
    int i = 0;
    char *anch_pt = NULL;		//pointer for first occurrence of anchor seq in hash_ipt
    char hash_seq[MAX_HASH_IPT+1];			//array for sequence to hash
    uintptr_t ntv_offset = 0;   //difference between anch_pt and hash_ipt[0] addresses

    //search for first hit of anchor seq in hash_ipt
    //DEF_ANCHOR_LEN is sufficiently long that any duplicates
    //will be aberrant fragments, so only checking for first hit is OK
    if ((anch_pt = strstr(hash_ipt, anchor.sq)) != NULL) {
        if ((ntv_offset = (uintptr_t)anch_pt - (uintptr_t)hash_ipt) < INT_MAX) { //test if cast to uint32_t is OK
            *offset = (int)ntv_offset;
        } else {
            printf("hash_rd2: ERROR unexpectedly large anchor point offset\n");
            abort();
        }
        
        for (i = 0; i < anchor.rnd_cnt && i < MAX_HASH_IPT; i++) {
            hash_seq[i] = *(anch_pt + anchor.rnd[i]);
        }
        hash_seq[i] = '\0';
        return (seq2bin_hash(hash_seq) % TABLE_SIZE);
    } else {
        return 0xFFFFFFFF; //failed to locate anchor sequence
    }
}

/* seq2bin_hash */
uint32_t seq2bin_hash(char *hash_seq) {
    
    int i = 0;
    int mshft = 0;				//maximum bit shift value
    int shft = 0;				//bit shift value for each iteration of loop
    uint32_t ctob = 0;			//value of nucleotide character to 2 bit conversion
    uint32_t hash_outpt = 0; 	//2 bit encoded nucleotide character sequence
    
    //perform hash
    mshft = (sizeof(hash_outpt) << 3) - 2; //initialize maximum shift to number of bits - 2
    for (i = 0, ctob = 0, shft = mshft; hash_seq[i] && shft >= 0; i++, shft -= 2) {
        ctob = (hash_seq[i] & 7) >> 1; //lowest 3 bit mask + 1 bit rshift produces unique int for ATGC
        hash_outpt |= ctob << shft; //shift to OR 2bit encoded nucleotide into place
    }
    return hash_outpt;
}

/* extend_h_bank */
void extend_h_bank(struct h_node_bank *crrnt_hn_rd1_bank)
{
    if ((crrnt_hn_rd1_bank->nxt = calloc(1, sizeof(*(crrnt_hn_rd1_bank->nxt)))) == NULL) {
        printf("srch_h_tbl: error - hash table node bank memory allocation failed\n");
        abort();
    }
    crrnt_hn_rd1_bank = crrnt_hn_rd1_bank->nxt;
    if ((crrnt_hn_rd1_bank->hn = calloc(BLOCK_SIZE, sizeof(*(crrnt_hn_rd1_bank->hn)))) == NULL) {
        printf("srch_h_tbl: error - hash table node bank memory allocation failed\n");
        abort();
    }
}

/* map_reads: use hash tables and mismatch search to map reads 1 and 2 to a mutant and end */
int map_reads(FILE *fp_rd1, FILE *fp_rd2, struct h_node **htbl_rd1, struct h_node **htbl_rd2, struct target *trgts)
{
    printf("mapping reads\n");
    char line_rd1[MAX_LINE], line_rd2[MAX_LINE]; 	//arrays to fill with fastq lines that are currently ignored
    char rd1[MAX_LINE], rd2[MAX_LINE];				//arrays for fastq read lines
    char *trim_rd1 = NULL;							//pointer for trimming from read 1
    
    int end = 0;					//value of mapped end
    int rd1_len = 0;				//length of read 1
    int rd2_len = 0;				//length of read 2
    
    struct h_node *p_rd1nd = NULL; 	//pointer to mapped rd1 targets
    struct h_node *p_rd2nd = NULL;	//pointer to mapped rd2 targets
    
	int mapped_ends = 0;			//read 2 and read 1 mapped
	int unmapped_ends = 0;			//only read 2 mapped
    int read_count = 0;				//total reads
    
    //the calls to fgets() below will not fail unless the user
    //has edited the fastq file (which should NEVER happen)
    //or has changed the value of MAX_LINE
    while (fgets(line_rd1, MAX_LINE, fp_rd1) != NULL && fgets(line_rd2, MAX_LINE, fp_rd2) != NULL) { //bypass fastq line1
        rd1_len = get_line(rd1, fp_rd1); 	//get read 1
        rd2_len = get_line(rd2, fp_rd2);	//get read 2
        trim_rd1 = &rd1[READ1_TRIM]; 		//set pointer to ignore trimmed sequence

        //map read
        if ((p_rd2nd = *srch_htbl(rd2, htbl_rd2, READ2)) != NULL) {	//hash table search for read 2
            //printf("%s\n", rd2);
            p_rd2nd->trg->val++; //track total read 2 alignments
            if ((p_rd1nd = *srch_htbl(trim_rd1, htbl_rd1, READ1)) != NULL) { //hash table search for read 1
                end = p_rd1nd->trg->val;
            } else {
                end = lnkrMM_srch(trim_rd1, trgts); //mismatch search for read 1
            }
            if (end) {	//end was mapped
                mapped_ends++;
                if (end >= MIN_TERM_LEN && end <= MAX_TERM_LEN) { //determine if mapped end is terminated
                    p_rd2nd->trg->term++;
                } else if (end > MAX_TERM_LEN){	//determine if mapped end is antiterminated
                    p_rd2nd->trg->antiterm++;
                }
            } else {
                unmapped_ends++;
            }
        }
		//bypass fastq lines 3 and 4
        fgets(line_rd1, MAX_LINE, fp_rd1);
        fgets(line_rd2, MAX_LINE, fp_rd2);
        fgets(line_rd1, MAX_LINE, fp_rd1);
        fgets(line_rd2, MAX_LINE, fp_rd2);
        
        read_count++;
        if ((read_count % 1000000) == 0) {
            printf("time:%d reads %f\n", read_count, ((float)clock())/CLOCKS_PER_SEC);
        }
    }
    printf("\n%d mapped ends\n%d unmapped ends\n", mapped_ends, unmapped_ends);
    return 1;
}

/* lnkrMM_srch */
int lnkrMM_srch(char *read1, struct target *trgts)
{
    //the search below maps reads that did not map using the hash table by testing
    //whether substitutions, insertions, or deletions permits alignment
    //reads are only considered for this alignment if they have a perfect linker
    //and no insertions preceeding the linker that shorten the mappable sequence
    //in each alignment only one mismatch is allowed
    //reads mapped this way will frequently map to several targets
    //multi-mapping is allowed if categorization as terminated or antiterminated is
    //unambiguous and mapping clusters across three targets with an end difference of
    //no more than two
    
    static const char lnkr[21] = "GTCCTTGGTGCCCGAGTCAG"; //linker sequence, constant to each read
    static const int lnkr_len = 20;	//length of linker sequence
    
    int i = 0;
    int j = 0;
    int trgt = 0;					//index for target
    int mtch_cnt = 0;				//number of matches for a given read, tracks multi-maps
    int full_len = 0;				//full length of read
    int shrt_len = 0;				//length of sequence after linker sequence
    int mxMM = 1; 					//maximum number of mismatches allowed TODO: could make this an option
    char *lnkr_strt = NULL; 		//pointer to start of linker in read1
    char *map_ptr = NULL;			//pointer for read being mapped
    struct target *p_trgmtch[64]; 	//pointer to target matches TODO: cap number of multi-maps so it doesn't exceed array size
    
    //variables for substitution, insertions, and deletion searches, respectively
    int MM_sub, rsub, tsub;
    int MM_ins, rins, tins;
    int MM_del, rdel, tdel, last_mtch;
    MM_sub = rsub = tsub = 0;
    MM_ins = rins = tins = 0;
    MM_del = rdel = tdel = last_mtch = 0;
    int last_trgt = -1;
    
    //set lnkr_strt to first linker sequence nt in read 1
    lnkr_strt = strstr(read1, lnkr);
    full_len = strlen(read1); //length of read
    if (lnkr_strt != NULL) {
        map_ptr = &lnkr_strt[lnkr_len];
        shrt_len = strlen(map_ptr);
        if ((shrt_len+lnkr_len) == full_len) { //test if read segment plus linker matches full read length
            //printf("%s\n", read1);
            //TODO: confirm that shrt_len is always a good stand in for target length measurement
            
            /*	right now these three searches are really inefficient but don't meaningfully
             impact performance, so leaving as is */
            
            /*	search for substitutions - simply moves past mismatch position to test
             	remaining sequence */
            for (trgt = 0; trgt < trgt_cnt[READ1]; trgt++) {
                for (rsub =0, tsub = lnkr_len, MM_sub = 0; rsub < shrt_len && MM_sub <= mxMM; rsub++, tsub++) {
                    if (map_ptr[rsub] != trgts[trgt].sq[tsub]) {
                        MM_sub++;
                    }
                }
                
                if (rsub == shrt_len && MM_sub <= mxMM) {
                    //printf("%s\t%d\tsub\n", trgts[trgt].sq, trgts[trgt].val);
                    p_trgmtch[mtch_cnt++] = &trgts[trgt];
                }
            }
            
            /*	search for single insertions - moves past mismatch position in read and compares
             	to remaining target sequence */
            for (trgt = 0; trgt < trgt_cnt[READ1]; trgt++) {
            	for (rins =0, tins = lnkr_len, MM_ins = 0; rins < shrt_len && MM_ins <= mxMM; rins++) {
                    if (map_ptr[rins] != trgts[trgt].sq[tins]) {
                        MM_ins++;
                    } else {
                        tins++; //only increment tins if no mismatch
                    }
            	}
                
                if (rins == shrt_len && MM_ins <= mxMM) {
                    //printf("%s\t%d\tins\n", trgts[trgt].sq, trgts[trgt].val);
                    p_trgmtch[mtch_cnt++] = &trgts[trgt];
                }
            }
        
            /* 	search for single deletions - moves past mismatch position in target and compares
             	to remaining read sequence, then checks last read character against last char of
             	preceeding target to confirm match */
            //last_trgt = -1;
            for (trgt = 0; trgt < trgt_cnt[READ1]; trgt++) {
                last_mtch = 0;
                for (rdel = 0, tdel = lnkr_len, MM_del = 0; trgts[trgt].sq[tdel] && MM_del <= mxMM; tdel++) {
                    if (map_ptr[rdel] != trgts[trgt].sq[tdel]) {
                        MM_del++;
                    } else {
                        rdel++; //only increment rdel if no mismatch
                    }
                }
                
                if (tdel == shrt_len+lnkr_len && MM_del <= mxMM) {
                    p_trgmtch[mtch_cnt++] = &trgts[trgt];
                }
                
            }
        }
    }
    
    /* 	test if mismatched read exclusively multi-maps
    	to terminated or antiterminated positions
    	return only reads that are unambigously
    	terminated or antiterminated and cluster across
    	a maximum of three targets */
    int t_mtch =0;			//counts multi-maps to terminated ends
    int at_mtch =0;			//counts multi-maps to antiterminated ends
    int mx_end = 0; 		//max mapped end value
    int mn_end = INT_MAX; 	//min mapped end value
    
    
    //printf("%s\n", read1);
    for (i = 0; i < mtch_cnt; i++) {
        //printf("%s\t%d\n", p_trgmtch[i]->sq, p_trgmtch[i]->val);
        if (p_trgmtch[i]->val > mx_end) {
            mx_end = p_trgmtch[i]->val;
        }
        if (p_trgmtch[i]->val < mn_end) {
            mn_end = p_trgmtch[i]->val;
        }
        if (p_trgmtch[i]->val >= MIN_TERM_LEN && p_trgmtch[i]->val <= MAX_TERM_LEN) {
            t_mtch++;
        } else if (p_trgmtch[i]->val > MAX_TERM_LEN) {
            at_mtch++;
        }
    }
    //printf("\n");
    
    if (mtch_cnt && ((t_mtch == mtch_cnt) || (at_mtch == mtch_cnt)) && (mx_end - mn_end <= 2)) { //TODO: decide how best to replace hardcoded 2
        return p_trgmtch[0]->val; //unambiguous match, return first observed target match
    } else {
        return 0;
    }
}

void print_outpt(struct target *trgts_rd2)
{
    FILE *tmp_fp;	//pointer for output file
    
    tmp_fp = fopen("out.txt", "w");
    
    int i;
    fprintf(tmp_fp, "mutant\trd2obs\trd2map\tterm\tantiterm\tfrac_readthrough\n");
    for (i = 0; i < trgt_cnt[READ2]; i++) {
        fprintf(tmp_fp, "%s\t%d\t%d\t%d\t%d\t%f\n",trgts_rd2[i].id, trgts_rd2[i].val, (int)(trgts_rd2[i].term + trgts_rd2[i].antiterm),(int)(trgts_rd2[i].term), (int)(trgts_rd2[i].antiterm), (trgts_rd2[i].antiterm /(trgts_rd2[i].term + trgts_rd2[i].antiterm)));
    }
    
    if (fclose(tmp_fp) != 0) {
        printf("file not closed\n");
        abort();
    }
    
}

