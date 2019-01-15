//
//  mk_3Pend_trgts.c
//  
//
//  Created by Eric Strobel on 4/6/18.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#define MAXLEN 128
#define MAXTARGET 256
#define N 0
#define R 1
#define Y 2

static struct varbase {	//stores information for variable positions
    int crnt;			//current index
    char not;			//IUPAC notation
    char alt[5];		//alternate bases at variable position
} vb_tbl[3] = {
    0, 'N', "ATGC",    //0
    0, 'R', "AG",      //1
    0, 'Y', "TC"       //2
};

static struct basemap {				//map of constant and variable positions in reference sequence
    char nts[MAXLEN];				//sequence, * denotes variable base
    struct varbase *p_vb[MAXLEN];	//pointer to vb_tbl table entry
} vbases;

struct trgt_vb {	//track variable base identity through recursion paths
    int indx;		//current index (== variable bases previouslyobserved)
    char bs[17];	//variable base identity
    int pos[17];	//variable base position
};

int isbase(char c);						//test that character matches IUPAC base notation
int get_ref(char * ipt, char * ref);	//get user supplied reference sequence input
int get_start_pos(char * ipt);			//get user supplied start nucleotide number
int get_match_len(char * ipt);			//get sequence length to append for target matching
int rev_comp(char * ipt, char * rc);	//reverse complement

int get_v_bases(char input[500], struct basemap vb[500]); //parses input string and points variable bases to vbase table entry
void mk_trgts_file(int nxt, char outpt[MAXLEN], struct trgt_vb lcl_bases); //recursively construct targets file

const char lnkr[21] = "GTCCTTGGTGCCCGAGTCAG";	//linker sequence for SHAPE-Seq v2.1

int trgt_cnt = 0;	//total targets
int start = 0;		//start nucleotide number
FILE * out = NULL;	//output file pointer

int nt = 0;
int mtchlen	= 0;

int main(int argc, char *argv[])
{
    int i = 0;
    
    char outpt[MAXLEN] = {0};
    char ref[MAXLEN] = {0};
    char out_name[256] = {0};
    
    /* parse options using getopt_long */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"ref",       	required_argument,  0,  'r'},
            {"start-pos",   required_argument,  0,  's'},
            {"match-len",   required_argument,  0,  'm'},
            {"output-name",	required_argument,	0,	'o'},
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "r:s:o:m:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0: /*printf("long option\n");*/ break;
            case 'r': get_ref(argv[optind-1], ref); break;
            case 's': get_start_pos(argv[optind-1]); break;
            case 'o': strcpy(out_name, argv[optind-1]); break;
            case 'm': get_match_len(argv[optind-1]); break;
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
    
    if (mtchlen < 1) { //TODO: add check that ref is longer than match length
        printf("main: error - must specify match length > 0\n");
        abort();
    }

    char ref_rc[MAXLEN] = {0};
    int len = 0;
    len = rev_comp(ref, ref_rc); //reverse complement input sequence
    
    printf("ref = %s\n", ref);
    printf("rc = %s\n", ref_rc);
    printf("start = %d\n", start);
    
    struct trgt_vb lcl_bases = {0};
    
    if (out_name[0]) {
        out = fopen(out_name, "w");
    } else {
        out = fopen("out.fasta", "w");
    }
    
    int j = 0;
    int k = 0;
    char tmp_trgt[14] = {0};
    int more_trgts = 1;
    for (i = len-mtchlen, nt = start; i >= 0; i--, nt++) {
        for (j = i, k = 0; k < mtchlen; j++, k++) {
            tmp_trgt[k] = ref_rc[j];
        }
        if (!get_v_bases(tmp_trgt, &vbases)) { //no variable bases, print tmptrgt
            fprintf(out, ">cbe_ZTP_3p%d\n%s%s\n", nt, lnkr, tmp_trgt);
        } else { //sequence contains variable base, construct all possible target variations
            mk_trgts_file(0, outpt, lcl_bases);
        }
    }
    
    fclose(out);
}

/* get_ref: get reference sequence */
int get_ref(char * ipt, char * ref)
{
    int i = 0;
    
    for (i = 0; ipt[i] && i <= MAXLEN; i++) {
        if (isbase(ipt[i])) {
            ref[i] = ipt[i];
        } else {
            printf("get_ref: error - unexpected character %c in reference sequence\n", ipt[i]);
            abort();
        }
    }
    ref[i] = '\0';
    
    return 1;
}

/* get_start_pos: get user supplied start position value */
int get_start_pos(char * ipt)
{
    extern int start;
    
    int i = 0;
    char s[32] = {0};
    
    strcpy(s, ipt);
    start = atoi(s);
    
    return 1;
}

/* get_match_len: get user supplied start position value */
int get_match_len(char * ipt)
{
    extern int mtchlen;
    
    int i = 0;
    char s[32] = {0};
    
    strcpy(s, ipt);
    mtchlen = atoi(s);
    
    return 1;
}

/* rev_comp: reverse complement input DNA sequence string */
int rev_comp(char * ipt, char * rc)
{
    int i = 0;
    int j = 0;
    int len = 0;
    
    len = strlen(ipt);
    
    printf("len = %d\n", len);
    
    for (i = len-1, j = 0; i >= 0; i--, j++) {
        switch (ipt[i]) {
            printf("ipt[%d] = %d\n", i, ipt[i]);
            case 'A': rc[j] = 'T'; break;
            case 'T': rc[j] = 'A'; break;
            case 'G': rc[j] = 'C'; break;
            case 'C': rc[j] = 'G'; break;
            case 'R': rc[j] = 'Y'; break;
            case 'Y': rc[j] = 'R'; break;
            case 'N': rc[j] = 'N'; break;
            case 'a': rc[j] = 'T'; break;
            case 't': rc[j] = 'A'; break;
            case 'g': rc[j] = 'C'; break;
            case 'c': rc[j] = 'G'; break;
            case 'r': rc[j] = 'Y'; break;
            case 'y': rc[j] = 'R'; break;
            case 'n': rc[j] = 'N'; break;
            default:
                printf("rev_comp: error - unexpected character %c/%d\n", ipt[i], ipt[i]);
                break;
        }
    }
    rc[j] = '\0';
    
    return len;
}


/* get_v_bases: parses input string and points variable bases to vbase table entry*/
int get_v_bases(char * ref, struct basemap *p_vbases) //TODO: allow compatibility with other IUPAC base identifiers
{
    int i = 0;
    char c = 0;
    int vb_cnt = 0;
    
    for (i = 0; ref[i]; i++) {
        switch (ref[i]) {
            case 'A':
                p_vbases->nts[i] = 'A';
                break;
            case 'a':
                p_vbases->nts[i] = 'A';
                break;
            case 'T':
                p_vbases->nts[i] = 'T';
                break;
            case 't':
                p_vbases->nts[i] = 'T';
                break;
            case 'G':
                p_vbases->nts[i] = 'G';
                break;
            case 'g':
                p_vbases->nts[i] = 'G';
                break;
            case 'C':
                p_vbases->nts[i] = 'C';
                break;
            case 'c':
                p_vbases->nts[i] = 'C';
                break;
            case 'N':
                p_vbases->nts[i] = '*';
                p_vbases->p_vb[i] = &(vb_tbl[N]);
                vb_cnt++;
                break;
            case 'R':
                p_vbases->nts[i] = '*';
                p_vbases->p_vb[i] = &(vb_tbl[R]);
                vb_cnt++;
                break;
            case 'Y':
                p_vbases->nts[i] = '*';
                p_vbases->p_vb[i] = &(vb_tbl[Y]);
                vb_cnt++;
                break;
            default:
                break;
        }
    }
    p_vbases->nts[i] = '\0';
    
    return vb_cnt;
}

/* mk_trgts_file: recursively construct targets using identified variable bases */
void mk_trgts_file(int nxt, char outpt[MAXLEN],struct trgt_vb lcl_bases)
{
    int i = 0;
    int j = 0;
    struct varbase tmp_vb;	//temporary copy of vb_tbl entry

    if (lcl_bases.bs[0]) {	//if variable base was previously found
        lcl_bases.indx++;	//increment lcl_bases.indx
    }
    
    //copy vbases to outpt until variable base is found
    for (i = nxt; vbases.nts[i] != '*' && vbases.nts[i] && i < MAXLEN; i++) {
        outpt[i] = vbases.nts[i];
    }
    
    if (vbases.nts[i] == '*') {								//found variable base
        tmp_vb = *vbases.p_vb[i];							//copy vb_tbl entry to tmp vb
        while (tmp_vb.alt[tmp_vb.crnt]) {					//repeat until all vb_tbl possibilities have been applied
            outpt[i] = tmp_vb.alt[tmp_vb.crnt];				//set output[i] to current vb_tbl entry
            lcl_bases.bs[lcl_bases.indx] = outpt[i];		//record current variable base entry in this recursion path
            lcl_bases.pos[lcl_bases.indx] = i+1;			//record current variable base position in this recursion path
            mk_trgts_file(i+1, outpt, lcl_bases);			//recursively call mk_trgts_file to generate all variants
            tmp_vb.crnt++;
        }
    }
    
    if (!vbases.nts[i]) {	//reached the end of vbase string
        trgt_cnt++;
        fprintf(out, ">cbe_ZTP_3p%d", nt);
        for (j = 0; j < lcl_bases.indx; j++) {
            fprintf(out, "_%c", lcl_bases.bs[j]); //TODO: set up to include position
        }
        fprintf(out, "\n%s%s\n", lnkr, outpt);
    }
    return;
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

