//
//  split_P3.c
//  
//
//  Created by Eric Strobel on 10/11/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>

#define MAXLINE 2048	//maximum line length
#define MAXLEN 256
#define DEFAULT_LEN 19	//expected P3 length

//nucleotide definitions
#define A 0
#define U 1
#define G 2
#define C 3
#define TOT 4

/* function prototypes */
int get_file(FILE **ifp, char *ipt);	//open input multiple sequence alignment

int isbase(char c);				//test that character matches IUPAC base notation
int ispair(char nt1, char nt2);	//test for base pair

void count_variant(int base[5], char pos90);	//count nt 90 variant

void print_freq_cnt_to_screen();	//print all frequency counts to screen
void print_freq_cnt(int base[5]);	//print nt 90 frequency counts
void print_datagraph_org();			//print datagraph organization file

int test_P3_contig(char * P3seq); //test for contiguous P3 base pairs

int const_len  = 1; //TODO: consider making this an option

/* position 90 frequency counts */
int frag1_90[5] = {0};
int frag2_90[5] = {0};
int frag3_90[5] = {0};
int ns_90[5] = {0};
int ns1_90[5] = {0};
int ns2_90[5] = {0};
int ns3_90[5] = {0};
int tot_90[5] = {0};

int main(int argc, char *argv[])
{
    //file pointers
    FILE *ifp = NULL;
    FILE *P3_fragile_1_ofp = NULL;
    FILE *P3_fragile_2_ofp = NULL;
    FILE *P3_fragile_3_ofp = NULL;
    FILE *P3_ns_1_ofp = NULL;
    FILE *P3_ns_2_ofp = NULL;
    FILE *P3_ns_3_ofp = NULL;
    FILE *P3_nonslctd_ofp = NULL;
    FILE *source_data = NULL;
    
    //make output files
    if ((P3_nonslctd_ofp = fopen("P3_non.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_fragile_1_ofp = fopen("P3_fragile_3bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_fragile_2_ofp = fopen("P3_fragile_4bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_fragile_3_ofp = fopen("P3_fragile_5bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_ns_1_ofp = fopen("P3_ns_3bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_ns_2_ofp = fopen("P3_ns_4bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((P3_ns_3_ofp = fopen("P3_ns_5bp.fa", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    if ((source_data = fopen("nt90_freq_source.txt", "w")) == NULL) {
        printf("failed to open file\n");
    }
    
    /* parse options using getopt_long */
    int c = -1;
    int option_index = 0;
    
    while (1) {
        static struct option long_options[] =
        {
            {"sequences",       required_argument,  0,  's'},
            {0, 0, 0, 0}
        };
        
        c = getopt_long(argc, argv, "s:", long_options, &option_index);
        
        if (c == -1) {
            break;
        }
        switch (c) {
            case 's': get_file(&ifp, argv[optind-1]); break;
            default: printf("error: unrecognized option. Aborting program...\n"); abort();
        }
    }
    /* end of option parsing */
    
    int i = 0;
    int j = 0;
    
    int match_frag1 = 0;		//flag for fragile_1 match
    int fragile1_cnt = 0;		//number of sequences with match to fragile_1
    int P3_nonslctd_cnt = 0;	//number of sequences without match to fragile_1
        
    char line[MAXLINE] = {0};	//input line
    char name[MAXLEN] = {0};	//sequence name
    char P3seq[MAXLEN] = {0};	//P3 sequence
    
    int start = 0;	//start index of P3 sequence
    int end = 0;	//end indexof P3 sequence
    int len = 0;    //length of P3 segment
    int pos90 = 0;	//index corresponding to Cbe pfl motif position 90 in multiple sequence alignment
    int nt_cnt = 0; //variable for nucleotide counting

    int P3_contig = 0;	//number of contiguous P3 base pairs
    int incld_var = 0;	//flag to include or omit variant

    FILE *ofp = NULL;	//output file pointer

    //print source data header lines
    fprintf(source_data, "\t\tMisfold?:\tyes\t\t\t\t\t\t\t\t\t\t\t\tno\t\t\t\t\t\t\t\t\t\t\t\n");
    fprintf(source_data, "\t\tP3Len:\t3\t\t\t\t4\t\t\t\t5\t\t\t\t3\t\t\t\t4\t\t\t\t5\t\t\t\n");
    fprintf(source_data, "variant\tsequence\tnt\tA\tU\tG\tC\tA\tU\tG\tC\tA\tU\tG\tC\tA\tU\tG\tC\tA\tU\tG\tC\tA\tU\tG\tC\n");
    
    const char fragile_1[20] = {".....G.CC....GGGC.."};
    
    while (fgets(line, MAXLINE, ifp) != NULL) {
        i = 0;
        len = strlen(line);
        if (line[i] != '#' && line[i] != '/') {
            /* get sequence name */
            for (i = 0, j = 0; !isspace(line[i]); i++, j++) {
                name[j] = line[i];
            }
            name[j] = '\0';

            /* find P3 start, magic number len-27 sets array position to last nt in upstream P3 stem
               counting 5 bases upstream finds 5'-most P3 bound */
            for (i = len-27, nt_cnt = 0; nt_cnt <= 4; i--) {
                if (isbase(line[i])) {
                    nt_cnt++;
                }
            }
            start = i+1;
            
            /* find P3 end, magic number sets len-15 array position to first nt in downstream P3 stem
             counting 5 bases downstream finds 53'-most P3 bound */
            for (i = len-15, nt_cnt = 0; nt_cnt < 5; i++) {
                if (isbase(line[i])) {
                    nt_cnt++;
                }
            }
            end = i;
            
            /* get P3 sequence using 5' and 3' bounds, store value of position 90 */
            for (i = start, j = 0; i <= end; i++) {
                if (i == len-24) { //magic number len-24 is position 90 of Cbe pfl in multiple sequence alignment
                    pos90 = j;
                }
                if (isbase(line[i])) {
                    P3seq[j++] = line[i];
                }
            }
            P3seq[j] = '\0';
            
            /* determine number of contiguous P3 base pairs */
            P3_contig = test_P3_contig(P3seq);

			/* determine whether P3 sequence contains 'fragile' sequence by search for fragile_1 */
            for (j = 0, match_frag1 = 1; P3seq[j] && fragile_1[j] && match_frag1; j++) {
                if (fragile_1[j] != '.' && P3seq[j] != fragile_1[j]) {
                    match_frag1 = 0;
                }
            }
            
            /* determine whether P3 sequence matches the predominant conserved length (i.e. no position 90 ambiguity) */
            
            if (const_len) {
                incld_var = (strlen(P3seq) != DEFAULT_LEN) ? 0 : 1;
            }
            
        
            if (incld_var && P3_contig >= 3) {
                count_variant(tot_90, P3seq[pos90]); //add position 90 identity to total frequency count
                fprintf(source_data, "%s\t%s\t%c", name, P3seq, P3seq[pos90]);
                if (match_frag1) {
                    switch (P3_contig) {
                        case 3:
                            ofp = P3_fragile_1_ofp;						//set output to fragile/P3=3
                            count_variant(frag1_90, P3seq[pos90]);		//add pos90 id to fragile/P3=3 freq count
                            break;
                        case 4:
                            fprintf(source_data, "\t\t\t\t"); 			//indent to fragile/P3=4 source data column
                            ofp = P3_fragile_2_ofp;						//set output to fragile/P3=4
                            count_variant(frag2_90, P3seq[pos90]);		//add pos90 id to fragile/P3=4 freq count
                            break;
                        case 5:
                            fprintf(source_data, "\t\t\t\t\t\t\t\t");	//indent to fragile/P3=5 source data column
                            ofp = P3_fragile_3_ofp;						//set output to fragile/P3=5
                            count_variant(frag3_90, P3seq[pos90]);		//add pos90 id to fragile/P3=5 freq count
                            break;
                        default:
                            break;
                    }
                } else {
                    fprintf(source_data, "\t\t\t\t\t\t\t\t\t\t\t\t");	//indent to non-fragile column start
                    count_variant(ns_90, P3seq[pos90]);					//add position 90 id to non-frag freq count
                    switch (P3_contig) {
                        case 3:
                            ofp = P3_ns_1_ofp;							//set output to non-fragile/P3=3
                            count_variant(ns1_90, P3seq[pos90]);		//add pos90 id to non-fragile/P3=3 freq count
                            break;
                        case 4:
                            fprintf(source_data, "\t\t\t\t");			//indent to frag/P3=4 source data column
                            ofp = P3_ns_2_ofp;							//set output to non-fragile/P3=4
                            count_variant(ns2_90, P3seq[pos90]);		//add pos90 id to non-fragile/P3=4 freq count
                            break;
                        case 5:
                            fprintf(source_data, "\t\t\t\t\t\t\t\t");	//indent to non-frag/P3=5 source data column
                            ofp = P3_ns_3_ofp;							//set output to non-fragile/P3=5
                            count_variant(ns3_90, P3seq[pos90]);		//add pos90 id to non-fragile/P3=5 freq count
                            break;
                        default:
                            break;
                    }
                }
                switch (P3seq[pos90]) {	//print '1' to appropriate source data nt column
                    case 'A': fprintf(source_data, "\t1\t\t\t"); break;
                    case 'U': fprintf(source_data, "\t\t1\t\t"); break;
                    case 'G': fprintf(source_data, "\t\t\t1\t"); break;
                    case 'C': fprintf(source_data, "\t\t\t\t1"); break;
                    default:
                        break;
                }
                fprintf(source_data, "\n");
                
                if (P3_contig > 2) { //print to output fasta if P3 > 2
                	fprintf(ofp, ">%s_wt%c\n%s\n", name, P3seq[pos90], P3seq);
                }
            }
        }
    }
    
    //close file pointers
    fclose(P3_nonslctd_ofp);
    fclose(P3_fragile_1_ofp);
    fclose(P3_fragile_2_ofp);
    fclose(P3_fragile_3_ofp);
    fclose(P3_ns_1_ofp);
    fclose(P3_ns_2_ofp);
    fclose(P3_ns_3_ofp);
    fclose(source_data);
    
    //print outputs
    print_freq_cnt_to_screen();
    print_datagraph_org();
}

/* get_file: open file input and assign to file pointer */
int get_file(FILE **ifp, char *ipt)
{
    //printf("getting file\n");
    if ((*ifp = fopen(ipt, "r")) == NULL) {
        printf("error: could not open %s as read one file. Aborting program...\n", ipt);
        abort();
    }
    return 1;
}

/* count_variant: add position 90 to frequency count*/
void count_variant(int base[5], char pos90)
{
    base[TOT]++;
    switch (pos90) {
        case 'A': base[A]++; break;
        case 'U': base[U]++; break;
        case 'G': base[G]++; break;
        case 'C': base[C]++; break;
        default:
            break;
    }
}

/* test_P3_contig: determine number of contiguous P3 pairs */
int test_P3_contig(char * P3seq)
{
    int i = 0;
    int P3_contig = 0;
    int contig = 0;
    int len = strlen(P3seq);
    
    for (i = 4, contig = 1; i >= 0 && contig; i--) {
        if (ispair(P3seq[i], P3seq[len-i-1])) {
            P3_contig++;
        } else {
            contig = 0;
        }
    }
    return P3_contig;
}

/* */
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

/* print_freq_cnt_to_screen: print position 90 frequency information to screen */
void print_freq_cnt_to_screen ()
{
    extern int frag1_90[5];
    extern int frag2_90[5];
    extern int frag3_90[5];
    extern int ns_90[5];
    extern int ns1_90[5];
    extern int ns2_90[5];
    extern int ns3_90[5];
    extern int tot_90[5];
    
    printf("fragile 90 count (type 1: 3bp P3):\t%d\n", frag1_90[TOT]);
    print_freq_cnt(frag1_90);
    printf("fragile 90 count (type 2: 4bp P3):\t%d\n", frag2_90[TOT]);
    print_freq_cnt(frag2_90);
    printf("fragile 90 count (type 3: 5bp P3):\t%d\n", frag3_90[TOT]);
    print_freq_cnt(frag3_90);
    printf("non-selected 90 count (type 1: 3bp P3):\t%d\n", ns1_90[TOT]);
    print_freq_cnt(ns1_90);
    printf("non-selected 90 count (type 2: 4bp P3):\t%d\n", ns2_90[TOT]);
    print_freq_cnt(ns2_90);
    printf("non-selected 90 count (type 3: 5bp P3):\t%d\n", ns3_90[TOT]);
    print_freq_cnt(ns3_90);
    printf("non-selected 90 count:\t%d\n", ns_90[TOT]);
    print_freq_cnt(ns_90);
    printf("tot 90 count:\t%d\n", tot_90[TOT]);
    print_freq_cnt(tot_90);
}

/* print_freq_cnt: print position 90 frequency information */
void print_freq_cnt(int base_90[5])
{
    int i = 0;
    char base[5] = {"AUGC"};
    
    for (i = 0; i < 4; i++) {
        printf("%c:\t%d\t%f\n", base[i], base_90[i], (float)(base_90[i])/(float)(base_90[TOT]));
    }
    
    printf("R:\t%d\t%f\n", base_90[A]+base_90[G],
           (float)(base_90[A])/(float)(base_90[TOT])+(float)(base_90[G])/(float)(base_90[TOT]));
    
    printf("Y:\t%d\t%f\n", base_90[U]+base_90[C],
           (float)(base_90[U])/(float)(base_90[TOT])+(float)(base_90[C])/(float)(base_90[TOT]));
    putchar('\n');
}

/* print_datagraph_org: print position 90 frequency for datagraph */
void print_datagraph_org()
{
    FILE * DG_ofp = NULL;
    if ((DG_ofp = fopen("P3_fragile_datagraph1.txt", "w")) == NULL) {
        printf("failed to open file\n");
    }
    char base[5] = {"AUGC"};
    
    fprintf(DG_ofp, "fA\tfU\tfG\tfC\tnsA\tnsU\tnsG\tnsC\n");
    
    fprintf(DG_ofp, "%f\t\t\t\t%f\t\t\t\n", (float)frag1_90[A]/(float)frag1_90[TOT], (float)ns1_90[A]/(float)ns1_90[TOT]);
    fprintf(DG_ofp, "\t%f\t\t\t\t%f\t\t\n", (float)frag1_90[U]/(float)frag1_90[TOT], (float)ns1_90[U]/(float)ns1_90[TOT]);
    fprintf(DG_ofp, "\t\t%f\t\t\t\t%f\t\n", (float)frag1_90[G]/(float)frag1_90[TOT], (float)ns1_90[G]/(float)ns1_90[TOT]);
    fprintf(DG_ofp, "\t\t\t%f\t\t\t\t%f\n", (float)frag1_90[C]/(float)frag1_90[TOT], (float)ns1_90[C]/(float)ns1_90[TOT]);
    fprintf(DG_ofp, "\n\n");
    fprintf(DG_ofp, "%f\t\t\t\t%f\t\t\t\n", (float)frag2_90[A]/(float)frag2_90[TOT], (float)ns2_90[A]/(float)ns2_90[TOT]);
    fprintf(DG_ofp, "\t%f\t\t\t\t%f\t\t\n", (float)frag2_90[U]/(float)frag2_90[TOT], (float)ns2_90[U]/(float)ns2_90[TOT]);
    fprintf(DG_ofp, "\t\t%f\t\t\t\t%f\t\n", (float)frag2_90[G]/(float)frag2_90[TOT], (float)ns2_90[G]/(float)ns2_90[TOT]);
    fprintf(DG_ofp, "\t\t\t%f\t\t\t\t%f\n", (float)frag2_90[C]/(float)frag2_90[TOT], (float)ns2_90[C]/(float)ns2_90[TOT]);
    fprintf(DG_ofp, "\n\n");
    fprintf(DG_ofp, "%f\t\t\t\t%f\t\t\t\n", (float)frag3_90[A]/(float)frag3_90[TOT], (float)ns3_90[A]/(float)ns3_90[TOT]);
    fprintf(DG_ofp, "\t%f\t\t\t\t%f\t\t\n", (float)frag3_90[U]/(float)frag3_90[TOT], (float)ns3_90[U]/(float)ns3_90[TOT]);
    fprintf(DG_ofp, "\t\t%f\t\t\t\t%f\t\n", (float)frag3_90[G]/(float)frag3_90[TOT], (float)ns3_90[G]/(float)ns3_90[TOT]);
    fprintf(DG_ofp, "\t\t\t%f\t\t\t\t%f\n", (float)frag3_90[C]/(float)frag3_90[TOT], (float)ns3_90[C]/(float)ns3_90[TOT]);
    fprintf(DG_ofp, "\n\n");
    
    fclose(DG_ofp);
}
