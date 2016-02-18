/*
 * =====================================================================================
 *
 *       Filename:  getSubSeq.c
 *
 *    Description:  get sequences from a fasta file according to coordinates
 *
 *        Version:  1.0
 *        Created:  2015/11/09 22时07分39秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JIN Xiaoyang (jinxiaoyang@jinxiaoyang.cn)
 *   Organization:  IBP - CAS
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXLEN 255
#define MAX_LINE_LEN 1024
/* #define DEBUG */

typedef struct faidx{
    unsigned total_len;
    unsigned start_pos;
    unsigned valid_char_per_line;
    unsigned char_per_line;
}faidx;

void usage(const char* arg)
{
    printf("Usage: %s [-q Query_String] [-f Fasta_File]\n", arg);
    printf("\te.g. %s -q chr1:1001-2000 -f genome.fa\n", arg);
}

/* the main function to read fasta index file */
faidx read_fa_idx(FILE *fp, const char* chr)
{
    char tmp_line[MAX_LINE_LEN];
    char tmp_chr[MAXLEN];

    faidx idx={0, 0, 0, 0};

    /* add a <tab> so as to distinct `xx` from `xxx` */
    sprintf(tmp_chr, "%s\t", chr);

    while(!feof(fp))
    {
	fgets(tmp_line, MAX_LINE_LEN, fp);

	if(strncmp(tmp_chr, tmp_line, strlen(tmp_chr))==0)
	{
#ifdef DEBUG
	    printf("idx=%s", tmp_line);
#endif
	    sscanf(tmp_line, "%s\t%u\t%u\t%u\t%u", tmp_chr, &idx.total_len, &idx.start_pos, &idx.valid_char_per_line, &idx.char_per_line);
	    break;
	}
    }
    return idx;
}

unsigned get_line_length(const char* str)
{
    char *tmp_p;
    tmp_p = strchr(str, '\n');
    return (tmp_p - str);
}

/* the main function to read fasta file */
int read_seq(FILE *fp, faidx idx, unsigned pos1, unsigned pos2)
{
#ifdef DEBUG
    printf("faidx.total_len=%u\n", idx.total_len);
    printf("faidx.start_pos=%u\n", idx.start_pos);
    printf("faidx.valid_char_per_line=%u\n", idx.valid_char_per_line);
    printf("faidx.char_per_line=%u\n", idx.char_per_line);
#endif
    char tmp_line[MAX_LINE_LEN];
    unsigned output_len=0;
    unsigned max_inc=idx.valid_char_per_line;

    if(pos1 > idx.total_len)
    {
	fprintf(stderr, "Error: invalid start position, must < %u.\n", idx.total_len);
	return 1;
    }

    if(pos2 > idx.total_len)
    {
	fprintf(stderr, "Warning: invalid end position, should < %u.\n", idx.total_len);
	pos2 = idx.total_len;
    }

    unsigned offset = idx.start_pos + pos1 / idx.valid_char_per_line * idx.char_per_line + pos1 % idx.valid_char_per_line - 1;
    if(pos1 % idx.valid_char_per_line == 0)
    {
	offset -= 1;
    }

#ifdef DEBUG
    printf("offset=%u\n", offset);
#endif

    fseek(fp, offset, SEEK_SET);
    while(!feof(fp))
    {
	fgets(tmp_line, MAX_LINE_LEN, fp);

	max_inc = get_line_length(tmp_line) < max_inc ? get_line_length(tmp_line) : idx.valid_char_per_line;

#ifdef DEBUG
    printf("\nmax_inc=%u\n", max_inc);
#endif

	if(output_len + max_inc >= pos2 - pos1 + 1)
	{
	    tmp_line[pos2-pos1+1-output_len]='\0';
	    printf("%s\n", tmp_line);
	    break;
	}
	else
	{
	    output_len += max_inc;
	    tmp_line[max_inc]='\0';
	    printf("%s", tmp_line);
	}
    }

    return 0;
}

int main(int argc, char * argv[])
{
    int idx, query_n=-1, fafile_n=-1;
    char fafile[MAXLEN];
    char query[MAXLEN];

    char faidxfile[MAXLEN];
    char chr[MAXLEN];
    char tmp_p0[MAXLEN];
    char *tmp_p1, *tmp_p2;
    unsigned pos1, pos2;


    if(argc==1)
    {
	usage(argv[0]);
	return 0;
    }

    for(idx=1; idx<argc; idx++)
    {
	if(strcmp(argv[idx], "-q")==0)
	{
	    query_n = ++idx;

	    /* if `-q` is the last parament */
	    if(query_n >= argc)
	    {
		usage(argv[0]);
		return 1;
	    }
	}

	else if(strcmp(argv[idx], "-f")==0)
	{
	    fafile_n = ++idx;

	    /* if `-f` is the last parament */
	    if(fafile_n >= argc)
	    {
		usage(argv[0]);
		return 1;
	    }
	}

	/* if `query_n` isn't set, then use the remained parament */
	else if(query_n == 0)
	{
	    query_n = idx;
	}
    }

    if(query_n < 0 || fafile_n < 0)
    {
	usage(argv[0]);
	return 1;
    }
    else
    {
	/* `.fa` file */
	strcpy(fafile, argv[fafile_n]);

	/* `.fa.fai` file */
	strcpy(faidxfile, argv[fafile_n]);
	strcat(faidxfile, ".fai");

	/* get `chr:pos1-pos2` */
	strcpy(query, argv[query_n]);
	if ((tmp_p1 = strchr(query, ':')) && (tmp_p2 = strchr(tmp_p1, '-')))
	{
	    strncpy(chr, query, tmp_p1 - query);
	    strncpy(tmp_p0, tmp_p1 + 1, tmp_p2 - tmp_p1 - 1);
	    pos1=atoi(tmp_p0);
	    pos2=atoi(tmp_p2 + 1);
	}
	else
	{
	    usage(argv[0]);
	    return 1;
	}

#ifdef DEBUG
	printf("chr=%s\n", chr);
	printf("pos1=%d\n", pos1);
	printf("pos2=%d\n", pos2);
	printf("fafile=%s\n", fafile);
	printf("faidxfile=%s\n", faidxfile);
#endif
    }

    FILE *fp_fa, *fp_faidx;
    if((fp_faidx = fopen(faidxfile, "r")))
    {
	faidx idx=read_fa_idx(fp_faidx, chr);

	if (idx.total_len == 0)
	{
	    fprintf(stderr, "Error: invalid CHR name - %s\n", chr);
	    return 1;
	}

	if((fp_fa = fopen(fafile, "r")))
	{
	    if(read_seq(fp_fa, idx, pos1, pos2))
	    {
		fprintf(stderr, "FAIL\n");
	    }
	}
	else
	{
	    fprintf(stderr, "Error: cannot open file %s\n", fafile);
	    return 1;
	}
    }
    else
    {
	fprintf(stderr, "Error: cannot open file %s\n", faidxfile);
	fprintf(stderr, "You may use command 'samtools faidx %s' first\n", fafile);
	return 1;
    }
    return 0;
}
