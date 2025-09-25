#ifndef _MM_H_
#define _MM_H_

#define MM_MAXLINE 1000

#include "HPNBM/util.h"

struct mmdata {

	int N, M, NNZ;
	int *x;
	int *y;
	double *v;

	int symmetricity;
	int binary;
	int ndiagonal;
	int realnnz;

	char *header;
};


int initialize_mm(char *infile, struct mmdata *mm);
int initialize_mm_compact(char *infile, struct mmdata *mm);
void printmm(char *infile, struct mmdata *mm, char *outfile);
void csr2mm(size_t* ia, size_t* ja, size_t rows, size_t cols, char* out_file);
void freemm(struct mmdata *mm);
void reordermm(struct mmdata *mm, int *rowordervec, int *colordervec, struct mmdata *mmordered);

void mm2csr(struct mmdata *mm, size_t* &ia, size_t* &ja, size_t &rows, size_t &cols);

#endif 