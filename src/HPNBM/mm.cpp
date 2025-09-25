#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "HPNBM/mm.h"


#define STRINGSIZE 1000

struct point {

	int x, y;
};


struct graphdata {

	int nvtxs;
	int nedges;
	int nconst;

	int *xadj;
	int *adjncy;
	int *vwgt;
	int *gids;

};

int comp_point(const void *pb, const void *pc) {

	struct point *b = (struct point *)pb;
	struct point *c = (struct point *)pc;
	
	return (b->x < c->x)? -1: (b->x > c->x)? 1: (b->y <= c->y)? -1: 1;
}

int initialize_mm_compact(char *infile, struct mmdata *mm) {

	char filename[STRINGSIZE];
	char line[MM_MAXLINE];
	char *saveptr, *token;
	sprintf(filename, "%s", infile);
	
	FILE *f = fopen(filename, "r");
	
	mm -> symmetricity = 0;
	
	fgets(line, MM_MAXLINE, f);
	if(line[0] == '%') {
	
		token = strtok_r(line, " \n", &saveptr);
		while(token != NULL) {
			if(!strcmp(token, "symmetric")) {
				mm -> symmetricity = 1;
				break;
			}
			token = strtok_r(NULL, " \n", &saveptr);
		}

		fgets(line, MM_MAXLINE, f);	
		while(line[0] == '%') {
			fgets(line, MM_MAXLINE, f);	
		}
	}
		
	sscanf(line, " %d %d %d", &mm->N, &mm->M, &mm->NNZ);
	
	mm -> ndiagonal = 0;
	mm -> x = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> y = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> v = NULL;	

	fgets(line, MM_MAXLINE, f);
	double v;
	if(2 == sscanf(line, " %d %d %lf", &mm -> x[0], &mm -> y[0], &v)) 
		mm -> binary = 1;
	else {
		mm->binary = 0;
		mm -> v = (double *)malloc(mm->NNZ * sizeof(double));	
		mm->v[0] = v;
	}
	mm -> x[0] --;
	mm -> y[0] --;
	if(mm->symmetricity && mm -> x[0] == mm -> y[0])
		mm -> ndiagonal ++;

	int *usedx = (int *)calloc(mm->N, sizeof(int));
	int *usedy = (int *)calloc(mm->M, sizeof(int));

	usedx[mm->x[0]] = usedy[mm->y[0]] = 1;
	
	int i;
	for(i=1; i<mm->NNZ; i++) {

		fgets(line, MM_MAXLINE, f);
		if(! mm -> binary) 
			sscanf(line, " %d %d %lf", &mm -> x[i], &mm -> y[i], &mm -> v[i]);
		else	
			sscanf(line, " %d %d", &mm -> x[i], &mm -> y[i]);

		mm -> x[i] --;
		mm -> y[i] --;
	
		usedx[mm->x[i]] = usedy[mm->y[i]] = 1;
	
		if(mm->symmetricity)
			usedy[mm->x[i]] = usedx[mm->y[i]] = 1;
				

		if(mm->symmetricity && mm -> x[i] == mm -> y[i])
			mm -> ndiagonal ++;
	}

	int cnt;
	for(cnt=i=0; i<mm->N; i++) 
		if(usedx[i])
			usedx[i] = cnt++;
	mm->N = cnt;
	for(cnt=i=0; i<mm->M; i++)
		if(usedy[i])
			usedy[i] = cnt++;
	mm->M = cnt;

	for(i=0; i<mm->NNZ; i++) 
		mm->x[i] = usedx[mm->x[i]];
	for(i=0; i<mm->NNZ; i++) 
		mm->y[i] = usedy[mm->y[i]];
	
	free(usedx);
	free(usedy);
	mm -> realnnz = (!mm -> symmetricity)? mm->NNZ: 2 * mm->NNZ - mm->ndiagonal;

	return 0;
}

int initialize_mm(char *infile, struct mmdata *mm) {

	char filename[STRINGSIZE];
	
	char line[MM_MAXLINE];
	char *saveptr, *token;
	sprintf(filename, "%s", infile);
	
	FILE *f = fopen(filename, "r");
	
	mm -> symmetricity = 0;

	mm->header = (char *)malloc(sizeof(char) * (1<<20));
	memset(mm->header, '\0', sizeof(char) * (1<<20));
	char *hptr = mm->header;
	
	fgets(line, MM_MAXLINE, f);
	if(line[0] == '%') {

		int tmp = strlen(line);
		strcpy(hptr, line);
		hptr += tmp;
	
		token = strtok_r(line, " \n", &saveptr);
		while(token != NULL) {
			if(!strcmp(token, "symmetric")) {
				mm -> symmetricity = 1;
				break;
			}
			token = strtok_r(NULL, " \n", &saveptr);
		}

		fgets(line, MM_MAXLINE, f);	
		while(line[0] == '%') {			
			tmp = strlen(line);
			strcpy(hptr, line);
			hptr += tmp;
			fgets(line, MM_MAXLINE, f);
		}
	}
		
	sscanf(line, " %d %d %d", &mm->N, &mm->M, &mm->NNZ);
	
	mm -> ndiagonal = 0;
	mm -> x = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> y = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> v = NULL;	

	fgets(line, MM_MAXLINE, f);
	double v;
	if(2 == sscanf(line, " %d %d %lf", &mm -> x[0], &mm -> y[0], &v)) 
		mm -> binary = 1;
	else {
		mm->binary = 0;
		mm -> v = (double *)malloc(mm->NNZ * sizeof(double));	
		mm->v[0] = v;
	}
	mm -> x[0] --;
	mm -> y[0] --;

	if(mm->symmetricity && mm -> x[0] == mm -> y[0])
		mm -> ndiagonal ++;
	
	int i;
	for(i=1; i<mm->NNZ; i++) {

		fgets(line, MM_MAXLINE, f);
		if(! mm -> binary) 
			sscanf(line, " %d %d %lf", &mm -> x[i], &mm -> y[i], &mm -> v[i]);
		else	
			sscanf(line, " %d %d", &mm -> x[i], &mm -> y[i]);

		mm -> x[i] --;
		mm -> y[i] --;

		if(mm->symmetricity && mm -> x[i] == mm -> y[i])
			mm -> ndiagonal ++;
	}

	mm -> realnnz = (!mm -> symmetricity)? mm->NNZ: 2 * mm->NNZ - mm->ndiagonal;

	return 0;
}

void printmm(char *infile, struct mmdata *mm, char *outfile) {
	
	FILE *fout = fopen(outfile, "w");	
	FILE *fin = fopen(infile, "r");
	char line[MM_MAXLINE];

	fgets(line, MM_MAXLINE, fin);
	while(line[0] == '%')
	{
		fprintf(fout, "%s", line);
		fgets(line, MM_MAXLINE, fin);
	}	
	fclose(fin);	

	fprintf(fout, "%d %d %d\n", mm->N, mm->M, mm->NNZ);
	int i;
	for(i=0; i<mm->NNZ; i++)
		if(mm->binary)
			fprintf(fout, "%d %d\n", mm->x[i]+1, mm->y[i]+1);
		else
			fprintf(fout, "%d %d %lf\n", mm->x[i]+1, mm->y[i]+1, mm->v[i]);
	fclose(fout);
}

void csr2mm(size_t* ia, size_t* ja, size_t rows, size_t cols, char* out_file) {
	struct mmdata *mm = (struct mmdata *)malloc(sizeof(struct mmdata));

	mm -> N = rows;
	mm -> M = cols;
	mm -> NNZ = ia[rows];

	mm -> ndiagonal = 0;
	mm -> x = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> y = (int *)malloc(mm->NNZ * sizeof(int));
	mm -> v = NULL;	
	mm -> binary = 1;
	mm -> symmetricity = 0;

	int idx = 0;
	for(int i=0; i<rows; i++) {
		for(int j=ia[i]; j<ia[i+1]; j++) {
			mm->x[idx] = i;
			mm->y[idx] = ja[j];
			if(i == ja[j])
				mm->ndiagonal++;
			idx++;
		}
	}

	mm -> realnnz = (!mm -> symmetricity)? mm->NNZ: 2 * mm->NNZ - mm->ndiagonal;
	FILE *fout = fopen(out_file, "w");
	char line[MM_MAXLINE];

	fprintf(fout, "%d %d %d\n", mm->N, mm->M, mm->NNZ);
	int i;
	for(i=0; i<mm->NNZ; i++)
		if(mm->binary)
			fprintf(fout, "%d %d\n", mm->x[i]+1, mm->y[i]+1);
		else
			fprintf(fout, "%d %d %lf\n", mm->x[i]+1, mm->y[i]+1, mm->v[i]);
	fclose(fout);

}

void reordermm(struct mmdata *mm, int *rowordervec, int *colordervec, struct mmdata *mmordered){

	mmordered -> N = mm -> N;
	mmordered -> M = mm -> M;
	mmordered -> NNZ = mm -> NNZ;

	mmordered -> x = (int *)malloc(sizeof(int) * mmordered->NNZ);
	mmordered -> y = (int *)malloc(sizeof(int) * mmordered->NNZ);
	mmordered -> v = (double *)malloc(sizeof(double) * mmordered->NNZ);

	mmordered -> binary = mm -> binary;
	mmordered -> symmetricity = 0;
	mmordered -> ndiagonal = 0;
	mmordered -> realnnz = mm -> realnnz;

	int i, nnz=0;
	for(i=0; i<mm->NNZ; i++) {

		mmordered -> x[nnz] = rowordervec[ mm -> x[i] ];
		mmordered -> y[nnz] = colordervec[ mm -> y[i] ];
		if(!mmordered -> binary)
			mmordered -> v[nnz] = mm -> v[i];
		nnz ++;

		/* if(mm -> symmetricity && mm -> x[i] != mm -> y[i]) { */

		/* 	mmordered -> x[nnz] = rowordervec[ mm -> y[i] ]; */
		/* 	mmordered -> y[nnz] = colordervec[ mm -> x[i] ]; */
		/* 	if(!mmordered -> binary) */
		/* 		mmordered -> v[nnz] = mm -> v[i]; */
		/* 	nnz ++; */
		/* } */
	}
}

void mm2csr(struct mmdata *mm, size_t* &ia, size_t* &ja, size_t &rows, size_t &cols) {

	extern int wflag;

	int i,x,y;
	
	struct graphdata *g = (struct graphdata *)malloc(sizeof(struct graphdata));

	g -> nvtxs = mm -> N;
	g -> xadj = (int *)calloc(g->nvtxs+2, sizeof(int));
	
	// struct point *points = (struct point *)calloc(mm -> NNZ * 2, sizeof(struct point));
	struct point *points = (struct point *)calloc((uint64_t)(mm -> NNZ) * 2, sizeof(struct point));
	int npoints = 0;
	for(i=0; i<mm->NNZ; i++) {
	
		x = mm->x[i];
		y = mm->y[i];
		
		if(x != y) {
			points[npoints].x = x<y? x: y;
			points[npoints].y = x<y? y: x;
			npoints ++;
		}
	}
	
	qsort(points, npoints, sizeof(struct point), comp_point);
	
	int xprev = -1;
	int yprev = -1;
	for(i=0; i<npoints; i++) {

		x = points[i].x;
		y = points[i].y;
		if(x == xprev && y == yprev)
			continue;
	
		g -> xadj[x+2] ++;
		g -> xadj[y+2] ++;
		
		xprev = x;
		yprev = y;
	}

	for(i=2; i<g->nvtxs+2; i++)
		g -> xadj[i] += g -> xadj[i-1];

	g -> nedges = g->xadj[g->nvtxs+1]/2;
	g -> adjncy = (int *)malloc(g->xadj[g->nvtxs+1] * sizeof(int));
	g -> vwgt = (int *)calloc(g -> nvtxs * g -> nconst, sizeof(int));	

	xprev = yprev = -1;
	for(i=0; i<npoints; i++) {

		x = points[i].x;
		y = points[i].y;
		if(x == xprev && y == yprev)
			continue;
	
		g -> adjncy[g -> xadj[x+1] ++] = y;
		g -> adjncy[g -> xadj[y+1] ++] = x;
		g -> vwgt[x*g->nconst]++;
		g -> vwgt[y*g->nconst]++;
		
		xprev = x;
		yprev = y;
	}

	free(points);

	g -> gids = (int *)malloc(sizeof(int) * g -> nvtxs);
	for(i=0; i<g->nvtxs; i++)
		g -> gids[i] = i;

	rows = g -> nvtxs;
	cols = mm -> M;
	ia = (size_t *)malloc(sizeof(size_t) * (rows + 1));
	ja = (size_t *)malloc(sizeof(size_t) * (g -> xadj[rows+1]));

	ia[0] = 0;
	for(i=0; i<rows; i++)
		ia[i+1] = g -> xadj[i+1];
	for(i=0; i<g->xadj[rows+1]; i++)
		ja[i] = g -> adjncy[i];

}


void freemm(struct mmdata *mm) {

	free(mm->x);
	free(mm->y);
	free(mm->v);
	free(mm);
}