#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "zlib.h"
#include "../LP/gencode.h"
#include "../LP/latticelib.h"
#include "general.h"
#include "structurelib.h"
#include "../PPI/bindinglib.h"

void Openfiles();
void Printinfo();
void Printout();

FILE *out0, *out1;

int structid[GENES];
int nucseq[GENES][NUCSEQLEN];
int aaseq[GENES][AASEQLEN];


int ii, jj, kk;

double pnat[GENES];
double pint;
double oldfitness;
double fixation;

int pop_size, criteria; //0 for monomers

double Tenv = 1.0;
int sim_time;

//change seq here
//char tempnucseq[GENES][NUCSEQLEN] = {"UCCGACAGGGAAAAAAAGAUAAAGGAGGCUUUGGCCGCUGCCUUAAAGGGAGGAAUCACAGCCUAUCGCCUCGGAUGGGAA", "UUCCGCGGUCACCCGCAAUCUAGGGAAGUACAUUCGCUAGAGCUGCUAUUGGGUAGCUAUGCGCUUAAAAGCUACGCCCUC"};
char tempnucseq[GENES][NUCSEQLEN] = {"UCAGUCUCGCUACCGGUCAUUCGUGCGCGGCUGGAGGAGGAAGAUUCUUCGGACCGUCCAAAGUGGUUGCUUGCCGCUGAU", "CUUAAGCAGGUCCCACCUAGUACGGAGAAAGGCGGGCUCGGUGUAUUAUACCGUAGUGCGAGAUAUAAAUCAUUCAUGGUG"};	//original
//char tempnucseq[GENES][NUCSEQLEN] = {"UCACGUGUAGACGCCGAAUCCCGGUUAUCGCAAACGAGGUUGAUGAUUAAGAUGGGAUGCGAUAACCUAAAGUUUGUUGAA", "AUCGUAGAGAAACCGUACCCAAUUGACUGCCUCCCAGCGUCGAGAAAGGUCCUCAUAUGUGCGUUCAAAGAGUGCUCGCUU"};	//for nat1 and nat2

char tempaaseq[GENES][AASEQLEN];

int bmode = 127;
int best_bmode;

int label;

void Openfiles(){
	char fopbuf[100];
	char rootdir[100];

	sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
 
	sprintf(fopbuf, "%s/v%d.dat", rootdir, label);
    out1 = fopen(fopbuf, "w");
}

void Printinfo(){
	int surface[GENES][AASURFACELEN];
	int s_i;
	char fopbuf[100];
	char rootdir[100];

	sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
 
	sprintf(fopbuf, "%s/info.txt", rootdir);
	out0 = fopen(fopbuf, "w");


//print out aa
	for (jj=0; jj<GENES; jj++){
		PrintAASequence(tempaaseq[jj], aaseq[jj], AASEQLEN);
		fprintf(out0, "%s\n", tempaaseq[jj]);
	}

	fprintf(out0, "Temp: %.3f (kT)\n", Tenv);
	fprintf(out0, "N: %d \n", pop_size);
	fprintf(out0, "criteria: %d \n", criteria);
	for (jj=0; jj<GENES; jj++) fprintf(out0, "Conform%d: %d, %.3E\n", jj, structid[jj], pnat[jj]);
	fprintf(out0, "Bmode: %d, %.3E\n", bmode, pint);
	GetSurfaceAAPositions(structid[0], structid[1], bmode, surface[0], surface[1]);
	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[0][s_i]);
	fprintf(out0, "\n");
	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[1][s_i]);

	fflush(out0);	
	fclose(out0);
}

void Printout(){
	fprintf(out1, "%d", ii);
	for (jj=0; jj<GENES; jj++){
		fprintf(out1, " %.3E", pnat[jj]);
	}	
	fprintf(out1, " %.3E", pint);	
	fprintf(out1, " %.3E", oldfitness);	
	fprintf(out1, " %.3E", fixation);	
	for (jj=0; jj<GENES; jj++){
	//	PrintNucCodeSequence(tempnucseq[jj], nucseq[jj], NUCSEQLEN);	
		fprintf(out1, " ");
		for (kk=0; kk<NUCSEQLEN; kk++){
		fprintf(out1, "%d", nucseq[jj][kk]);
		}
	}	
	fprintf(out1, "\n");
	fflush(out1);
}

//Generates stable AA sequence of conf specified in argv[1]
int main(int argc, char *argv[]){
	int oldnucseq[GENES][NUCSEQLEN];
	double fitness, select;
	int status;
	time_t current_time;

	if (argc != 5){
		printf("time | pop_size | criteria | label\n");
		printf("something is missing, bro\n");
		exit(-1);
	}

	ReadCommondata();
	sim_time = atoi(argv[1]);	//order of mag
	pop_size = atof(argv[2]);	
	criteria = atoi(argv[3]);
	label = atoi(argv[4]);
	current_time = time(NULL);
	fprintf(stderr, ctime(&current_time));
	srand(current_time);

	for (jj=0; jj<GENES; jj++){
		LetterToNucCodeSeq(tempnucseq[jj], nucseq[jj], NUCSEQLEN);
		NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		pnat[jj] = GetSequencePnat(aaseq[jj], Tenv, structid+jj);
		CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN);
	}
	pint = GetBindingP(aaseq[0], structid[0], aaseq[1], structid[1], &(best_bmode), Tenv);
	printf("best initial bmode: %d", best_bmode);
	printf("\nvalue: %f", pint);
	bmode=best_bmode;
	pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
	if (criteria==0){oldfitness = pnat[0]*pnat[1];}
	else {oldfitness = pnat[0]*pnat[1]*pint;}
	Openfiles();
	if (label==0) Printinfo();
	Printout(0);
	 // Find stabilizing mutations 
	for(ii=1; ii<pow(10, sim_time); ii++){
		jj=RandomBit();		//randomly choose which protein to mutate
		PointMutateNucSequence(nucseq[jj], NUCSEQLEN);
		status = NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		if (status==1) {CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}// Reject mutations that introduce a stop codon
		else{
			pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
			pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
			
			if (criteria==0){fitness = pnat[0]*pnat[1];}
			else {fitness = pnat[0]*pnat[1]*pint;}
			select = (fitness - oldfitness)/oldfitness;
			fixation = (1-exp(-2*select))/(1-exp(-2*pop_size*select));
			if (fixation > (double)rand()/RAND_MAX){
				CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldfitness=fitness; Printout();
			}
			else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
		}
	}
	fclose(out1);
	return 0;
}
