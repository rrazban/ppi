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

int ii, jj;

double pnat[GENES];
double pint;

double Tenv = 0.75, targetPint = 1.0;

//change seq here
char tempnucseq[GENES][NUCSEQLEN] = {"UCAGUCUCGCUACCGGUCAUUCGUGCGCGGCUGGAGGAGGAAGAUUCUUCGGACCGUCCAAAGUGGUUGCUUGCCGCUGAU", 
									 "CUUAAGCAGGUCCCACCUAGUACGGAGAAAGGCGGGCUCGGUGUAUUAUACCGUAGUGCGAGAUAUAAAUCAUUCAUGGUG"}; 
int bmode = 1;

int seed;

void Openfiles(){
	char fopbuf[100];
	char rootdir[100];

	sprintf(rootdir,"/n/regal/shakhnovich_lab/rrazban/");
 
	sprintf(fopbuf, "%s/v%d.dat", rootdir, seed);
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

	fprintf(out0, "Temp: %.3f (kT)\n", Tenv);
	for (jj=0; jj<GENES; jj++) fprintf(out0, "Conform%d: %d\n", jj, structid[jj]);
	fprintf(out0, "Bmode: %d\n", bmode);
	GetSurfaceAAPositions(structid[0], structid[1], bmode, surface[0], surface[1]);
	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[0][s_i]);
	fprintf(out0, "\n");
	for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out0, " %02d", surface[1][s_i]);

	fflush(out0);	
	fclose(out0);
}

void Printout(){
	fprintf(out1, "%d\t", ii);
	for (jj=0; jj<GENES; jj++){
		NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
		fprintf(out1, " %.3E", pnat[jj]);
		PrintNucCodeSequence(tempnucseq[jj], nucseq[jj], NUCSEQLEN); // Convert numeric nuc sequence to letter nuc sequence
		fprintf(out1, " %s", tempnucseq[jj]);
	}	
	pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
	fprintf(out1, " %.3E\n", pint);	
	fflush(out1);
}

//Generates stable AA sequence of conf specified in argv[1]
int main(int argc, char *argv[]){
	int oldnucseq[GENES][NUCSEQLEN];
	double oldpint;
	double oldpnat[GENES];
	int status;
	int accept, accept1;

	if (argc != 2){
		printf("seed is missing, bro\n");
		exit(-1);
	}

	ReadCommondata();
	seed = atoi(argv[1]);
	srand(seed);

	for (jj=0; jj<GENES; jj++){
		LetterToNucCodeSeq(tempnucseq[jj], nucseq[jj], NUCSEQLEN);
		NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		oldpnat[jj] = GetSequencePnat(aaseq[jj], Tenv, structid+jj);
		CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN);
	}	
	oldpint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);

	Openfiles();
	if (seed==0) Printinfo();
	Printout();
	 // Find stabilizing mutations 
	for(ii=1; ii<pow(10, 4); ii++){
		jj=RandomBit();		//randomly choose which protein to mutate
		PointMutateNucSequence(nucseq[jj], NUCSEQLEN);
		status = NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		if (status==1) {CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}// Reject mutations that introduce a stop codon
		else{
			pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
			pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);

			if(pnat[jj] >= oldpnat[jj] && pint >= oldpint) {CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout();} // Accept stabilizing mutation
			else if(pint >= oldpint){
				accept = AcceptOrRejectAttempt((oldpnat[jj]-pnat[jj])/pnat[jj], Tenv);
				if (accept==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout();} // Accept stabilizing mutation
				else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
			}
			else if(pnat[jj] >= oldpnat[jj]){
				accept1 = AcceptOrRejectAttempt((oldpint-pint)/pint, Tenv);
				if (accept1==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout();} // Accept stabilizing mutation
				else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
			}
			else{
				accept = AcceptOrRejectAttempt((oldpnat[jj]-pnat[jj])/pnat[jj], Tenv);
				accept1 = AcceptOrRejectAttempt((oldpint-pint)/pint, Tenv);
				if (accept==1 && accept1==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout();} // Accept stabilizing mutation
				else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
			}		
		}
		if(pint > targetPint) {break;} // Break if target is reached
	}
	fclose(out1);
	return 0;
}
