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
void Printout(int criteria);

FILE *out0, *out1;

int structid[GENES];
int nucseq[GENES][NUCSEQLEN];
int aaseq[GENES][AASEQLEN];


int ii, jj;

double pnat[GENES], Enat[GENES];
double pint, Eint;

int pop_size;
double Tenv = 1.0, targetPint = 1.0;
int sim_time;

//change seq here
char tempnucseq[GENES][NUCSEQLEN] = {"UCAGUCUCGCUACCGGUCAUUCGUGCGCGGCUGGAGGAGGAAGAUUCUUCGGACCGUCCAAAGUGGUUGCUUGCCGCUGAU", "CUUAAGCAGGUCCCACCUAGUACGGAGAAAGGCGGGCUCGGUGUAUUAUACCGUAGUGCGAGAUAUAAAUCAUUCAUGGUG"};

int bmode = 127;
int best_bmode;

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

void Printout(int criteria){
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
	int selection=0; //0:fixation 1:Pint/Pnat; 2:Eint/Enat
	int oldnucseq[GENES][NUCSEQLEN];
	double oldpint, oldEint;
	double oldpnat[GENES], oldEnat[GENES];
	double fitness, oldfitness, select, fixation;
	int status;
	int accept, accept1;
	double Tsel=1.0;

	if (argc != 4){
		printf("time | pop_size | seed\n");
		printf("something is missing, bro\n");
		exit(-1);
	}

	ReadCommondata();
	sim_time = atoi(argv[1]);	//order of mag
	pop_size = atof(argv[2]);	
	seed = atoi(argv[3]);
	srand(seed);

	for (jj=0; jj<GENES; jj++){
		LetterToNucCodeSeq(tempnucseq[jj], nucseq[jj], NUCSEQLEN);
		NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		oldpnat[jj] = GetSequencePnat(aaseq[jj], Tenv, structid+jj);
		oldEnat[jj] = SequenceEnergy(aaseq[jj], structid[jj]);
		CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN);
	}
	oldpint = GetBindingP(aaseq[0], structid[0], aaseq[1], structid[1], &(best_bmode), Tenv);
	printf("best initial bmode: %d", best_bmode);
	printf("\nvalue: %f", oldpint);
//	bmode=best_bmode;
	oldpint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
	oldEint = GetBindingEnergy2(aaseq[0], structid[0], aaseq[1], structid[1], bmode);
	oldfitness = oldpnat[0]*oldpnat[1]*oldpint;
	Openfiles();
	if (seed==0) Printinfo();
	Printout(0);
	 // Find stabilizing mutations 
	for(ii=1; ii<pow(10, sim_time); ii++){
		jj=RandomBit();		//randomly choose which protein to mutate
		PointMutateNucSequence(nucseq[jj], NUCSEQLEN);
		status = NucSeqToAASeq(nucseq[jj], NUCSEQLEN, aaseq[jj]);
		if (status==1) {CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}// Reject mutations that introduce a stop codon
		else{
			if (selection==0){
				pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
				pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);
				
				fitness = pnat[0]*pnat[1]*pint;
				select = (fitness - oldfitness)/oldfitness;
				fixation = (1-exp(-2*select))/(1-exp(-2*pop_size*select));
				if (fixation > (double)rand()/RAND_MAX){
					CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldfitness=fitness; Printout(0);
				}
				else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
			}
			else if (selection==1){
				pnat[jj] = GetStructurePnat(aaseq[jj], Tenv, structid[jj]);
				pint = GetBindingP2(aaseq[0], structid[0], aaseq[1], structid[1], bmode, Tenv);

				if(pnat[jj] >= oldpnat[jj] && pint >= oldpint) {CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout(1);} // Accept stabilizing mutation
				else if(pint >= oldpint){
					accept = AcceptOrRejectAttempt((oldpnat[jj]-pnat[jj])/pnat[jj], Tsel);
					if (accept==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout(2);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}
				else if(pnat[jj] >= oldpnat[jj]){
					accept1 = AcceptOrRejectAttempt((oldpint-pint)/pint, Tsel);
					if (accept1==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout(3);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}
				else{
					accept = AcceptOrRejectAttempt2((oldpnat[jj]-pnat[jj])/pnat[jj], (oldpint-pint)/pint, Tsel);
				//	accept = AcceptOrRejectAttempt((oldpnat[jj]-pnat[jj])/pnat[jj], Tenv);
				//	accept1 = AcceptOrRejectAttempt((oldpint-pint)/pint, Tenv);
					if (accept==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldpnat[jj]=pnat[jj]; oldpint=pint; Printout(4);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}		
			
		//		if(pint > targetPint) {break;} // Break if target is reached
			}
			else{
				Enat[jj] = SequenceEnergy(aaseq[jj], structid[jj]);
				Eint = GetBindingEnergy2(aaseq[0], structid[0], aaseq[1], structid[1], bmode);

				if(Enat[jj] <= oldEnat[jj] && Eint <= oldEint) {CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldEnat[jj]=Enat[jj]; oldEint=Eint; Printout(1);} // Accept stabilizing mutation
				else if(Eint <= oldEint){
					accept = AcceptOrRejectAttempt(Enat[jj]-oldEnat[jj], Tenv);
					if (accept==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldEnat[jj]=Enat[jj]; oldEint=Eint; Printout(2);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}
				else if(Enat[jj] <= oldEnat[jj]){
					accept1 = AcceptOrRejectAttempt(Eint-oldEint, Tenv);
					if (accept1==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldEnat[jj]=Enat[jj]; oldEint=Eint; Printout(3);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}
				else{
					accept = AcceptOrRejectAttempt(Enat[jj]-oldEnat[jj], Tenv);
					accept1 = AcceptOrRejectAttempt(Eint-oldEint, Tenv);
					if (accept==1 && accept1==1){CopySeq(oldnucseq[jj], nucseq[jj], NUCSEQLEN); oldEnat[jj]=Enat[jj]; oldEint=Eint; Printout(4);}
					else{CopySeq(nucseq[jj], oldnucseq[jj], NUCSEQLEN);}
				}		
			}
		}
	}
	fclose(out1);
	return 0;
}
