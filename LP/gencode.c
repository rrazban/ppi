/*
 * gencode.c
 *
 *  Genetic code and nucleotide sequences routines
 *  (c) 2005 K.Zeldovich
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#include<stdio.h>
#include<stdlib.h>

#define A_STOP -1

// code as in Miyazawa,Jernigan JMB 1996
#define A_Cys 0
#define A_Met 1
#define A_Phe 2
#define A_Ile 3
#define A_Leu 4
#define A_Val 5
#define A_Trp 6
#define A_Tyr 7
#define A_Ala 8
#define A_Gly 9
#define A_Thr 10
#define A_Ser 11
#define A_Asn 12
#define A_Gln 13
#define A_Asp 14
#define A_Glu 15
#define A_His 16
#define A_Arg 17
#define A_Lys 18
#define A_Pro 19

#define N_UUU 0x0000
#define N_UUC 0x0001
#define N_UUA 0x0002
#define N_UUG 0x0003

#define N_UCU 0x0010
#define N_UCC 0x0011
#define N_UCA 0x0012
#define N_UCG 0x0013

#define N_UAU 0x0020
#define N_UAC 0x0021
#define N_UAA 0x0022
#define N_UAG 0x0023

#define N_UGU 0x0030
#define N_UGC 0x0031
#define N_UGA 0x0032
#define N_UGG 0x0033

#define N_CUU 0x0100
#define N_CUC 0x0101
#define N_CUA 0x0102
#define N_CUG 0x0103

#define N_CCU 0x0110
#define N_CCC 0x0111
#define N_CCA 0x0112
#define N_CCG 0x0113

#define N_CAU 0x0120
#define N_CAC 0x0121
#define N_CAA 0x0122
#define N_CAG 0x0123

#define N_CGU 0x0130
#define N_CGC 0x0131
#define N_CGA 0x0132
#define N_CGG 0x0133


#define N_AUU 0x0200
#define N_AUC 0x0201
#define N_AUA 0x0202
#define N_AUG 0x0203

#define N_ACU 0x0210
#define N_ACC 0x0211
#define N_ACA 0x0212
#define N_ACG 0x0213

#define N_AAU 0x0220
#define N_AAC 0x0221
#define N_AAA 0x0222
#define N_AAG 0x0223

#define N_AGU 0x0230
#define N_AGC 0x0231
#define N_AGA 0x0232
#define N_AGG 0x0233

#define N_GUU 0x0300
#define N_GUC 0x0301
#define N_GUA 0x0302
#define N_GUG 0x0303

#define N_GCU 0x0310
#define N_GCC 0x0311
#define N_GCA 0x0312
#define N_GCG 0x0313

#define N_GAU 0x0320
#define N_GAC 0x0321
#define N_GAA 0x0322
#define N_GAG 0x0323

#define N_GGU 0x0330
#define N_GGC 0x0331
#define N_GGA 0x0332
#define N_GGG 0x0333


// U C A G -> 0 1 2 3


int TripletToAA(unsigned short int triplet)
{
int out;

out = -2;

switch(triplet)
{

case N_UUU: { out =  A_Phe; break; }
case N_UUC: { out =  A_Phe; break; }
case N_UUA: { out =  A_Leu; break; }
case N_UUG: { out =  A_Leu; break; }

case N_UCU: { out =  A_Ser; break; }
case N_UCC: { out =  A_Ser; break; }
case N_UCA: { out =  A_Ser; break; }
case N_UCG: { out =  A_Ser; break; }

case N_UAU: { out =  A_Tyr; break; }
case N_UAC: { out =  A_Tyr; break; }
case N_UAA: { out =  A_STOP; break; }
case N_UAG: { out =  A_STOP; break; }

case N_UGU: { out =  A_Cys; break; }
case N_UGC: { out =  A_Cys; break; }
case N_UGA: { out =  A_STOP; break; }
case N_UGG: { out =  A_Trp; break; }

//-----------

case N_CUU: { out =  A_Leu; break; }
case N_CUC: { out =  A_Leu; break; }
case N_CUA: { out =  A_Leu; break; }
case N_CUG: { out =  A_Leu; break; }

case N_CCU: { out =  A_Pro; break; }
case N_CCC: { out =  A_Pro; break; }
case N_CCA: { out =  A_Pro; break; }
case N_CCG: { out =  A_Pro; break; }

case N_CAU: { out =  A_His; break; }
case N_CAC: { out =  A_His; break; }
case N_CAA: { out =  A_Gln; break; }
case N_CAG: { out =  A_Gln; break; }

case N_CGU: { out =  A_Arg; break; }
case N_CGC: { out =  A_Arg; break; }
case N_CGA: { out =  A_Arg; break; }
case N_CGG: { out =  A_Arg; break; }

//-----------

case N_AUU: { out =  A_Ile; break; }
case N_AUC: { out =  A_Ile; break; }
case N_AUA: { out =  A_Ile; break; }
case N_AUG: { out =  A_Met; break; }

case N_ACU: { out =  A_Thr; break; }
case N_ACC: { out =  A_Thr; break; }
case N_ACA: { out =  A_Thr; break; }
case N_ACG: { out =  A_Thr; break; }

case N_AAU: { out =  A_Asn; break; }
case N_AAC: { out =  A_Asn; break; }
case N_AAA: { out =  A_Lys; break; }
case N_AAG: { out =  A_Lys; break; }

case N_AGU: { out =  A_Ser; break; }
case N_AGC: { out =  A_Ser; break; }
case N_AGA: { out =  A_Arg; break; }
case N_AGG: { out =  A_Arg; break; }

//-----------

case N_GUU: { out =  A_Val; break; }
case N_GUC: { out =  A_Val; break; }
case N_GUA: { out =  A_Val; break; }
case N_GUG: { out =  A_Val; break; }

case N_GCU: { out =  A_Ala; break; }
case N_GCC: { out =  A_Ala; break; }
case N_GCA: { out =  A_Ala; break; }
case N_GCG: { out =  A_Ala; break; }

case N_GAU: { out =  A_Asp; break; }
case N_GAC: { out =  A_Asp; break; }
case N_GAA: { out =  A_Glu; break; }
case N_GAG: { out =  A_Glu; break; }

case N_GGU: { out =  A_Gly; break; }
case N_GGC: { out =  A_Gly; break; }
case N_GGA: { out =  A_Gly; break; }
case N_GGG: { out =  A_Gly; break; }

}

//printf("dec %.3x to %d\n",triplet,out);

if (out==-2){
fprintf(stderr,"genetic code error\n");
fprintf(stderr,"triplet: %.3x\n",triplet);
}

return out;
}
int NucSeqToAASeq(int *NucSeq, int N, int *AASeq)
{
unsigned short int a;
unsigned char a1, a2, a3;
int k=0,i,j=0;

for(i=0;i<N;i+=3)
{
a1 = NucSeq[i];
a2 = NucSeq[i+1];
a3 = NucSeq[i+2];

a = 0;
a = (a1 << 8 ) | (a2 << 4) | a3;

//printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

AASeq[j] = TripletToAA(a);
if (AASeq[j]==-1) { //fprintf(stderr,"%.3x to stop\n",a);*/
 k=1; }
if (AASeq[j]==-2) {
  fprintf(stderr,"NucSeqToAASeq() : Genetic Code Error!!!\n");
  exit(0);
}
j++;
}

return k;
}



int CharNucSeqToAASeq(char *NucSeq, int N, int *AASeq)
{
unsigned short int a;
unsigned char a1, a2, a3;
int k=0,i,j=0;

for(i=0;i<N;i+=3)
{
a1 = NucSeq[i];
a2 = NucSeq[i+1];
a3 = NucSeq[i+2];

a = 0;
a = (a1 << 8 ) | (a2 << 4) | a3;

printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);

AASeq[j] = TripletToAA(a);
if (AASeq[j]==-1) { //fprintf(stderr,"%.3x to stop\n",a);*/
 k=1; }
if (AASeq[j]==-2) {
  fprintf(stderr,"CharNucSeqToAASeq() : Genetic Code Error!!!\n");
   //exit(0);
}
 j++;
}

return k;
}

void CreateRandomNucSequence2(int *Seq, int Len)
{
  unsigned short int a;
  unsigned char a1, a2, a3;
  int i,j, aa;
  for(i=0;i<Len/3;i++) {
    do {
      for(j=0;j<3;j++){
        do{aa = (int)(((double)rand()/RAND_MAX) * 4 );} while(aa==4);
        Seq[3*i+j]=aa;
      }
      a1 = Seq[3*i];
      a2 = Seq[3*i+1];
      a3 = Seq[3*i+2];
      a = 0;
      a = (a1 << 8 ) | (a2 << 4) | a3;
      aa=TripletToAA(a);
      //printf(" %x %x %x -> 0x%.4x\n",a1,a2,a3, a);
    } while(aa==A_STOP);
  }
  return;
}

void CreateRandomNucSequence(int *Seq, int Len)
{
int i,j;
for(i=0;i<Len;i++) {
do{
j = (int)( ((double)rand()/RAND_MAX) * 4 );
} while(j==4);
Seq[i]=j;
}
return;
}

int PointMutateNucSequence(int *Seq, int Len)
{
int i,j,k;
int AAS1[10000], AAS2[10000];

NucSeqToAASeq(Seq, Len, AAS1);

do{
j = (int)( ((double)rand()/RAND_MAX) * Len );
} while(j==Len);

do
  {
  do{
  i = (int)( ((double)rand()/RAND_MAX) * 4 );
  } while(i==4);
  }
while(Seq[j]==i);
Seq[j]=i;

k = NucSeqToAASeq(Seq, Len, AAS2);

if (k==1) return -1;

j=0;
for(i=0; i<Len/3; i++)
if (AAS1[i]!=AAS2[i]) {j=1; break;}

return j;
}



int PointMutateCharNucSequence(char *Seq, int Len)
{
int i,j,k;
int AAS1[10000], AAS2[10000];

CharNucSeqToAASeq(Seq, Len, AAS1);

do{
j = (int)( ((double)rand()/RAND_MAX) * Len );
} while(j==Len);

do
  {
    do{
      i = (int)( ((double)rand()/RAND_MAX) * 4 );
      } while(i==4);
  }
while(Seq[j]==i);
Seq[j]=i;

k = CharNucSeqToAASeq(Seq, Len, AAS2);

if (k==1) return -1;

j=0;
for(i=0; i<Len/3; i++)
if (AAS1[i]!=AAS2[i]) {j=1; break;}
return j;
}



void LetterToNucCodeSeq(char *buf, int *Seq, int Len)
{
int i;
int c;
for(i=0;i<Len;i++)
{
c = -1;
switch (buf[i])
{
case  'U': { c = 0; break; }
case  'C': { c = 1; break; }
case  'A': { c = 2; break; }
case  'G': { c = 3; break; }
}
Seq[i]=c;
}
return;
}


void PrintNucCodeSequence(char *buf, int *Seq, int Len)
{
int i;
char c;
for(i=0;i<Len;i++)
{
c = 'Z';
switch (Seq[i])
{
case  0: { c = 'U'; break; }
case  1: { c = 'C'; break; }
case  2: { c = 'A'; break; }
case  3: { c = 'G'; break; }
}
buf[i]=c;
}
buf[i]=0;
return;
}

void PrintCharNucCodeSequence(char *buf, char *Seq, int Len)
{
int i;
char c;
for(i=0;i<Len;i++)
{
c = 'Z';
switch (Seq[i])
{
case  0: { c = 'U'; break; }
case  1: { c = 'C'; break; }
case  2: { c = 'A'; break; }
case  3: { c = 'G'; break; }
}
buf[i]=c;
}
buf[i]=0;
return;
}


void CopyIntToCharSeq(char *dest, int *src, int Len)
{
int i;
for(i=0; i<Len; i++) dest[i]=src[i];
return;
}
