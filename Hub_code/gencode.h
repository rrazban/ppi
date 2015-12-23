/*
 * gencode.h
 *
 *  Created on: Dec 25, 2009
 *      Author: mheo
 */

#ifndef GENCODE_H_
#define GENCODE_H_

/*
Routines for manipulating nucleotide sequences
*/

void CreateRandomNucSequence(int *Seq, int Len);
void CreateRandomNucSequence2(int *Seq, int Len);
/* creates a random sequence of int{0,1,2,3} of length Len */

int PointMutateCharNucSequence(char *Seq, int Len);
int PointMutateNucSequence(int *Seq, int Len);
/*
int PointMutateCharNucSequence(char *Seq, int Len);
int PointMutateNucSequence(int *Seq, int Len);
char/int -> always 0,1,2,3 for nucleotides, not UCAG for char!
return values:
0:  synonymous mutation (no changes in aminoacid sequence)
1:  nonsynonymous mutation
-1: mutation to STOP codon
*/

int NucSeqToAASeq(int *NucSeq, int N, int *AASeq);
int CharNucSeqToAASeq(char *NucSeq, int N, int *AASeq);
/*
translates nucleotide sequence to amino acid sequence
nucleotides: 0123 -> UCAG
amino acids: 01234... -> CMFIL... MJ96 order, see latticelib.c
N is the length of nucleotide sequence
return values:
0: conversion successful
1: stop codon encoutered. STOP codon is -1 in AASeq.
*/

void PrintCharNucCodeSequence(char *buf, char *Seq, int Len);
void PrintNucCodeSequence(char *buf, int *Seq, int Len);
/*input: int* or char* Seq {0,1,2,3}, output char *buf UCAG */

void LetterToNucCodeSeq(char *buf, int *Seq, int Len);
/*input: char *buf UCAG, output int *Seq 0123 */

void CopyIntToCharSeq(char *dest, int *src, int Len);
/*copies Len values from src to dest*/

#endif /* GENCODE_H_ */
