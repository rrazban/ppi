/*************************
 Modified by Orit Peleg
 opeleg@fas.harvard.edu
 on 2/26/13
 ************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/times.h>
#include <time.h>
#include <errno.h>
#include <stddef.h>

#include "zlib.h"
#include "define.h"
#include "../LP/gencode.h"
#include "../LP/latticelib.h"
#include "bindinglib.h"
#include "cell.h"
#include "compare.h"


char aacode[] = {
    'C', 'M', 'F', 'I', 'L',
    'V', 'W', 'Y', 'A', 'G',
    'T', 'S', 'N', 'Q', 'D',
    'E', 'H', 'R', 'K', 'P'
};

int aamat[AASEQLEN][20];
int aaseq[AASEQLEN], aaseq2[AASEQLEN];
int domi_species;

int curr_MAXGENES;
int curr_MAXPPIS;
int curr_MAXSTATES;
int allow_fold_change;
int allow_gene_exp;
int allow_unfolded_states;
int allow_chaps;
float chaps_x;
int hub_ID;
int singlish;
double state_3rd;
int POST_Proc;
int sequenceversion;
int solo;
int homo;

//int DumpProteome(char *filename); //un-used
//int RecountOrgDB(int divisioncycle); //un-used
//int DumpDB(parameter *myParam, int divisioncycle); //un-used

int ResetOrgDB(parameter *myParam, int divisioncycle);
void SetMeanOrg(organism *mean);
void AddMeanOrg(organism *mean, int who);
void GetMeanOrg(organism *mean, int orgcount);
int seq_entropy(char **seq, int N, double *entropy);


/*time evaluation functions:*/
void start_clock(void);
int end_clock();

/* time evaluation:*/
static clock_t st_time;
static clock_t en_time;
static struct tms st_cpu;
static struct tms en_cpu;

int time2;
int runTime;

clock_t HDP_time_start, HDP_time_end;
double HDP_cpu_time_used;
clock_t HDP_time_start2, HDP_time_end2;
double HDP_cpu_time_used2;

void   WriteConfig();
void   ReadConfig();


//output variables:

void PrintOutput();
void Flushfiles();
void Closefiles();
void Openfiles();
void PrepareOutput();
void PostProcessing();


FILE *error_op;
FILE *it_solver;
FILE *out1, *out2, *out3, *out5, *out6, *out7, *out8, *out10, *out11, *finitvals, *out12, *out15, *out16, *out40, *out41;
FILE *out42,*out43,*out44,*out45;


FILE *out50;


//FILE *zout4;
FILE *out4;

parameter myParam;
organism mean;

/* parameters for sequence entropy calculation */
double entropy[MAXGENES][AASEQLEN], entropy_sum[MAXGENES];
char *seq[MAXGENES][MAXORGANISMS];
double entropy_sumtot;


int divisioncycle, lastdecimtime;
int mutatororigin[7];
int orgcount;
int mutatorcount, count;

char seqbuf[800];
int sizeRank[MAXORGANISMS] = {0};
int speciesSizeTotal;

int generation;
double TIME;

int RUN_BASED_ON_CONFIG=0;

int start_divisioncycle;

/***********************************************************************************************************
 **********************************************************************************************************/
int main(int argc, char *argv[]){
    
    int kk;
    int ii, jj, kia;
    int who, status;
    int count1, count2, prev_orgcount, newborncount;
    int fate, overflowflag;
    int deathcount, death2count, maxorgind;
    double dC;
    double dt,cum_br[MAXORGANISMS],max_br,r1,r2;
    
    
    
    
    /*initialize time evaluation:*/
    start_clock();
    ReadConfig();
    
    printf("OPEN...\n"); 
    
    SetupParameter(argc, argv, &myParam, &orgcount);
    
    
    if (POST_Proc==1){
        fprintf(stderr,"Start post-proccessing ...\n"); fflush(stdout);
        PostProcessing();
        fprintf(stderr,"Done post-proccessing, exit(0)\n"); fflush(stdout);
        exit(0); 
    }
    
    Openfiles();
    if (RUN_BASED_ON_CONFIG==0){start_divisioncycle = 0;} else{start_divisioncycle=divisioncycle;}
    
    fprintf(stderr,"Start simulations...\n");
    
    divisioncycle=start_divisioncycle;
    ResetOrgDB(&myParam, divisioncycle);
    RankSpeciesSizeDB(sizeRank, nOrgDB);
    
    
    printf("Starting Output at zero step...\n");
    
    PrepareOutput(); PrintOutput();
    
    
    printf("Starting Main Loop...\n");
    
    
        
    
    
    // MAIN LOOP
    /***********************************************************************************************************/
    for (divisioncycle=start_divisioncycle+1; divisioncycle <= myParam.maxdivcycle; divisioncycle++) {
        if (end_clock() < runTime){
            
      //      if (divisioncycle>1000) myParam.printoutcycle = 250;
    //        if (divisioncycle>10000) myParam.printoutcycle = 500;
  //          if (divisioncycle>1000000) myParam.printoutcycle = 1000;
            
            //divisioncycle=1;while (end_clock() < runTime && divisioncycle <= myParam.maxdivcycle){
            switch (ALGORITHM) {
                case 0: //Muyoung
                    overflowflag = 0; count1=count2=0; prev_orgcount=orgcount;
                    //if(divisioncycle == myParam.timeLow+1) myParam.birthrate /= 10;
                    for(who=0; who<MAXORGANISMS; who++) {
                        if (myOrgstatus[who]==S_DEAD) continue;
                        if (myOrgstatus[who]==S_NEWBORN) continue;
                        if (myOrgstatus[who]==S_TODIE) continue;
                        if (myOrgstatus[who]==S_TODIE2) continue;
                        if (myOrgstatus[who]!=S_ALIVE) {fprintf(stderr,"bad org status\n"); exit(1);}
                        // only living ones are processed
                        
                        if (allow_gene_exp==1){
                            do { fate = (int) ((double) rand() / RAND_MAX * 3); }
                            while(fate == 3);
                            // randomly assigns a fate: 0, 1, or 2
                        }
                        else{
                            do { fate = (int) ((double) rand() / RAND_MAX * 2); }
                            while(fate == 2);
                            // randomly assigns a fate: 0, or 1
                        }
                        
                        
                        if (fate == 0) {
                            status = OrgChildBirth(&myParam, who, divisioncycle, 3);
                            if (status > 0) { count1 += status;
                            } else if(status == -9) {
                                overflowflag = 1;}
                        } else if (fate == 1) {status=OrgDeath(&myParam, who, 3);
                        } else if (fate == 2) {status=GeneExpress(&myParam, who, 3, &dC);
                            //myOrg[who].action_count[3] += sqrt(dC*dC);
                        }
                    }
                    
                    deathcount=death2count=maxorgind=newborncount=0;
                    count=0;
                    for(who=0;who<MAXORGANISMS;who++) {
                        if(myOrgstatus[who] == S_NEWBORN){ myOrgstatus[who] = S_ALIVE;newborncount++;orgcount++;}
                        if(myOrgstatus[who] == S_TODIE) { deathcount++; orgcount--; KillOrganism(who);}
                        if(myOrgstatus[who] == S_TODIE2) {death2count++;KillOrganism(who);}
                        if(myOrgstatus[who] == S_ALIVE) {count++;maxorgind=who;}
                    }
                    if(count!=orgcount) {
                        fprintf(stderr, "t=%d inconsistent orgcount %d count %d\n", divisioncycle, orgcount, count);
                        fprintf(stderr, "# of new born cells : %d %d\n", count1, newborncount);
                        fprintf(stderr, "# of killed cells : %d, %d\n", deathcount, death2count);
                        exit(1);
                    }
                    if(orgcount<=0) {fprintf(stdout, "died out at gen %d\n", divisioncycle); fprintf(out1, "died out at gen %d\n", divisioncycle);
                        //fclose(zout4);
                        exit(1);}
                    if(overflowflag) {fprintf(stdout, "organism overflow %d\n", orgcount); break;}
                    
                    /* decimation */
                    kia=0;
                    if(orgcount>myParam.decimthresh) {
                        count1=0;
                        do{
                            ii = (int)(((double)rand()/RAND_MAX)*(maxorgind+1));
                            if(myOrgstatus[ii]==S_ALIVE){KillOrganism(ii);myOrgstatus[ii]=S_DEAD;kia++;count--;}
                        } while(count>myParam.decimto);
                        orgcount=count;
                        lastdecimtime=divisioncycle;
                    }
                    
                    TIME = divisioncycle;
                    break;
                    
                    
                case 1: //Gillespie
                    generation = divisioncycle;
                    for (orgcount=0; orgcount < POPSIZE; orgcount++) {
                        cum_br[0]=myOrg[0].birthrate;
                        for(who=1; who<MAXORGANISMS; who++) {
                            cum_br[who] = cum_br[who-1] + myOrg[who].birthrate;
                        }
                        max_br=cum_br[MAXORGANISMS-1];
                        r1=( (double) rand()/RAND_MAX );
                        r2=( (double) rand()/RAND_MAX );
                        if (r1 < 1.0e-50){ r1 = 1.0e-50;}
                        dt=-(1.0e0/max_br)*log(r1);
                        TIME=TIME+dt;
                        r2=r2*max_br;
                        for(who=0; who<MAXORGANISMS; who++) {
                            if ( r2 < cum_br[who]) { break ; }
                        }
                        if ( r2 >= max_br ){ who = MAXORGANISMS-1 ; }
                        
                        //find org to kill:
                        do{ ii = (int)( ((double) rand()/RAND_MAX) * MAXORGANISMS );
                        } while( ii==MAXORGANISMS || ii==who );
                        
                        //kill it:
                        memmove(&myOrg[ii], &myOrg[who], sizeof(organism));
                        
                        if (allow_gene_exp==1){
                            do{ kk = (int)( ((double) rand()/RAND_MAX) * 2 );
                            } while( kk==2 );
                            // randomly assigns a fate: 0, 1
                        }
                        else{kk = 0;} // randomly assigns a fate: 0 ...
                        
                        if(kk==0){
                            //printf("before OrgGeneMutate ... \n");
                            OrgGeneMutate(&myParam,who,2);
                        }
                        if(kk==1){
                            GeneExpress(&myParam,who,2, &dC);
                        }
                        
                        if (allow_gene_exp==1){
                            do{ jj = (int)( ((double) rand()/RAND_MAX) * 2 );
                            } while( jj==2 );
                            // randomly assigns a fate: 0, 1
                        }
                        else{jj = 0;} // randomly assigns a fate: 0 ...
                        
                        if(jj==0){
                            OrgGeneMutate(&myParam,ii,2);
                        }
                        if(jj==1){
                            GeneExpress(&myParam,ii,2, &dC);
                        }
                    }
                    
                    break;
                default:
                    break;
            }
            
            
            
            //output
            if ((divisioncycle % myParam.printoutcycle == 0) || (divisioncycle%myParam.dumpcycle == 0)){PrepareOutput();}
            if ((divisioncycle % myParam.printoutcycle == 0)) {PrintOutput(); WriteConfig();}
            
            Flushfiles();
        }//end if time condition
        else{
            printf("reached time limit %f min\n", (runTime/6000.0));
            break;
        }
    } //end divisioncycle
    
    WriteConfig();
    fprintf(error_op, "1\n");
    Closefiles();
    return 0;
}


/***********************************************************************************************************
 **********************************************************************************************************/
void PostProcessing(){
    int t_pp;
    int generation_pp, Nclon_pp, PopSize_pp;//; , nOrgDB_pp; //, diff_pp;
    int index_pp,ClonSize_pp;
    int ii,j;
    int tmp_pnat,br_pp;
    
    char fopbuf[200];
    char seq_pp[1000];

    int aaseq_pp[AASEQLEN], aaseq2_pp[AASEQLEN],aaseq_face[AASURFACELEN],aaseq_face_hub[AASURFACELEN],aaseq_face_par[AASURFACELEN];
    int chargeseq_face_hub[AASURFACELEN],chargeseq_face_par[AASURFACELEN];

    int fi,gi,int_i;
    
    int int_genome[MAXGENES*NUCSEQLEN];
    
    int eofReached;
    
    /*properties to calculate: */
    
    double surface_frac_hydro[curr_MAXGENES][6], surface_net_charge[curr_MAXGENES][6],
        surface_pos_charge[curr_MAXGENES][6], surface_neg_charge[curr_MAXGENES][6], surface_frac_charge[curr_MAXGENES][6];
    
    double int_surface_frac_hydro[curr_MAXGENES], int_surface_net_charge[curr_MAXGENES],
        int_surface_pos_charge[curr_MAXGENES], int_surface_neg_charge[curr_MAXGENES], int_surface_frac_charge[curr_MAXGENES];
    
    int int_surface_counter[curr_MAXGENES];
    
    double non_int_surface_frac_hydro[curr_MAXGENES], non_int_surface_net_charge[curr_MAXGENES],
        non_int_surface_pos_charge[curr_MAXGENES], non_int_surface_neg_charge[curr_MAXGENES], non_int_surface_frac_charge[curr_MAXGENES];
    
    int interaction_partners[curr_MAXPPIS][2];
    int interaction_surfaces[curr_MAXPPIS][2];
    int interaction_bmode[curr_MAXPPIS];
    
    int protein_structid[curr_MAXGENES];
    int protein_interacting_faced[curr_MAXGENES][6];
    
    double interaction_energy[curr_MAXPPIS];
    double interaction_energy_hydro[curr_MAXPPIS];
    double interaction_energy_charge[curr_MAXPPIS];
    //double non_interaction_energy[curr_MAXGENES*];
    
    
    int inter_surface_charge[curr_MAXPPIS][2][AASURFACELEN];
    int inter_surface_hydro[curr_MAXPPIS][2][AASURFACELEN];
    int inter_surface_hydro_avil[curr_MAXPPIS][2][AASURFACELEN];
    
    int iii;
    
    
    
    //printf("\nin postproc...\n");fflush(stdout);
    //exit(0);
    
    for(gi=0;gi<curr_MAXGENES;gi++){
        protein_structid[gi] = myOrg[1].structid[gi];
    }
    
    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
        interaction_partners[int_i][0]=myOrg[1].ppi_pair[int_i][0];
        interaction_partners[int_i][1]=myOrg[1].ppi_pair[int_i][1];
        interaction_bmode[int_i] = myOrg[1].bmode[int_i];
        interaction_surfaces[int_i][0] = interaction_bmode[int_i]/24; 
        interaction_surfaces[int_i][1] = (interaction_bmode[int_i]%24)/4;
        printf("interaction_surfaces[int_i][0] = %d, interaction_surfaces[int_i][1]=%d\n",interaction_surfaces[int_i][0],interaction_surfaces[int_i][1]);
        
    }
    
    for(gi=0;gi<curr_MAXGENES;gi++){
        for(fi=0;fi<6;fi++){
            protein_interacting_faced[gi][fi]=0;
        }
    }
    
    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
        protein_interacting_faced[interaction_partners[int_i][0]][interaction_surfaces[int_i][0]]=1;
        protein_interacting_faced[interaction_partners[int_i][1]][interaction_surfaces[int_i][1]]=1;
    }
    
    //printf("beep1\n");fflush(stdout);
    
    sprintf(fopbuf,"seqlog-%s.dat", myParam.targetname);
    printf("%s\n", fopbuf);fflush(stdout);
    out4=fopen(fopbuf,"r");
    
    //printf("beep1.1\n");fflush(stdout);
    
    sprintf(fopbuf,"faceslog-%s.dat", myParam.targetname);
    out40=fopen(fopbuf,"w");
    
    //printf("beep1.2\n");fflush(stdout);
    
    sprintf(fopbuf,"bindingenergies-%s.dat", myParam.targetname);
    out41=fopen(fopbuf,"w");
    
    //printf("beep1.3\n");fflush(stdout);
    
    sprintf(fopbuf,"surfacemaps_charge-%s.dat", myParam.targetname);
    out42=fopen(fopbuf,"w");
    
    //printf("beep1.35\n");fflush(stdout);
    
    sprintf(fopbuf,"surfacemaps_hydro_avil-%s.dat", myParam.targetname);
    out43=fopen(fopbuf,"w");
    
    //printf("beep1.4\n");fflush(stdout);
    
    sprintf(fopbuf,"surfacemaps_hydro-%s.dat", myParam.targetname);
    out44=fopen(fopbuf,"w");
    
    //printf("beep1.5\n");fflush(stdout);
    
    sprintf(fopbuf,"non_funct_energies-%s.dat", myParam.targetname);
    out45=fopen(fopbuf,"w");
    
    sprintf(fopbuf,"all_surfaces_stats-%s.dat", myParam.targetname);
    out50=fopen(fopbuf,"w");
    
    //printf("beep2\n");fflush(stdout);
    
    for (t_pp=0; t_pp < divisioncycle ; t_pp++){
    //for (t_pp=0; t_pp < 500 ; t_pp++){
        //printf("time=%d\n",t_pp);
        eofReached = (fscanf(out4, "%d %4d %4d\n",&generation_pp, &Nclon_pp, &PopSize_pp)==EOF);
        //printf("generation_pp %d, Nclon_pp %d, PopSize_pp %d\n", generation_pp, Nclon_pp, PopSize_pp);
        //printf("eofReached = %d\n",eofReached);
        
        
        if (eofReached) {
            printf("Post Processing DONE \n\n");
            fclose(out40); fclose(out41);
            exit(0);
            
            
        }
        //printf("Read first line \n\n");
        
        for(ii=0;ii<Nclon_pp;ii++) {
            fscanf(out4, "%4d %4d %1d ",&index_pp,&ClonSize_pp,&br_pp);
            fscanf(out4, "%s ",seq_pp);
            
            //fprintf(stdout,"SEQ : %s\n", seq_pp);fflush(stdout);

            
            LetterToNucCodeSeq(seq_pp, int_genome, 7*NUCSEQLEN);
            //PrintNucCodeSequence(buf, int_genome, 7*NUCSEQLEN); //just checking ..
            //fprintf(stdout,"buf : %s\n", buf);fflush(stdout);
            
            
            for(j=0;j<curr_MAXGENES;j++) {
                fscanf(out4, "%5d ", &tmp_pnat);
            }
            fscanf(out4, "\n");
            
            if (ii==0){ //first clone

                
                for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                    NucSeqToAASeq(int_genome+interaction_partners[int_i][0]*NUCSEQLEN,NUCSEQLEN,aaseq_pp);
                    NucSeqToAASeq(int_genome+interaction_partners[int_i][1]*NUCSEQLEN,NUCSEQLEN,aaseq2_pp);
                    
                    interaction_energy[int_i] = GetBindingEnergy2 (aaseq_pp, protein_structid[interaction_partners[int_i][0]], aaseq2_pp, protein_structid[interaction_partners[int_i][1]], interaction_bmode[int_i]);
                    interaction_energy_hydro[int_i] = GetBindingEnergyHydro (aaseq_pp, protein_structid[interaction_partners[int_i][0]], aaseq2_pp, protein_structid[interaction_partners[int_i][1]], interaction_bmode[int_i]);
                    interaction_energy_charge[int_i] = GetBindingEnergyCharge (aaseq_pp, protein_structid[interaction_partners[int_i][0]], aaseq2_pp, protein_structid[interaction_partners[int_i][1]], interaction_bmode[int_i]);
                }
                
                fprintf(out41, "%d %d", generation_pp, generation_pp);
                for(int_i=0;int_i<curr_MAXPPIS;int_i++) fprintf(out41, " %E", interaction_energy[int_i]);
                fprintf(out41,"\n");
                
                
                
                
                if (generation_pp > (divisioncycle - 10000)){
                    
                    //printf("beep1\n");
                    //exit(0);
                    
                    
                    ////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////////////////////
                    
                    //print all non functional binding energies:
                    
                    int gi_1, gi_2;
                    int face1, face2, rotate;
                    int k, surfacetmp[9],surface1[9], surface2[9];
                    double e, e_hydro, e_charge, e_hydro_neg, e_charge_neg;
                    
                    double min_energy = -50; //-49;;
                    double max_energy = 1;
                    double energy_interval = 0.05;
                    int energy_histogram[3000];
                    int energy_histogram_non_funct[3000];
                    int energy_histogram_funct[3000];
                    
                    int energy_histogram_funct_hydro[3000];
                    int energy_histogram_funct_charge[3000];
                    int energy_histogram_non_funct_hydro[3000];
                    int energy_histogram_non_funct_charge[3000];
                    int energy_histogram_non_funct_hydro_neg[3000];
                    int energy_histogram_non_funct_charge_neg[3000];
                    
                    double energy_cum;
                    int energy_hist_i;
                    
                    int curr_bmode;
                    int is_funct_int;
                    
                    int surface1_t_charge[9], surface2_t_charge[9];
                    int surface1_t_hydro[9], surface2_t_hydro[9];
                    
                    
                    energy_hist_i=0;
                    
                    
                    energy_hist_i=0;
                    energy_cum = min_energy;
                    while (energy_cum<=max_energy){
                        energy_histogram[energy_hist_i] = 0;
                        energy_histogram_non_funct[energy_hist_i] = 0;
                        energy_histogram_funct[energy_hist_i] = 0;
                        energy_histogram_funct_hydro[energy_hist_i] = 0;
                        energy_histogram_funct_charge[energy_hist_i] = 0;
                        energy_histogram_non_funct_hydro[energy_hist_i] = 0;
                        energy_histogram_non_funct_charge[energy_hist_i] = 0;
                        energy_histogram_non_funct_hydro_neg[energy_hist_i] = 0;
                        energy_histogram_non_funct_charge_neg[energy_hist_i] = 0;
                        
                        energy_cum += energy_interval;
                        energy_hist_i++;
                    }
                    
                    
//                    interaction_partners[int_i][0]=myOrg[1].ppi_pair[int_i][0];
//                    interaction_partners[int_i][1]=myOrg[1].ppi_pair[int_i][1];
//                    interaction_bmode[int_i]
                    
                    
                    for(gi_1=0; gi_1<curr_MAXGENES; gi_1++){
                        for(gi_2=gi_1; gi_2<curr_MAXGENES; gi_2++){
                            NucSeqToAASeq(int_genome+gi_1*NUCSEQLEN,NUCSEQLEN,aaseq_pp);
                            NucSeqToAASeq(int_genome+gi_2*NUCSEQLEN,NUCSEQLEN,aaseq2_pp);
                            
                            
                            
                            
                            for(face1=0;face1<6;face1++){
                                for(k=0; k<9; k++) surfacetmp[k] = aaseq_pp[(int) AllFaces[24*protein_structid[gi_1]+4*face1+0][k]];
                                MirrorWall(surface1, surfacetmp);
                                for(face2=0;face2<6;face2++){
                                    for(rotate=0;rotate<4;rotate++){
                                        
                                        
                                        for(k=0; k<9; k++) surface2[k] = aaseq2_pp[(int) AllFaces[24*protein_structid[gi_2]+4*face2+rotate][k]];
                                        
                                        ConvertAAtoCharge(surface1, surface1_t_charge, AASURFACELEN);
                                        ConvertAAtoCharge(surface2, surface2_t_charge, AASURFACELEN);
                                        ConvertAAtoHydro(surface1, surface1_t_hydro, AASURFACELEN);
                                        ConvertAAtoHydro(surface2, surface2_t_hydro, AASURFACELEN);
                                        
                                        e=0.0; e_charge=0.0; e_hydro=0.0;e_charge_neg=0.0; e_hydro_neg=0.0;
                                        for(k=0; k<9; k++) {
                                            e+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;;
                                            //if (((surface1_t_charge[k]==1) & (surface2_t_charge[k]==-1)) | ((surface1_t_charge[k]==-1) & (surface2_t_charge[k]==1))){
                                            if (((surface1_t_charge[k]==1) & (surface2_t_charge[k]==-1)) |
                                                    ((surface1_t_charge[k]==-1) & (surface2_t_charge[k]==1))){
                                                e_charge+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;;
                                            }
                                            if ((surface1_t_hydro[k]==1) & (surface2_t_hydro[k]==1)){
                                                e_hydro+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;;
                                                
                                            }
                                            
                                            if (((surface1_t_charge[k]==1) & (surface2_t_charge[k]==1)) |
                                                ((surface1_t_charge[k]==-1) & (surface2_t_charge[k]==-1))){
                                                e_charge_neg+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;;
                                            }
                                            if (((surface1_t_hydro[k]==1) & (surface2_t_hydro[k]==0)) |
                                                ((surface1_t_hydro[k]==0) & (surface2_t_hydro[k]==1))){
                                                e_hydro_neg+=EnergyMatrix[surface1[k]][surface2[k]]-3.16;;
                                            }
                                            
                                            //printf ("e_charge_neg=%f e_hydro_neg=%f\n",e_charge_neg,e_hydro_neg);

                                            
                                        }
                                        //exit(0);
                                        
                                        // check if interaction is functional:
                                        curr_bmode=face1*24+face2*4+rotate;
                                        is_funct_int = 0;
                                        for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                                            if (hub_ID!=10){
                                                if((interaction_partners[int_i][0]==gi_1)&(interaction_partners[int_i][1]==gi_2)&(interaction_bmode[int_i]==curr_bmode)){
                                                    is_funct_int = 1;
                                                }
                                            }
                                            else{

                                                    //take into account only first pair
                                                    if (int_i==0){
                                                        if((interaction_partners[int_i][0]==gi_1)&(interaction_partners[int_i][1]==gi_2)&(interaction_bmode[int_i]==curr_bmode)){
                                                            is_funct_int = 1;
                                                        }
                                                    }

                                                
                                                
                                            }
                                                
                                        }
                                        
                                        
                                        //total energy:                                        
                                        e-=min_energy;
                                        e_hydro-=min_energy;
                                        e_charge-=min_energy;
                                        e_hydro_neg-=min_energy;
                                        e_charge_neg-=min_energy;
                                        
                                        
                                        
//                                        energy_histogram[energy_hist_i] = 0;
//                                        energy_histogram_non_funct[energy_hist_i] = 0;
//                                        energy_histogram_funct[energy_hist_i] = 0;
//                                        energy_histogram_funct_hydro[energy_hist_i] = 0;
//                                        energy_histogram_funct_charge[energy_hist_i] = 0;
//                                        energy_histogram_non_funct_hydro[energy_hist_i] = 0;
//                                        energy_histogram_non_funct_charge[energy_hist_i] = 0;
                                        
                                        energy_histogram[(int) floor(e/energy_interval)]++;
                                        
                                        if (is_funct_int == 1){
                                            energy_histogram_funct[(int) floor(e/energy_interval)]++;
//                                            e_hydro /=e;
//                                            e_charge /=e;
                                            energy_histogram_funct_hydro[(int) floor(e_hydro/energy_interval)]++;
                                            energy_histogram_funct_charge[(int) floor(e_charge/energy_interval)]++;
                                        }
                                        else{
//                                            e_hydro /=e;
//                                            e_charge /=e;
//                                            e_hydro_neg /=e;
//                                            e_charge_neg /=e;
                                            energy_histogram_non_funct[(int) floor(e/energy_interval)]++;
                                            energy_histogram_non_funct_hydro[(int) floor(e_hydro/energy_interval)]++;
                                            energy_histogram_non_funct_charge[(int) floor(e_charge/energy_interval)]++;
                                            energy_histogram_non_funct_hydro_neg[(int) floor(e_hydro_neg/energy_interval)]++;
                                            energy_histogram_non_funct_charge_neg[(int) floor(e_charge_neg/energy_interval)]++;
                                        }
                                    }
                                }
                            }
                            
                            
                            
                            
                        }
                    }
                    
                    energy_hist_i=0;
                    energy_cum = min_energy;
                    while (energy_cum<=max_energy){
                        fprintf(out45, "%f %d %d %d %d %d %d %d %d %d\n", energy_cum,
                                energy_histogram[energy_hist_i],
                                energy_histogram_funct[energy_hist_i],energy_histogram_funct_hydro[energy_hist_i],energy_histogram_funct_charge[energy_hist_i],
                                energy_histogram_non_funct[energy_hist_i],energy_histogram_non_funct_hydro[energy_hist_i],energy_histogram_non_funct_charge[energy_hist_i],
                                energy_histogram_non_funct_hydro_neg[energy_hist_i],energy_histogram_non_funct_charge_neg[energy_hist_i]);
                        energy_cum += energy_interval;
                        energy_hist_i++;
                    }
                    
                    ////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////////////////////

                
                    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                        // bmode mirror rotate
                        GetSingleSurfaceAA(aaseq_face_hub, aaseq_pp, protein_structid[interaction_partners[int_i][0]], interaction_surfaces[int_i][0],interaction_bmode[int_i],1,0);
                        GetSingleSurfaceAA(aaseq_face_par, aaseq_pp, protein_structid[interaction_partners[int_i][1]], interaction_surfaces[int_i][1],interaction_bmode[int_i],0,1);
                        
                        
                        NucSeqToAASeq(int_genome+interaction_partners[int_i][0]*NUCSEQLEN,NUCSEQLEN,aaseq_pp);
                        NucSeqToAASeq(int_genome+interaction_partners[int_i][1]*NUCSEQLEN,NUCSEQLEN,aaseq2_pp);

                        
                        //GetDoubleSurfaceAA(int *seq1, int struct1, int *seq2, int struct2, int bmode, int *surface1, int *surface2);
                        GetDoubleSurfaceAA(aaseq_pp, protein_structid[interaction_partners[int_i][0]], aaseq2_pp, protein_structid[interaction_partners[int_i][1]], interaction_bmode[int_i], aaseq_face_hub, aaseq_face_par);
                        
                        
                        
                        ConvertAAtoCharge(aaseq_face_hub, chargeseq_face_hub, AASURFACELEN);
                        ConvertAAtoCharge(aaseq_face_par, chargeseq_face_par, AASURFACELEN);
                        
//                        //control:
//                        printf("hub seq\n");
//                        for (iii=0;iii<AASURFACELEN;iii++){
//                            printf("%d %d\n",aaseq_face_hub[iii],chargeseq_face_hub[iii]);
//                        }
                        
                        for (iii=0;iii<AASURFACELEN;iii++){
                            inter_surface_charge[int_i][0][iii] = chargeseq_face_hub[iii];
                            inter_surface_charge[int_i][1][iii] = chargeseq_face_par[iii];
                        }
                        
                        ConvertAAtoAVIL(aaseq_face_hub, chargeseq_face_hub, AASURFACELEN);
                        ConvertAAtoAVIL(aaseq_face_par, chargeseq_face_par, AASURFACELEN);

                        for (iii=0;iii<AASURFACELEN;iii++){
                            inter_surface_hydro_avil[int_i][0][iii] = chargeseq_face_hub[iii];
                            inter_surface_hydro_avil[int_i][1][iii] = chargeseq_face_par[iii];
                        }

                        
                        ConvertAAtoHydro(aaseq_face_hub, chargeseq_face_hub, AASURFACELEN);
                        ConvertAAtoHydro(aaseq_face_par, chargeseq_face_par, AASURFACELEN);
                        
                        for (iii=0;iii<AASURFACELEN;iii++){
                            inter_surface_hydro[int_i][0][iii] = chargeseq_face_hub[iii];
                            inter_surface_hydro[int_i][1][iii] = chargeseq_face_par[iii];
                        }

                        
                    }
                    
                    //fprintf(out42, "%d %d\n", generation_pp, generation_pp);
                    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                        fprintf(out42, "%f %f ", interaction_energy[int_i], interaction_energy_charge[int_i]);
                        
                        
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out42, "%d ", inter_surface_charge[int_i][0][iii]);
                        }
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out42, "%d ", inter_surface_charge[int_i][1][iii]);
                        }
                        
                        fprintf(out42,"\n");
                    }
                    fprintf(out42,"\n");
                    
                    //fprintf(out43, "%d %d\n", generation_pp, generation_pp);
                    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                        fprintf(out43, "%f %f ", interaction_energy[int_i],interaction_energy_hydro[int_i]);
                        
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out43, "%d ", inter_surface_hydro_avil[int_i][0][iii]);
                        }
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out43, "%d ", inter_surface_hydro_avil[int_i][1][iii]);
                        }
                        
                        fprintf(out43,"\n");
                    }
                    fprintf(out43,"\n");
                    
                    //fprintf(out43, "%d %d\n", generation_pp, generation_pp);
                    for(int_i=0;int_i<curr_MAXPPIS;int_i++){
                        fprintf(out44, "%f ", interaction_energy[int_i]);
                        
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out44, "%d ", inter_surface_hydro[int_i][0][iii]);
                        }
                        for (iii=0;iii<AASURFACELEN;iii++){
                            fprintf(out44, "%d ", inter_surface_hydro[int_i][1][iii]);
                        }
                        
                        fprintf(out44,"\n");
                    }
                    fprintf(out44,"\n");
                    
                    
                    exit(0);
                }
                
                                
                
                
                
//                //OP TO DO, define a new function with pre defined bmode?
//                double GetBindingK(int *seq1, int struct1, int *seq2, int struct2, double T){
//                    int face1, face2, rotate;
//                    int k, surfacetmp[9],surface1[9], surface2[9];
//                    double e, z=0, emin=1e10;
//                    
//                    for(face1=0;face1<6;face1++)
//                    {
//                        for(k=0; k<9; k++) surfacetmp[k] = seq1[(int) AllFaces[24*struct1+4*face1+0][k]];
//                        MirrorWall(surface1, surfacetmp);
//                        for(face2=0;face2<6;face2++)
//                            for(rotate=0;rotate<4;rotate++)
//                            {
//                                e=0;
//                                for(k=0; k<9; k++) surface2[k] = seq2[(int) AllFaces[24*struct2+4*face2+rotate][k]];
//                                for(k=0; k<9; k++) e+=EnergyMatrix[surface1[k]][surface2[k]];
//                                z+=exp(-e/T);
//                                if (e<emin) emin=e;
//                            }
//                    }
//                    
//                    //if(z<1.0e-12) { fprintf(stderr,"Error!!! Partition Function z = %8.3lf\n", z); exit(1); }
//                    
//                    return z;
//                }

                
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                
                
                for(gi=0;gi<curr_MAXGENES;gi++){
                    NucSeqToAASeq(int_genome+gi*NUCSEQLEN,NUCSEQLEN,aaseq_pp);
                    
//                    for (pppi=0;pppi<AASEQLEN;pppi++){
//                        printf("%d ", aaseq_pp[pppi]);
//                    }
//                    printf("\n\n");
                    
                    for(fi=0;fi<6;fi++){
                        GetSingleSurfaceAA(aaseq_face, aaseq_pp, protein_structid[gi], fi,0,0,0);
                        
//                        for (pppi=0;pppi<AASURFACELEN;pppi++){
//                            printf("%d ", aaseq_face[pppi]);
//                        }
//                        printf("\n\n");
                        
                        
                        //fprintf(stdout,"aaseq_face SEQ : %d\n", aaseq_face[0]);fflush(stdout);
                        surface_frac_hydro[gi][fi] = GetFracHydrophobicity_avil(aaseq_face, AASURFACELEN);
                        surface_frac_charge[gi][fi] = GetFracCharge(aaseq_face, AASURFACELEN);
                        surface_net_charge[gi][fi] = GetNetCharge(aaseq_face, AASURFACELEN);
                        surface_pos_charge[gi][fi] = GetPosCharge(aaseq_face, AASURFACELEN);
                        surface_neg_charge[gi][fi] = GetNegCharge(aaseq_face, AASURFACELEN);
                    }
                }
                
                //exit(0);
                
                for(gi=0;gi<curr_MAXGENES;gi++){
                    
                    int_surface_frac_hydro[gi] = 0.0;
                    int_surface_frac_charge[gi] = 0.0;
                    int_surface_net_charge[gi] = 0.0;
                    int_surface_pos_charge[gi] = 0.0;
                    int_surface_neg_charge[gi] = 0.0;
                    
                    non_int_surface_frac_hydro[gi] = 0.0;
                    non_int_surface_frac_charge[gi] = 0.0;
                    non_int_surface_net_charge[gi] = 0.0;
                    non_int_surface_pos_charge[gi] = 0.0;
                    non_int_surface_neg_charge[gi] = 0.0;

                    
                    
                    int_surface_counter[gi]=0;
                    for(fi=0;fi<6;fi++){
                        if (protein_interacting_faced[gi][fi]==1){
                            int_surface_frac_hydro[gi] += surface_frac_hydro[gi][fi];
                            int_surface_frac_charge[gi] += surface_frac_charge[gi][fi];
                            int_surface_net_charge[gi] += surface_net_charge[gi][fi];
                            int_surface_pos_charge[gi] += surface_pos_charge[gi][fi];
                            int_surface_neg_charge[gi] += surface_neg_charge[gi][fi];
                            int_surface_counter[gi] ++; 
                        }
                        else {
                            non_int_surface_frac_hydro[gi] += surface_frac_hydro[gi][fi];
                            non_int_surface_frac_charge[gi] += surface_frac_charge[gi][fi];
                            non_int_surface_net_charge[gi] += surface_net_charge[gi][fi];
                            non_int_surface_pos_charge[gi] += surface_pos_charge[gi][fi];
                            non_int_surface_neg_charge[gi] += surface_neg_charge[gi][fi];
                        }
                    }
                    
                    if (int_surface_counter[gi]>0) {
                        
                        int_surface_frac_hydro[gi] /=(double) int_surface_counter[gi];
                        int_surface_frac_charge[gi] /=(double) int_surface_counter[gi];
                        int_surface_net_charge[gi] /=(double) int_surface_counter[gi];
                        int_surface_pos_charge[gi] /=(double) int_surface_counter[gi];
                        int_surface_neg_charge[gi] /=(double) int_surface_counter[gi];
                        
                    }
                    
                    if ((6-int_surface_counter[gi])>0) {
                        
                        non_int_surface_frac_hydro[gi] /=(double) (6-int_surface_counter[gi]);
                        non_int_surface_frac_charge[gi] /=(double) (6-int_surface_counter[gi]);
                        non_int_surface_net_charge[gi] /=(double) (6-int_surface_counter[gi]);
                        non_int_surface_pos_charge[gi] /=(double) (6-int_surface_counter[gi]);
                        non_int_surface_neg_charge[gi] /=(double) (6-int_surface_counter[gi]);
                        
                    }
           
                }//gi
                
                fprintf(out40, "%d %d", generation_pp, generation_pp);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", int_surface_frac_hydro[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", int_surface_frac_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", int_surface_net_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", int_surface_pos_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", int_surface_neg_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", non_int_surface_frac_hydro[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", non_int_surface_frac_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", non_int_surface_net_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", non_int_surface_pos_charge[gi]);
                for (gi=0; gi<curr_MAXGENES; gi++) fprintf(out40, " %E", non_int_surface_neg_charge[gi]);
                fprintf(out40,"\n");
                
                fprintf(out50, "%d %d 0 0 0 0\n", generation_pp, generation_pp);
                for(gi=0;gi<curr_MAXGENES;gi++){
                    for(fi=0;fi<6;fi++){
                        fprintf(out50, "%d ", protein_interacting_faced[gi][fi]);
                    }
                    fprintf(out50,"\n");
                }
                
                for(gi=0;gi<curr_MAXGENES;gi++){
                    for(fi=0;fi<6;fi++){fprintf(out50, "%E ", surface_frac_hydro[gi][fi]);}fprintf(out50,"\n");
                }

                for(gi=0;gi<curr_MAXGENES;gi++){
                    for(fi=0;fi<6;fi++){fprintf(out50, "%E ", surface_frac_charge[gi][fi]);}fprintf(out50,"\n");
                }

                for(gi=0;gi<curr_MAXGENES;gi++){
                    for(fi=0;fi<6;fi++){fprintf(out50, "%E ", surface_net_charge[gi][fi]);}fprintf(out50,"\n");
                }


                
            }//end first clone
            
        }
        fscanf(out4, ";\n\n");
        
        //eofReached = (fscanf(out4, ";\n\n")==EOF);
        
                 
        //printf("beep END\n");fflush(stdout);
        
       // exit(0);
    }
    
    

}





/***********************************************************************************************************
 **********************************************************************************************************/
void PrepareOutput(){
    int ii, who, m, jj;
    
    ResetOrgDB(&myParam, divisioncycle);
    RankSpeciesSizeDB(sizeRank, nOrgDB);
    
    /* statistics for organisms */
    m=count=mutatorcount=0;
    for(ii=0;ii<curr_MAXGENES;ii++) mutatororigin[ii]=0;
    SetMeanOrg(&mean);
    for(who=0;who<POPSIZE;who++){
        if(myOrgstatus[who]!=S_ALIVE) continue;
        AddMeanOrg(&mean, who);
        if(myOrg[who].mutlevel == 1) mutatorcount++;
        mutatororigin[myOrg[who].mutorigin]++;
        for (ii=0; ii<curr_MAXGENES; ii++) seq[ii][m] = myOrg[who].genome+ii*NUCSEQLEN;
        
        m++;
    }
    GetMeanOrg(&mean,m);
    
    entropy_sumtot = 0.0;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        seq_entropy(seq[ii], m, entropy[ii]);
        
        entropy_sum[ii]=0.0;
        for(jj=0;jj<AASEQLEN;jj++) entropy_sum[ii] += entropy[ii][jj];
        entropy_sum[ii] /= (double)AASEQLEN;
        entropy_sumtot += entropy_sum[ii];
    }
    entropy_sumtot /= (double)curr_MAXGENES;
}

/***********************************************************************************************************
 **********************************************************************************************************/
void WriteConfig(){
    int i, j, who;
    
    FILE *fp=fopen(FILE_CONFIG, "w");
    
    // variables from evo-cell
    fprintf(fp,"%d %d %d %d %d %d %d %d\n", divisioncycle, lastdecimtime, orgcount, mutatorcount, generation, speciesSizeTotal, domi_species, nOrgDB);
    fprintf(fp,"%E\n", TIME);
    
    
    // variables from cell.h: myParam, myOrg
    fprintf(fp,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",myParam.seed, myParam.startcode, myParam.orgcount, myParam.maxdivcycle, myParam.decimthresh, myParam.decimto, myParam.initpop, myParam.dumpcycle, myParam.printoutcycle, myParam.plotoutcycle, myParam.screenoutcycle, myParam.seqlogcycle, myParam.timeLow, myParam.timeHigh);
    fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", myParam.Tenv, myParam.tol, myParam.TLow, myParam.THigh, myParam.birthrate, myParam.deathrate, myParam.expressrate, myParam.alpha, myParam.pnatcutoff, myParam.speciessizefactor, myParam.mutrate[0], myParam.mutrate[1], myParam.mutrate0, myParam.x0, myParam.b0, myParam.fixed_mutrate, myParam.mutthreshfac, myParam.mutthresh, myParam.bindindex); //19
    //fprintf(fp,"%s\n", myParam.input);
    fprintf(fp,"%s\n", myParam.targetname);
    
    for (i=0; i<7; i++){
        fprintf(fp,"%d ", mutatororigin[i]);
    }
    fprintf(fp,"\n");
    
    for (who=0; who<MAXORGANISMS; who++){
        fprintf(fp,"%d %d\n", sizeRank[who], myOrgstatus[who]);
    }
    
    // variables from cell.h: myOrg
    for (who=0; who<MAXORGANISMS; who++){
        fprintf(fp,"%d %d %d %d %d %d %d\n", myOrg[who].genecount, myOrg[who].mutlevel, myOrg[who].mutorigin, myOrg[who].ppicount, myOrg[who].dob, myOrg[who].numkids, myOrg[who].generation);
        fprintf(fp,"%lf %lf %lf %lf\n", myOrg[who].minpnat, myOrg[who].meanpnat, myOrg[who].G0, myOrg[who].birthrate);
        fprintf(fp,"%f %f\n",myOrg[who].mutcount, myOrg[who].mutrate);
        
        for (i=0; i<(MAXGENES*NUCSEQLEN); i++){
            fprintf(fp,"%d ",myOrg[who].genome[i]); //MAXGENES*NUCSEQLEN;
        }
        fprintf(fp,"\n");
        
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%d ",myOrg[who].structid[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES*(MAXGENES-1)/2); i++){
            fprintf(fp,"%d ",myOrg[who].SeqID[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].ppi_pair[i][0]);
        }
        fprintf(fp,"\n");
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].ppi_pair[i][1]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<MAXPPIS; i++){
            fprintf(fp,"%d ",myOrg[who].bmode[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXSTATES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].F[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXSTATES+1); i++){
            for (j=0; j<(MAXSTATES+1); j++){
                fprintf(fp,"%lf ",myOrg[who].K[i][j]);
            }
            fprintf(fp,"\n");
        }
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].C[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].Kc[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%lf ",myOrg[who].Nch[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXGENES); i++){
            fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f\n",
                    myOrg[who].pnat[i],
                    myOrg[who].hydro[i],
                    myOrg[who].netcharge[i],
                    myOrg[who].poscharge[i],
                    myOrg[who].negcharge[i],
                    myOrg[who].fraccharge[i],
                    myOrg[who].frachydro[i],
                    myOrg[who].frachydro_avil[i],
                    myOrg[who].synmut[i],
                    myOrg[who].nonsynmut[i],
                    myOrg[who].tot_mut[i]);
        }

        for (i=0; i<(MAXPPIS); i++){
            fprintf(fp,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                    myOrg[who].hydro_s_hub[i],
                    myOrg[who].netcharge_s_hub[i],
                    myOrg[who].poscharge_s_hub[i],
                    myOrg[who].negcharge_s_hub[i],
                    myOrg[who].fraccharge_s_hub[i],
                    myOrg[who].frachydro_s_hub[i],
                    myOrg[who].frachydro_avil_s_hub[i],
                    
                    myOrg[who].hydro_s_partner[i],
                    myOrg[who].netcharge_s_partner[i],
                    myOrg[who].poscharge_s_partner[i],
                    myOrg[who].negcharge_s_partner[i],
                    myOrg[who].fraccharge_s_partner[i],
                    myOrg[who].frachydro_s_partner[i],
                    myOrg[who].frachydro_avil_s_partner[i]);
        }

        
        
        for (i=0; i<(MAXGENES+1); i++){
            fprintf(fp,"%f %f %f %f\n",
                    myOrg[who].nsi[i],
                    myOrg[who].si[i],
                    myOrg[who].nsi2[i],
                    myOrg[who].si2[i]);
        }
        
        
        for (i=0; i<(MAXPPIS+1); i++){
            fprintf(fp,"%f ",myOrg[who].pint[i]);
        }
        fprintf(fp,"\n");
        
        for (i=0; i<(MAXPPIS+1); i++){
            fprintf(fp,"%f ",myOrg[who].Gij[i]);
        }
        fprintf(fp,"\n");
        fprintf(fp,"\n");
    }
    
    
    fclose(fp);
}

/***********************************************************************************************************
 fscanf(fp,"%s %s %s %s %s %s\n",t_string, step_count_string, stretching_string,
 output_count_string,initial_L_string,strain_stoped_once_string );
 t = atof(t_string);
 output_count=atoi(output_count_string);
 step_count=atoi(step_count_string);
 stretching=atoi(stretching_string);
 initial_L = atof(initial_L_string);
 strain_stoped_once = atoi(strain_stoped_once_string);
 **********************************************************************************************************/
void ReadConfig(){
    int i, j, who;
    FILE *fp;
    char  divisioncycle_s[100], lastdecimtime_s[100], orgcount_s[100], mutatorcount_s[100], generation_s[100],
    speciesSizeTotal_s[100], domi_species_s[100], nOrgDB_s[100], TIME_s[100];
    char  seed_s[100], startcode_s[100], mp_orgcount_s[100], maxdivcycle_s[100], decimthresh_s[100], decimto_s[100], initpop_s[100], dumpcycle_s[100], printoutcycle_s[100], plotoutcycle_s[100], screenoutcycle_s[100], seqlogcycle_s[100], timeLow_s[100], timeHigh_s[100];
    char Tenv_s[100], tol_s[100], TLow_s[100], THigh_s[100], birthrate_s[100], deathrate_s[100], expressrate_s[100], alpha_s[100], pnatcutoff_s[100], speciessizefactor_s[100], mutrate0_s[100], mutrate01_s[100], mutrate00_s[100], x0_s[100], b0_s[100], fixed_mutrate_s[100], mutthreshfac_s[100], mutthresh_s[100], bindindex_s[100];
    //char targetname_s[500];
    char mutatororigin_i_s[100];
    char sizeRank_who_s[100], myOrgstatus_who_s[100];
    char mo_genecount_s[100], mo_mutlevel[100], mo_mutorigin[100], mo_ppicount[100], mo_dob[100], mo_numkids[100], mo_generation[100], mo_minpnat[100], mo_meanpnat[100], mo_G0[100], mo_birthrate[100], mo_mutcount[100], mo_mutrate[100], mo_genome_i[100];
    char mo_structid_i_s[100], mo_SeqID_i_s[100], mo_ppi_pair0_i_s[100], mo_ppi_pair1_i_s[100], mo_bmode_i_s[100], mo_F_i_s[100], mo_K_ij_s[100], mo_C_i_s[100], mo_Kc_i_s[100], mo_Nch_i_s[100];
    char mo_pnat_i_s[100], mo_hydro_i_s[100], mo_netcharge_i_s[100], mo_poscharge_i_s[100], mo_negcharge_i_s[100], mo_fraccharge_i_s[100], mo_frachydro_i_s[100], mo_frachydro_avil_i_s[100], mo_synmut_i_s[100], mo_nonsynmut_i_s[100], mo_tot_mut_i_s[100], mo_nsi_i_s[100], mo_si_i_s[100], mo_nsi2_i_s[100], mo_si2_i_s[100], mo_pint_i_s[100], mo_Gij_ij_s[100];
    char mo_hydro_i_s2[100], mo_netcharge_i_s2[100], mo_poscharge_i_s2[100], mo_negcharge_i_s2[100], mo_fraccharge_i_s2[100], mo_frachydro_i_s2[100], mo_frachydro_avil_i_s2[100]; 
    
    
    
    
    fp=fopen(FILE_CONFIG, "rb+");
    if(fp == NULL) { //if file does not exist
        RUN_BASED_ON_CONFIG = 0;
    }else {
        RUN_BASED_ON_CONFIG = 1;
    }
    printf("\nRUN_BASED_ON_CONFIG = %d\n",RUN_BASED_ON_CONFIG);
    
    if (RUN_BASED_ON_CONFIG){
        printf("Reading Config file ... ");
        
        // variables from evo-cell
        fscanf(fp,"%s %s %s %s %s %s %s %s\n",divisioncycle_s, lastdecimtime_s, orgcount_s, mutatorcount_s, generation_s,
               speciesSizeTotal_s, domi_species_s, nOrgDB_s);
        divisioncycle = atoi(divisioncycle_s);
        lastdecimtime = atoi(lastdecimtime_s);
        orgcount = atoi(orgcount_s);
        mutatorcount = atoi(mutatorcount_s);
        generation = atoi(generation_s);
        speciesSizeTotal = atoi(speciesSizeTotal_s);
        domi_species = atoi(domi_species_s);
        nOrgDB = atoi(nOrgDB_s);
        
        fscanf(fp,"%s\n",TIME_s);
        TIME = atof(TIME_s);
        
        
        printf("divisioncycle = %d\n",divisioncycle);
        printf("TIME = %E\n",TIME);
        //exit(0);
        
        // variables from cell.h: myParam, myOrg
        
        
        fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",seed_s,  startcode_s,  mp_orgcount_s,  maxdivcycle_s,  decimthresh_s,  decimto_s,  initpop_s,  dumpcycle_s,  printoutcycle_s,  plotoutcycle_s,  screenoutcycle_s,  seqlogcycle_s,  timeLow_s,  timeHigh_s);
        
        myParam.seed = atoi(seed_s);
        myParam.startcode = atoi(startcode_s);
        myParam.orgcount = atoi(orgcount_s);
        myParam.maxdivcycle = atoi(maxdivcycle_s);
        myParam.decimthresh = atoi(decimthresh_s);
        myParam.decimto = atoi(decimto_s);
        myParam.initpop = atoi(initpop_s);
        myParam.dumpcycle = atoi(dumpcycle_s);
        myParam.printoutcycle = atoi(printoutcycle_s);
        myParam.plotoutcycle = atoi(plotoutcycle_s);
        myParam.screenoutcycle = atoi(screenoutcycle_s);
        myParam.seqlogcycle = atoi(seqlogcycle_s);
        myParam.timeLow = atoi(timeLow_s);
        myParam.timeHigh = atoi(timeHigh_s);
        
        printf("myParam.timeHigh = %d\n",myParam.timeHigh);
        //exit(0);
        
        fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n",Tenv_s,  tol_s,  TLow_s,  THigh_s,  birthrate_s,  deathrate_s,  expressrate_s,  alpha_s,  pnatcutoff_s,  speciessizefactor_s,  mutrate0_s,  mutrate01_s,  mutrate00_s,  x0_s,  b0_s,  fixed_mutrate_s,  mutthreshfac_s,  mutthresh_s,  bindindex_s);
        
        myParam.Tenv = atof(Tenv_s);
        myParam.tol = atof(tol_s);
        myParam.TLow = atof(TLow_s);
        myParam.THigh = atof(THigh_s);
        myParam.birthrate = atof(birthrate_s);
        myParam.deathrate = atof(deathrate_s);
        myParam.expressrate = atof(expressrate_s);
        myParam.alpha = atof(alpha_s);
        myParam.pnatcutoff = atof(pnatcutoff_s);
        myParam.speciessizefactor = atof(speciessizefactor_s);
        myParam.mutrate[0] = atof(mutrate00_s);
        myParam.mutrate[1] = atof(mutrate01_s);
        myParam.mutrate0 = atof(mutrate0_s);
        myParam.x0 = atof(x0_s);
        myParam.b0 = atof(b0_s);
        myParam.fixed_mutrate = atof(fixed_mutrate_s);
        myParam.mutthreshfac = atof(mutthreshfac_s);
        myParam.mutthresh = atof(mutthresh_s);
        myParam.bindindex  = atof(bindindex_s);
        
        
        //printf("myParam.b0  = %lf\n",myParam.b0 );
        //exit(0);
        
        fscanf(fp,"%s\n", myParam.targetname);
        //printf("myParam.targetname  = %s\n",myParam.targetname);
        //exit(0);
        
        
        for (i=0; i<7; i++){
            fscanf(fp,"%s ", mutatororigin_i_s);
            mutatororigin[i] = atoi(mutatororigin_i_s);
            //printf("mutatororigin[i]   = %d\n",mutatororigin[i] );
            
        }
        //exit(0);
        fscanf(fp,"\n");
        
        for (who=0; who<MAXORGANISMS; who++){
            fscanf(fp,"%s %s\n", sizeRank_who_s, myOrgstatus_who_s);
            sizeRank[who] = atoi(sizeRank_who_s);
            myOrgstatus[who] = atoi(myOrgstatus_who_s);
            //printf("sizeRank  = %d myOrgstatus=%d\n",sizeRank[who] , myOrgstatus[who] );
        }
        
        //exit(0);
        
        // variables from cell.h: myOrg
        for (who=0; who<MAXORGANISMS; who++){
            fscanf(fp,"%s %s %s %s %s %s %s\n", mo_genecount_s, mo_mutlevel, mo_mutorigin, mo_ppicount, mo_dob, mo_numkids, mo_generation);
            myOrg[who].genecount = atoi(mo_genecount_s);
            myOrg[who].mutlevel = atoi(mo_mutlevel);
            myOrg[who].mutorigin = atoi(mo_mutorigin);
            myOrg[who].ppicount = atoi(mo_ppicount);
            myOrg[who].dob = atoi(mo_dob);
            myOrg[who].numkids = atoi(mo_numkids);
            myOrg[who].generation = atoi(mo_generation);
            
            //printf("generation  = %d\n",generation);
            //exit(0);
            
            fscanf(fp,"%s %s %s %s\n", mo_minpnat, mo_meanpnat, mo_G0, mo_birthrate);
            myOrg[who].minpnat = atof(mo_minpnat);
            myOrg[who].meanpnat = atof(mo_meanpnat);
            myOrg[who].G0 = atof(mo_G0);
            myOrg[who].birthrate = atof(mo_birthrate);
            //printf("myOrg[who].meanpnat  = %lf\n",myOrg[who].meanpnat);
            //printf("myOrg[who].birthrate  = %lf\n",myOrg[who].birthrate);
            //exit(0);
            
            
            fscanf(fp,"%s %s\n",mo_mutcount, mo_mutrate);
            myOrg[who].mutcount = atof(mo_mutcount);
            myOrg[who].mutrate = atof(mo_mutrate);
            
            //printf("myOrg[who].mutrate  = %f\n",myOrg[who].mutrate);
            //exit(0);
            
            for (i=0; i<(MAXGENES*NUCSEQLEN); i++){
                fscanf(fp,"%s ",mo_genome_i);
                myOrg[who].genome[i] = atoi(mo_genome_i);
                //printf("%s", mo_genome_i); //empty.....
                //printf("%d ", myOrg[who].genome[i]);
            }
            //printf("\n");
            fscanf(fp,"\n");
            //exit(0);
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_structid_i_s);
                myOrg[who].structid[i] = atoi(mo_structid_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES*(MAXGENES-1)/2); i++){
                fscanf(fp,"%s ",mo_SeqID_i_s);
                myOrg[who].SeqID[i] = atoi(mo_SeqID_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_ppi_pair0_i_s);
                myOrg[who].ppi_pair[i][0] = atoi(mo_ppi_pair0_i_s);
            }
            fscanf(fp,"\n");
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_ppi_pair1_i_s);
                myOrg[who].ppi_pair[i][1] = atoi(mo_ppi_pair1_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<MAXPPIS; i++){
                fscanf(fp,"%s ",mo_bmode_i_s);
                myOrg[who].bmode[i] = atoi(mo_bmode_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXSTATES+1); i++){
                fscanf(fp,"%s ",mo_F_i_s);
                myOrg[who].F[i] = atof(mo_F_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXSTATES+1); i++){
                for (j=0; j<(MAXSTATES+1); j++){
                    fscanf(fp,"%s ",mo_K_ij_s);
                    myOrg[who].K[i][j] = atof(mo_K_ij_s);
                }
                fscanf(fp,"\n");
            }
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_C_i_s);
                myOrg[who].C[i] = atof(mo_C_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_Kc_i_s);
                myOrg[who].Kc[i] = atof(mo_Kc_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s ",mo_Nch_i_s);
                myOrg[who].Nch[i] = atof(mo_Nch_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXGENES); i++){
                fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s\n",mo_pnat_i_s, mo_hydro_i_s, mo_netcharge_i_s, mo_poscharge_i_s, mo_negcharge_i_s, mo_fraccharge_i_s, mo_frachydro_i_s, mo_frachydro_avil_i_s, mo_synmut_i_s, mo_nonsynmut_i_s, mo_tot_mut_i_s);
                
                myOrg[who].pnat[i] = atof(mo_pnat_i_s);
                myOrg[who].hydro[i]= atof(mo_hydro_i_s);
                myOrg[who].netcharge[i]= atof(mo_netcharge_i_s);
                myOrg[who].poscharge[i]= atof(mo_poscharge_i_s);
                myOrg[who].negcharge[i]= atof(mo_negcharge_i_s);
                myOrg[who].fraccharge[i]= atof(mo_fraccharge_i_s);
                myOrg[who].frachydro[i]= atof(mo_frachydro_i_s);
                myOrg[who].frachydro_avil[i]= atof(mo_frachydro_avil_i_s);
                myOrg[who].synmut[i]= atof(mo_synmut_i_s);
                myOrg[who].nonsynmut[i]= atof(mo_nonsynmut_i_s);
                myOrg[who].tot_mut[i]= atof(mo_tot_mut_i_s);
            }
            
            for (i=0; i<(MAXPPIS); i++){
                fscanf(fp,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", mo_hydro_i_s, mo_netcharge_i_s, mo_poscharge_i_s, mo_negcharge_i_s, mo_fraccharge_i_s, mo_frachydro_i_s, mo_frachydro_avil_i_s, mo_hydro_i_s2, mo_netcharge_i_s2, mo_poscharge_i_s2, mo_negcharge_i_s2, mo_fraccharge_i_s2, mo_frachydro_i_s2, mo_frachydro_avil_i_s2);
                
                myOrg[who].hydro_s_hub[i]= atof(mo_hydro_i_s);
                myOrg[who].netcharge_s_hub[i]= atof(mo_netcharge_i_s);
                myOrg[who].poscharge_s_hub[i]= atof(mo_poscharge_i_s);
                myOrg[who].negcharge_s_hub[i]= atof(mo_negcharge_i_s);
                myOrg[who].fraccharge_s_hub[i]= atof(mo_fraccharge_i_s);
                myOrg[who].frachydro_s_hub[i]= atof(mo_frachydro_i_s);
                myOrg[who].frachydro_avil_s_hub[i]= atof(mo_frachydro_avil_i_s);
                
                myOrg[who].hydro_s_partner[i]= atof(mo_hydro_i_s2);
                myOrg[who].netcharge_s_partner[i]= atof(mo_netcharge_i_s2);
                myOrg[who].poscharge_s_partner[i]= atof(mo_poscharge_i_s2);
                myOrg[who].negcharge_s_partner[i]= atof(mo_negcharge_i_s2);
                myOrg[who].fraccharge_s_partner[i]= atof(mo_fraccharge_i_s2);
                myOrg[who].frachydro_s_partner[i]= atof(mo_frachydro_i_s2);
                myOrg[who].frachydro_avil_s_partner[i]= atof(mo_frachydro_avil_i_s2);
            }

            
            
            for (i=0; i<(MAXGENES+1); i++){
                fscanf(fp,"%s %s %s %s\n",mo_nsi_i_s, mo_si_i_s, mo_nsi2_i_s, mo_si2_i_s);
                myOrg[who].nsi[i] = atof(mo_nsi_i_s);
                myOrg[who].si[i] = atof(mo_si_i_s);
                myOrg[who].nsi2[i] = atof(mo_nsi2_i_s);
                myOrg[who].si2[i] = atof(mo_si2_i_s);
            }
            
            
            for (i=0; i<(MAXPPIS+1); i++){
                fscanf(fp,"%s ",mo_pint_i_s);
                myOrg[who].pint[i] = atof(mo_pint_i_s);
            }
            fscanf(fp,"\n");
            
            for (i=0; i<(MAXPPIS+1); i++){
                fscanf(fp,"%s ",mo_Gij_ij_s);
                myOrg[who].Gij[i] = atof(mo_Gij_ij_s);
            }
            fscanf(fp,"\n");
            fscanf(fp,"\n");
        }
        
        
        printf("Done\n\n");
    }
    //exit(0);
}


/***********************************************************************************************************
 **********************************************************************************************************/
void PrintOutput(){
    int ii, jj;

    
    //int aaseq[AASEQLEN];
    
    //char fopbuf[100];
    
    printf("gen %d orgcount %d\n", divisioncycle, orgcount);
    
    //printf("test\n");
    
    fprintf(out3, "%6d %E %6.3lf %5d %5d %9.6lf %8.5lf ", divisioncycle, TIME, myParam.Tenv, orgcount, nOrgDB, mean.birthrate, (double) domi_species/orgcount);
    fprintf(out3, "%10.7lf %8.5lf %8.5lf", mean.Gij[4], (double) mutatorcount/orgcount, mean.mutrate);
    fprintf(out3, "\n");
    
    //printf("beep1\n");
    
    fprintf(out5, "%6d %E", divisioncycle, TIME);
    for(ii=0;ii<curr_MAXGENES;ii++) fprintf(out5," %8.5f", (float) mutatororigin[ii]/orgcount);
    fprintf(out5, "\n");
    
    //printf("beep2\n");
    
    
    //fprintf(out1,"%d %E %lf %d",divisioncycle, TIME, myParam.Tenv, orgcount);
    //fprintf(out1," %lf %lf", mean.mutcount, mean.birthrate);
    //fprintf(out1,"\n");
    fprintf(out2,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out2, " %E", mean.C[ii]);
    for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out2, " %E", mean.F[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out2, " %.10E", mean.pnat[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out2, " %E", mean.hydro[ii]);
    fprintf(out2,"\n");
    
    
    //printf("beep3\n");
    
    //chaps
    fprintf(out12,"%d %E",divisioncycle, TIME);
    fprintf(out12, " %E", mean.C[curr_MAXGENES]);
    fprintf(out12, " %E", mean.F[2*curr_MAXGENES]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out12, " %E", mean.Nch[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out12, " %E", mean.Kc[ii]);
    fprintf(out12,"\n");
    
    //printf("beep4\n");
    
    
    fprintf(out6,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out6, " %E", mean.pint[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out6, " %E", mean.Gij[ii]);
    for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out6, " %E", mean.K[0][ii+1]); //OP, think about
    fprintf(out6,"\n");
    
    //printf("beep5\n");
    
    fprintf(out7,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out7, " %E", mean.nonsynmut[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out7, " %E", mean.synmut[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out7, " %E", mean.tot_mut[ii]);
    fprintf(out7,"\n");
    
    //printf("beep6\n");
    
    
    //nsilog:
    fprintf(out10,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out10, " %E", mean.nsi[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out10, " %E", mean.si[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out10, " %E", mean.nsi2[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out10, " %E", mean.si2[ii]);
    fprintf(out10,"\n");
    
    //printf("beep7\n");
    
    //seq prop:
    fprintf(out11,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.hydro[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.frachydro[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.frachydro_avil[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.netcharge[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.poscharge[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.negcharge[ii]);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out11, " %E", mean.fraccharge[ii]);
    fprintf(out11,"\n");

    
    //seq surfaces prop:
    fprintf(out15,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.hydro_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.frachydro_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.frachydro_avil_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.netcharge_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.poscharge_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.negcharge_s_hub[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out15, " %E", mean.fraccharge_s_hub[ii]);
    fprintf(out15,"\n");
    
    fprintf(out16,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.hydro_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.frachydro_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.frachydro_avil_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.netcharge_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.poscharge_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.negcharge_s_partner[ii]);
    for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out16, " %E", mean.fraccharge_s_partner[ii]);
    fprintf(out16,"\n");
    
    // print out entropy
    fprintf(out8,"%d %E",divisioncycle, TIME);
    for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out8, " %.10lE", entropy_sum[ii]);
    fprintf(out8, " %.10lE", entropy_sumtot);
    fprintf(out8,"\n");
    
    //printf("beep8\n");
    
    
    // seqlog printout
    //printf("myParam->seqlogcycle=%d\n",myParam.seqlogcycle);
    
    if(divisioncycle%myParam.seqlogcycle == 0) {
        //printf("beep8.-1\n");
        //sprintf(fopbuf,"seqlog-%s.dat.gz", myParam.targetname);
        //printf("beep8.0\n");
        //zout4=gzopen(fopbuf,"a");
        RankSpeciesSizeDB(sizeRank, nOrgDB);
        //gzprintf(zout4, "%6d %4d %4d\n", divisioncycle, nOrgDB, orgcount);
        fprintf(out4, "%6d %4d %4d\n", divisioncycle, nOrgDB, orgcount);
        speciesSizeTotal=0;
        //printf("beep8.1\n");
        
        for(ii=0;ii<nOrgDB;ii++) {
            //gzprintf(zout4, "%4d %4d %1d ", sizeRank[ii], myOrgDB[sizeRank[ii]].count, myOrgDB[sizeRank[ii]].genecount);
            fprintf(out4, "%4d %4d %1d ", sizeRank[ii], myOrgDB[sizeRank[ii]].count, myOrgDB[sizeRank[ii]].genecount);
            PrintCharNucCodeSequence(seqbuf, myOrgDB[sizeRank[ii]].genome, myOrgDB[sizeRank[ii]].genecount*NUCSEQLEN);
            //gzprintf(zout4, "%s ", seqbuf);
            fprintf(out4, "%s ", seqbuf);
            for(jj=0;jj<curr_MAXGENES;jj++){
                //gzprintf(zout4, "%5d ", myOrgDB[sizeRank[ii]].structid[jj]);
                fprintf(out4, "%5d ", myOrgDB[sizeRank[ii]].structid[jj]);

            }
            speciesSizeTotal += myOrgDB[sizeRank[ii]].count;
            //gzprintf(zout4, "\n");
            fprintf(out4, "\n");
        }
        //gzprintf(zout4, ";\n\n");
        //gzclose(zout4);
        fprintf(out4, ";\n\n");
        //fclose(zout4);
        //printf("beep8.2\n");
        
        
        fprintf(out1,"%08d %d",divisioncycle, orgcount);
        for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %E", mean.C[ii]);
        for (ii=0; ii<curr_MAXSTATES; ii++) fprintf(out1, " %E", mean.F[ii]);
        for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %.10E", mean.pnat[ii]);
        for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %E", mean.pint[ii]);
        for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %E", mean.Gij[ii]);
        
        
        for (ii=0; ii<curr_MAXGENES; ii++) fprintf(out1, " %06d", myOrgDB[sizeRank[0]].structid[ii]);
        for (ii=0; ii<curr_MAXPPIS; ii++) fprintf(out1, " %03d", myOrgDB[sizeRank[0]].bmode[ii]);
        
        
        int aaseq_surface_hub[AASURFACELEN], aaseq_surface_partner[AASURFACELEN], s_i, hub_i, par_i;
        
        for (ii=0;ii<curr_MAXPPIS; ii++){
            hub_i=0; par_i=ii+1;
            GetSurfaceAAPositions(myOrgDB[sizeRank[0]].structid[hub_i], myOrgDB[sizeRank[0]].structid[par_i], myOrgDB[sizeRank[0]].bmode[ii],
                                  aaseq_surface_hub, aaseq_surface_partner);
            
            for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out1, " %02d", aaseq_surface_hub[s_i]);
            for (s_i=0;s_i<AASURFACELEN; s_i++) fprintf(out1, " %02d", aaseq_surface_partner[s_i]);
            
        }
        
        
        PrintCharNucCodeSequence(seqbuf, myOrgDB[sizeRank[0]].genome, myOrgDB[sizeRank[0]].genecount*NUCSEQLEN);
        fprintf(out1, " %s ", seqbuf);
        fprintf(out1,"\n");
        
        
        
        
        
//        for (ii=0;ii<curr_MAXPPIS; ii++){
//            i=myOrgDB[sizeRank[0]].ppi_pair[ii][0], j=myOrgDB[sizeRank[0]].ppi_pair[ii][1];
//            //printf("i=%d j=%d\n", i,j);
//            CharNucSeqToAASeq(myOrgDB[sizeRank[0]].genome+i*NUCSEQLEN,NUCSEQLEN,aaseq_hub);
//            CharNucSeqToAASeq(myOrgDB[sizeRank[0]].genome+j*NUCSEQLEN,NUCSEQLEN,aaseq_partner);
//            
//            GetSurfaceAA(aaseq_hub, aaseq_partner, aaseq_surface_hub, aaseq_surface_partner, myOrg[who].structid[i], myOrg[who].structid[j], myOrg[who].bmode[ii]);
//            
//        }
        
        //CharNucSeqToAASeq(seqbuf,NUCSEQLEN,aaseq);
        
        
        //for (ii=0; ii<AASEQLEN; ii++) fprintf(out1, "%d", aaseq[ii]);
        
        
        

    }
    
    //printf("beep9\n");
    
}

/***********************************************************************************************************
 **********************************************************************************************************/
void Openfiles(){
    char fopbuf[100];
    char filetype_buf[100];
    
    if (RUN_BASED_ON_CONFIG == 0){
        sprintf(filetype_buf,"w");
    }
    else {
        sprintf(filetype_buf,"a");
    }
    printf("in Openfiles, filetype_buf=%s\n",filetype_buf);
    
    sprintf(fopbuf,"init_vals-%s.dat", myParam.targetname);
    printf("in Openfiles, fopbuf=%s\n",fopbuf);
    
    
    finitvals=fopen(fopbuf,filetype_buf);
    fprintf(finitvals, "%s ", myParam.targetname);
    fprintf(finitvals, "%d ", myParam.seed);
    fprintf(finitvals, "%d ", myParam.maxdivcycle);
    fprintf(finitvals, "%d ", curr_MAXGENES);
    fprintf(finitvals, "%d ", hub_ID);
    fprintf(finitvals, "%d ", allow_fold_change);
    fprintf(finitvals, "%d ", allow_gene_exp);
    fprintf(finitvals, "%d ", allow_unfolded_states);
    fclose(finitvals);
    
    //printf("beep 1\n"); exit(0);
    
    sprintf(fopbuf,"initial-%s.dat", myParam.targetname);
    out1=fopen(fopbuf,filetype_buf);
    PrintInitialCondition(out1,&myParam);
    fclose(out1);
    
    sprintf(fopbuf,"basiclog-%s.dat", myParam.targetname);
    out1=fopen(fopbuf,filetype_buf);
    sprintf(fopbuf,"statlog-%s.dat", myParam.targetname);
    out2=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"statlog_chaps-%s.dat", myParam.targetname);
    out12=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"plotlog-%s.dat", myParam.targetname);
    out3=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"mutlog-%s.dat", myParam.targetname);
    out5=fopen(fopbuf,filetype_buf);
    sprintf(fopbuf,"problog-%s.dat", myParam.targetname);
    out6=fopen(fopbuf,filetype_buf);
    sprintf(fopbuf,"ratelog-%s.dat", myParam.targetname);
    out7=fopen(fopbuf,filetype_buf);
    sprintf(fopbuf,"entlog-%s.dat", myParam.targetname);
    out8=fopen(fopbuf,filetype_buf);
    sprintf(fopbuf,"nsilog-%s.dat", myParam.targetname);
    out10=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"seqproplog-%s.dat", myParam.targetname);
    out11=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"seqproplog_surfaces_hub-%s.dat", myParam.targetname);
    out15=fopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"seqproplog_surfaces_par-%s.dat", myParam.targetname);
    out16=fopen(fopbuf,filetype_buf);
    
    //sprintf(fopbuf,"seqlog-%s.dat.gz", myParam.targetname);
    //zout4=gzopen(fopbuf,filetype_buf);
    
    sprintf(fopbuf,"seqlog-%s.dat", myParam.targetname);
    out4=fopen(fopbuf,filetype_buf);
    
}

/***********************************************************************************************************
 **********************************************************************************************************/
void Flushfiles(){
    fflush(out1); fflush(out2); fflush(out3); fflush(out5); fflush(out6); fflush(out7); fflush(out8);
    fflush(out10); fflush(out11); fflush(out12); fflush(error_op); fflush(stdout); fflush(it_solver);
    fflush(out15); fflush(out16);
    fflush(out4);
}

/***********************************************************************************************************
 **********************************************************************************************************/
void Closefiles(){
    fclose(out1); fclose(out2); fclose(out3); fclose(out5); fclose(out6); fclose(out10); fclose(out11); fclose(out12);
    fclose(out7);fclose(error_op);
    fclose(it_solver); fclose(out8);
    fclose(out15); fclose(out16);
    fclose(out4);
}


/***********************************************************************************************************
 // calculates entropy at each position of an amino acid
 // based on distribution among organisms
 **********************************************************************************************************/
int seq_entropy(char **seq, int N, double *entropy){
    int i, j;
    double f;
    
    for(j=0;j<AASEQLEN;j++)
        for(i=0;i<20;i++)
            aamat[j][i]=0;
    
    for(i=0;i<N;i++){   // N: number of copies of amino acid
        CharNucSeqToAASeq(seq[i], NUCSEQLEN, aaseq);
        for(j=0;j<AASEQLEN;j++){
            aamat[j][aaseq[j]]++;
            // counts number of each amino acid at each position
        }
    }
    
    for(j=0;j<AASEQLEN;j++){
        entropy[j]=0.0;
        //fprintf(stdout, "%4d :", j);
        for(i=0;i<20;i++){
            if(aamat[j][i] > 0) {
                f = (double) aamat[j][i]/N;
                entropy[j]-= f*logf(f);
            } else {
                f =0.0;
            }
            //fprintf(stdout, " %8.3f", f);
        }
        //fprintf(stdout,"\n");
    }
    return 0;
}


/***********************************************************************************************************
 **********************************************************************************************************/

// Updates database of all organism types (i.e. distinct genomes)
// in which each distinct genome gets its own space
int ResetOrgDB(parameter *myParam, int divisioncycle){
    int i, j, genediff, who;
    //int face1, face2, rotate;
    //int surfacetmp[9],surface1[9], surface2[9];
    //double e, z, emin=1e10;
    //double T = myParam->Tenv;
    
    nOrgDB=0;
    domi_species = 0;
    
    for(who=0;who<MAXORGANISMS;who++) {
        if (myOrgstatus[who]!=S_ALIVE) continue;
        
        for(i=0;i<nOrgDB;i++) {
            if(myOrgDB[i].genecount != myOrg[who].genecount) continue;
            genediff=0;
            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) {
                // comparing if genome of who matches any genome in the database
                genediff+=(myOrgDB[i].genome[j]-myOrg[who].genome[j])*(myOrg[i].genome[j]-myOrg[who].genome[j]);
            }
            if(!genediff) break;
        }
        if(i == nOrgDB) {
            // if the genome does not match any previous entry, a new space is
            // made for it in OrgDB
            myOrgDB[i].genecount=myOrg[who].genecount;
            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) myOrgDB[i].genome[j]=myOrg[who].genome[j];
            for(j=0;j<myOrg[who].genecount;j++) {
                myOrgDB[i].structid[j]=myOrg[who].structid[j];
                myOrgDB[i].pnat[j]=myOrg[who].pnat[j];
                myOrgDB[i].hydro[j]=myOrg[who].hydro[j];
            }
            
            for(j=0;j<myOrg[who].genecount-1  ;j++) {
                myOrgDB[i].bmode[j]=myOrg[who].bmode[j];
            }

            
            myOrgDB[i].reporg = who;
            myOrgDB[i].count = 1;
            //      myOrgDB[i].meanmutrate = myOrg[who].mutrate; //mod myOrgDBMut
            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
            nOrgDB++;
        } else {
            myOrgDB[i].count++;
            //      myOrgDB[i].meanmutrate += myOrg[who].mutrate; //mod myOrgDBMut
            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
        }
    } // for(who=0;
    //  for(i=0;i<nOrgDB;i++) myOrgDB[i].meanmutrate /= myOrgDB[i].count; //mod myOrgDBMut
    
    //fprintf(stderr, "ResetOrgDB : %d types of organisms...\n", nOrgDB);
    
    return 0;
}


void SetMeanOrg(organism *mean){
    int ii, jj;
    //mean->genecount = mean->ppicount = 0;
    mean->dob = mean->generation = mean->numkids = 0;
    mean->G0 = mean->birthrate = mean->meanpnat = mean->minpnat = mean->mutrate = mean->mutcount = 0.0e0;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->si2[ii] = mean->nsi2[ii] = mean->si[ii] = mean->nsi[ii] = mean->C[ii] = mean->pnat[ii] = mean-> hydro[ii] = mean->tot_mut[ii] = mean->nonsynmut[ii] = mean->synmut[ii] = 0.0e0;
        
        mean->netcharge[ii] = 0.0e0;
        mean->poscharge[ii] = 0.0e0;
        mean->negcharge[ii] = 0.0e0;
        mean->fraccharge[ii] = 0.0e0;
        mean->frachydro[ii] = 0.0e0;
        mean->frachydro_avil[ii] = 0.0e0;
        
        mean->Nch[ii] = 0.0e0;
        mean->Kc[ii] = 0.0e0;
        
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] =  0.0e0;
        for(jj=0;jj<curr_MAXSTATES;jj++) mean->K[ii][jj] = 0.0e0;
    }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] =  0.0e0;
        mean->C[curr_MAXGENES] =  0.0e0;
    }
    
    
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] = mean->pint[ii] = 0.0e0;
        mean->netcharge_s_partner[ii] = 0.0e0;
        mean->poscharge_s_partner[ii] = 0.0e0;
        mean->negcharge_s_partner[ii] = 0.0e0;
        mean->fraccharge_s_partner[ii] = 0.0e0;
        mean->frachydro_s_partner[ii] = 0.0e0;
        mean->frachydro_avil_s_partner[ii] = 0.0e0;
        
        mean->netcharge_s_hub[ii] = 0.0e0;
        mean->poscharge_s_hub[ii] = 0.0e0;
        mean->negcharge_s_hub[ii] = 0.0e0;
        mean->fraccharge_s_hub[ii] = 0.0e0;
        mean->frachydro_s_hub[ii] = 0.0e0;
        mean->frachydro_avil_s_hub[ii] = 0.0e0;
    }
    return;
}

void AddMeanOrg(organism *mean, int who){
    int ii, jj;
    mean->birthrate += myOrg[who].birthrate;
    mean->mutrate += myOrg[who].mutrate;
    mean->mutcount += myOrg[who].mutcount;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->C[ii] += myOrg[who].C[ii];
        mean->nsi[ii] += myOrg[who].nsi[ii];
        mean->si[ii] += myOrg[who].si[ii];
        mean->nsi2[ii] += myOrg[who].nsi2[ii];
        mean->si2[ii] += myOrg[who].si2[ii];
        //mean->F[ii] += myOrg[who].F[ii];
        mean->pnat[ii] += myOrg[who].pnat[ii];
        
        mean->hydro[ii] += myOrg[who].hydro[ii];
        mean->netcharge[ii] += myOrg[who].netcharge[ii];
        mean->poscharge[ii] += myOrg[who].poscharge[ii];
        mean->negcharge[ii] += myOrg[who].negcharge[ii];
        mean->fraccharge[ii] += myOrg[who].fraccharge[ii];
        mean->frachydro[ii] += myOrg[who].frachydro[ii];
        mean->frachydro_avil[ii] += myOrg[who].frachydro_avil[ii];
        
        mean->nonsynmut[ii] += myOrg[who].nonsynmut[ii];
        mean->tot_mut[ii] += myOrg[who].tot_mut[ii];
        mean->synmut[ii] += myOrg[who].synmut[ii];
        
        mean->Nch[ii] += myOrg[who].Nch[ii];
        mean->Kc[ii] += myOrg[who].Kc[ii];
        
        //for(jj=0;jj<curr_MAXGENES;jj++) mean->K[ii][jj] += myOrg[who].K[ii][jj];
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] += myOrg[who].F[ii];
        for(jj=0;jj<curr_MAXSTATES;jj++) mean->K[ii][jj] += myOrg[who].K[ii][jj];
    }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] += myOrg[who].F[2*curr_MAXGENES];
    }
    
    if (allow_chaps==1){
        mean->C[curr_MAXGENES] += myOrg[who].C[curr_MAXGENES];
    }
    
    for(ii=0;ii<curr_MAXGENES*(curr_MAXGENES-1)/2;ii++) mean->SeqID[ii] +=myOrg[who].SeqID[ii];
    mean->G0 += myOrg[who].G0;
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] += myOrg[who].Gij[ii];
        mean->pint[ii] += myOrg[who].pint[ii];
        
        mean->hydro_s_partner[ii] += myOrg[who].hydro_s_partner[ii];
        mean->netcharge_s_partner[ii] += myOrg[who].netcharge_s_partner[ii];
        mean->poscharge_s_partner[ii] += myOrg[who].poscharge_s_partner[ii];
        mean->negcharge_s_partner[ii] += myOrg[who].negcharge_s_partner[ii];
        mean->fraccharge_s_partner[ii] += myOrg[who].fraccharge_s_partner[ii];
        mean->frachydro_s_partner[ii] += myOrg[who].frachydro_s_partner[ii];
        mean->frachydro_avil_s_partner[ii] += myOrg[who].frachydro_avil_s_partner[ii];

        
        mean->hydro_s_hub[ii] += myOrg[who].hydro_s_hub[ii];
        mean->netcharge_s_hub[ii] += myOrg[who].netcharge_s_hub[ii];
        mean->poscharge_s_hub[ii] += myOrg[who].poscharge_s_hub[ii];
        mean->negcharge_s_hub[ii] += myOrg[who].negcharge_s_hub[ii];
        mean->fraccharge_s_hub[ii] += myOrg[who].fraccharge_s_hub[ii];
        mean->frachydro_s_hub[ii] += myOrg[who].frachydro_s_hub[ii];
        mean->frachydro_avil_s_hub[ii] += myOrg[who].frachydro_avil_s_hub[ii];

    }
    return;
}

void GetMeanOrg(organism *mean, int orgcount){
    int ii, jj;
    mean->birthrate /= (double) orgcount;
    mean->mutrate /= (double)orgcount;
    mean->mutcount /= (double)orgcount;
    for(ii=0;ii<curr_MAXGENES;ii++) {
        mean->C[ii] /= (double)orgcount;
        mean->nsi[ii] /= (double)orgcount;
        mean->si[ii] /= (double)orgcount;
        mean->nsi2[ii] /= (double)orgcount;
        mean->si2[ii] /= (double)orgcount;
        
        mean->Kc[ii] /= (double)orgcount;
        mean->Nch[ii] /= (double)orgcount;
        
        //mean->F[ii] /= (double)orgcount;
        mean->pnat[ii] /= (double)orgcount;
        mean->hydro[ii] /= (double)orgcount;
        mean->netcharge[ii] /= (double)orgcount;
        mean->poscharge[ii] /= (double)orgcount;
        mean->negcharge[ii] /= (double)orgcount;
        mean->fraccharge[ii] /= (double)orgcount;
        mean->frachydro[ii] /= (double)orgcount;
        mean->frachydro_avil[ii] /= (double)orgcount;
        
        
        mean->nonsynmut[ii] /= (double)orgcount;
        mean->tot_mut[ii] /= (double)orgcount;
        mean->synmut[ii] /= (double)orgcount;
        
        //for(jj=0;jj<curr_MAXGENES;jj++) mean->K[ii][jj] /= (double)orgcount;
    }
    
    for(ii=0;ii<curr_MAXSTATES;ii++) {
        mean->F[ii] /= (double)orgcount;
        for(jj=0;jj<curr_MAXSTATES;jj++) mean->K[ii][jj] /= (double)orgcount;
    }
    
    if (allow_chaps==1){
        mean->F[2*curr_MAXGENES] /= (double)orgcount;
        mean->C[curr_MAXGENES] /= (double)orgcount;
        
    }
    
    
    for(ii=0;ii<curr_MAXGENES*(curr_MAXGENES-1)/2;ii++) mean->SeqID[ii] /= (double)orgcount;
    mean->G0 /= (double)orgcount;
    for(ii=0;ii<curr_MAXPPIS;ii++) {
        mean->Gij[ii] /= (double)orgcount;
        mean->pint[ii] /= (double)orgcount;
        
        mean->hydro_s_partner[ii] /= (double)orgcount;
        mean->netcharge_s_partner[ii] /= (double)orgcount;
        mean->poscharge_s_partner[ii] /= (double)orgcount;
        mean->negcharge_s_partner[ii] /= (double)orgcount;
        mean->fraccharge_s_partner[ii] /= (double)orgcount;
        mean->frachydro_s_partner[ii] /= (double)orgcount;
        mean->frachydro_avil_s_partner[ii] /= (double)orgcount;

        mean->hydro_s_hub[ii] /= (double)orgcount;
        mean->netcharge_s_hub[ii] /= (double)orgcount;
        mean->poscharge_s_hub[ii] /= (double)orgcount;
        mean->negcharge_s_hub[ii] /= (double)orgcount;
        mean->fraccharge_s_hub[ii] /= (double)orgcount;
        mean->frachydro_s_hub[ii] /= (double)orgcount;
        mean->frachydro_avil_s_hub[ii] /= (double)orgcount;

    }
    return;
}



/************************************************************************
 void start_clock(): time evaluation funtion.
 This example assumes that the result of each subtraction
 is within the range of values that can be represented in
 an integer type.
 ************************************************************************/

void start_clock()
{
    st_time = times(&st_cpu);
}


/************************************************************************
 void end_clock(): time evaluation funtion.
 This example assumes that the result of each subtraction
 is within the range of values that can be represented in
 an integer type.
 
 printf("Real Time: %jd, User Time %jd, System Time %jd\n")
 CPUtime = user time + system time.
 fprintf(fpTime,"%d\t%jd\t%jd\t%jd\n", stepCount,
 (intmax_t)(en_time - st_time),
 (intmax_t)(en_cpu.tms_utime - st_cpu.tms_utime),
 (intmax_t)(en_cpu.tms_stime - st_cpu.tms_stime));
 fclose(fpTime);
 ************************************************************************/

int end_clock()
{
    en_time = times(&en_cpu);
    
    return ((intmax_t)(en_cpu.tms_utime - st_cpu.tms_utime)+(intmax_t)(en_cpu.tms_stime - st_cpu.tms_stime));
}








/***********************************************************************************************************
 UNUSED!!
 **********************************************************************************************************/

//int seq_entropy(char **seq, int N, float *entropy)
//// calculates entropy at each position of an amino acid
//// based on distribution among organisms
//{
//    int i, j;
//    float f;
//
//    for(j=0;j<AASEQLEN;j++)
//        for(i=0;i<20;i++)
//            aamat[j][i]=0;
//
//    for(i=0;i<N;i++){   // N: number of copies of amino acid
//        //PrintCharNucCodeSequence(pcharbuf, seq[i], NUCSEQLEN);
//        //fprintf(stdout,"%6d %s\n", i, pcharbuf);
//        CharNucSeqToAASeq(seq[i], NUCSEQLEN, aaseq);
//        //PrintAACodeSequence(pcharbuf2, aaseq, AASEQLEN);
//        //fprintf(stdout,"%6d %s\n", i, pcharbuf2);
//        for(j=0;j<AASEQLEN;j++){
//            aamat[j][aaseq[j]]++;
//            // counts number of each amino acid at each position
//        }
//    }
//
//    for(j=0;j<AASEQLEN;j++){
//        entropy[j]=0.0;
//        //fprintf(stdout, "%4d :", j);
//        for(i=0;i<20;i++){
//            if(aamat[j][i] > 0) {
//                f = (float) aamat[j][i]/N;
//                entropy[j]-= f*logf(f);
//            } else {
//                f =0.0;
//            }
//            //fprintf(stdout, " %8.3f", f);
//        }
//        //fprintf(stdout,"\n");
//    }
//    return 0;
//}
//
//
//int DumpDB(parameter *myParam, int divisioncycle) {
//
//    int i, j, k, ii;
//    char filename[200], buf[200];
//    FILE *fp;
//
//    sprintf(filename, "BindingDB-%s-t%d.dat", myParam->targetname, divisioncycle);
//    fp=fopen(filename, "w");
//    for(i=0;i<nOrgDB;i++) {
//        fprintf(fp, "OrgDB %5d RepOrg %5d Count %5d pint12 = %8.3f \n", i, myOrgDB[i].reporg, myOrgDB[i].count, myOrgDB[i].pint12);
//        CharNucSeqToAASeq(myOrgDB[i].genome+NUCSEQLEN, NUCSEQLEN, aaseq);
//        CharNucSeqToAASeq(myOrgDB[i].genome+NUCSEQLEN*2, NUCSEQLEN, aaseq2);
//        fprintf(fp, "%6d", myOrgDB[i].structid[1]);
//        for(k=0;k<9;k++) fprintf(fp, " %3d (%1c)", myOrgDB[i].binding[0][k], aacode[aaseq[myOrgDB[i].binding[0][k]]]);
//        fprintf(fp, "\n");
//        fprintf(fp, "%6d", myOrgDB[i].structid[2]);
//        for(k=0;k<9;k++) fprintf(fp, " %3d (%1c)", myOrgDB[i].binding[1][k], aacode[aaseq2[myOrgDB[i].binding[1][k]]]);
//        fprintf(fp, "\n");
//
//        /*
//         CharNucSeqToAASeq(myOrgDB[i].genome+NUCSEQLEN*3, NUCSEQLEN, aaseq);
//         CharNucSeqToAASeq(myOrgDB[i].genome+NUCSEQLEN*4, NUCSEQLEN, aaseq2);
//         fprintf(fp, "%6d", myOrgDB[i].structid[3]);
//         for(k=0;k<9;k++) fprintf(fp, " %3d (%1c)", myOrgDB[i].binding[2][k], aacode[aaseq[myOrgDB[i].binding[2][k]]]);
//         fprintf(fp, "\n");
//         fprintf(fp, "%6d", myOrgDB[i].structid[4]);
//         for(k=0;k<9;k++) fprintf(fp, " %3d (%1c)", myOrgDB[i].binding[3][k], aacode[aaseq2[myOrgDB[i].binding[3][k]]]);
//         fprintf(fp, "\n");
//         */
//    }
//    fclose(fp);
//
//
//    for(ii=0;ii<MAXGENES;ii++) {
//        k=0;
//        for(i=0;i<nOrgDB;i++) {
//            if(myOrgDB[i].count > 10) k++;
//        }
//        sprintf(filename,"SeqDB-%s-%d-t%d.dat", myParam->targetname, ii, divisioncycle);
//        fp=fopen(filename, "w");
//        fprintf(fp, "%d %d I\n\n", k, NUCSEQLEN);
//        for(i=0;i<nOrgDB;i++) {
//            if(myOrgDB[i].count > 10) fprintf(fp, "s%d.%d\n", i, myOrgDB[i].count);
//        }
//
//        for(k=0;k<NUCSEQLEN/60;k++) {
//            fprintf(fp,"%d\n",k*60+1);
//            for(i=0;i<nOrgDB;i++) {
//                if(myOrgDB[i].count > 10) {
//                    PrintCharNucCodeSequence(buf, myOrgDB[i].genome+ii*NUCSEQLEN, NUCSEQLEN);
//                    for(j=k*60;j<(k+1)*60;j++) {
//                        fprintf(fp, "%c",buf[j]);
//                        if(j%3==2) fprintf(fp, " ");
//                    }
//                    fprintf(fp,"\n");
//                }
//            }
//        }
//        if(NUCSEQLEN%60 != 0 ) {
//            fprintf(fp,"%d\n",k*60+1);
//            for(i=0;i<nOrgDB;i++) {
//                if(myOrgDB[i].count > 10) {
//                    PrintCharNucCodeSequence(buf, myOrgDB[i].genome+ii*NUCSEQLEN, NUCSEQLEN);
//                    for(j=k*60;j<NUCSEQLEN;j++) {
//                        fprintf(fp, "%c",buf[j]);
//                        if(j%3==2) fprintf(fp, " ");
//                    }
//                    fprintf(fp,"\n");
//                }
//            }
//        }
//        fclose(fp);
//    } // for(ii=0;
//    return 0;
//}
//
//int DumpProteome(char *filename){
//    int who, ii;
//    gzFile gzout1p;
//    char buf[200];
//    gzout1p = gzopen(filename, "w");
//    for(who=0;who<MAXORGANISMS;who++) {
//        if(myOrgstatus[who]!=S_ALIVE) continue;
//        gzprintf(gzout1p,"ORG %d\t%d\t%d\t%d\t%f\t%f", who, myOrg[who].dob, myOrg[who].numkids, myOrg[who].generation, myOrg[who].mutcount, myOrg[who].pint[0]);
//        for(ii=0;ii<MAXGENES*(MAXGENES-1)/2;ii++)  gzprintf(gzout1p,"\t%lf", myOrg[who].SeqID[ii]); // count 3
//        gzprintf(gzout1p,"\n");
//
//        for(ii=0;ii<myOrg[who].genecount;ii++) {
//            PrintCharNucCodeSequence(buf, myOrg[who].genome+ii*NUCSEQLEN, NUCSEQLEN);
//            gzprintf(gzout1p,"%d\t%d\t%d\t%f\t%f\t%f\t%f\t%s\n",who,ii,myOrg[who].structid[ii], myOrg[who].pnat[ii],myOrg[who].C[ii], myOrg[who].F[ii], myOrg[who].hydro[ii], buf);
//        }
//    }
//    gzclose(gzout1p);
//    //ResetGeneDB();
//
//    return 0;
//}
//
//// This is a trimmed-down version of ResetOrgDB used to recalculate dominant species and species number information
//int RecountOrgDB(int divisioncycle){
//    int i, j, genediff, who;
//
//    nOrgDB=0;
//    domi_species=0;
//
//    for(who=0;who<MAXORGANISMS;who++) {
//        if (myOrgstatus[who]!=S_ALIVE) continue;
//
//        for(i=0;i<nOrgDB;i++) {
//            if(myOrgDB[i].genecount != myOrg[who].genecount) continue;
//
//            genediff=0;
//            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) {
//                // comparing if genome of who matches any genome in the database
//                genediff+=(myOrgDB[i].genome[j]-myOrg[who].genome[j])*(myOrg[i].genome[j]-myOrg[who].genome[j]);
//            }
//            if(!genediff) break;
//        }
//        if(i == nOrgDB) {
//            // if the genome does not match any previous entry, a new space is
//            // made for it in OrgDB
//            myOrgDB[i].genecount = myOrg[who].genecount;
//            for(j=0;j<myOrg[who].genecount*NUCSEQLEN;j++) myOrgDB[i].genome[j]=myOrg[who].genome[j];
//
//            myOrgDB[i].reporg = who;
//            myOrgDB[i].count = 1;
//            //      myOrgDB[i].meanmutrate = myOrg[who].mutrate; //mod myOrgDBMut
//            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
//            nOrgDB++;
//        } else {
//            myOrgDB[i].count++;
//            //      myOrgDB[i].meanmutrate += myOrg[who].mutrate; //mod myOrgDBMut
//            if(myOrgDB[i].count > domi_species) domi_species = myOrgDB[i].count;
//        }
//    } // for(who=0;
//    //  for(i=0;i<nOrgDB;i++) myOrgDB[i].meanmutrate /= myOrgDB[i].count; //mod myOrgDBMut
//
//    return 0;
//}
//
//
