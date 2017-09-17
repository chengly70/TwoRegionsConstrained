#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mex.h"

/* Box-Muller; faster than using sin,cos */
 void z_randn(double *normal_rv, double *nrm_rv2)
 {
	 double u, v, s=2.0, z, z2;
	 while(s>=1 || s==0){
		 u=(double)(rand())/RAND_MAX*2-1;
		 v=(double)(rand())/RAND_MAX*2-1;
		 s=u*u+v*v;
	 }
	 z=sqrt(-2.0*log(s)/s)*u;
	 z2=z/u*v;
	 
	 *normal_rv=z;
	 *nrm_rv2=z2;
 }

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
/* 
 
[nuPC,nuOB,mn_PC,mn_OB,var_PC,var_OB,icov_pp,icov_oo]=mx_jCorv_LIF2Pop_delay(W_oo,W_pp,W_op,W_po,g_vec,id1_ob,id2_ob,id1_pc,id2_pc,ThresOB,ThresPC,[crOB;crPC],[sigO;sigP],Iap_OB,Iap_PC,N);
called in Matlab script
must compile first: >> mex mx_jCorv_LIF2Pop_delay.c
ALL TIMES IN MILLISECONDS
 (Allows) Het. thres.
 With separate Background corr for OB and PC
 ASSSUMING: N[ob/pc]*0.8 and N[ob/pc]*0.2 are both integers, and first 80% of N[pc] are excitatory, BUT first 20% of N[ob] are excitatory (last 80% INHIB,granule)
 
 Have synaptic delay from OB->PC and PC->OB, can change in tdel_frOB, tdel_frPC
*/

/*  Check for proper number of arguments */
if (nrhs !=  16) {
	mexErrMsgTxt("16 input arguments required.");
} else if (nlhs > 8) {
	mexErrMsgTxt("Too many output arguments.");
}

int Nob, Npc, Nobe, Npce, szRpc, szRob, i, j, k, ind1, N;
double taum=0.02, tRef=0.002, t_re=0.001, t_de=0.005, t_ri=0.002, t_di=0.01, Esyn=6.5, Isyn=-2.5, tdel_frOB=0.01, tdel_frPC=0.005; /* intrins neuro-params */
double ampE=1.0, ampI=2.0, sqtn; /* intrins neuro-params */

int l_woo, l_wpp, l_wop, l_wpo;
double *Woo, *Wpp, *Wop, *Wpo, *g_vec, *id1_obd, *id2_obd, *id1_pcd, *id2_pcd, *ThresOB, *ThresPC, *cr_pi, *sig_pi, *Iap_OB, *Iap_PC, *N_inp; /* passed-in */
double g_IO, g_ieo, g_eeo, g_iio, g_IP, g_iep, g_eep, g_iip, gCrs_oep, g_EP, gCrs_peo, g_EO; /* passed-in */
double GIO, Gieo, Geeo, Giio, GIP, Giep, Geep, Giip, GcrsOep, GEP, GcrsPeo, GEO; /* aux, used in time-loop to add-up all synaptic inputs */

double sigP, sigO, crOB, crPC; /* passed-in, not directly */
double c1, cob, c2, cpc, etac1, etac2, randc1, randc2=0.0, dt=0.0001; /* used in for-loop; etac is common-noise */
int Lt=20000, totTime=2, nrnSpace, Ltrans, LdelFrOB, LdelFrPC; /* 2 s for each realiz, store in totTime [sec] */

/* correl window params (ms); hard-coding in; COULD pass these pieces in as argument(s) */
int nmT=8, totWn=0, indx_wn, subIndx;
int numWins_rlz[nmT];
numWins_rlz[0]=40; /* 2s/50ms=20*2  */
numWins_rlz[1]=20; /* 2s/100ms=10*2 */
numWins_rlz[2]=10; /* 2s/200ms=5*2  */
numWins_rlz[3]=8;  /* 2s/250ms=4*2  */
numWins_rlz[4]=5;  /* 2s/400ms=5  */
numWins_rlz[5]=4;  /* 2s/500ms=4  */
numWins_rlz[6]=2;  /* 2s/1s=2  */
numWins_rlz[7]=1;  /* 2s/2s=1  */
for (k=0; k<nmT; k++) {
    totWn+=numWins_rlz[k];
}
   
/* outputs */
double *nuPC, *nuOB, *mn_PC, *mn_OB, *var_PC, *var_OB, *icov_pp, *icov_oo;
	
/* Using mxGetScalar to retreive input arguments */
Woo=mxGetPr(prhs[0]);		  /* Nob x 'Nob', (l_woo, smaller) vector */
Wpp=mxGetPr(prhs[1]);		  /* Npc x 'Npc' l_wpp vector */
Wop=mxGetPr(prhs[2]);         /* Nob x 'Npc' l_wop vector */
Wpo=mxGetPr(prhs[3]);         /* Npc x 'Nob' l_wpo vector */
g_vec=mxGetPr(prhs[4]);
id1_obd=mxGetPr(prhs[5]);    /* vector szRob x 1, 1st index of OB; CONVERT all 4 into INTs in initializ */
id2_obd=mxGetPr(prhs[6]);    /* vector szRob x 1, 2nd index of OB; */
id1_pcd=mxGetPr(prhs[7]);    /* vector szRpc x 1, 1st index of PC; */
id2_pcd=mxGetPr(prhs[8]);    /* vector szRpc x 1, 2nd index of PC; */
ThresOB=mxGetPr(prhs[9]);
ThresPC=mxGetPr(prhs[10]);
cr_pi=mxGetPr(prhs[11]);
sig_pi=mxGetPr(prhs[12]);
Iap_OB=mxGetPr(prhs[13]);
Iap_PC=mxGetPr(prhs[14]);
N_inp=mxGetPr(prhs[15]); /* Number of realizations */
    
Nob=mxGetM(prhs[9]);
Npc=mxGetM(prhs[10]);    /* # rows in connect matrix */
 Npce=(int)(0.8*Npc); /* making E 80% of Npc */
 /* Npci=(int)(0.2*Npc); not explicitly used (hard-code .8,.2), so remove */
 Nobe=(int)(0.2*Nob); /* making E 20% of Nob */
 /* Nobi=(int)(0.8*Nob); not explicitly used (hard-code .8,.2), so remove */
    
l_woo=mxGetN(prhs[0]); /* # cols in conn matrix; smaller than Npc and Nob to save time */
l_wpp=mxGetN(prhs[1]);
l_wop=mxGetN(prhs[2]);
l_wpo=mxGetN(prhs[3]);	
szRob=mxGetM(prhs[5]); /* #rows in id1_obd */
szRpc=mxGetM(prhs[7]);
    
g_IO=g_vec[0];
g_ieo=g_vec[1];
g_eeo=g_vec[2];
g_iio=g_vec[3];
g_IP=g_vec[4];
g_iep=g_vec[5];
g_eep=g_vec[6];
g_iip=g_vec[7];
gCrs_oep=g_vec[8];
g_EP=g_vec[9];
gCrs_peo=g_vec[10];
g_EO=g_vec[11];
    
crOB=cr_pi[0];  /*crOB is 1st entry, crPC is 2nd entry*/
crPC=cr_pi[1];
sigO=sig_pi[0]; /*sigO is 1st entry, sigP is 2nd entry*/
sigP=sig_pi[1];

N=(int)N_inp[0]; /* Number of realizations */
    
 Ltrans=(int)(0.25*Lt); /* quick run to get rid of transients... maybe not even necessary */
    
/* number of time steps before out of refractory */
nrnSpace=(int)(tRef/dt);
/* Size (number of time steps) delay until OB/PC pre input across */
LdelFrOB=(int)(tdel_frOB/dt);
LdelFrPC=(int)(tdel_frPC/dt);

	
/* OUTPUTs ... */
plhs[0]=mxCreateDoubleMatrix( Npc, 1, mxREAL); /* only mxREAL part */
nuPC = mxGetPr(plhs[0]);
plhs[1]=mxCreateDoubleMatrix( Nob, 1, mxREAL); /* only mxREAL part */
nuOB = mxGetPr(plhs[1]);
plhs[2]=mxCreateDoubleMatrix( Npc, nmT, mxREAL); /* only mxREAL part */
mn_PC = mxGetPr(plhs[2]);
plhs[3]=mxCreateDoubleMatrix( Nob, nmT, mxREAL); /* only mxREAL part */
mn_OB = mxGetPr(plhs[3]);
plhs[4]=mxCreateDoubleMatrix( Npc, nmT, mxREAL); /* only mxREAL part */
var_PC = mxGetPr(plhs[4]);
plhs[5]=mxCreateDoubleMatrix( Nob, nmT, mxREAL); /* only mxREAL part */
var_OB = mxGetPr(plhs[5]);
plhs[6]=mxCreateDoubleMatrix( szRpc, nmT, mxREAL); /* only mxREAL part */
icov_pp = mxGetPr(plhs[6]);
plhs[7]=mxCreateDoubleMatrix( szRob, nmT, mxREAL); /* only mxREAL part */
icov_oo = mxGetPr(plhs[7]);

/* store spkCnt in vector/array for each window size */
/* correl and counts for running sum */
double spkC_ob[Nob*totWn];  /* spike counts all FS(I) realiz; (Nob x totWn ) */
double spkC_pc[Npc*totWn];  /* spike counts all RS(E) realiz; (Npc x totWn ) */
double muPC1[szRpc*nmT]; /* mean of spike counts for PC-PC */
double muPC2[szRpc*nmT]; /* mean of spike counts for PC-PC */
double muOB1[szRob*nmT];  /* mean of spike counts for OB-OB (FS) */
double muOB2[szRob*nmT];  /* mean of spike counts for OB-OB (RS) */

    
/* used in for-loop; ODEs */
double vlt_ob[Nob];    /* voltage */
double synOB[Nob];  /* synapses */
double asynOB[Nob]; /* aux syn */
double vlt_pc[Npc];     
double synPC[Npc];
double asynPC[Npc];
double etaV[Nob+Npc]; /* indiv noise */
int TmPCspk[Npc];   /* for ref period, time last spike */
int TmOBspk[Nob];
int id1_pc[szRpc];
int id2_pc[szRpc];
int id1_ob[szRob];
int id2_ob[szRob];
int W_oo[Nob*l_woo];
int W_pp[Npc*l_wpp];
int W_op[Npc*l_wop];
int W_po[Nob*l_wpo];
 double delSynOB[Nob*LdelFrOB]; /*use in time-delay for cross syn input */
 double delSynPC[Npc*LdelFrPC]; /*use in time-delay for cross syn input */
	
/*  initialization */
srand ( time(NULL) ); /* seeing random # gen. */
sqtn=sqrt(1/dt);
c1 = sqrt((1-crOB));
cob = sqrt(crOB);
c2 = sqrt((1-crPC));
cpc = sqrt(crPC);
    

for (j=0; j<(Nob+Npc); j++) {
    if(j<Nob){
        nuOB[j]=0.0;
    }
    else{
        nuPC[j-Nob]=0.0;
    }
}
/* more initialization; set corr pieces to 0 */
for (j=0; j<nmT; j++) {
	for (k=0; k<szRpc; k++) {
		muPC1[j*szRpc+k]=0.0; 
		muPC2[j*szRpc+k]=0.0;
		icov_pp[j*szRpc+k]=0.0;
        if (j==0) { /* Only need to do it once */
            id1_pc[k]=(int)id1_pcd[k];
            id2_pc[k]=(int)id2_pcd[k];
        }
	}
	for (k=0; k<szRob; k++) {
		muOB1[j*szRob+k]=0.0; 
		muOB2[j*szRob+k]=0.0;
		icov_oo[j*szRob+k]=0.0;
        if (j==0) { /* Only need to do it once */
            id1_ob[k]=(int)id1_obd[k];
            id2_ob[k]=(int)id2_obd[k];
        }
	}
    for (k=0; k<Nob; k++) {
        mn_OB[j*Nob+k]=0.0;
        var_OB[j*Nob+k]=0.0;
    }
    for (k=0; k<Npc; k++) {
        mn_PC[j*Npc+k]=0.0;
        var_PC[j*Npc+k]=0.0;
    }
}

/* typecast all of the double inputs to ints (coupling matrices) */
for (j=0; j<Nob; j++) {
	for (k=0; k<l_woo; k++)
		W_oo[k*Nob+j]=(int)Woo[k*Nob+j];
	for (k=0; k<l_wop; k++)
		W_op[k*Nob+j]=(int)Wop[k*Nob+j];
}
for (j=0; j<Npc; j++) {
	for (k=0; k<l_wpp; k++)
		W_pp[k*Npc+j]=(int)Wpp[k*Npc+j];
	for (k=0; k<l_wpo; k++)
		W_po[k*Npc+j]=(int)Wpo[k*Npc+j];
}
	

/* set for FIRST run only; then Initi Cond. determined from t(end) of previous run */
for (j=0; j<(Nob+Npc); j++) {
	etaV[j] = ((double)(rand())/RAND_MAX-0.5)*2; /* rand unif i.c. [-1,1] */
	if(j<Nob){
		vlt_ob[j]=(double)(rand())/RAND_MAX; /* rand unif i.c. */
		synOB[j]=0.0;
		asynOB[j]=0.0;
        TmOBspk[j]=-nrnSpace;
	}
	else{
		vlt_pc[j-Nob]=(double)(rand())/RAND_MAX; /* rand unif i.c. */
		synPC[j-Nob]=0.0;
		asynPC[j-Nob]=0.0;
        TmPCspk[j-Nob]=-nrnSpace;
	}
}
 for (k=0; k<LdelFrOB; k++) { /* Initialize delayed syn to 0 */
    for (j=0; j<Nob; j++) {
        delSynOB[k*Nob+j]=0.0;
    }
}
 for (k=0; k<LdelFrPC; k++) { /* Initialize delayed syn to 0 */
    for (j=0; j<Npc; j++) {
        delSynPC[k*Npc+j]=0.0;
    }
}
etac1=((double)(rand())/RAND_MAX-0.5)*2;
etac2=((double)(rand())/RAND_MAX-0.5)*2;
		

/* --- Run it once to get rid of transients -- */
    for (j=0; j<Ltrans; j++){
        
        for (k=0; k<Nob; k++){
            if((j-TmOBspk[k]) >= nrnSpace ){ /* v change only if not in refrac */
                if(k<Nobe){
                    Geeo=0.0; /* recalc */
                    ind1=0;
                    while (W_oo[ind1*Nob+k]<Nobe && ind1<l_woo && W_oo[ind1*Nob+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one is >= Nobe (I-to-E in OB) */
                        Geeo+=(g_eeo*synOB[W_oo[ind1*Nob+k]]);
                        ind1++;
                    } /* no need to reset ind1; advances along the kth row */
                    GIO=0.0; /* recalc */
                    while (ind1<l_woo && W_oo[ind1*Nob+k]!=-1) {
                        GIO+=(g_IO*synOB[W_oo[ind1*Nob+k]]);
                        ind1++;
                    }
                    /* cross coupling, PC->OB */
                    GcrsOep=0.0; /* recalc */
                    ind1=0;
                    while (W_op[ind1*Nob+k]<Npce && W_op[ind1*Nob+k]!=-1 && ind1<l_wop){ /* ASSuming indices of coupling sorted (ascending) */
		      GcrsOep+=(gCrs_oep*delSynPC[W_op[ind1*Nob+k]]); /* using delayed synapse */
                        ind1++;
                    }/* no need to reset ind1; advances along the kth row */
                     /*  GcrsOeip=0.0; no I(pc)->OB; I's too short, so MUST have !=-1 check in loop */
                    
                    vlt_ob[k]=vlt_ob[k]+dt/taum*(Iap_OB[k]-vlt_ob[k]-(Geeo+GcrsOep)*(vlt_ob[k]-Esyn)-GIO*(vlt_ob[k]-Isyn)+sigO*sqtn*(c1*etaV[k] +cob*etac1));
                    /* if (vlt_ob[k]<Isyn) {  lower barrier 
                        vlt_ob[k]=2*Isyn-vlt_ob[k];
			} */
                    
                }
                else{ /*inhibitory OB cells */
                    Gieo=0.0; /*recalc */
                    ind1=0;
                    while (W_oo[ind1*Nob+k]<Nobe && ind1<l_woo && W_oo[ind1*Nob+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one is >= Nobe (I-to-I in OB) */
                        Gieo+=(g_ieo*synOB[W_oo[ind1*Nob+k]]);
                        ind1++;
                    } /* no need to reset ind1; advances along the kth row */
                    Giio=0.0; /* recalc */
                    while (ind1<l_woo && W_oo[ind1*Nob+k]!=-1) {
                        Giio+=(g_iio*synOB[W_oo[ind1*Nob+k]]);
                        ind1++;
                    }
                    /* cross coupling, PC->OB */
                    GEP=0.0; /* recalc */
                    ind1=0;
                    while (W_op[ind1*Nob+k]<Npce && W_op[ind1*Nob+k]!=-1 && ind1<l_wop){ /* ASSuming indices of coupling sorted (ascending) */
		      GEP+=(g_EP*delSynPC[W_op[ind1*Nob+k]]); /*using delayed synapse */
                        ind1++;
                    }/*no need to reset ind1; advances along the kth row */
		       /*  GcrsOeip=0.0; no I(pc)->OB; I's too short, so MUST have !=-1 check in loop */

                    vlt_ob[k]=vlt_ob[k]+dt/taum*(Iap_OB[k]-vlt_ob[k]-(Gieo+GEP)*(vlt_ob[k]-Esyn)-Giio*(vlt_ob[k]-Isyn)+sigO*sqtn*(c1*etaV[k] +cob*etac1));
                }
            }
            if (k<Nobe) {
                synOB[k]=synOB[k]+dt/t_de*(-synOB[k]+asynOB[k]);
                asynOB[k]=asynOB[k]+dt/t_re*(-asynOB[k]);
                /* spiked */
                if(vlt_ob[k]>=ThresOB[k]){
                    asynOB[k]+=ampE; /* update synapse */
                    vlt_ob[k]=0.0;		 /* reset voltage */
                    TmOBspk[k]=j;
                }
            }
            else{
                synOB[k]=synOB[k]+dt/t_di*(-synOB[k]+asynOB[k]);
                asynOB[k]=asynOB[k]+dt/t_ri*(-asynOB[k]);
                /* spiked */
                if(vlt_ob[k]>=ThresOB[k]){
                    asynOB[k]+=ampI; /* update synapse */
                    vlt_ob[k]=0.0;		 /* reset voltage */
                    TmOBspk[k]=j;
                }
            }
            /* Update delayed synapses */
            for (ind1=0; ind1<(LdelFrOB-1); ind1++) { /* move columns to left, except last column; only do kth row here! */
                delSynOB[ind1*Nob+k]=delSynOB[(ind1+1)*Nob+k];
            }
            delSynOB[(LdelFrOB-1)*Nob+k]=synOB[k]; /* set the last column to the current synOB value */
            
            if(k%2 == 0){
                z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                etaV[k]=randc1;
            }
            else
                etaV[k]=randc2;
        } /* end of for k=0 to Nob-1 */
        
        for (k=0; k<Npc; k++){
            if((j-TmPCspk[k]) >= nrnSpace ){ /* v change only if not in refrac */
                if(k<Npce){
                    Geep=0.0; /* recalc */
                    ind1=0;
                    while (W_pp[ind1*Npc+k]<Npce && ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one >=Npce (I-to-E in PC) */
                        Geep+=(g_eep*synPC[W_pp[ind1*Npc+k]]);
                        ind1++;
                    } /*no need to reset ind1; advances along the kth row */
                    GIP=0.0; /* recalc */
                    while (ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) {
                        GIP+=(g_IP*synPC[W_pp[ind1*Npc+k]]);
                        ind1++;
                    }
                    /* cross coupling, OB->PC */
                    GcrsPeo=0.0; /* recalc */
                    ind1=0;
                    while (W_po[ind1*Npc+k]<Nobe && W_po[ind1*Npc+k]!=-1 && ind1<l_wpo){ /* ASSuming indices of coupling sorted (ascending) */
		      GcrsPeo+=(gCrs_peo*delSynOB[W_po[ind1*Npc+k]]); /*using delayed synapse */
                        ind1++;
                    }/* no need to reset ind1; advances along the kth row */
                    /* GcrsPeio=0.0; no I(ob)->PC; I's too short, so MUST have !=-1 check in loop */

                    vlt_pc[k]=vlt_pc[k]+dt/taum*(Iap_PC[k]-vlt_pc[k]-(Geep+GcrsPeo)*(vlt_pc[k]-Esyn)-GIP*(vlt_pc[k]-Isyn)+sigP*sqtn*(c2*etaV[Nob+k] +cpc*etac2));
                }
                else{ /* inhibitory PC cells */
                    Giep=0.0; /*recalc */
                    ind1=0;
                    while (W_pp[ind1*Npc+k]<Npce && ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one >=Npce (I-to-I in PC) */
                        Giep+=(g_iep*synPC[W_pp[ind1*Npc+k]]);
                        ind1++;
                    } /* no need to reset ind1; advances along the kth row */
                    Giip=0.0; /* recalc */
                    while (ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) {
                        Giip+=(g_iip*synPC[W_pp[ind1*Npc+k]]);
                        ind1++;
                    }
                    /* cross coupling, OB->PC */
                    GEO=0.0; /* recalc */
                    ind1=0;
                    while (W_po[ind1*Npc+k]<Nobe && W_po[ind1*Npc+k]!=-1 && ind1<l_wpo){ /* ASSuming indices of coupling sorted (ascending) */
		      GEO+=(g_EO*delSynOB[W_po[ind1*Npc+k]]); /* using delayed synapse */
                        ind1++;
                    }/*no need to reset ind1; advances along the kth row */
                    /*GcrsPeio=0.0; no I(ob)->PC; I's too short, so MUST have !=-1 check in loop */
                
                    vlt_pc[k]=vlt_pc[k]+dt/taum*(Iap_PC[k]-vlt_pc[k]-(Giep+GEO)*(vlt_pc[k]-Esyn)-Giip*(vlt_pc[k]-Isyn)+sigP*sqtn*(c2*etaV[Nob+k] +cpc*etac2 ));
                }
            }
            if (k<Npce) {
                synPC[k]=synPC[k]+dt/t_de*(-synPC[k]+asynPC[k]);
                asynPC[k]=asynPC[k]+dt/t_re*(-asynPC[k]);
                /* spiked */
                if(vlt_pc[k]>=ThresPC[k]){
                    asynPC[k]+=ampE; /* update synapse */
                    vlt_pc[k]=0.0;		 /* reset voltage */
                    TmPCspk[k]=j;
                }
            }
            else{
                synPC[k]=synPC[k]+dt/t_di*(-synPC[k]+asynPC[k]);
                asynPC[k]=asynPC[k]+dt/t_ri*(-asynPC[k]);
                /* spiked */
                if(vlt_pc[k]>=ThresPC[k]){
                    asynPC[k]+=ampI; /* update synapse */
                    vlt_pc[k]=0.0;		 /* reset voltage */
                    TmPCspk[k]=j;
                }
            }
            /* Update delayed synapses */
            for (ind1=0; ind1<(LdelFrPC-1); ind1++) { /*move columns to left, except last column; only do kth row here! */
                delSynPC[ind1*Npc+k]=delSynPC[(ind1+1)*Npc+k];
            }
            delSynPC[(LdelFrPC-1)*Npc+k]=synPC[k]; /*set the last column to the current synPC value*/
            
            if(k%2 == 0){
                z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                etaV[Nob+k]=randc1;
            }
            else
                etaV[Nob+k]=randc2;
        } /*end of for k=0 to Npc-1 */
        
        z_randn(&etac1,&etac2); /* common noise for both OB and PC, independ */
        
}/* ending j-loop (transient time) */

    
    /*  MAIN....N realizations. */
    for (i=0; i<N; i++){
            for (j=0; j<(Nob+Npc); j++) { /* SET INITIAL CONDITIONS */
                if (j<Nob) {
                    for (k=0; k<totWn; k++) {
                        spkC_ob[j*totWn+k]=0.0; /* sloppy */
                    }
                    if (i==0) {
                        TmOBspk[j]=TmOBspk[j]-Ltrans;
                        /*    TmOBspk[j]=-nrnSpace; use if no transient loop Ltrans, on 1st run */
                    }
                    else{
                        TmOBspk[j]=TmOBspk[j]-Lt;
                    }
                }
                else {
                    for (k=0; k<totWn; k++) {
                        spkC_pc[(j-Nob)*totWn+k]=0.0; /* sloppy */
                    }
                    if (i==0) {
                        TmPCspk[j-Nob]=TmPCspk[j-Nob]-Ltrans;
                        /*    TmPCspk[j-Nob]=-nrnSpace; use if no transient loop Ltrans, on 1st run */
                    }
                    else{
                        TmPCspk[j-Nob]=TmPCspk[j-Nob]-Lt;
                    }
                        
                }
            }

        /* start of time-loop */
        for (j=0; j<Lt; j++){
                
            for (k=0; k<Nob; k++){
                if((j-TmOBspk[k]) >= nrnSpace ){ /* v change only if not in refrac */
                    if(k<Nobe){
                        Geeo=0.0; /* recalc */
                        ind1=0;
                        while (W_oo[ind1*Nob+k]<Nobe && ind1<l_woo && W_oo[ind1*Nob+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one is >= Nobe (I-to-E in OB) */
                            Geeo+=(g_eeo*synOB[W_oo[ind1*Nob+k]]);
                            ind1++;
                        } /* no need to reset ind1; advances along the kth row */
                        GIO=0.0; /* recalc */
                        while (ind1<l_woo && W_oo[ind1*Nob+k]!=-1) {
                            GIO+=(g_IO*synOB[W_oo[ind1*Nob+k]]);
                            ind1++;
                        }
                        /* cross coupling, PC->OB */
                        GcrsOep=0.0; /* recalc */
                        ind1=0;
                        while (W_op[ind1*Nob+k]<Npce && W_op[ind1*Nob+k]!=-1 && ind1<l_wop){ /* ASSuming indices of coupling sorted (ascending) */
			  GcrsOep+=(gCrs_oep*delSynPC[W_op[ind1*Nob+k]]); /* Using delayed synapse */
                            ind1++;
                        }/* no need to reset ind1; advances along the kth row */
                        /* GcrsOeip=0.0; no I(pc)->OB; I's too short, so MUST have !=-1 check in loop */
                        
                        vlt_ob[k]=vlt_ob[k]+dt/taum*(Iap_OB[k]-vlt_ob[k]-(Geeo+GcrsOep)*(vlt_ob[k]-Esyn)-GIO*(vlt_ob[k]-Isyn)+sigO*sqtn*(c1*etaV[k] +cob*etac1 ));
                        /* if (vlt_ob[k]<Isyn) {  lower barrier 
                            vlt_ob[k]=2*Isyn-vlt_ob[k];
                        } */
                        
                    }
                    else{ /* inhibitory OB cells */
                        Gieo=0.0; /*recalc */
                        ind1=0;
                        while (W_oo[ind1*Nob+k]<Nobe && ind1<l_woo && W_oo[ind1*Nob+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one is >= Nobe (I-to-I in OB) */
                            Gieo+=(g_ieo*synOB[W_oo[ind1*Nob+k]]);
                            ind1++;
                        } /* no need to reset ind1; advances along the kth row */
                        Giio=0.0; /* recalc */
                        while (ind1<l_woo && W_oo[ind1*Nob+k]!=-1) {
                            Giio+=(g_iio*synOB[W_oo[ind1*Nob+k]]);
                            ind1++;
                        }
                        /* cross coupling, PC->OB */
                        GEP=0.0; /* recalc */
                        ind1=0;
                        while (W_op[ind1*Nob+k]<Npce && W_op[ind1*Nob+k]!=-1 && ind1<l_wop){ /* ASSuming indices of coupling sorted (ascending) */
			  GEP+=(g_EP*delSynPC[W_op[ind1*Nob+k]]); /* Using delayed synapse */
                            ind1++;
                        }/*no need to reset ind1; advances along the kth row */
                        /*GcrsOeip=0.0; no I(pc)->OB; I's too short, so MUST have !=-1 check in loop */
                        
                        vlt_ob[k]=vlt_ob[k]+dt/taum*(Iap_OB[k]-vlt_ob[k]-(Gieo+GEP)*(vlt_ob[k]-Esyn)-Giio*(vlt_ob[k]-Isyn)+sigO*sqtn*(c1*etaV[k] +cob*etac1 ));
                    }
                }
                if (k<Nobe) {
                    synOB[k]=synOB[k]+dt/t_de*(-synOB[k]+asynOB[k]);
                    asynOB[k]=asynOB[k]+dt/t_re*(-asynOB[k]);
                    /* spiked */
                    if(vlt_ob[k]>=ThresOB[k]){
                        asynOB[k]+=ampE; /* update synapse */
                        vlt_ob[k]=0.0;		 /* reset voltage */
                        TmOBspk[k]=j;
                        
                        nuOB[k]+=(double)(1/(dt*Lt)); /* record spike */
                        for (ind1=0; ind1<nmT; ind1++) {/* add to window count */
                            indx_wn=(int)(j*dt*numWins_rlz[ind1]/totTime);
                            for (subIndx=ind1-1; subIndx>=0; subIndx--) { /*adding up previous windows so in the correct column (spkC_ob) */
                                indx_wn+=numWins_rlz[subIndx];
                            }
                            spkC_ob[indx_wn*Nob+k]+=1;
                        }
                    }
                }
                else{
                    synOB[k]=synOB[k]+dt/t_di*(-synOB[k]+asynOB[k]);
                    asynOB[k]=asynOB[k]+dt/t_ri*(-asynOB[k]);
                    /* spiked */
                    if(vlt_ob[k]>=ThresOB[k]){
                        asynOB[k]+=ampI; /* update synapse */
                        vlt_ob[k]=0.0;		 /* reset voltage */
                        TmOBspk[k]=j;
                        
                        nuOB[k]+=(double)(1/(dt*Lt)); /* record spike */
                        for (ind1=0; ind1<nmT; ind1++) {/* add to window count */
                            indx_wn=(int)(j*dt*numWins_rlz[ind1]/totTime);
                            for (subIndx=ind1-1; subIndx>=0; subIndx--) { /*adding up previous windows so in the correct column (spkC_ob) */
                                indx_wn+=numWins_rlz[subIndx];
                            }
                            spkC_ob[indx_wn*Nob+k]+=1;
                        }
                    }
                }
                /* Update delayed synapses */
                for (ind1=0; ind1<(LdelFrOB-1); ind1++) { /*move columns to left, except last column; only do kth row here! */
                    delSynOB[ind1*Nob+k]=delSynOB[(ind1+1)*Nob+k];
                }
                delSynOB[(LdelFrOB-1)*Nob+k]=synOB[k]; /*set the last column to the current synOB value */
                
                if(k%2 == 0){
                    z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                    etaV[k]=randc1;
                }
                else
                    etaV[k]=randc2;
            }/* end of for k=0 to Nob-1 */
            
            for (k=0; k<Npc; k++){
                if((j-TmPCspk[k]) >= nrnSpace ){ /* v change only if not in refrac */
                    if(k<Npce){
                        Geep=0.0; /* recalc */
                        ind1=0;
                        while (W_pp[ind1*Npc+k]<Npce && ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one >=Npce (I-to-E in PC) */
                            Geep+=(g_eep*synPC[W_pp[ind1*Npc+k]]);
                            ind1++;
                        } /*no need to reset ind1; advances along the kth row */
                        GIP=0.0; /* recalc */
                        while (ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) {
                            GIP+=(g_IP*synPC[W_pp[ind1*Npc+k]]);
                            ind1++;
                        }
                        /* cross coupling, OB->PC */
                        GcrsPeo=0.0; /* recalc */
                        ind1=0;
                        while (W_po[ind1*Npc+k]<Nobe && W_po[ind1*Npc+k]!=-1 && ind1<l_wpo){ /* ASSuming indices of coupling sorted (ascending) */
			  GcrsPeo+=(gCrs_peo*delSynOB[W_po[ind1*Npc+k]]); /*using delayed synapse */
                            ind1++;
                        }/* no need to reset ind1; advances along the kth row */
                        /* GcrsPeio=0.0; no I(ob)->PC; I's too short, so MUST have !=-1 check in loop */
                        
                        vlt_pc[k]=vlt_pc[k]+dt/taum*(Iap_PC[k]-vlt_pc[k]-(Geep+GcrsPeo)*(vlt_pc[k]-Esyn)-GIP*(vlt_pc[k]-Isyn)+sigP*sqtn*(c2*etaV[Nob+k] +cpc*etac2 ));
                    }
                    else{ /* inhibitory PC cells */
                        Giep=0.0; /*recalc */
                        ind1=0;
                        while (W_pp[ind1*Npc+k]<Npce && ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) { /* ASSuming indices of coupling sorted (ascending) AND at least one >=Npce (I-to-I in PC) */
                            Giep+=(g_iep*synPC[W_pp[ind1*Npc+k]]);
                            ind1++;
                        } /* no need to reset ind1; advances along the kth row */
                        Giip=0.0; /* recalc */
                        while (ind1<l_wpp && W_pp[ind1*Npc+k]!=-1) {
                            Giip+=(g_iip*synPC[W_pp[ind1*Npc+k]]);
                            ind1++;
                        }
                        /* cross coupling, OB->PC */
                        GEO=0.0; /* recalc */
                        ind1=0;
                        while (W_po[ind1*Npc+k]<Nobe && W_po[ind1*Npc+k]!=-1 && ind1<l_wpo){ /* ASSuming indices of coupling sorted (ascending) */
			  GEO+=(g_EO*delSynOB[W_po[ind1*Npc+k]]); /*using delayed synapse */
                            ind1++;
                        }/* no need to reset ind1; advances along the kth row */
                        /* GcrsPeio=0.0; no I(ob)->PC; I's too short, so MUST have !=-1 check in loop */
                        
                        vlt_pc[k]=vlt_pc[k]+dt/taum*(Iap_PC[k]-vlt_pc[k]-(Giep+GEO)*(vlt_pc[k]-Esyn)-Giip*(vlt_pc[k]-Isyn)+sigP*sqtn*(c2*etaV[Nob+k] +cpc*etac2 ));
                    }
                }
                if (k<Npce) {
                    synPC[k]=synPC[k]+dt/t_de*(-synPC[k]+asynPC[k]);
                    asynPC[k]=asynPC[k]+dt/t_re*(-asynPC[k]);
                    /* spiked */
                    if(vlt_pc[k]>=ThresPC[k]){
                        asynPC[k]+=ampE; /* update synapse */
                        vlt_pc[k]=0.0;		 /* reset voltage */
                        TmPCspk[k]=j;
                        
                        nuPC[k]+=(double)(1/(dt*Lt));         /* record spike */
                        for (ind1=0; ind1<nmT; ind1++) {/* add to window count */
                            indx_wn=(int)(j*dt*numWins_rlz[ind1]/totTime);
                            for (subIndx=ind1-1; subIndx>=0; subIndx--) {   /* adding up previous windows so in the correct column (spkC_pc) */
                                indx_wn+=numWins_rlz[subIndx];
                            }
                            spkC_pc[indx_wn*Npc+k]+=1;
                        }
                    }
                }
                else{
                    synPC[k]=synPC[k]+dt/t_di*(-synPC[k]+asynPC[k]);
                    asynPC[k]=asynPC[k]+dt/t_ri*(-asynPC[k]);
                    /* spiked */
                    if(vlt_pc[k]>=ThresPC[k]){
                        asynPC[k]+=ampI; /* update synapse */
                        vlt_pc[k]=0.0;		 /* reset voltage */
                        TmPCspk[k]=j;
                        
                        nuPC[k]+=(double)(1/(dt*Lt));         /* record spike */
                        for (ind1=0; ind1<nmT; ind1++) {/* add to window count */
                            indx_wn=(int)(j*dt*numWins_rlz[ind1]/totTime);
                            for (subIndx=ind1-1; subIndx>=0; subIndx--) {   /* adding up previous windows so in the correct column (spkC_pc) */
                                indx_wn+=numWins_rlz[subIndx];
                            }
                            spkC_pc[indx_wn*Npc+k]+=1;
                        }
                    }
                }
                /* Update delayed synapses */
                for (ind1=0; ind1<(LdelFrPC-1); ind1++) { /* move columns to left, except last column; only do kth row here! */
                    delSynPC[ind1*Npc+k]=delSynPC[(ind1+1)*Npc+k];
                }
                delSynPC[(LdelFrPC-1)*Npc+k]=synPC[k]; /* set the last column to the current synPC value */
                
                if(k%2 == 0){
                    z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                    etaV[Nob+k]=randc1;
                }
                else
                    etaV[Nob+k]=randc2;
            } /* end of for k=0 to Npc-1 */
            
            z_randn(&etac1,&etac2); /* common noise for both OB and PC, independ */
            
        } /* ending time (j=1:Lt) loop */
        
        /* UPDATE running sum so unbiased-estim of var,cov, etc */
        for (j=0; j<nmT; j++) {
	  /* starting index, depending on j (nmT); to loop over all windows */
            subIndx=0;
            for (k=j-1; k>=0; k--) {  /* calc correct starting index for each windowSize */
                subIndx+=numWins_rlz[k];
            }
            for (k=0; k<szRpc; k++) {
	      for (ind1=subIndx; ind1<(subIndx+numWins_rlz[j]); ind1++) { /* loop through all windows */
                    muPC1[j*szRpc+k]+=spkC_pc[ind1*Npc+id1_pc[k]];
                    muPC2[j*szRpc+k]+=spkC_pc[ind1*Npc+id2_pc[k]];
                    icov_pp[j*szRpc+k]+=spkC_pc[ind1*Npc+id1_pc[k]]*spkC_pc[ind1*Npc+id2_pc[k]]; /*2nd moment for now */
                }
            }
            for (k=0; k<szRob; k++) {
	      for (ind1=subIndx; ind1<(subIndx+numWins_rlz[j]); ind1++) { /*loop through all windows */
                    muOB1[j*szRob+k]+=spkC_ob[ind1*Nob+id1_ob[k]];
                    muOB2[j*szRob+k]+=spkC_ob[ind1*Nob+id2_ob[k]];
                    icov_oo[j*szRob+k]+=spkC_ob[ind1*Nob+id1_ob[k]]*spkC_ob[ind1*Nob+id2_ob[k]]; /*2nd moment for now */
                }
            }
            for (k=0; k<(Nob+Npc); k++) {
                if (k<Nob) {
		  for (ind1=subIndx; ind1<(subIndx+numWins_rlz[j]); ind1++) { /* loop through all windows */
                        mn_OB[j*Nob+k]+=spkC_ob[ind1*Nob+k];
                        var_OB[j*Nob+k]+=spkC_ob[ind1*Nob+k]*spkC_ob[ind1*Nob+k];  /*2nd moment for now */
                    }
                }
                else{
		  for (ind1=subIndx; ind1<(subIndx+numWins_rlz[j]); ind1++) { /* loop through all windows */
                        mn_PC[j*Npc+k-Nob]+=spkC_pc[ind1*Npc+k-Nob];
                        var_PC[j*Npc+k-Nob]+=spkC_pc[ind1*Npc+k-Nob]*spkC_pc[ind1*Npc+k-Nob];  /*2nd moment for now */
                    }
                }
            }
        } /*end of j=0 to nmT */
        
        
    } /* ending i-loop (realizations) */

/* Normalize firing rates; diving by N */
for (k=0; k<Npc; k++) {
    nuPC[k]/=(N);
}
for (k=0; k<Nob; k++) {
    nuOB[k]/=(N);
}
    
            /*  Normalization by N (mu,v,cv).. unraveling var,cov, etc */
            for (j=0; j<nmT; j++) {
                for (k=0; k<szRpc; k++) {
                    muPC1[j*szRpc+k]/=(N*numWins_rlz[j]);
                    muPC2[j*szRpc+k]/=(N*numWins_rlz[j]);
                    icov_pp[j*szRpc+k] = ( icov_pp[j*szRpc+k] - (N*numWins_rlz[j])*muPC1[j*szRpc+k]*muPC2[j*szRpc+k] )/((N*numWins_rlz[j])-1);
                }
                for (k=0; k<szRob; k++) {
                    muOB1[j*szRob+k]/=(N*numWins_rlz[j]);
                    muOB2[j*szRob+k]/=(N*numWins_rlz[j]);
                    icov_oo[j*szRob+k] = ( icov_oo[j*szRob+k] - (N*numWins_rlz[j])*muOB1[j*szRob+k]*muOB2[j*szRob+k] )/((N*numWins_rlz[j])-1);
                }
                for (k=0; k<(Nob+Npc); k++) {
                    if (k<Nob) {
                        mn_OB[j*Nob+k]/=(N*numWins_rlz[j]);
                        var_OB[j*Nob+k] = ( var_OB[j*Nob+k] - (N*numWins_rlz[j])*mn_OB[j*Nob+k]*mn_OB[j*Nob+k] )/((N*numWins_rlz[j])-1);
                    }
                    else{
                        mn_PC[j*Npc+k-Nob]/=(N*numWins_rlz[j]);
                        var_PC[j*Npc+k-Nob] = ( var_PC[j*Npc+k-Nob] - (N*numWins_rlz[j])*mn_PC[j*Npc+k-Nob]*mn_PC[j*Npc+k-Nob] )/((N*numWins_rlz[j])-1);
                    }
                }
            }
    


return;
                
}
