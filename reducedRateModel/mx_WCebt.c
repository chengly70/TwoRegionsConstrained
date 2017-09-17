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
/* the transfer function F */
double F_io(double xinp, double s_rv, double s_sp, double nuMax)
{
    return nuMax/2*(1+tanh( (xinp-s_rv)/s_sp ));
}

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ])
{
/* 
 
[cov_pcFx,cov_obFx,vr_Fs,mn_Fs,cov_pcX,cov_obX,vr_Xs,mn_Xs,EBcov_pcFx,EBcov_obFx,EBvr_Fs,EBmn_Fs,EBcov_pcX,EBcov_obX,EBvr_Xs,EBmn_Xs]=mx_WCebt(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,gMat);
called in Matlab script
must compile first: >> mex mx_WCebt.c, very similar to mx_WC.c but NOW have error bars
!!!! NOTE: MUST take sqrt of EB[] to get the std. dev; this program returns the VARIANCE of deviations (error bars) !!!!
Time is dimensionless

*/  

/*  Check for proper number of arguments */
if (nrhs !=  8) {
	mexErrMsgTxt("8 input arguments required.");
} else if (nlhs > 16) {
	mexErrMsgTxt("Too many output arguments.");
}

int i, j, k, ind1, N_relz=3000, nstd=10, numStr=300;

double *mu_vec, *sig_vec, *cr_obp, *cr_pcp, *nu_maxp, *s_rvp, *s_spp, *gMat; /* passed-in */
double crOB, crPC, nuMax, s_rv, s_sp; /*  passed-in */

/* used in for-loop; etac is common-noise */
double sqtn, dt=0.01;
double c1, cob, c2, cpc, etac1, etac2, randc1, randc2=0.0, G_ob, G_pc;
    
int Lt=500000; /* 500[time] for N_relz realiz */
  
/* outputs */
double *cov_pcFx, *cov_obFx, *vr_Fs, *mn_Fs, *cov_pcX, *cov_obX, *vr_Xs, *mn_Xs, *EBcov_pcFx, *EBcov_obFx, *EBvr_Fs, *EBmn_Fs, *EBcov_pcX, *EBcov_obX, *EBvr_Xs, *EBmn_Xs;
	
/* Using mxGetScalar to retreive input arguments */
mu_vec=mxGetPr(prhs[0]);      /* 6 x 1 vector */
sig_vec=mxGetPr(prhs[1]);      /* 6 x 1 vector */
cr_obp=mxGetPr(prhs[2]);	   /* scalar */
cr_pcp=mxGetPr(prhs[3]);	   /* scalar */
nu_maxp=mxGetPr(prhs[4]);	   /* scalar */
s_rvp=mxGetPr(prhs[5]);        /* scalar */
s_spp=mxGetPr(prhs[6]);        /* scalar */
gMat=mxGetPr(prhs[7]);		  /* 6 x 6 vector */

/* set input scalar values to correct var */
crOB=cr_obp[0];
crPC=cr_pcp[0];
nuMax=nu_maxp[0];
s_rv=s_rvp[0];
s_sp=s_spp[0];
	
	
/* OUTPUTs ... */
plhs[0]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
cov_pcFx = mxGetPr(plhs[0]);
plhs[1]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
cov_obFx= mxGetPr(plhs[1]);
plhs[2]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
vr_Fs = mxGetPr(plhs[2]);
plhs[3]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
mn_Fs = mxGetPr(plhs[3]);
    plhs[4]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
    cov_pcX = mxGetPr(plhs[4]);
    plhs[5]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
    cov_obX = mxGetPr(plhs[5]);
    plhs[6]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
    vr_Xs = mxGetPr(plhs[6]);
    plhs[7]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
    mn_Xs = mxGetPr(plhs[7]);
plhs[8]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
EBcov_pcFx = mxGetPr(plhs[8]);
plhs[9]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
EBcov_obFx= mxGetPr(plhs[9]);
plhs[10]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
EBvr_Fs = mxGetPr(plhs[10]);
plhs[11]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
EBmn_Fs = mxGetPr(plhs[11]);
    plhs[12]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
    EBcov_pcX = mxGetPr(plhs[12]);
    plhs[13]=mxCreateDoubleMatrix( 3, 1, mxREAL); /* only mxREAL part */
    EBcov_obX = mxGetPr(plhs[13]);
    plhs[14]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
    EBvr_Xs = mxGetPr(plhs[14]);
    plhs[15]=mxCreateDoubleMatrix( 6, 1, mxREAL); /* only mxREAL part */
    EBmn_Xs = mxGetPr(plhs[15]);
    
/* mean, var, cov for running sum */
double mu_x_rs[6];         /* for running sum of mean Xs */
double var_x_rs[6];
double cov_pcx_rs[3];
double cov_obx_rs[3];
double mu_F_rs[6];         /* for running sum of mean F(x)s */
double var_F_rs[6];
double cov_pcF_rs[3];
double cov_obF_rs[3];
    double mu_x_EB[6];         /* for Error Bars, tmp */
    double var_x_EB[6];
    double cov_pcx_EB[3];
    double cov_obx_EB[3];
    double mu_F_EB[6];
    double var_F_EB[6];
    double cov_pcF_EB[3];
    double cov_obF_EB[3];
double tmp_eb1, tmp_eb2;
    
/* used in for-loop; ODEs */
double x_ob[3];    /* voltage or activity */
double x_pc[3];    /* voltage or activity */
double etaV[6]; /* indiv noise */
	
/*  initialization */
srand ( time(NULL) ); /* seeing random # gen. */
sqtn=sqrt(1/dt);
c1 = sqrt((1-crOB));
cob = sqrt(crOB);
c2 = sqrt((1-crPC));
cpc = sqrt(crPC);
    /* initializing 16 outputs */
for (j=0; j<6; j++) {
    vr_Fs[j]=0.0;
    mn_Fs[j]=0.0;
    vr_Xs[j]=0.0;
    mn_Xs[j]=0.0;
    EBvr_Fs[j]=0.0;
    EBmn_Fs[j]=0.0;
    EBvr_Xs[j]=0.0;
    EBmn_Xs[j]=0.0;
    if (j<3) {
        cov_obFx[j]=0.0;
        cov_obX[j]=0.0;
        EBcov_obFx[j]=0.0;
        EBcov_obX[j]=0.0;
    }
    else{
        cov_pcFx[j-3]=0.0;
        cov_pcX[j-3]=0.0;
        EBcov_pcFx[j-3]=0.0;
        EBcov_pcX[j-3]=0.0;
    }
}

	
srand ( time(NULL) ); /* seeing random # gen. */

/* set for FIRST run only; then Initi Cond. determined from t(end) of previous run */
for (j=0; j<6; j++) {
	etaV[j] = ((double)(rand())/RAND_MAX-0.5)*2; /* rand unif i.c. [-1,1] */
}
etac1=((double)(rand())/RAND_MAX-0.5)*2;
etac2=((double)(rand())/RAND_MAX-0.5)*2;
		

/* --- Run it once to get rid of transients -- */
    for (j=0; j<1000; j++){
        
        for (k=0; k<3; k++){
                G_ob=0.0; /* recalc */
            for (ind1=0; ind1<6; ind1++) {
                if (ind1<3) {
                    G_ob+=( gMat[ind1*6+k]*F_io(x_ob[ind1],s_rv,s_sp,nuMax) );
                }
                else{
                    G_ob+=( gMat[ind1*6+k]*F_io(x_pc[ind1-3],s_rv,s_sp,nuMax) );
                }
            }

            x_ob[k]=x_ob[k]+dt*( -x_ob[k]+mu_vec[k]+sig_vec[k]*sqtn*(c1*etaV[k]+cob*etac1) + G_ob );
            

            if(k%2 == 0){
                z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                etaV[k]=randc1;
            }
            else
                etaV[k]=randc2;
        }
        for (k=0; k<3; k++){
                G_pc=0.0; /* recalc */
            for (ind1=0; ind1<6; ind1++) {
                if (ind1<3) {
                    G_pc+=( gMat[ind1*6+(k+3)]*F_io(x_ob[ind1],s_rv,s_sp,nuMax) );
                }
                else{
                    G_pc+=( gMat[ind1*6+(k+3)]*F_io(x_pc[ind1-3],s_rv,s_sp,nuMax) );
                }
            }
            
            x_pc[k]=x_pc[k]+dt*(-x_pc[k]+mu_vec[k+3]+sig_vec[k+3]*sqtn*(c2*etaV[k+3] +cpc*etac2) + G_pc);
            

            if(k%2 == 0){
                z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
                etaV[k+3]=randc1;
            }
            else
                etaV[k+3]=randc2;
        }
        z_randn(&etac1,&etac2); /* common noise for both OB and PC, independ */
        
    }/* ending j-loop (transient time) */
    


/*  MAIN....N realizations. */
for (i=0; i<N_relz; i++){
    
    /* initaliz; restart running sums */
    for (j=0; j<6; j++) {
        mu_x_rs[j]=0.0;
        var_x_rs[j]=0.0;
        mu_F_rs[j]=0.0;
        var_F_rs[j]=0.0;
        if (j<3) {
            cov_obF_rs[j]=0.0;
            cov_obx_rs[j]=0.0;
        }
        else{
            cov_pcF_rs[j-3]=0.0;
            cov_pcx_rs[j-3]=0.0;
        }
    }
    /* for error bars: only update every numStr; set all tmp pieces []_EB to 0 */
    if (i%numStr == 0) {
        for (j=0; j<6; j++) {
            mu_x_EB[j]=0.0;
            var_x_EB[j]=0.0;
            mu_F_EB[j]=0.0;
            var_F_EB[j]=0.0;
            if (j<3) {
                cov_obF_EB[j]=0.0;
                cov_obx_EB[j]=0.0;
            }
            else{
                cov_pcF_EB[j-3]=0.0;
                cov_pcx_EB[j-3]=0.0;
            }
        }
    }

    
/* ------ start of time-loop ------ */
for (j=0; j<Lt; j++){
    for (k=0; k<3; k++){
        G_ob=0.0; /* recalc */
        for (ind1=0; ind1<6; ind1++) {
            if (ind1<3) {
                G_ob+=( gMat[ind1*6+k]*F_io(x_ob[ind1],s_rv,s_sp,nuMax) );
            }
            else{
                G_ob+=( gMat[ind1*6+k]*F_io(x_pc[ind1-3],s_rv,s_sp,nuMax) );
            }
        }
        
        x_ob[k]=x_ob[k]+dt*( -x_ob[k]+mu_vec[k]+sig_vec[k]*sqtn*(c1*etaV[k]+cob*etac1) + G_ob );
        
        
        if(k%2 == 0){
            z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
            etaV[k]=randc1;
        }
        else
            etaV[k]=randc2;
    }
    for (k=0; k<3; k++){
        G_pc=0.0; /* recalc */
        for (ind1=0; ind1<6; ind1++) {
            if (ind1<3) {
                G_pc+=( gMat[ind1*6+(k+3)]*F_io(x_ob[ind1],s_rv,s_sp,nuMax) );
            }
            else{
                G_pc+=( gMat[ind1*6+(k+3)]*F_io(x_pc[ind1-3],s_rv,s_sp,nuMax) );
            }
        }
        
        x_pc[k]=x_pc[k]+dt*(-x_pc[k]+mu_vec[k+3]+sig_vec[k+3]*sqtn*(c2*etaV[k+3] +cpc*etac2) + G_pc);
        
        
        if(k%2 == 0){
            z_randn(&randc1,&randc2); /*k starts at 0, so randc2 will be assigned */
            etaV[k+3]=randc1;
        }
        else
            etaV[k+3]=randc2;
    }
    z_randn(&etac1,&etac2); /* common noise for both OB and PC, independ */
    
    /* Store running sums; hardcode the cov parts */
    cov_obx_rs[0]+= (x_ob[0]*x_ob[1]); //2nd moment for now
    cov_obx_rs[1]+= (x_ob[0]*x_ob[2]); //2nd moment for now
    cov_obx_rs[2]+= (x_ob[1]*x_ob[2]); //2nd moment for now
    cov_obF_rs[0]+= F_io(x_ob[0],s_rv,s_sp,nuMax)*F_io(x_ob[1],s_rv,s_sp,nuMax); //2nd moment for now
    cov_obF_rs[1]+= F_io(x_ob[0],s_rv,s_sp,nuMax)*F_io(x_ob[2],s_rv,s_sp,nuMax); //2nd moment for now
    cov_obF_rs[2]+= F_io(x_ob[1],s_rv,s_sp,nuMax)*F_io(x_ob[2],s_rv,s_sp,nuMax); //2nd moment for now
    cov_pcx_rs[0]+= (x_pc[0]*x_pc[1]); //2nd moment for now
    cov_pcx_rs[1]+= (x_pc[0]*x_pc[2]); //2nd moment for now
    cov_pcx_rs[2]+= (x_pc[1]*x_pc[2]); //2nd moment for now
    cov_pcF_rs[0]+= F_io(x_pc[0],s_rv,s_sp,nuMax)*F_io(x_pc[1],s_rv,s_sp,nuMax); //2nd moment for now
    cov_pcF_rs[1]+= F_io(x_pc[0],s_rv,s_sp,nuMax)*F_io(x_pc[2],s_rv,s_sp,nuMax); //2nd moment for now
    cov_pcF_rs[2]+= F_io(x_pc[1],s_rv,s_sp,nuMax)*F_io(x_pc[2],s_rv,s_sp,nuMax); //2nd moment for now
    for (k=0; k<6; k++) {
        if (k<3) {
            mu_x_rs[k]+=x_ob[k];
            var_x_rs[k]+=(x_ob[k]*x_ob[k]); //2nd moment for now
            mu_F_rs[k]+= F_io(x_ob[k],s_rv,s_sp,nuMax);
            var_F_rs[k]+= F_io(x_ob[k],s_rv,s_sp,nuMax)*F_io(x_ob[k],s_rv,s_sp,nuMax); //2nd moment for now
        }
        else{
            mu_x_rs[k]+=x_pc[k-3];
            var_x_rs[k]+=(x_pc[k-3]*x_pc[k-3]); //2nd moment for now
            mu_F_rs[k]+= F_io(x_pc[k-3],s_rv,s_sp,nuMax);
            var_F_rs[k]+= F_io(x_pc[k-3],s_rv,s_sp,nuMax)*F_io(x_pc[k-3],s_rv,s_sp,nuMax); //2nd moment for now
        }
    }
    /* repeat store commands above but now for ERROR BARS */
    cov_obx_EB[0]+= (x_ob[0]*x_ob[1]); //2nd moment for now
    cov_obx_EB[1]+= (x_ob[0]*x_ob[2]); //2nd moment for now
    cov_obx_EB[2]+= (x_ob[1]*x_ob[2]); //2nd moment for now
    cov_obF_EB[0]+= (F_io(x_ob[0],s_rv,s_sp,nuMax)*F_io(x_ob[1],s_rv,s_sp,nuMax)); //2nd moment for now
    cov_obF_EB[1]+= (F_io(x_ob[0],s_rv,s_sp,nuMax)*F_io(x_ob[2],s_rv,s_sp,nuMax)); //2nd moment for now
    cov_obF_EB[2]+= (F_io(x_ob[1],s_rv,s_sp,nuMax)*F_io(x_ob[2],s_rv,s_sp,nuMax)); //2nd moment for now
    cov_pcx_EB[0]+= (x_pc[0]*x_pc[1]); //2nd moment for now
    cov_pcx_EB[1]+= (x_pc[0]*x_pc[2]); //2nd moment for now
    cov_pcx_EB[2]+= (x_pc[1]*x_pc[2]); //2nd moment for now
    cov_pcF_EB[0]+= (F_io(x_pc[0],s_rv,s_sp,nuMax)*F_io(x_pc[1],s_rv,s_sp,nuMax)); //2nd moment for now
    cov_pcF_EB[1]+= (F_io(x_pc[0],s_rv,s_sp,nuMax)*F_io(x_pc[2],s_rv,s_sp,nuMax)); //2nd moment for now
    cov_pcF_EB[2]+= (F_io(x_pc[1],s_rv,s_sp,nuMax)*F_io(x_pc[2],s_rv,s_sp,nuMax)); //2nd moment for now
    for (k=0; k<6; k++) {
        if (k<3) {
            mu_x_EB[k]+=x_ob[k];
            var_x_EB[k]+=(x_ob[k]*x_ob[k]); //2nd moment for now
            mu_F_EB[k]+= F_io(x_ob[k],s_rv,s_sp,nuMax);
            var_F_EB[k]+= (F_io(x_ob[k],s_rv,s_sp,nuMax)*F_io(x_ob[k],s_rv,s_sp,nuMax)); //2nd moment for now
        }
        else{
            mu_x_EB[k]+=x_pc[k-3];
            var_x_EB[k]+=(x_pc[k-3]*x_pc[k-3]); //2nd moment for now
            mu_F_EB[k]+= F_io(x_pc[k-3],s_rv,s_sp,nuMax);
            var_F_EB[k]+= (F_io(x_pc[k-3],s_rv,s_sp,nuMax)*F_io(x_pc[k-3],s_rv,s_sp,nuMax)); //2nd moment for now
        }
    }
    
}/* ending j-loop (time) */

    /* Update the actual output variables; not scaling by N_relz yet */
    for (k=0; k<6; k++) {
        mn_Xs[k] += mu_x_rs[k]/Lt;
        mn_Fs[k] += mu_F_rs[k]/Lt;
        vr_Xs[k] += var_x_rs[k]/Lt;
        vr_Fs[k] += var_F_rs[k]/Lt;
    }
    cov_obX[0] += cov_obx_rs[0]/Lt;
    cov_obX[1] += cov_obx_rs[1]/Lt;
    cov_obX[2] += cov_obx_rs[2]/Lt;
    cov_obFx[0] += cov_obF_rs[0]/Lt;
    cov_obFx[1] += cov_obF_rs[1]/Lt;
    cov_obFx[2] += cov_obF_rs[2]/Lt;
    cov_pcX[0] += cov_pcx_rs[0]/Lt;
    cov_pcX[1] += cov_pcx_rs[1]/Lt;
    cov_pcX[2] += cov_pcx_rs[2]/Lt;
    cov_pcFx[0] += cov_pcF_rs[0]/Lt;
    cov_pcFx[1] += cov_pcF_rs[1]/Lt;
    cov_pcFx[2] += cov_pcF_rs[2]/Lt;

    
    /* Only update the running sum every numStr, store resutls from []_EB in EB[]; squaring stats b/c computing std. error below */
    /* after this, EB[] is the 2nd moment of the std. err bar */
    if (i%numStr == numStr-1) {
        for (k=0; k<6; k++) {
            EBmn_Xs[k] += (mu_x_EB[k]/(Lt*numStr)*mu_x_EB[k]/(Lt*numStr));
            EBmn_Fs[k] += (mu_F_EB[k]/(Lt*numStr)*mu_F_EB[k]/(Lt*numStr));
                tmp_eb1=(var_x_EB[k]-mu_x_EB[k]*mu_x_EB[k]/(Lt*numStr))/(Lt*numStr-1);
                tmp_eb2=(var_x_EB[k]-mu_x_EB[k]*mu_x_EB[k]/(Lt*numStr))/(Lt*numStr-1);
            EBvr_Xs[k] += tmp_eb1*tmp_eb2;
                tmp_eb1=(var_F_EB[k]-mu_F_EB[k]*mu_F_EB[k]/(Lt*numStr))/(Lt*numStr-1);
                tmp_eb2=(var_F_EB[k]-mu_F_EB[k]*mu_F_EB[k]/(Lt*numStr))/(Lt*numStr-1);
            EBvr_Fs[k] += tmp_eb1*tmp_eb2;

        }
            tmp_eb1=(cov_obx_EB[0]-mu_x_EB[0]*mu_x_EB[1]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obx_EB[0]-mu_x_EB[0]*mu_x_EB[1]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obX[0] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_obx_EB[1]-mu_x_EB[0]*mu_x_EB[2]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obx_EB[1]-mu_x_EB[0]*mu_x_EB[2]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obX[1] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_obx_EB[2]-mu_x_EB[1]*mu_x_EB[2]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obx_EB[2]-mu_x_EB[1]*mu_x_EB[2]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obX[2] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_obF_EB[0]-mu_F_EB[0]*mu_F_EB[1]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obF_EB[0]-mu_F_EB[0]*mu_F_EB[1]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obFx[0] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_obF_EB[1]-mu_F_EB[0]*mu_F_EB[2]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obF_EB[1]-mu_F_EB[0]*mu_F_EB[2]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obFx[1] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_obF_EB[2]-mu_F_EB[1]*mu_F_EB[2]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_obF_EB[2]-mu_F_EB[1]*mu_F_EB[2]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_obFx[2] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcx_EB[0]-mu_x_EB[3]*mu_x_EB[4]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcx_EB[0]-mu_x_EB[3]*mu_x_EB[4]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcX[0] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcx_EB[1]-mu_x_EB[3]*mu_x_EB[5]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcx_EB[1]-mu_x_EB[3]*mu_x_EB[5]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcX[1] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcx_EB[2]-mu_x_EB[4]*mu_x_EB[5]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcx_EB[2]-mu_x_EB[4]*mu_x_EB[5]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcX[2] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcF_EB[0]-mu_F_EB[3]*mu_F_EB[4]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcF_EB[0]-mu_F_EB[3]*mu_F_EB[4]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcFx[0] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcF_EB[1]-mu_F_EB[3]*mu_F_EB[5]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcF_EB[1]-mu_F_EB[3]*mu_F_EB[5]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcFx[1] += tmp_eb1*tmp_eb2;
            tmp_eb1=(cov_pcF_EB[2]-mu_F_EB[4]*mu_F_EB[5]/(Lt*numStr))/(Lt*numStr-1);
            tmp_eb2=(cov_pcF_EB[2]-mu_F_EB[4]*mu_F_EB[5]/(Lt*numStr))/(Lt*numStr-1);
        EBcov_pcFx[2] += tmp_eb1*tmp_eb2;
    }
	
} /* ending i-loop (realizations) */


    /* normaliz be N_relz and 2ndMom->cov/var */
    for (k=0; k<6; k++) {
        mn_Xs[k] /= N_relz;
        mn_Fs[k] /= N_relz;
        vr_Xs[k] = vr_Xs[k]/(N_relz-1)-mn_Xs[k]*mn_Xs[k];
        vr_Fs[k] = vr_Fs[k]/(N_relz-1)-mn_Fs[k]*mn_Fs[k];
    }
    cov_obX[0]=cov_obX[0]/(N_relz-1)-mn_Xs[0]*mn_Xs[1];
    cov_obX[1]=cov_obX[1]/(N_relz-1)-mn_Xs[0]*mn_Xs[2];
    cov_obX[2]=cov_obX[2]/(N_relz-1)-mn_Xs[1]*mn_Xs[2];
    cov_obFx[0]=cov_obFx[0]/(N_relz-1)-mn_Fs[0]*mn_Fs[1];
    cov_obFx[1]=cov_obFx[1]/(N_relz-1)-mn_Fs[0]*mn_Fs[2];
    cov_obFx[2]=cov_obFx[2]/(N_relz-1)-mn_Fs[1]*mn_Fs[2];
    cov_pcX[0]=cov_pcX[0]/(N_relz-1)-mn_Xs[3]*mn_Xs[4];
    cov_pcX[1]=cov_pcX[1]/(N_relz-1)-mn_Xs[3]*mn_Xs[5];
    cov_pcX[2]=cov_pcX[2]/(N_relz-1)-mn_Xs[4]*mn_Xs[5];
    cov_pcFx[0]=cov_pcFx[0]/(N_relz-1)-mn_Fs[3]*mn_Fs[4];
    cov_pcFx[1]=cov_pcFx[1]/(N_relz-1)-mn_Fs[3]*mn_Fs[5];
    cov_pcFx[2]=cov_pcFx[2]/(N_relz-1)-mn_Fs[4]*mn_Fs[5];

    
    /* Unravel the std err Bars; have a total of nstd instances, measure std dev cf. whole sim */
    /* overwrite the EB[], which previously saved the 2nd moments in the i-loop*/
    for (k=0; k<6; k++) {
        EBmn_Xs[k] = ( EBmn_Xs[k]-nstd*mn_Xs[k]*mn_Xs[k] )/(nstd-1);
        EBmn_Fs[k] = ( EBmn_Fs[k]-nstd*mn_Fs[k]*mn_Fs[k] )/(nstd-1);
        EBvr_Xs[k] = ( EBvr_Xs[k]-nstd*vr_Xs[k]*vr_Xs[k] )/(nstd-1);
        EBvr_Fs[k] = ( EBvr_Fs[k]-nstd*vr_Fs[k]*vr_Fs[k] )/(nstd-1);
        if (k<3) {
            EBcov_obX[k] = ( EBcov_obX[k]-nstd*cov_obX[k]*cov_obX[k] )/(nstd-1);
            EBcov_obFx[k] = ( EBcov_obFx[k]-nstd*cov_obFx[k]*cov_obFx[k] )/(nstd-1);
        }
        else{
            EBcov_pcX[k-3] = ( EBcov_pcX[k-3]-nstd*cov_pcX[k-3]*cov_pcX[k-3] )/(nstd-1);
            EBcov_pcFx[k-3] = ( EBcov_pcFx[k-3]-nstd*cov_pcFx[k-3]*cov_pcFx[k-3] )/(nstd-1);
        }
    }
    //end of errorbar EB calculations
    
return;
                
}
