/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/* NMO: Optimized Version */

#include "suhdr.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" Stolt NMO - stolt NMO for migrating common scatterpoint gathers            ",
"									     ",
" susnmo3 <stdin >stdout [optional parameters]				     ",
"									     ",
" Required parameters:							     ",
" f=			filename of the nmo file (This is a su segy file)    ",
"									     ",
"									     ",
" The seismic traces and velocity traces are matched to each other	     ",
" using the cdp header word in both					     ",
"									     ",
" Optional Parameters:							     ",
" fac=1.0		a factor to multiply velocities			     ",
" w=0.6			Stolt stretch factor                                 ",
"									     ",
NULL};

/**************** end self doc *******************************************/

void getvel (int nt, int ncdp, float *cdp, 
	float **ovv, int cdpt, 
	float *ovvt);
void susnmo(segy **rec_o, float *vrms, int nt, float dt, int ntr, segy *outtrace, float w);
void makeut (float vstolt, float fmax, float *vrms, int nt, float dt, float **ut, int *nu, float *du, float **tu);
void makeu (float vstolt, float *v, int nt, float dt, float *u);
void stolt1k (float k, float v, float s, float fmax, int nt, float dt, complex *p, complex *q);

segy tr, trv, outtr;

int
main(int argc, char **argv)
{
	float ft;	/* time of first sample */
	int it;		/* time sample index */
	int ncdp;	/* number of cdps specified */
	float *cdp = NULL;	/* array[ncdp] of cdps */
	int icdp;	/* index into cdp array */
	int nvnmo;	/* number of vnmos specified */
	float *vnmo = NULL;	/* array[nvnmo] of vnmos */
	float *tnmo = NULL;	/* array[ntnmo] of tnmos */
	float **ovv = NULL;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
	float *ovvt = NULL;	/* array[nt] of sloth for a particular trace */
 	float v;
	float w;
	
        cwp_String key;         /* header key word from segy.h          */
        cwp_String type;        /* ... its type                         */
        Value val;
        segy **rec_o;           /* trace header+data matrix */
        int first=0;            /* true when we passed the first gather */
        float dt;
        int nt;
        int ntr;
	
	cwp_String f="";	/* su velocity file  */
	FILE *fp;		/* file pointer */
	float fac=1.0;
	
	int verbose;
	
	/* hook up getpar */
	initargs(argc, argv);
	requestdoc(1);

	MUSTGETPARSTRING("f",&f);
	if (!getparfloat("fac",&fac)) fac = 1.0;
	if (!getparfloat("w",&w)) w = 0.6;
        if (!getparstring("key", &key)) key = "cdp";
	if (!getparint("verbose",&verbose)) verbose = 1;
 	
	/* get information from the first header */
        rec_o = get_gather(&key, &type, &val, &nt, &ntr, &dt, &first);
	ft = (*rec_o[0]).delrt/1000.0;
	
	/* Read cdps */
	if((fp=fopen(f,"r"))==NULL)
                        err("cannot open file=%s\n",f);
	
	ncdp = fgettra(fp, &trv, 0);
	nvnmo = trv.ns;

	/* Allocate memory only once, outside the loop */
	ovv = ealloc2float(nt, ncdp);
	cdp = ealloc1float(ncdp);
	vnmo = ealloc1float(nvnmo);
	tnmo = ealloc1float(nvnmo);
	ovvt = ealloc1float(nt);  /* Allocate workspace outside the loop */
	
	icdp = 0;
	
	/* Precompute tnmo values outside the loop */
	for(it=0; it<nvnmo; it++)
		tnmo[it] = ft + it * (double)trv.dt/1000000.0;
	
	/* Read velocity data */
        do {
		cdp[icdp] = trv.cdp;
		for(it=0; it<nvnmo; it++) 
			vnmo[it] = fac * trv.data[it];
			
		float tn;
		for (it=0, tn=ft; it<nt; ++it, tn+=dt) {
			intlin(nvnmo, tnmo, vnmo, vnmo[0], vnmo[nvnmo-1], 1, &tn, &v);
			ovv[icdp][it] = (double)v;
		}
		icdp++;
	} while(fgettr(fp, &trv));

	fclose(fp);
        
	/* loop over gathers */
	do {
		if(verbose) fprintf(stderr," %d\n", vtoi(type, val));
		
		getvel(nt, ncdp, cdp, ovv, vtoi(type, val), ovvt);
		
		susnmo(rec_o, ovvt, nt, dt, ntr, &outtr, w);
		
		/* write output trace */
		puttr(rec_o[0]);
		
		/* Free previous gather */
		for(int i=0; i<ntr; i++)
			free1((void *)rec_o[i]);
		
		/* Get next gather */
		rec_o = get_gather(&key, &type, &val, &nt, &ntr, &dt, &first);
	} while (rec_o);

	/* Clean up memory */
	free1float(vnmo);
	free1float(tnmo);
	free1float(ovvt);
	free1float(cdp);
	free2float(ovv);

	return EXIT_SUCCESS;
}


/* match cdp in trace to cdp in velocity file */
void getvel (int nt, int ncdp, float *cdp, float **ovv, 
                      int cdpt, float *ovvt)
{
	int indx=0;

	/* if before first cdp, err */
	if (cdpt < (int)cdp[0]) 
		err(" trace CDP # %d is smaller than one in velocity file %d \n", cdpt, (int)cdp[0]);
	
	/* else if beyond last cdp, error*/
	if (cdpt > (int)cdp[ncdp-1]) 
		err(" trace CDP # %d is larger than largets in velocity file %d ", cdpt, (int)cdp[ncdp-1]);
	
	/* else, copy */
	else {
		xindex(ncdp, cdp, cdpt, &indx);
		if((int)(cdp[indx]) != cdpt) 
			err(" CDP # %d not in velocity file exiting \n", cdpt);
		memcpy((void *) ovvt, (const void *) ovv[indx], nt*FSIZE);
	}
}

void susnmo(segy **rec_o, float *vrms, int nt, float dt, int ntr, segy *outtrace, float w)
{
	float **px = NULL, *ut = NULL, *tu = NULL;
	complex **pk = NULL;
	float du, vstolt;
	int itr, ntr2;
	float nxpad;
	int nxfft;
	int nk, ik, iu, nu = 0;
	float dx, dk, scale, k;
	int ntm;
	float offmax, offmin;
	int ntr0;	/* number of traces in the gather padded up to zero offset */
	int ftr;	/* location of the first trace of the csp gather in px array */
	float fmax;
	float vmin=99999.9, vmax=-99999.9;
	
	/* Find velocity min/max in a single pass */
	for(itr=0; itr<nt; itr++) {
		if(vmin > vrms[itr]) 
			vmin = vrms[itr];
		if(vmax < vrms[itr]) 
			vmax = vrms[itr];
	}
	
	vstolt = vmin;	
	fmax = 0.5/dt;
	offmax = (*rec_o[ntr-1]).offset;
	offmin = (*rec_o[0]).offset;
	dx = (offmax-offmin)/(ntr-1);
	ntr0 = NINT(offmax/dx);
	ftr = NINT(offmin/dx);
	ntr2 = 2*ntr0;
	
	/* make u(t) and t(u) for Stolt stretch */
	makeut(vstolt, fmax, vrms, nt, dt, &ut, &nu, &du, &tu);
	
	/* wavenumber (k) sampling */
	nxpad = 0.5*vmax*nt*dt/dx;
	nxfft = npfar(ntr2+nxpad);
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);
	scale = 1.0/nxfft;
	
	/* Allocate gathers for the migration */
	ntm = MAX(nu, nt);
	px = ealloc2float(ntm, nxfft);
	pk = ealloc2complex(ntm, nk);
	
	/* setup the gather for migration */
	float *temp = NULL;
	
	/* If we need to stretch, allocate temp only once */
	if (nu != nt && tu != NULL) {
		temp = ealloc1float(ntm);
	}
	
	for(itr=0; itr<ntr; itr++) {
		/* if necessary, stretch */
		if (nu != nt && tu != NULL) {
			ints8r(nt, dt, 0.0, (*rec_o[itr]).data, 0.0, 0.0, nu, tu, temp);
			memcpy((void *) px[ntr-1-itr], (const void *) temp, nu*FSIZE);
			memcpy((void *) px[ntr+itr], (const void *) temp, nu*FSIZE);
		} else {
			memcpy((void *) px[ntr-1-itr], (const void *) (*rec_o[itr]).data, nt*FSIZE);
			memcpy((void *) px[ntr+itr], (const void *) (*rec_o[itr]).data, nt*FSIZE);
		}
	}
	
	/* Free temp if allocated */
	if (temp) free1float(temp);
	
	/* Zero from ntr2 to nxfft */
	for(itr=ntr2; itr<nxfft; itr++) 
		memset((void *) px[itr], (int) '\0', nt*FSIZE);	

	/* Fourier transform p(u,x) to p(u,k) */
	pfa2rc(-1, 2, ntm, nxfft, px[0], pk[0]);

	/* migrate each wavenumber */
	for (ik=1, k=dk; ik<nk; ++ik, k+=dk)
		stolt1k(k, vstolt, w, fmax, nu, du, pk[ik], pk[ik]);

	/* Fourier transform p(u,k) to p(u,x) and scale */
	pfa2cr(1, 2, ntm, nxfft, pk[0], px[0]);
	
	for (itr=0; itr<ntr2; ++itr)
		for (iu=0; iu<nu; ++iu)
			px[itr][iu] *= scale;

	/* if necessary, unstretch */
	if (nu != nt && ut != NULL) {
		temp = ealloc1float(ntm);
		for (itr=0; itr<ntr2; ++itr) {
			ints8r(nu, du, 0.0, px[itr], 0.0, 0.0, nt, ut, temp);
			memcpy((void *) px[itr], (const void *) temp, nu*FSIZE);
		}
		free1float(temp);
	}
			
	/* Copy result to output */
	memcpy((void *) (*rec_o[0]).data, (const void *) px[ntr-1], nt*FSIZE);
	
	/* Clean up memory */
	if (ut) free1float(ut);
	if (tu) free1float(tu);
	free2float(px);
	free2complex(pk);
}

void makeut(float vstolt, float fmax, float *vrms, int nt, float dt, float **ut, int *nu, float *du, float **tu)
{
	int it;
	
	/* check for constant velocity */
	for (it=1; it<nt; ++it)
		if (vrms[it] != vrms[0]) break;
		
	/* if constant velocity */
	if (it==nt) {
		*ut = NULL;
		*tu = NULL;
		*nu = nt;
		*du = dt;

	/* else if velocity not constant */
	} else {
		int it, nuu;
		float duu, delu, umax, *u, *t;
		
		/* u(t) */
		u = alloc1float(nt);
		makeu(vstolt, vrms, nt, dt, u);

		/* smallest du and maximum u */
		duu = FLT_MAX;
		for (it=1; it<nt; ++it) {
			delu = u[it]-u[it-1];
			if (delu < duu) duu = delu;
		}
		umax = u[nt-1];

		/* u sampling */
		duu = duu/(2.0*fmax*dt);
		nuu = 1+NINT(umax/duu);

		/* t(u) */
		t = alloc1float(nuu);
		yxtoxy(nt, dt, 0.0, u, nuu, duu, 0.0, 0.0, (nt-1)*dt, t);

		/* set output parameters before returning */
		*ut = u;
		*tu = t;
		*nu = nuu;
		*du = duu;
	}
}

void makeu(float vstolt, float *v, int nt, float dt, float *u)
{
	int it;
	float t, scale, sum;

	scale = 2.0/(vstolt*vstolt);
	u[0] = sum = 0.0;
	for (it=1, t=dt; it<nt; ++it, t+=dt) {
		sum += 0.5*dt*(t*v[it]*v[it]+(t-dt)*v[it-1]*v[it-1]);
		u[it] = sqrt(scale*sum);
	}
}

void stolt1k(float k, float v, float s, float fmax, int nt, float dt, complex *p, complex *q)
{
	int nw, it, nwtau, iwtau, ntau, itau, iwtaul, iwtauh;
	float vko2s, wmax, dw, fw, dwtau, fwtau, wtau, dtau,
		wtauh, wtaul, scale, fftscl, a, b, *wwtau;
	complex czero=cmplx(0.0, 0.0), *pp, *qq;

	/* modify stolt stretch factor to simplify calculations below */
	if (s!=1.0) s = 2.0-s;

	/* (v*k/2)^2 */
	vko2s = 0.25*v*v*k*k;

	/* maximum frequency to migrate in radians per unit time */
	wmax = 2.0*PI*MIN(fmax, 0.5/dt);

	/* frequency sampling - must pad to avoid interpolation error;
	 * pad by factor of 2 because time axis is not centered;
	 * pad by factor of 1/0.6 because 8-point sinc is valid
	 * only to about 0.6 Nyquist
	 */
	nw = nt*2/0.6;
	nw = npfao(nw, nw*2);
	dw = 2.0*PI/(nw*dt);
	fw = -PI/dt;

	/* migrated time */
	ntau = nt;
	dtau = dt;

	/* migrated frequency - no need to pad since no interpolation */
	nwtau = npfao(ntau, ntau*2);
	dwtau = 2.0*PI/(nwtau*dtau);
	fwtau = -PI/dtau;

	/* tweak first migrated frequency to avoid wtau==0.0 below */
	fwtau += 0.001*dwtau;

	/* high and low migrated frequencies - don't migrate evanescent */
	wtauh = sqrt(MAX(0.0, wmax*wmax-s*vko2s));
	iwtauh = MAX(0, MIN(nwtau-1, NINT((wtauh-fwtau)/dwtau)));
	iwtaul = MAX(0, MIN(nwtau-1, NINT((-wtauh-fwtau)/dwtau)));
	wtauh = fwtau+iwtauh*dwtau;
	wtaul = fwtau+iwtaul*dwtau;

	/* workspace */
	pp = alloc1complex(nw);
	qq = alloc1complex(nwtau);
	wwtau = alloc1float(nwtau);

	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it+=2)
		pp[it] = p[it];
	for (it=1; it<nt; it+=2) {
		pp[it].r = -p[it].r;
		pp[it].i = -p[it].i;
	}
	for (it=nt; it<nw; ++it)
		pp[it].r = pp[it].i = 0.0;
	pfacc(1, nw, pp);

	/* zero -Nyquist frequency for symmetry */
	pp[0] = czero;

	/* frequencies at which to interpolate */
	if (s==1.0) {
		for (iwtau=iwtaul, wtau=wtaul; wtau<0.0; ++iwtau, wtau+=dwtau)
			wwtau[iwtau] = -sqrt(wtau*wtau+vko2s);
		for (; iwtau<=iwtauh; ++iwtau, wtau+=dwtau)
			wwtau[iwtau] = sqrt(wtau*wtau+vko2s);
	} else {
		a = 1.0/s;
		b = 1.0-a;
		for (iwtau=iwtaul, wtau=wtaul; wtau<0.0; ++iwtau, wtau+=dwtau)
			wwtau[iwtau] = b*wtau-a*sqrt(wtau*wtau+s*vko2s);
		for (; iwtau<=iwtauh; ++iwtau, wtau+=dwtau)
			wwtau[iwtau] = b*wtau+a*sqrt(wtau*wtau+s*vko2s);
	}
	
	/* interpolate */
	ints8c(nw, dw, fw, pp, czero, czero,
		iwtauh-iwtaul+1, wwtau+iwtaul, qq+iwtaul);

	/* fft scaling and obliquity factor */
	fftscl = 1.0/nwtau;
	if (s==1.0) {
		for (iwtau=iwtaul, wtau=wtaul; iwtau<=iwtauh; 
			++iwtau, wtau+=dwtau) {
			scale = fftscl*wtau/wwtau[iwtau];
			qq[iwtau].r *= scale;
			qq[iwtau].i *= scale;
		}
	} else {
		a = 1.0/(s*s);
		b = 1.0-1.0/s;
		for (iwtau=iwtaul, wtau=wtaul; iwtau<=iwtauh; 
			++iwtau, wtau+=dwtau) {
			scale = fftscl*(b+a*wtau/(wwtau[iwtau]-b*wtau));
			qq[iwtau].r *= scale;
			qq[iwtau].i *= scale;
		}
	}

	/* zero evanescent frequencies */
	for (iwtau=0; iwtau<iwtaul; ++iwtau)
		qq[iwtau] = czero;
	for (iwtau=iwtauh+1; iwtau<nwtau; ++iwtau)
		qq[iwtau] = czero;

	/* Fourier transform wtau to tau, accounting for centered wtau */
	pfacc(-1, nwtau, qq);
	for (itau=0; itau<ntau; itau+=2)
		q[itau] = qq[itau];
	for (itau=1; itau<ntau; itau+=2) {
		q[itau].r = -qq[itau].r;
		q[itau].i = -qq[itau].i;
	}
	
	/* free workspace */
	free1complex(pp);
	free1complex(qq);
	free1float(wwtau);
}