#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "dsp.h"
#include "funcs.h"

/****************************************************************
 *	FFT
 ****************************************************************/
int nn;
cplx *ff,*zz,*fh;

static void rfft (LONG m0, LONG p0, LONG q0, LONG n)
/***** rfft 
	make a fft on x[m],x[m+q0],...,x[m+(p0-1)*q0] (p points).
	one has xi_p0 = xi_n^n = zz[n] ; i.e., p0*n=nn.
*****/
{	LONG p,q,m,l;
	LONG mh,ml;
	int found=0;
	cplx sum,h;
	if (p0==1) return;
//	if (test_key()==escape) cc_error("fft interrupted")
	if (p0%2==0) {
		p=p0/2; q=2;
	} else {
		q=3;
		while (q*q<=p0) {
			if (p0%q==0) {
				found=1; break;
			}
			q+=2;
		}
		if (found) p=p0/q;
		else { q=p0; p=1; }
	}
	if (p>1) for (m=0; m<q; m++) 
		rfft((m0+m*q0)%nn,p,q*q0,nn/p);
	mh=m0;
	for (l=0; l<p0; l++) {
		ml=l%p;
		c_copy(ff[(m0+ml*q*q0)%nn],sum);
		for (m=1; m<q; m++) {
			c_mul(ff[(m0+(m+ml*q)*q0)%nn],zz[(n*l*m)%nn],h);
			c_add(sum,h,sum);
		}
		sum[0]/=q; sum[1]/=q;
		c_copy(sum,fh[mh]);
		mh+=q0; if (mh>=nn) mh-=nn;
	}
	for (l=0; l<p0; l++) {
		c_copy(fh[m0],ff[m0]);
		m0+=q0; if (m0>=nn) mh-=nn;
	}
}

static void fft (Calc *cc, real *a, int n, int signum)
{	cplx z;
	real h;
	int i;
	
	if (n==0) return;
	nn=n;

	char *ram=cc->newram;
	ff=(cplx *)a;
	zz=(cplx *)ram;
	ram+=n*sizeof(cplx);
	fh=(cplx *)ram;
	ram+=n*sizeof(cplx);
	if (ram>cc->udfstart)  cc_error(cc,"Memory overflow!");
	
	/* compute zz[k]=e^{-2*pi*i*k/n}, k=0,...,n-1 */
	h=2*M_PI/n; z[0]=cos(h); z[1]=signum*sin(h);
	zz[0][0]=1; zz[0][1]=0;
	for (i=1; i<n; i++) {
		if (i%16) { zz[i][0]=cos(i*h); zz[i][1]=signum*sin(i*h); }
		else c_mul(zz[i-1],z,zz[i]);
	}
	rfft(0,n,1,1);
	if (signum==1)
		for (i=0; i<n; i++) {
			ff[i][0]*=n; ff[i][1]*=n;
		}
}

header* mfft (Calc *cc, header *hd)
{	header *st=hd,*result;
	real *m,*mr;
	int r,c;
	hd=getvalue(cc,hd);
	if (hd->type==s_real || hd->type==s_matrix)	{
		make_complex(cc,st); hd=st;
	}
	getmatrix(hd,&r,&c,&m);
	if (r!=1) cc_error(cc,"row vector expected");
	result=new_cmatrix(cc,1,c,"");
	mr=matrixof(result);
    memmove((char *)mr,(char *)m,(ULONG)2*c*sizeof(real));
	fft(cc,mr,c,-1);
	return pushresults(cc,result);
}

header* mifft (Calc *cc, header *hd)
{	header *st=hd,*result;
	real *m,*mr;
	int r,c;
	hd=getvalue(cc,hd);
	if (hd->type==s_real || hd->type==s_matrix) {
		make_complex(cc,st); hd=st;
	}
	getmatrix(hd,&r,&c,&m);
	if (r!=1) cc_error(cc,"row vector expected");
	result=new_cmatrix(cc,1,c,"");
	mr=matrixof(result);
    memmove((char *)mr,(char *)m,(ULONG)2*c*sizeof(real));
	fft(cc,mr,c,1);
	return pushresults(cc,result);
}

/* accelerometer */
header* maccel (Calc *cc, header *hd)
{
	return NULL;
}

header* mpqcos(Calc* cc, header* hd) {
	return NULL;
}

static int32_t fft_out[1024]; // 512 points
static int32_t fft_in[1024];

header* mpqfft (Calc *cc, header *hd)
{
	return NULL;
}

header* mpqifft (Calc* cc, header* hd)
{
	return NULL;
}

