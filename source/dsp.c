#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "dsp.h"
#include "funcs.h"
#include "fsl_powerquad.h"


extern uint32_t timerCounter; 


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


static real random_real(real min, real max) {
    return min + ((real)rand() / RAND_MAX) * (max - min);
}

/* accelerometer */
header* maccel(Calc *cc, header *hd)
{

    // Define the default dimensions
    int rows = 1;
    int cols = 3; // We want a matrix of 1 row and 3 columns for x, y, and z

    header *result = new_matrix(cc, rows, cols, "");

    // Fill the matrix with random values for x, y, and z
    real *m = matrixof(result);
    
	int32_t data[3];
	mma8652_read_xyz(data);
	m[0] = (real)data[0]/1000.0;
	m[1] = (real)data[1]/1000.0;
	m[2] = (real)data[2]/1000.0;

    return pushresults(cc, result);
}


header* mpqcos(Calc* cc, header* hd) {
    int r, c; // rows and columns
    real *m, *mr; // pointers to input and output matrices
    header *result; // Output header
    // PRINTF("$s\r\n", __func__ ):

	if (hd->type == s_reference) { // if the variable is a reference 
        hd = getvalue(cc, hd);
    }

    // Check input type
    if (hd->type != s_real && hd->type != s_matrix) {
        cc_error(cc, "Invalid input type for mpqcos");
        return NULL;
    }
    
    // Get the input matrix
    getmatrix(hd, &r, &c, &m);
    
    // Allocate a new output matrix (or scalar)
    if (hd->type == s_real) {
        result = new_real(cc, 0, "");
        mr = realof(result);
        PQ_CosF32(&m[0], &mr[0]); // Update to use the function correctly
    } else {
        result = new_matrix(cc, r, c, "");
        mr = matrixof(result);
        PQ_VectorCosF32(m, mr, r * c); // Using vectorized function for matrices
    }

    // Push result to stack and return
    return pushresults(cc, result);
}

void real_cos(real *in, real *out) {
    PQ_CosF32(in, out);
}

header* mpqcos2(Calc* cc, header* hd) {
	if (hd->type == s_reference) { // if the variable is a reference 
        hd = getvalue(cc, hd);
    }
	
    if (hd->type != s_real && hd->type != s_matrix) {
        cc_error(cc, "Invalid input type for mpqcos");
        return NULL;
    }

    // Use the map1 function to handle real or matrix inputs.
    // real_cos is passed as the function to be applied elementwise.
    return map1(cc, real_cos, NULL, hd);
}

static int32_t fft_out[1024]; // 512 points
static int32_t fft_in[1024];

header* mpqfft (Calc *cc, header *hd)
{
    int r, c;
    real *m, *mr;
    header *result;



	if (hd->type == s_reference) { // if the variable is a reference 
        hd = getvalue(cc, hd);
    }

    // Vérifier le type d'entrée
    if (hd->type != s_real && hd->type != s_complex && hd->type != s_matrix) {
        cc_error(cc, "Type d'entrée non valide");
        return NULL;
    }

    getmatrix(hd, &r, &c, &m);
    if (r != 1) cc_error(cc, "Vecteur ligne attendu");

    // Allouer un nouveau vecteur pour le résultat
    result = new_cmatrix(cc, 1, c, "");
    mr = matrixof(result);

	PQ_Init(POWERQUAD);

	// Configurer les régions de mémoire pour PowerQuad
	pq_config_t pgConfig;
	pgConfig.inputAFormat = kPQ_16Bit; // Format fixed-point for input A
	pgConfig.inputAPrescale = 9; // Prescale for input A
	pgConfig.tmpFormat = kPQ_32Bit; 
	pgConfig.tmpPrescale = 0; 
	pgConfig.outputFormat = kPQ_32Bit; 
	pgConfig.outputPrescale = 0; 
	pgConfig.machineFormat = kPQ_32Bit;

	PQ_SetConfig(POWERQUAD, &pgConfig);

	// Convert floating point input to fixed point
	
	for (int i = 0; i < c; i++) {
		m[i] = (int16_t)(m[i] * (1 << 9)); // Shift by 5 bits to the left
		fft_in[i] = (int32_t)m[i];
	}
    // // Appel à la fonction PowerQuad pour le calcul FFT
	// for (int i=0;i<r;i++) {
	//     PQ_TransformRFFT(POWERQUAD,c,&m[i*c],&mr[i*(c/2+1)]);
	// }

    PQ_TransformRFFT(POWERQUAD, 512, fft_in, fft_out); 
	for (int i = 0; i < c; i++) {
		mr[i] = fft_out[i];
	}
	PQ_WaitDone(POWERQUAD);
    // Renvoyer le résultat
    return pushresults(cc, result);
}

// Exemple de mpqifft, le cas IFFT avec PowerQuad
header* mpqifft (Calc* cc, header* hd)
{
    int r, c;
    real *m, *mr;
    header *result;

	if (hd->type == s_reference) { // if the variable is a reference 
        hd = getvalue(cc, hd);
    }
    // Vérifier le type d'entrée
    if (hd->type != s_real && hd->type != s_complex && hd->type != s_matrix) {
        cc_error(cc, "Type d'entrée non valide pour mpqifft");
        return NULL;
    }

    getmatrix(hd, &r, &c, &m);
    if (r != 1) cc_error(cc, "Vecteur ligne attendu");

    // Allouer un nouveau vecteur pour le résultat
    result = new_cmatrix(cc, 1, c, "");
    mr = matrixof(result);

    // Appel à la fonction PowerQuad pour le calcul IFFT
    PQ_Init(POWERQUAD);

	// Configurer les régions de mémoire pour PowerQuad
	pq_config_t pgConfig;
	pgConfig.inputAFormat = kPQ_16Bit; // Format fixed-point for input A
	pgConfig.inputAPrescale = 9; // Prescale for input A
	pgConfig.tmpFormat = kPQ_32Bit; 
	pgConfig.tmpPrescale = 0; 
	pgConfig.outputFormat = kPQ_32Bit; 
	pgConfig.outputPrescale = 0; 
	pgConfig.machineFormat = kPQ_32Bit;

	PQ_SetConfig(POWERQUAD, &pgConfig);

	// Convert floating point input to fixed point
	
	for (int i = 0; i < c; i++) {
		m[i] = (int16_t)(m[i] * (1 << 9)); // Shift by 5 bits to the left
		fft_in[i] = (int32_t)m[i];
	}
    // // Appel à la fonction PowerQuad pour le calcul FFT
	// for (int i=0;i<r;i++) {
	//     PQ_TransformRFFT(POWERQUAD,c,&m[i*c],&mr[i*(c/2+1)]);
	// }

    PQ_TransformIFFT(POWERQUAD, 512, fft_in, fft_out); 
	for (int i = 0; i < c; i++) {
		mr[i] = fft_out[i];
	}

	PQ_WaitDone(POWERQUAD);

    // Renvoyer le résultat
    return pushresults(cc, result);
}

