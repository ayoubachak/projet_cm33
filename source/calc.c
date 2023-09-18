/****************************************************************
 * calc
 *  (C) 1993-2021 R. Grothmann
 *  (C) 2021-2022 E. Bouchare
 *
 * calc.c
 *
 ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "sysdep.h"
#include "calc.h"
#include "stack.h"
#include "edit.h"
#include "funcs.h"
#include "graphics.h"
#include "dsp.h"

/***************************************************************************
 * dumping output to the screen or to a file 
 ***************************************************************************/
void output (Calc *cc, char *s)
{	text_mode();
	if (CC_ISSET(cc,CC_OUTPUTING)) sys_print(s);
	if (cc->outfile) {
		fprintf(cc->outfile,"%s",s);
		if (ferror(cc->outfile)) {	
			fclose(cc->outfile); cc->outfile=NULL;
			sys_out_mode(CC_ERROR);
			output(cc,"Error on dump file (disk full?).\n");
		}
	}
}

void outputf (Calc *cc, char *fmt, ...)
{	char text [256];
	va_list v;
	text_mode();
	va_start(v,fmt);
	vsnprintf(text,256,fmt,v);
	va_end(v);
	if (CC_ISSET(cc,CC_OUTPUTING)) sys_print(text);
	if (cc->outfile) {
		fprintf(cc->outfile,"%s",text);
		if (ferror(cc->outfile)) {	
			fclose(cc->outfile); cc->outfile=NULL;
			sys_out_mode(CC_ERROR);
			output(cc,"Error on dump file (disk full?).\n");
		}
	}
}

int outputfhold (Calc *cc, int f, char *fmt, ...)
{	static char text [256];
	int len;
	va_list v;
	if (f==0) text[0]=0;
	len=strlen(text);
	text_mode();
	va_start(v,fmt);
	vsnprintf(text+len,256-len,fmt,v);
	va_end(v);
	if (f<=0) return strlen(text);
	len=strlen(text);
	if (len<f) {
		memmove(text+(f-len),text,len+1);
		memset(text,' ',f-len);
	}
	if (CC_ISSET(cc,CC_OUTPUTING)) sys_print(text);
	if (cc->outfile) {
		fprintf(cc->outfile,"%s",text);
		if (ferror(cc->outfile)) {
			fclose(cc->outfile); cc->outfile=NULL;
			sys_out_mode(CC_ERROR);
			output(cc,"Error on dump file (disk full?).\n");
		}
	}
	return strlen(text);
}

#define PREFIX_START (-24)
static const char *eng_prefix[] = {
  "y", "z", "a", "f", "p", "n", "u", "m", "",
  "k", "M", "G", "T", "P", "E", "Z", "Y"
}; 
#define PREFIX_END (PREFIX_START+\
(int)((sizeof(eng_prefix)/sizeof(char *)-1)*3))

static void eng_out(Calc *cc, real value, int digits, int numeric, int hold)
{
	int expof10;
	int is_signed = signbit(value);
	char* sign = is_signed ? "-" : "";

	if (is_signed) value = -value;
	
	switch(fpclassify(value)) {
	case FP_NORMAL:
		expof10 = (int) log10(value);
		if(expof10 > 0)
			expof10 = (expof10/3)*3;
		else
			expof10 = (-expof10+3)/3*(-3); 
	 
		value *= pow(10,-expof10);
	
		if (value >= 1000.) { value /= 1000.0; expof10 += 3; }
	
		if(numeric || (expof10 < PREFIX_START) || (expof10 > PREFIX_END))
			outputfhold(cc, hold, "%s%.*GE%+.2d", sign, digits, value, expof10); 
		else
			outputfhold(cc, hold, "%s%.*G%s", sign, digits, value, 
	          eng_prefix[(expof10-PREFIX_START)/3]);
		break;
	case FP_INFINITE:
		outputfhold(cc, hold, "%sINF", sign);
		break;
	case FP_NAN:
		outputfhold(cc, hold, "%sNAN", sign);
		break;
	case FP_SUBNORMAL:
	case FP_ZERO:
	default:
		if(numeric) {
			outputfhold(cc, hold, "%s%.*GE+00", sign, digits, 0.0);
		} else {
			outputfhold(cc, hold, "%s%.*G", sign, digits, 0.0);
		}
		break;
	}
}

/* now part of the global structure
int disp_mode=0;
int disp_digits=6;
int disp_fieldw=14;
int disp_eng_sym=0;

real maxexpo=1.0e6,minexpo=1.0e-4;
char expoformat[16]="%0.4E";
char fixedformat[16]="%0.5G";
*/

void real_out (Calc *cc, real x)
/***** real_out
	print a real number.
*****/
{
	switch (cc->disp_mode) {
	case 0:		/* smart STD */
		if ((fabs(x)>cc->maxexpo || fabs(x)<cc->minexpo) && x!=0.0) 
			outputfhold(cc,0,cc->expoformat,x);
		else outputfhold(cc,0,cc->fixedformat,x);
		break;
	case 1:		/* ENG1 */
	case 2:		/* ENG2 */
		eng_out(cc,x,cc->disp_digits,!cc->disp_eng_sym,0);
		break;
	case 3:		/* SCI */
		outputfhold(cc,0,cc->expoformat,x);
		break;
	case 4:		/* FIXED */
		outputfhold(cc,0,cc->fixedformat,x);
		break;
	case 5:		/* FRAC */
		break;
	default:	/* never used */
		break;
	}
	outputfhold(cc,cc->disp_fieldw,"");
}

void out_matrix (Calc *cc, header *hd)
/***** out_matrix
   print a matrix.
*****/
{	int c,r,i,j,c0,cend;
	real *m,*x;
	
	int linew=cc->termwidth/cc->disp_fieldw;

	getmatrix(hd,&r,&c,&m);
	for (c0=0; c0<c; c0+=linew) {
		cend=c0+linew-1; 
		if (cend>=c) cend=c-1;
		if (c>linew) outputf(cc,"Column %d to %d:\n",c0+1,cend+1);
		for (i=0; i<r; i++) {
			x=mat(m,c,i,c0);
			for (j=c0; j<=cend; j++) real_out(cc,*x++);
			output(cc,"\n");
			if (sys_test_key()==escape) return;
		}
	}
}

void complex_out (Calc *cc, real x, real y)
/***** complex_out
	print a complex number.
*****/
{
	real m=sqrt(x*x+y*y);
	x=fabs(x)/m>EPSILON ? x : 0.0;
	y=fabs(y)/m>EPSILON ? y : 0.0;
	switch (cc->disp_mode) {
	case 0:		/* smart STD */
		if ((fabs(x)>cc->maxexpo || fabs(x)<cc->minexpo) && x!=0.0) 
			outputfhold(cc,0,cc->expoformat,x);
		else outputfhold(cc,0,cc->fixedformat,x);
		if (y>=0) outputfhold(cc,-1,"+");
		else outputfhold(cc,-1,"-");
		y=fabs(y);
		if ((y>cc->maxexpo || y<cc->minexpo) && y!=0.0)
			outputfhold(cc,-1,cc->expoformat,y);
		else outputfhold(cc,-1,cc->fixedformat,y);
		break;
	case 1:		/* ENG1 */
	case 2:		/* ENG2 */
		eng_out(cc,x,cc->disp_digits,!cc->disp_eng_sym,0);
		if (y>=0) outputfhold(cc,-1,"+");
		else outputfhold(cc,-1,"-");
		y=fabs(y);
		eng_out(cc,y,cc->disp_digits,!cc->disp_eng_sym,-1);
		break;
	case 3:		/* SCI */
		outputfhold(cc,0,cc->expoformat,x);
		if (y>=0) outputfhold(cc,-1,"+");
		else outputfhold(cc,-1,"-");
		y=fabs(y);
		outputfhold(cc,-1,cc->expoformat,y);
		break;
	case 4:		/* FIXED */
		outputfhold(cc,0,cc->fixedformat,x);
		if (y>=0) outputfhold(cc,-1,"+");
		else outputfhold(cc,-1,"-");
		y=fabs(y);
		outputfhold(cc,-1,cc->fixedformat,y);
		break;
	case 5:		/* FRAC */
		break;
	default:	/* never used */
		break;
	}
	outputfhold(cc,2*cc->disp_fieldw,"i");
}

void out_cmatrix (Calc *cc, header *hd)
/***** out_matrix
   print a complex matrix.
*****/
{	int c,r,i,j,c0,cend;
	real *m,*x;

	int linew=cc->termwidth/(2*cc->disp_fieldw);

	getmatrix(hd,&r,&c,&m);
	for (c0=0; c0<c; c0+=linew) {
		cend=c0+linew-1; 
		if (cend>=c) cend=c-1;
		if (c>linew) outputf(cc,"Column %d to %d:\n",c0+1,cend+1);
		for (i=0; i<r; i++) {
			x=cmat(m,c,i,c0);
			for (j=c0; j<=cend; j++) {
				complex_out(cc,*x,*(x+1)); 
				x+=2;
			}
			output(cc,"\n");
			if (sys_test_key()==escape) return;
		}
	}
}

void give_out (Calc *cc, header *hd)
/***** give_out
	print a value.
*****/
{
	sys_out_mode(CC_OUTPUT);
	if (hd->type==s_reference) hd=getvalue(cc,hd);
	switch(hd->type) {
		case s_real: real_out(cc,*realof(hd)); output(cc,"\n"); break;
		case s_complex: complex_out(cc,*realof(hd),*(realof(hd)+1));
			output(cc,"\n"); break;
		case s_matrix: out_matrix(cc,hd); break;
		case s_cmatrix: out_cmatrix(cc,hd); break;
		case s_string: output(cc,stringof(hd)); output(cc,"\n"); break;
		case s_funcref: {
			char* name=hd->name;
			if (hd->flags & FLAG_BINFUNC) {
				while (hd && hd->type==s_funcref) hd=referenceof(hd);
				outputf(cc,"binary function %s (%p)\n",name,(char*)hd);
			} else {
				while (hd && hd->type==s_funcref) hd=referenceof(hd);
				outputf(cc,"user function %s (0x%08x)\n",name,(char*)hd-cc->udfstart);
			}
			break;
		}
		default : output(cc,"?\n");
	}
}

/***************************************************************************
 *	binfunc function core
 ***************************************************************************/
#ifndef EMBED
static binfunc_t binfunc_list[] = {
	{"time",0,mtime},
	{"wait",1,mwait},
	
	{"index",0,mindex},
	{"argn",0,margn},
	{"args",1,margs},
	{"xargs",0,mxargs},

	{"format",2,mformat},
	{"printf",2,mprintf},

	{"name",1,mname},
	{"isreal",1,misreal},
	{"iscomplex",1,miscomplex},
	{"isstring",1,misstring},
	{"isfunction",1,misfunction},
	{"isvar",1,misvar},

	{"epsilon",0,mepsilon},
	{"epsilon",1,msetepsilon},

	{"input",1,minput},
	{"error",1,merror},

	{"complex",1,mcomplex},
	{"re",1,mre},
	{"im",1,mim},
	{"abs",1,mabs},
	{"arg",1,marg},
	{"conj",1,mconj},

	{"sin",1,msin},
	{"cos",1,mcos},
	{"tan",1,mtan},
	{"atan",1,matan},
	{"acos",1,macos},
	{"asin",1,masin},
	{"exp",1,mexp},
	{"log",1,mlog},
	{"sqrt",1,msqrt},
	{"mod",2,mmod},
	{"sign",1,msign},
	{"floor",1,mfloor},
	{"ceil",1,mceil},
	{"round",2,mround},
	
	{"fac",1,mfac},
	{"bin",2,mbin},
	{"logfac",1,mlogfac},
	{"logbin",2,mlogbin},

	{"size",-1,msize},
	{"rows",1,mrows},
	{"cols",1,mcols},
	{"extrema",1,mextrema},
	{"any",1,many},
	{"nonzeros",1,mnonzeros},

	{"zeros",1,mzeros},
	{"ones",1,mones},
	{"matrix",2,mmatrix},
	{"diag",3,mdiag},
	{"diag",2,mdiag2},
	{"band",3,mband},
	{"setdiag",3,msetdiag},

	{"bandmult",2,wmultiply},
	{"symmult",2,smultiply},
	
	{"dup",2,mdup},
	{"redim",2,mredim},
	{"sum",1,msum},
	{"prod",1,mprod},
	{"colsum",1,mcolsum},
	{"cumsum",1,mcumsum},
	{"cumprod",1,mcumprod},
	{"flipx",1,mflipx},
	{"flipy",1,mflipy},

	{"max",1,mmax1},
	{"min",1,mmin1},
	{"max",2,mmax},
	{"min",2,mmin},
	{"sort",1,msort},
	{"count",2,mstatistics},
	
	{"polyval",2,mpolyval},
	{"polyadd",2,mpolyadd},
	{"polymult",2,mpolymult},
	{"polydiv",2,mpolydiv},
	{"polytrunc",1,mpolytrunc},
	{"polycons",1,mpolycons},
	{"polysolve",1,mpolysolve},
	{"polyroot",2,mpolyroot},
	
	{"interp",2,dd},
	{"interpval",3,ddval},
	{"polytrans",2,polydd},
	
	{"hb",1,mtridiag},
	{"charpoly",1,mcharpoly},
	{"lu",1,mlu},
	{"lusolve",2,mlusolve},
	
	{"fft",1,mfft},
	{"ifft",1,mifft},
	
	{"seed",1,mseed},
	{"random",1,mrandom},
	{"shuffle",1,mshuffle},
	{"normal",1,mnormal},
	{"find",2,mfind},
	
	{"subplot",1,msubplot},
	{"setplot",1,msetplot},
	{"plotarea",2,mplotarea},
	{"plotstyle",1,mplotstyle},
	{"plot",3,mplot},
	{"plot",2,mplot1},
	{"xgrid",5,mxgrid},
	{"ygrid",5,mygrid},
	{"text",5,mtext},
	{"title",1,mtitle},
	{"xlabel",1,mxlabel},
	{"ylabel",1,mylabel},
	
	{"rmfir",6,mrmfir},

#ifdef PCM_ENABLE
	{"pcmvol",1,mpcmvol},
	{"pcmfreq",0,mpcmfreq},
	{"pcmplay",1,mpcmplay},
	{"pcmrec",1,mpcmrec},
	{"pcmloop",0,mpcmloop},
#endif

	{"pqcos",1,mpqcos},
	{"pqfft",1,mpqfft},
	{"pqifft",1,mpqifft},
};
#else
static const binfunc_t binfunc_list[] = {
	{"abs",1,mabs},
	{"accel",0,maccel},
	{"acos",1,macos},
	{"any",1,many},
	{"arg",1,marg},
	{"argn",0,margn},
	{"args",1,margs},
	{"asin",1,masin},
	{"atan",1,matan},
	{"band",3,mband},
	{"bandmult",2,wmultiply},
	{"bin",2,mbin},
#ifdef PCM_ENABLE
	{"captaudio",2,mcaptaudio},
#endif
	{"ceil",1,mceil},
	{"charpoly",1,mcharpoly},
	{"cols",1,mcols},
	{"colsum",1,mcolsum},
	{"complex",1,mcomplex},
	{"conj",1,mconj},
	{"cos",1,mcos},
	{"count",2,mstatistics},
	{"cumprod",1,mcumprod},
	{"cumsum",1,mcumsum},
	{"diag",2,mdiag2},
	{"diag",3,mdiag},
	{"dup",2,mdup},
	{"epsilon",0,mepsilon},
	{"epsilon",1,msetepsilon},
	{"error",1,merror},
	{"exp",1,mexp},
	{"extrema",1,mextrema},
	{"fac",1,mfac},
	{"fft",1,mfft},
	{"find",2,mfind},
	{"flipx",1,mflipx},
	{"flipy",1,mflipy},
	{"floor",1,mfloor},
	{"format",2,mformat},
	{"hb",1,mtridiag},
	{"ifft",1,mifft},
	{"im",1,mim},
	{"index",0,mindex},
	{"input",1,minput},
	{"interp",2,dd},
	{"interpval",3,ddval},
	{"iscomplex",1,miscomplex},
	{"isfunction",1,misfunction},
	{"isreal",1,misreal},
	{"isstring",1,misstring},
	{"isvar",1,misvar},
	{"log",1,mlog},
	{"logbin",2,mlogbin},
	{"logfac",1,mlogfac},
	{"lu",1,mlu},
	{"lusolve",2,mlusolve},
	{"matrix",2,mmatrix},
	{"max",1,mmax1},
	{"max",2,mmax},
	{"min",1,mmin1},
	{"min",2,mmin},
	{"mod",2,mmod},
	{"name",1,mname},
	{"nonzeros",1,mnonzeros},
	{"normal",1,mnormal},
	{"ones",1,mones},
#ifdef PCM_ENABLE
	{"pcmfreq",0,mpcmfreq},
	{"pcmloop",0,mpcmloop},
	{"pcmplay",1,mpcmplay},
	{"pcmrec",1,mpcmrec},
	{"pcmvol",1,mpcmvol},
#endif
	{"plot",2,mplot1},
	{"plot",3,mplot},
	{"plotarea",2,mplotarea},
	{"plotstyle",1,mplotstyle},
	{"polyadd",2,mpolyadd},
	{"polycons",1,mpolycons},
	{"polydiv",2,mpolydiv},
	{"polymult",2,mpolymult},
	{"polyroot",2,mpolyroot},
	{"polysolve",1,mpolysolve},
	{"polytrans",2,polydd},
	{"polytrunc",1,mpolytrunc},
	{"polyval",2,mpolyval},
	{"pqcos",1,mpqcos},
	{"pqfft",1,mpqfft},
	{"pqifft",1,mpqifft},
	{"printf",2,mprintf},
	{"prod",1,mprod},
	{"random",1,mrandom},
	{"re",1,mre},
	{"redim",2,mredim},
	{"round",2,mround},
	{"rows",1,mrows},
#ifdef PCM_ENABLE
	{"rtfilter",2,mrtfilter},
#endif
	{"setdiag",3,msetdiag},
	{"setplot",1,msetplot},
	{"shuffle",1,mshuffle},
	{"sign",1,msign},
	{"sin",1,msin},
	{"size",-1,msize},
	{"sort",1,msort},
	{"sqrt",1,msqrt},
	{"subplot",1,msubplot},
	{"sum",1,msum},
	{"symmult",2,smultiply},
	{"tan",1,mtan},
	{"text",5,mtext},
	{"time",0,mtime},
	{"title",1,mtitle},
	{"wait",1,mwait},
	{"xargs",0,mxargs},
	{"xgrid",5,mxgrid},
	{"xlabel",1,mxlabel},
	{"ygrid",5,mygrid},
	{"ylabel",1,mylabel},
	{"zeros",1,mzeros},
};
#endif
#if 0
	{"lineinput",1,mlineinput},
	{"interpret",1,minterpret},
	{"eval",-1,mdo},
	{"free",0,mfree},
	
	{"char",1,mchar},
	{"key",0,mkey},
	{"setkey",2,msetkey},
	
	
	{"jacobi",1,mjacobi},

	{"normaldis",1,mgauss},
	{"invnormaldis",1,minvgauss},
	{"tdis",2,mtd},
	{"invtdis",2,minvtd},
	{"chidis",2,mchi},
	{"fdis",3,mfdis},

	{"mesh",1,mmesh},
	{"view",1,mview},
	{"view",0,mview0},
	{"textsize",0,mtextsize},
	{"wire",3,mwire},
	{"solid",3,msolid},
	{"solid",4,msolid1},
	{"pixel",0,mpixel},
	{"contour",2,mcontour},
	{"color",1,mcolor},
	{"framecolor",1,mfcolor},
	{"wirecolor",1,mwcolor},
	{"textcolor",1,mtcolor},
	{"linewidth",1,mlinew},
	{"window",1,mwindow},
	{"window",0,mwindow0},
	{"scale",1,mscale},
	{"mouse",0,mmouse},
	{"project",3,mproject},
	{"scaling",1,mscaling},
	{"holding",1,mholding},
	{"holding",0,mholding0},
	{"logscale",1,mlogscale},
	{"logscale",0,mlogscale0},
	{"twosides",1,mtwosides},
	{"triangles",1,mtriangles},
	{"meshfactor",1,mmeshfactor},
	{"frame",0,mframe},
	{"density",1,mdensity},
	{"huecolor",1,mdcolor},
	{"huegrid",1,mdgrid},
	{"solidhue",4,msolidh},
	
	{"store",1,mstore},
	{"restore",1,mrestore},
	
	{"errorlevel",1,merrlevel},
	{"setepsilon",1,msetepsilon}, /* redefined to epsilon(val) */
#endif


#define BINFUNCS	((int)(sizeof(binfunc_list)/sizeof(binfunc_t)))

/* Format string for */
#define XSTR(x) #x
#define STR(x) XSTR(x)
#define OUTFMT "%-" STR(LABEL_LEN_MAX) "s"

static void binfunc_print (Calc *cc)
{	int i, c, cend, lw=cc->termwidth/LABEL_LEN_MAX;
	
	for (i=0; i<BINFUNCS; i+=lw) {
		cend = i+lw<=BINFUNCS ? i+lw : BINFUNCS;
		for (c=i; c<cend ; c++) {
			outputf(cc,OUTFMT,binfunc_list[c].name);
			if (sys_test_key()==escape) return;
		}
		output(cc,"\n");
	}
	output(cc,"\n");
}

static int binfunc_compare (const binfunc_t *p1, const binfunc_t *p2)
{	int h;
	h=strcmp(p1->name,p2->name);
	if (h) return h;
	else {
		if (p1->nargs==-1 || p2->nargs==-1) return 0;
		else if (p1->nargs<p2->nargs) return -1; 
		else if (p1->nargs>p2->nargs) return 1;
		else return 0;
	}
}

#ifndef EMBED
static void binfunc_init (void)
{
	qsort((void*)binfunc_list,BINFUNCS,sizeof(binfunc_t),
		(int (*) (const void *, const void *))binfunc_compare);
}
#endif

header* binfunc_exec (Calc *cc, char *name, int nargs, header *hd)
{
	binfunc_t *b, h;
	h.name= name[0]!='$' ? name : name+1;
	h.nargs=nargs;
	b=bsearch(&h,binfunc_list,BINFUNCS,sizeof(binfunc_t),
		(int (*) (const void *, const void *))binfunc_compare);
	if (b) {
		return b->f(cc,hd);
	}
	return NULL;
}

binfunc_t *binfunc_find (char *name)
{	binfunc_t h;
	h.name= name[0]!='$' ? name : name+1;
	h.nargs=-1;
	return (binfunc_t *)bsearch(&h,binfunc_list,BINFUNCS,sizeof(binfunc_t),
		(int (*) (const void *, const void *))binfunc_compare);
}

/***************************************************************************
 * commands and language structure
 ***************************************************************************/
static void do_quit (Calc *cc);
static void do_cls (Calc *cc);
static void do_list (Calc *cc);
static void do_listvar (Calc *cc);
static void do_help (Calc *cc);
static void do_mdump (Calc *cc);
static void do_hexdump (Calc *cc);
static void do_clear (Calc *cc);
static void do_forget (Calc *cc);
static void do_parse_udf (Calc *cc);
static void do_if (Calc *cc);
static void do_then (Calc *cc);
static void do_else (Calc *cc);
static void do_elseif (Calc *cc);
static void do_endif (Calc *cc);
static void do_end (Calc *cc);
static void do_return (Calc *cc);
static void do_endfunction (Calc *cc);
static void do_show (Calc *cc);
static void do_trace(Calc *cc);
static void do_repeat (Calc *cc);
static void do_loop (Calc *cc);
static void do_for (Calc *cc);
static void do_end (Calc *cc);
static void do_break (Calc *cc);
static void do_to (Calc *cc);
static void do_do (Calc *cc);
static void do_step (Calc *cc);
static void do_global (Calc *cc);
static void do_cd (Calc *cc);
static void do_ls (Calc *cc);
static void do_rm (Calc *cc);
static void do_mkdir (Calc *cc);
static void do_load (Calc *cc);
static void do_dump (Calc *cc);
static void do_const (Calc *cc);
#ifdef WITHGRAPHICS
void do_show_graphics (Calc *cc);
#endif

typedef struct {
	char*	name;
	comtyp	nr;
	void 	(*f)(Calc *cc);
} cmd_t;

#ifndef EMBED
/*sorted according to functionnality, needs alphabet sorting */
static cmd_t cmd_list[] = {
	{"quit",c_quit,do_quit},
	{"cls",c_cls,do_cls},
	{"help",c_global,do_help},
#ifdef WITHGRAPHICS
	{"shg",c_shg,do_show_graphics},
#endif

	{"cd",c_global,do_cd},
	{"ls",c_global,do_ls},
	{"rm",c_global,do_rm},
	{"mkdir",c_global,do_mkdir},
	{"dump",c_global,do_dump},

	{"list",c_global,do_list},
	{"listvar",c_global,do_listvar},
	{"memdump",c_global,do_mdump},
	{"hexdump",c_global,do_hexdump},
	{"clear",c_clear,do_clear},
	{"forget",c_forget,do_forget},

	{"load",c_load,do_load},

	{"if",c_if,do_if},
	{"then",c_then,do_then},
	{"else",c_else,do_else},
	{"elseif",c_elseif,do_elseif},
	{"endif",c_endif,do_endif},

	{"repeat",c_repeat,do_repeat},
	{"loop",c_loop,do_loop},
	{"for",c_for,do_for},
	{"end",c_end,do_end},
	{"break",c_break,do_break},
	{"to", c_to, do_to},
	{"do", c_do, do_do},
	{"step", c_step, do_step},
	
	{"const",c_global,do_const},
	{"global",c_global,do_global},
	
	{"function",c_udf,do_parse_udf},
	{"return",c_return,do_return},
	{"endfunction",c_endfunction,do_endfunction},
	
	{"show",c_global,do_show},
	{"trace",c_global,do_trace},
};
#else
	/* sorted according to alphabet order, can be sent to FLASH memory on
	   embedded platforms */
static const cmd_t cmd_list[] = {
	{"break",c_break,do_break},
	{"cd",c_global,do_cd},
	{"clear",c_clear,do_clear},
	{"cls",c_cls,do_cls},
	{"const",c_global,do_const},
	{"do", c_do, do_do},
	{"dump",c_global,do_dump},
	{"else",c_else,do_else},
	{"elseif",c_elseif,do_elseif},
	{"end",c_end,do_end},
	{"endfunction",c_endfunction,do_endfunction},
	{"endif",c_endif,do_endif},
	{"for",c_for,do_for},
	{"forget",c_forget,do_forget},
	{"function",c_udf,do_parse_udf},
	{"global",c_global,do_global},
	{"help",c_global,do_help},
	{"hexdump",c_global,do_hexdump},
	{"if",c_if,do_if},
	{"list",c_global,do_list},
	{"listvar",c_global,do_listvar},
	{"load",c_load,do_load},
	{"loop",c_loop,do_loop},
	{"ls",c_global,do_ls},
	{"memdump",c_global,do_mdump},
	{"mkdir",c_global,do_mkdir},
	{"quit",c_quit,do_quit},
	{"repeat",c_repeat,do_repeat},
	{"return",c_return,do_return},
	{"rm",c_global,do_rm},
#ifdef WITHGRAPHICS
	{"shg",c_shg,do_show_graphics},
#endif
	{"show",c_global,do_show},
	{"step", c_step, do_step},
	{"then",c_then,do_then},
	{"to", c_to, do_to},
	{"trace",c_global,do_trace},
};
#endif
#if 0
	{"do",c_global,do_do},
	
//	{"echo", c_global, do_echo},
	{"output",c_global,do_output},
	{"meta",c_global,do_meta},
	{"comment",c_global,do_comment},
	
	{"exec",c_exec,do_exec},
#endif

#define CMDS	((int)(sizeof(cmd_list)/sizeof(cmd_t)))

static void cmd_print (Calc *cc)
{	int i, c, cend, lw=cc->termwidth/LABEL_LEN_MAX;
	
	for (i=0; i<CMDS; i+=lw) {
		cend = i+lw<CMDS ? i+lw : CMDS;
		for (c=i; c<cend ; c++) {
			outputf(cc,OUTFMT,cmd_list[c].name);
			if (sys_test_key()==escape) return;
		}
		output(cc,"\n");
	}
	output(cc,"\n");
}

static int cmd_compare (const cmd_t *p1, const cmd_t *p2)
{	return strcmp(p1->name,p2->name);
}

#ifndef EMBED
static void cmd_init (void)
{
	qsort(cmd_list,CMDS,sizeof(cmd_t),
		(int (*)(const void *, const void *))cmd_compare);
}
#endif

static char cmd_find (Calc *cc, int *l)
{	cmd_t h;
	char name[LABEL_LEN_MAX],*a,*n, c;
	*l=0;
	/* parse the name of the command */
	a=cc->next; n=name;
	while (*l<LABEL_LEN_MAX && (((c=*a)>='A' && c<='Z') || (c>='a' && c<='z'))) { *n++=*a++; *l+=1; }
	*n++=0;
	/* name too long! */
	if ((c>='A' && c<='Z') || (c>='a' && c<='z')) return -1;
	/* look for name in the command_list table */
	h.name=name;
	a=bsearch(&h,cmd_list,CMDS,sizeof(cmd_t),
		(int (*)(const void *, const void *))cmd_compare);
	
	if (a) {
		cc->next+=*l;
		return (a-(char*)cmd_list)/sizeof(cmd_t);
	}
	return -1;
}

static int cmd_parse_and_exec (Calc *cc)
/***** command_run
	interpret a binfunc command, number no.
*****/
{
	signed char cmd=-1;
	int l;
	if (*cc->next==3) {		/* run a command from a user defined function */
		cmd=*(cc->next+1);
		cc->next+=2;
//	} else if (udfon) {
//		return 0;
	} else cmd=cmd_find(cc,&l);
	if (cmd!=-1) {
		cmd_list[(int)cmd].f(cc);
		return cmd_list[(int)cmd].nr;
	}
	return cmd;
}

token_t cmd2tok(int cmd)
{
	token_t tok;
	switch (cmd_list[(int)cmd].nr) {
	case c_then: tok=T_THEN; break;
	case c_endif: tok=T_ENDIF; break;
	case c_else: tok=T_ELSE; break;
	case c_elseif: tok=T_ELSEIF; break;
	case c_to: tok=T_TO; break;
	case c_do: tok=T_DO; break;
	case c_end: tok=T_END; break;
	case c_step:tok=T_STEP; break;
	default: tok=T_NONE; break;
	}
	return tok;
}

/*****************************************************************************/
static void do_quit (Calc *cc)
{
	cc->quit=1;
}

static void do_cls (Calc *cc)
{	text_mode(); sys_clear();
}

#ifdef WITHGRAPHICS
void do_show_graphics (Calc *cc)
{	scan_t scan;
	graphic_mode(); sys_wait_key(&scan); text_mode();
}
#endif

static void do_list (Calc *cc)
{	header *hd;
	int i=0, lw=cc->termwidth/LABEL_LEN_MAX;
	sys_out_mode(CC_OUTPUT);
	output(cc,"  *** Builtin functions:\n");
	binfunc_print(cc);
	output(cc,"  *** Commands:\n");
	cmd_print(cc);
	output(cc,"  *** Your functions:\n");
	hd=(header *)cc->udfstart;
	while ((char*)hd<cc->udfend)
	{	if (hd->type!=s_udf) break;
		outputf(cc,OUTFMT,hd->name);
		if (sys_test_key()==escape) return;
		hd=nextof(hd);
		i++;
		if (i==lw) {
			i=0;
			output(cc,"\n");
		}
	}
	if (i) output(cc,"\n");
}

static const char * sname[] = {
	"real", "complex", "real matrix", "complex matrix", "string", "user function",
	"reference", "real submatrix reference", "complex submatrix reference", "function reference", "command"
};

static void do_listvar (Calc *cc)
{	header *hd=(header *)cc->startlocal;
	sys_out_mode(CC_OUTPUT);
	while (hd<(header*)cc->endlocal) {
		switch (hd->type) {
		case s_real:
		case s_complex:
		case s_string:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s",hd->name,sname[hd->type]);
			break;
		case s_matrix:
		case s_cmatrix:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s (%dx%d)",hd->name,
				sname[hd->type],dimsof(hd)->r,dimsof(hd)->c);
			break;
		case s_reference:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s",hd->name,sname[hd->type]);
			break;
		case s_udf:
		case s_funcref:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s",hd->name,sname[hd->type]);
			break;
		case s_submatrixref:
		case s_csubmatrixref:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s (%dx%0d)\n",hd->name,
				sname[hd->type],submdimsof(hd)->r,submdimsof(hd)->c);
			break;
		default:
			outputf(cc,"%-" STR(LABEL_LEN_MAX) "s : %s",hd->name,"unknown type");
			break;
		}
		if (hd->flags & FLAG_CONST) output(cc," [const]\n");
		else output(cc,"\n");
		hd=nextof(hd);
		if (sys_test_key()==escape) break;
	}
}

static void do_help (Calc *cc)
{
	header *hd;
	int count=LABEL_LEN_MAX,i,defaults;
	char name[LABEL_LEN_MAX+1];
	char *p=name,*end,*pnote;
	char c;
#ifndef EMBED
	binfunc_t *b=binfunc_list;
#else
	binfunc_t const *b=binfunc_list;
#endif
//	token_t tok=scan(cc);
	/* skip space */
	sys_out_mode(CC_OUTPUT);
	while ((c=*cc->next)==' ' || c=='\t') cc->next++;
	while (*cc->next && count) {
		count--; *p++=*cc->next++;
	}
	*p=0;
	if ((c=*cc->next)!=0 && c!=';') cc_error(cc,"bad topic in help!");

//	if (tok!=T_LABEL) cc_error(cc,"name of variable or binfunc or user function expected!");
	b=binfunc_find(name);
	if ((b=binfunc_find(name))!=NULL) {
		if (b->nargs>=0) {
			outputf(cc,"%s is a binary function with %d argument(s).\n",name,b->nargs);
		} else {
			outputf(cc,"%s is a binary function.\n",name);
		}
	}
	if ((hd=searchudf(cc,name))!=NULL && hd->type==s_udf) {
	   if (b) outputf(cc,"%s is also a user defined function.\n",name);
		outputf(cc,"function %s (",name);
		end=udfof(hd);
		p=udfargsof(hd);
		/* get the number of arguments */
		count=*((int *)p); p+=sizeof(int);
		/* get the default value bitmap */
		defaults=*(unsigned int*)p; p+=sizeof(unsigned int);
		pnote=p;
		for (i=0; i<count; i++) {
			udf_arg* arg=(udf_arg*)p;
			outputf(cc,"%s",arg->name);
			if (defaults & (1<<i)) {
				output(cc,"=...");
			}
			if (i!=count-1) output(cc,",");
			p=udfnextarg(p, defaults & 1<<i);
		}
		output(cc,")\n");
		p=pnote;
		for (i=0; i<count; i++) {
			if (defaults & 1<<i) {
				header* arg=(header*)p;
				outputf(cc,"## Default for %s :\n",arg->name);
				give_out(cc, arg);
			}
			p=udfnextarg(p, defaults & 1<<i);
		}
		while (p<end) {
			output(cc,p); output(cc,"\n");
			p+=strlen(p); p++;
		}
	} else if ((hd=searchvar(cc,name))!=NULL) {
		switch (hd->type) {
		case s_real:
		case s_complex:
		case s_string:
			outputf(cc,"%s is a %s variable\n",hd->name,sname[hd->type]);
			break;
		case s_matrix:
		case s_cmatrix:
			outputf(cc,"%s is a %s (%dx%d) variable\n",hd->name,
				sname[hd->type],dimsof(hd)->r,dimsof(hd)->c);
			break;
		case s_reference:
			outputf(cc,"%s is a %s variable\n",hd->name,sname[hd->type]);
			break;
		case s_submatrixref:
		case s_csubmatrixref:
			outputf(cc,"%s is a %s (%dx%0d) variable",hd->name,
				sname[hd->type],submdimsof(hd)->r,submdimsof(hd)->c);
			break;
		default:
			outputf(cc,"%s is a %s\n variable",hd->name,"unknown type");
			break;
		}
	}
	char buf[256];
	FILE *file=fopen(HELPFILE,"r");
	if (file) {
		while (fgets(buf,sizeof(buf),file)) {
			if (buf[0]=='!') {
				char *s=buf+1, *e=buf+1;
				while (*s) {
					char c;
					while ((c=*e)!=0 && c!='\n' && c!=',') e++;
					*e=0;
					if (strcmp(name,s)==0) {	/* key found! */
						*e=c;
						//output(cc,buf+1);
						while (fgets(buf,sizeof(buf),file)) {
							if (buf[0]=='#') goto done;
							output(cc,buf);
						}
						goto done;
					}
					*e=c;
					if (c==',') s=e=e+1;
					else break;
				}
			}
		}
done:
		fclose(file);
	}
}

static void do_mdump (Calc *cc)
{	header *hd;
	sys_out_mode(CC_OUTPUT);
	output(cc,"ramstart   : 0x00000000\n");
	outputf(cc,"startlocal : 0x%08X\n",cc->startlocal-cc->ramstart);
	outputf(cc,"endlocal   : 0x%08X\n",cc->endlocal-cc->ramstart);
	outputf(cc,"newram     : 0x%08X\n",cc->newram-cc->ramstart);
	outputf(cc,"ramend     : 0x%08X\n\n",cc->ramend-cc->ramstart);
	hd=(header *)cc->ramstart;
	while ((char *)hd<cc->newram) {
		outputf(cc,"0x%08X : %-" STR(LABEL_LEN_MAX) "s, ",(char *)hd-cc->ramstart,hd->name);
		outputf(cc,"size %6ld, ",(long)hd->size);
		outputf(cc,"type %s",sname[hd->type]);
		if (hd->flags & FLAG_CONST) output(cc," [const]\n");
		else output(cc,"\n");
		hd=nextof(hd);
	}
	hd=(header *)cc->udfstart;
	while ((char *)hd<cc->udfend) {
		outputf(cc,"0x%08X : %-" STR(LABEL_LEN_MAX) "s, ",(char *)hd-cc->ramstart,hd->name);
		outputf(cc,"size %6ld, ",(long)hd->size);
		outputf(cc,"type %s\n",sname[hd->type]);
		hd=nextof(hd);
	}
}

static void string_out (Calc *cc, unsigned char *p)
{	int i;
	unsigned char a;
	for (i=0; i<16; i++) 
	{	a=*p++;
		outputf(cc,"%c",(a<' ')?'_':a);
	}
}

static void do_hexdump (Calc *cc)
{
	unsigned char *p,*end;
	int i=0,j;
	ULONG count=0;
	header *hd;
	token_t tok=scan(cc);
	if (tok!=T_LABEL) cc_error(cc,"name of variable or user function expected!");

	hd=searchvar(cc,cc->str);
	if (!hd) hd=searchudf(cc,cc->str);
	if (!hd) cc_error(cc,"%s not a variable or user function!", cc->str);
	p=(unsigned char *)hd; end=p+hd->size;
	sys_out_mode(CC_OUTPUT);
	outputf(cc,"\n%5lx: ",count);
	while (p<end) {
		outputf(cc,"%02X ",*p++); i++; count++;
		if (i>=16) {
			i=0; string_out(cc,p-16);
			outputf(cc,"\n%5lx: ",count);
			if (sys_test_key()==escape) break;
		}
	}
	for (j=i; j<16; j++) output(cc,"   ");
	string_out(cc,p-i);
	output(cc,"\n");
}

static void do_clear (Calc *cc)
/***** clear
	clear all the global variables from the stack or just named ones.
*****/
{
	token_t tok;
	if (CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"Cannot use \"clear\" in a function!");
	tok=scan(cc);
	if (tok==T_SEMICOL || tok==T_COMMA || tok==T_EOS) {
		/* clear all the variables, but const */
		header*hd;
		for (hd=(header*)cc->startlocal;(hd->flags & FLAG_CONST) && (hd!=(header*)cc->endlocal);hd=nextof(hd)) ;
		cc->endlocal=(char*)hd;
	} else while(tok==T_LABEL) {
		kill_local(cc,cc->str);
		if (scan(cc)!=T_COMMA) break;
		tok=scan(cc);
	}
}

static void do_forget (Calc *cc)
{
	if (CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"Cannot 'forget' functions in a function!");
	while (1) {
		token_t tok=scan(cc);
		if (tok==T_LABEL) {
			if (!kill_udf(cc,cc->str)) cc_error(cc,"Function %s not found!",cc->str);
			if ((tok=scan(cc))==T_COMMA) continue;
			if (tok==T_SEMICOL || tok==T_EOS) break;
		}
	}
}

static void do_cd (Calc *cc)
{
	token_t tok;
	char *path;
	
	tok=scan_path(cc);
	if (cc->result) {
		path=fs_cd(stringof(cc->result));
		if (tok!=T_SEMICOL) {
			sys_out_mode(CC_OUTPUT);
			outputf(cc,"%s\n",path);
		}
	}
}

static void do_ls (Calc *cc)
{
	char pattern[16]="*";
	int len, npl, i, j, k, imax, outputing=1;
	char **entries=NULL;
	int n_entries=0;
	token_t tok;
	char *path=NULL;
	
	tok=scan_path(cc);
	if (cc->result) path=stringof(cc->result);
	if (!path || (path && !path[0])) path=".";
	
	if (tok==T_SEMICOL) return;
	sys_out_mode(CC_OUTPUT);
	len = fs_dir(path,pattern,&entries,&n_entries);
	len +=2;
	npl = cc->termwidth/len;
	imax = n_entries/npl;
	if (n_entries % npl) imax++;
	for (i=0; i<imax;i++) {
		for (j=0; j<npl; j++) {
			if (outputing) {
				int l; 
				if (npl*i+j>=n_entries) break;
				outputf(cc,"%s", entries[npl*i+j]);
				l=strlen(entries[npl*i+j]);
				for (k=0;k<len-l;k++)
					output(cc," ");
				if (sys_test_key()==escape) outputing=0;
			}
			free(entries[npl*i+j]);
		}
		output(cc,"\n");
	}
	free(entries);
	if (!outputing) output(cc,"\n");
}

static void do_rm (Calc *cc)
{
	int res;
	
	scan_path(cc);
	if (cc->result && strlen(stringof(cc->result))!=0) {
		res=fs_rm(stringof(cc->result));
		if (res<0) cc_error(cc,"error: can't remove that file or directory");
	}
}

static void do_mkdir (Calc *cc)
{
	int res;
	
	scan_path(cc);
	if (cc->result && strlen(stringof(cc->result))!=0) {
		res=fs_mkdir(stringof(cc->result));
		if (res<0) cc_error(cc,"error: can't make that directory");
	}
}

static void do_dump (Calc *cc)
{
	if (cc->outfile) {
		fclose(cc->outfile);
		cc->outfile=NULL;
	}
	scan_path(cc);
	if (cc->result && strlen(stringof(cc->result))!=0) {
		cc->outfile=fopen(stringof(cc->result),"a");
		if (!cc->outfile) cc_error(cc,"Could not open %s.",stringof(cc->result));
	}
}

/******************************* module loading *****************************/
/***** load_file
	interpret a file.
*****/
static void do_load (Calc *cc)
{
	/* saved context */
	jmp_buf *oldenv;
	char *oldline,*oldnext;
	FILE *oldinfile;
	int oldlinenb;
	/* locals */
	jmp_buf env;
	char *filename;
	char input[LINEMAX]="";
	char file[LINEMAX];
	
	if (CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"Cannot load a file in a function!");
	
	/* get a file name or path*/
	scan_path(cc);
	if (!(cc->result && strlen(stringof(cc->result))!=0)) 
		cc_error(cc,"load <filename>");
	filename=stringof(cc->result);
	strncpy(file,filename,LINEMAX-1);file[LINEMAX-1]=0;
	
	/* try to open it */
	oldinfile=cc->infile;
	if (filename[0]==PATH_DELIM_CHAR) {	/* an absolute path, use it */
		cc->infile=fopen(filename,"r");
	} else {							/* use standard path */
		for (int k=0; k<npath; k++) {
			char fn[strlen(path[k])+strlen(filename)+strlen(EXTENSION)+strlen(PATH_DELIM_STR)+1];
			
			strcpy(fn,path[k]);strcat(fn,PATH_DELIM_STR);strcat(fn,filename);
			cc->infile=fopen(fn,"r");
			if (!cc->infile) {
				strcat(fn,EXTENSION);
				cc->infile=fopen(fn,"r");
				if (cc->infile) break;
			} else break;
		}
	}
	
	/* interpret the file if it exists */
	if (cc->infile) {
		/* set synchronisation point for error handling.
		   save the context, setup the new one before call setjmp to help
		   compiler optimization
		 */
		oldenv=cc->env;
		oldnext=cc->next;
		oldline=cc->line;
		oldlinenb=cc->linenb;
		cc->linenb=0;
		cc->env=&env;
		
		switch (setjmp(env)) {
		case 0:
			break;
		case 4:
			/* end of file */
			cc->next=oldnext;
			cc->line=oldline;
			cc->linenb=oldlinenb;
			cc->env=oldenv;
			fclose(cc->infile);
			cc->infile=oldinfile;
			return;
		default:
			/* error in a statement in the package or in a child package */
			outputf(cc,"Error when interpreting '%s', line %d\n", file, cc->linenb);
			cc->next=oldnext;
			cc->line=oldline;
			cc->env=oldenv;
			fclose(cc->infile);
			cc->infile=oldinfile;
			longjmp(*cc->env,2);	/* back to enclosing error handler */
			break;
		}

		cc->next=cc->line=input;*cc->line=0;
		while (1) {
			/* reset global context for commands evaluated in the 
	    	   lower level context (global scope) */
		    CC_UNSET(cc,CC_NOSUBMREF|CC_PARSE_INDEX|CC_PARSE_PARAM_LIST|CC_PARSE_UDF);
			cc->newram=cc->globalend=cc->endlocal;
			parse(cc);
		}
		/* we should never be here! */
		cc->next=oldnext;
		cc->line=oldline;
		cc->linenb=oldlinenb;
		cc->env=oldenv;
		fclose(cc->infile);
		cc->infile=oldinfile;
	} else {
		cc->infile=oldinfile;
		outputf(cc,"  %s\n",cc->line);
		cc_error(cc,"Could not open %s!",filename);
	}
}

/*********************** programming language structure **********************/

char *type_udfline (Calc *cc, char *start)
{	char outline[LINEMAX],*p=start,*q;
	real x;
	signed char cmd;
	q=outline;
	while (*p) {
		if (*p==2) {
			/* a constant in IEEE simple/double precision, convert it back */
			p++; memmove((char *)(&x),p,sizeof(real));
			p+=sizeof(real);
			sprintf(q,"%g",x);
			q+=strlen(q);
		} else if (*p==3) {
			/* a command */
			p++; cmd=*p++;
			sprintf(q,"%s",cmd_list[(int)cmd].name);
			q+=strlen(q);
			if (cmd_list[(int)cmd].f==do_endfunction) {
				/* end of function */
				*q++=0;
				output(cc,outline); output(cc,"\n");
				return NULL;
			}
		} else *q++=*p++;
		if (q>outline+LINEMAX-2) {
			q=outline+LINEMAX-1;
			break;
		}
	}
	*q++=0;
	output(cc,outline); output(cc,"\n");
	return p+1;
}


void trace_udfline (Calc* cc, char *next)
{	int oldtrace;
	scan_t scan;
	header *hd;
	sys_out_mode(CC_OUTPUT);
	outputf(cc,"%s: ",cc->running->name); type_udfline(cc,next);
again:
	sys_wait_key(&scan);
	switch (scan)
	{	case fk1 :
		case cursor_down :
			break;
		case fk2 :
		case cursor_up :
			cc->trace=2; break;
		case fk3 :
		case cursor_right :
			cc->trace=0; break;
		case fk4 :
		case help :
			hd=(header *)cc->newram;
			oldtrace=cc->trace; cc->trace=0;
			new_cstring(cc,"Expression",12,"");
			minput(cc,hd);
	 		cc->newram=(char *)hd;
			cc->trace=oldtrace;
			goto again;
		case fk9 :
		case escape :
			output(cc,"Trace interrupted\n"); break;
		case fk10 :
			cc->trace=-1; break;
		default :
			output(cc,
				"\nKeys:\n"
				"F1 (cursor_down)  Single step\n"
				"F2 (cursor_up)    Step over subroutines\n"
				"F3 (cursor_right) Go until return\n"
				"F4 (help)         Evaluate expression\n"
				"F9 (escape)       Abort execution\n"
				"F10               End trace\n\n");
			goto again;
	}
}

static void do_trace(Calc *cc)
/**** do_trace
	toggles tracing or sets the trace bit of a udf.
****/
{	header *f;
	char c;
	while ((c=*cc->next)=='\t' || c==' ') cc->next++;	
	if (!strncmp(cc->next,"off",3)) {
		cc->trace=0; cc->next+=3;
	} else if (!strncmp(cc->next,"alloff",6)) {
		cc->next+=6;
		f=(header *)cc->udfstart;
		while ((char *)f<cc->udfend && f->type==s_udf) {
			f->flags&=~1;
			f=nextof(f);
		}
		cc->trace=0;
	} else if (!strncmp(cc->next,"on",2)) {
		cc->trace=1; cc->next+=2;
	} else if (*cc->next==';' || *cc->next==',' || *cc->next==0) {
		cc->trace=!cc->trace;
	} else if (scan(cc)==T_LABEL) {
		f=searchudf(cc,cc->str);
		if (!f || f->type!=s_udf) cc_error(cc,"Function not found!");
		f->flags^=1;
		if (f->flags & 1) outputf(cc,"Tracing %s\n",cc->str);
		else outputf(cc,"No longer tracing %s\n",cc->str);
		while ((c=*cc->next)=='\t' || c==' ') cc->next++;	
	}
	if (*cc->next==';' || *cc->next==',') cc->next++;
}

static void do_show (Calc *cc)
{
	header *hd;
	char *p,*pnote;
	int i,count;
	unsigned int defaults;
	token_t tok=scan(cc);
	if (tok!=T_LABEL) cc_error(cc,"A UDF name is expected");
	hd=searchudf(cc,cc->str);
	if (hd && hd->type==s_udf) {
		sys_out_mode(CC_OUTPUT);
		outputf(cc,"function %s (",cc->str);
		p=udfargsof(hd);
		/* get the number of arguments */
		count=*((int *)p); p+=sizeof(int);
		/* get the default value bitmap */
		defaults=*(unsigned int*)p; p+=sizeof(unsigned int);
		pnote=p;
		for (i=0; i<count; i++) {
			udf_arg* arg=(udf_arg*)p;
			outputf(cc,"%s",arg->name);
			if (defaults & (1<<i)) {
				output(cc,"=...");
			}
			if (i!=count-1) output(cc,",");
			p=udfnextarg(p, defaults & 1<<i);
		}
		output(cc,")\n");
		p=pnote;
		for (i=0; i<count; i++) {
			if (defaults & 1<<i) {
				header* arg=(header*)p;
				outputf(cc,"## Default for %s :\n",arg->name);
				give_out(cc, arg);
			}
			p=udfnextarg(p, defaults & 1<<i);
		}
		p=udfof(hd);
		while (p && p<(char *)nextof(hd))
			p=type_udfline(cc,p);
	} else {
		cc_error(cc,"No such UDF function!");
	}
}

void scan_endif (Calc *cc)
/***** scan_endif
	look for 'endif' to step over instructions. If an inner if is found look for
	its corresponding 'endif'
*****/
{
	signed char cmd;
	void (*f)(Calc *cc);
	char *oldline=cc->line;
	while (1) {
		switch (*cc->next) {
		case 0 : cc->line=cc->next+1; cc->next++; break;
		case 2 : cc->next+=1+sizeof(real); break;
		case 3 :
			cc->next++;
			cmd=*cc->next++;
			f=cmd_list[(int)cmd].f;
			if (f==do_endif) {
				if (cc->trace>0) trace_udfline(cc,cc->line);
				return;
			} else if (f==do_if) {
				scan_endif(cc);
			} else if (f==do_endfunction) {
				cc->line=oldline;
				cc_error(cc,"'endif' missing!");
			}
			break;
		default :
			cc->next++;
		}
	}
}

static int scan_else (Calc *cc)
/***** scan_else
	the condition is false, so do execute the 'then' case and look for
	'else' or 'endif' or 'elseif'
	
	If an inner if is defined in the 'then' section, look over it (so look for
	its corresponding 'endif'
*****/
{
	signed char cmd;
	void (*f)(Calc *cc);
	char *oldline=cc->line;
	while (1) {
		switch (*cc->next) {
		case 0 : cc->line=cc->next+1; cc->next++; break;
		case 2 : cc->next+=1+sizeof(real); break;
		case 3 : cc->next++; 
			cmd=*cc->next++;
			f=cmd_list[(int)cmd].f;
			if (f==do_endif || f==do_else) {
				if (cc->trace>0) trace_udfline(cc,cc->line);
				return 0;
			} else if (f==do_elseif) {	/* else if : go back and evaluate the new condition */
				return 1;
			} else if (f==do_if) {		/* inner if */
				scan_endif(cc);
			} else if (f==do_endfunction) {
				cc->line=oldline;
				cc_error(cc,"'endif' missing!");
			}
			break;
		default :
			cc->next++;
		}
	}
}

static void scan_end (Calc *cc)
/***** scan_end
	scan for "end".
*****/
{
	signed char cmd;
	void (*f)(Calc *cc);
	char *oldline=cc->line;

	while (1) {
		switch (*cc->next) {
		case 0 : cc->line=cc->next+1; cc->next++; break;
		case 2 : cc->next+=1+sizeof(real); break;
		case 3 : cc->next++;
			cmd=*cc->next++;
			f=cmd_list[(int)cmd].f;
			if (f==do_end) {
				if (cc->trace>0) trace_udfline(cc,cc->line);
				return;
			} else if (f==do_repeat /*|| f==do_loop || f==do_for*/) {
				scan_end(cc); return;
			} else if (f==do_endfunction) {
				cc->line=oldline;
				cc_error(cc,"'end' missing!");
			}
			break;
		default : cc->next++;
		}
	}
}

static void do_parse_udf (Calc *cc)
/***** parse_udf
	define a user defined function.

   user defined function on the stack
   
   udf header			size (of the udf),name[LABEL_LEN_MAX+1],xor,s_udf,flags=0
   ULONG offset			offset to the beginning of the code
   int paramcnt			formal parameter count
   unsigned def_bitmap	bitmap showing parameter with default value (32 max)
   -------------- parameter list ----------------
   udf_arg				just a formal parameter (name+xor)
     or
   header+data			includes the formal parameter name+xor, but also, the 
   						complete header and data for the default value.
   						
   						parameter list can be run through with the next_arg()
   						macro
   ...
   ------------------- body ---------------------
   function body
   with constants precompiled: char=2+double value
   and pointers to basic statement handlers: char=3+handler pointer
   byte = 1 			endfunction encountered
*****/
{
	char name[LABEL_LEN_MAX+1],*p,*firstchar,*startp;
	int *pcount, count=0;				/* argument counter */
	unsigned int *pdefmap,defmap=0;		/* pointer to the default value bitmap */ 
	int l;
	header *var,*result,*hd;
	signed char cmd;
	real x;
	token_t tok;
	
	if (CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"Cannot define a function in a function!");
	if ((tok=scan(cc))==T_FUNCREF || tok==T_LABEL) {
		/* we got function fname( ... or function fname */
		kill_udf(cc,cc->str);			/* kill any udf already defined with his name */
		var=new_reference(cc,0,cc->str);	/* push an empty reference to the function name (for assignment) */
		result=new_udf(cc,"");				/* create a new function */
		CC_SET(cc,CC_PARSE_UDF);		/* switch to udf input prompt! */
		p=(char*)(result+1)+sizeof(ULONG);	/* start of parameter section */
		pcount=(int *)p; p+=sizeof(int);/* leave room for parameter count, keep the pointer ph to this place */
		pdefmap=(unsigned int*)p; p+=sizeof(unsigned int); /* leave room for defaut value bitmap */
		cc->newram=p;
		
		/* parse parameter list () ou (arg1     ,
		 *                             arg2=val )
		 */
		if (tok==T_FUNCREF) {			/* parse parameter list */
			CC_SET(cc,CC_PARSE_PARAM_LIST);
			while(1) {
				tok=scan(cc);
				if (tok==T_LABEL) {
					strcpy(name,cc->str);
					tok=scan(cc);
					if (tok==T_ASSIGN) {
						/* get header + data of the default parameter value */
						cc->newram=p; hd=(header *)p;
						tok=parse_expr(cc);
						if (cc->result) cc->result=getvalue(cc,cc->result);
						/* allow only binary function references as default? */
			// ?			if (cc->result->type==s_funcref && !(cc->result->flags & FLAG_BINFUNC)) cc_error(cc,"references to UDF not allowed in default parameters");
						if (cc->result!=hd) {	/* result was a reference */
							memmove(hd,cc->result,cc->result->size);
							cc->newram=(char*)hd+cc->result->size;
						}
						strcpy(hd->name,name); hd->xor=xor(name);
						p=cc->newram;			/* update pointer for the next parameter */
						defmap |= 1<<count;		/* update the default value bitmap */
					} else {
						/* parameter without default value, juste store the name */
						udf_arg* arg=(udf_arg*)p;
						strcpy(arg->name,name);
						arg->xor=xor(name);
						p+=sizeof(udf_arg);		/* update pointers */
						cc->newram=p;
					}
					count++;
					if (tok==T_COMMA) continue;
				}
				if (tok==T_RPAR) break;
				else cc_error(cc,"Error in parameter list!");
			}
			CC_UNSET(cc,CC_PARSE_PARAM_LIST);
		}
		*pcount=count;
		*pdefmap=defmap;
		
		/* parse for the help section of the udf */
		if ((tok=scan(cc))==T_EOS) {
			next_line(cc);
			while (1) {
				if (*cc->next=='#' && *(cc->next+1)=='#') {
					while (*cc->next) {
						*p++=*cc->next++;
						if (p>=cc->udfstart) cc_error(cc,"Memory overflow!");
					}
					*p++=0; next_line(cc);
				} else break;
			}
		}
		
		/* parse the body of the function */
		*udfstartof(result)=(p-(char *)result);
		startp=p;
		firstchar=cc->next;
		while (1) {
			if (*cc->next=='#' && *(cc->next+1)=='#') {
				*p++=0; next_line(cc); firstchar=cc->next;
			} else if (*cc->next) {	
				if (*cc->next=='"') {	
#ifndef MULTILINE_STRING
					*p++=*cc->next++;
					while (*cc->next!='"' && *cc->next) *p++=*cc->next++;
					while(1) {
						if ((*cc->next=='\"' && *(cc->next-1)!='\\') || *cc->next==0) break;
						*p++=*cc->next++;
					}
					if (*cc->next=='"') *p++=*cc->next++;
					else {
						cc_error(cc,"\" missing while parsing user defined function %s.",var->name);
					}
#else
					int len=0;		/* verify the string is terminated by 
									   a '"' and its size<MAXLINE, allow multiline */
					*p++=*cc->next++;
					while (len<MAXLINE) {
						if (*cc->next==0) {
							len++;		// for '\n'
							*p++=0; next_line(cc); firstchar=cc->next;
						} else {
							if (*next=='"' && *(cc->next-1)!='\\') {
								break;
							}
							len++;
							*p++=*cc->next++;
						}
					}
					if (*cc->next=='"') *p++=*cc->next++;
					else {
						cc_error(cc,"\" missing while parsing user defined function %s.",var->name);
					}
#endif
				} else if ( ISDIGIT(*cc->next) || 
					(*cc->next=='.' && ISDIGIT(*(cc->next+1))) ) {
					/* precompile numbers */
					if (cc->next!=firstchar && ISALPHA(*(cc->next-1))) {
						*p++=*cc->next++;
						while (ISDIGIT(*cc->next)) *p++=*cc->next++;
					} else {
						/* write byte=0x02 to signal a precompiled float */
						*p++=2;
						scan(cc); x=cc->val;
						// push the number to the function body
			   			memmove(p,(char *)(&x),sizeof(real));
			   			p+=sizeof(real);
			   			if (*(cc->next-1)=='i' || *(cc->next-1)=='j') *p++=*(cc->next-1);
				   	}
				} else if (ISALPHA(*cc->next) &&
					(cc->next==firstchar || !ISALPHA(*(cc->next-1))) &&
					(cmd=cmd_find(cc,&l))!=-1) {
					/* new command found: cmd = cmd index in command_list table
					   cast cmd to char ==> no more than 255 commands can be 
					   handled (It leaves quite some room) */
					/* push byte=0x03 to signal a precompiled command */
					*p++=3;*p++=cmd;
					if (cmd_list[(int)cmd].nr==c_endfunction) {
						if (p==startp || *(p-1)) *p++=0;
						break;
					}
				} else if (*cc->next=='.' && *(cc->next+1)=='.') {
					*p++=' '; next_line(cc); firstchar=cc->next;
				} else *p++=*cc->next++;
			} else { *p++=0; next_line(cc); firstchar=cc->next; }
			if (p>=cc->udfstart-80) cc_error(cc,"Memory overflow!");
		}
		CC_UNSET(cc,CC_PARSE_UDF);
		result->size=ALIGN(p-(char *)result);
		cc->newram=(char *)result+result->size;
		assign(cc, var,result);
	} else cc_error(cc,"No name defined!");
}

static void do_global (Calc *cc)
{
	token_t tok;
	char r;
	header *hd, *hd1;
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'global' only allowed in functions!");

	/* parse 'global *' */
	tok=scan(cc);
	if (tok==T_MUL && ((tok=scan(cc))==T_SEMICOL || tok==T_EOS)) {
		CC_SET(cc,CC_SEARCH_GLOBALS);
	} else {
		/* parse list of identifiers */
		while (tok==T_LABEL) {
			r=xor(cc->str);
			hd=(header *)cc->globalstart;
			if (hd==(header *)cc->globalend) break;
			while ((char *)hd<cc->globalend) {
				if (r==hd->xor && !strcmp(hd->name,cc->str)) break;
				hd=nextof(hd);
			}
			if ((char *)hd>=cc->startlocal) cc_error(cc,"Variable %s not found!",cc->str);
			// test if variable already exists as a local variable
			for (hd1=(header*)cc->startlocal ; (char*)hd1<cc->endlocal ; hd1=nextof(hd1)) {
				if (r==hd1->xor && !strcmp(hd1->name,cc->str)) cc_error(cc,"Local variable with same name");
			}
			cc->newram=cc->endlocal;
			hd=new_reference(cc,hd,cc->str);
			cc->newram=cc->endlocal=(char *)nextof(hd);
			if ((tok=scan(cc))!=T_COMMA) break;
			tok=scan(cc);
		}
		if (!(tok==T_SEMICOL || tok==T_EOS)) cc_error(cc,"Bad separator");
	}
}

static void do_const (Calc *cc)
{
	token_t tok;
	tok=scan(cc);
	if (tok==T_LABEL) {
		header *hd=searchvar(cc,cc->str);
		if (hd) {
			tok=scan(cc);
			if (tok==T_ASSIGN) cc_error(cc,"const value already defined");
			else {
				hd->flags |= FLAG_CONST;
				if (tok!=T_SEMICOL) give_out(cc,hd);
			}
		} else {
			/* defining a new const, push a reference to the new name,
			   and parse the value;
			 */
			hd=new_reference(cc,NULL,cc->str);
			tok=scan(cc);
			if (tok!=T_ASSIGN) cc_error(cc,"new constant needs a value");
			CC_SET(cc,CC_NOSUBMREF); tok=parse_expr(cc); CC_UNSET(cc,CC_NOSUBMREF);
			cc->result=assign(cc,hd,cc->result);
			cc->result->flags |= FLAG_CONST;
			if (tok!=T_SEMICOL) give_out(cc,cc->result);
		}
		if (!(tok==T_SEMICOL || tok==T_EOS)) cc_error(cc,"Bad separator");
	} else cc_error(cc,"variable name expected");
}

static void do_return (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"No user defined function active!");
	CC_SET(cc, CC_EXEC_RETURN);
}

static void do_endfunction (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"No user defined function active!");
	CC_UNSET(cc, CC_EXEC_UDF);
}


static int ctest (Calc *cc, header *hd)
/**** ctest
	test, if a matrix contains nonzero elements.
****/
{	real *m;
	LONG n,i;
	hd=getvalue(cc,hd);
	if (hd->type==s_string) return (*stringof(hd)!=0);
	if (hd->type==s_real) return (*realof(hd)!=0.0);
	if (hd->type==s_complex) return (*realof(hd)!=0.0 ||
		*imagof(hd)!=0.0);
	if (hd->type==s_matrix) {
		n=(LONG)(dimsof(hd)->r)*dimsof(hd)->c;
		m=matrixof(hd);
		for (i=0; i<n; i++) if (*m++!=0.0) return 1;
		return 0;
	}
	if (hd->type==s_cmatrix) {
		n=(LONG)(dimsof(hd)->r)*dimsof(hd)->c;
		m=matrixof(hd);
		for (i=0; i<n; i++) {
			if (*m!=0.0 || *(m+1)!=0.0) return 1;
			m+=2;
		}
		return 0;
	}
	return 0;
}

static void do_if (Calc *cc)
{
	int flag=0;
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'if' only allowed in functions!");
	do {
		token_t tok=parse_expr(cc);
		if (!(tok==T_THEN || tok==T_SEMICOL || tok==T_COMMA || tok==T_EOS))
			cc_error(cc,"'then' or ';' or ',' or '\\n' expected!");
		if (cc->result) flag=ctest(cc,cc->result);
		else cc_error(cc,"Bad condition");
	} while (!flag && scan_else(cc));
}

static void do_then (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'then' only allowed in functions!");
}

static void do_else (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'else' only allowed in functions!");
	scan_endif(cc);
}

static void do_elseif (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'elseif' only allowed in functions!");
	scan_endif(cc);
}

static void do_endif (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'endif' only allowed in functions!");
}

static void do_for (Calc *cc)
/***** do_for
	do a for command in a UDF.
	for i=value to value step value do .... end
*****/
{
	token_t tok;
	int h,signum;
	char name[LABEL_LEN_MAX+1],*jump;
	header *hd,*init=NULL,*end=NULL,*step=NULL;
	real vend,vstep;
	struct { header hd; real value; } rv;
	
	rv.hd.type=s_real; *rv.hd.name=0;
	rv.hd.size=sizeof(header)+sizeof(real); rv.value=0.0;

	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'for' only allowed in functions!");
	/* parse loop variable name */
	tok=scan(cc);
	if (tok!=T_LABEL) cc_error(cc,"'for' var= ...");
	strcpy(name,cc->str);
	tok=scan(cc);
	if (tok!=T_ASSIGN) cc_error(cc,"'for' var=<start> 'to' ...");

	/* parse the init value */
	tok=parse_expr(cc);
	if (tok!=T_TO) goto err;
	if (cc->result) init=getvalue(cc,cc->result);
	if (!(init && init->type==s_real))
		cc_error(cc,"<init> value must be real!");
	rv.value=*realof(init);
	
	/* parse the end value */
	tok=parse_expr(cc);
	if (!(tok==T_DO || tok==T_STEP || tok==T_SEMICOL || tok==T_COMMA || tok==T_EOS)) goto err;
	if (cc->result) end=getvalue(cc,cc->result);
	if (!(end && end->type==s_real))
		cc_error(cc,"<end> value must be real!");
	vend=*realof(end);
	
	/* parse the optional step */
	if (tok==T_STEP) {
		tok=parse_expr(cc);
		if (!(tok==T_DO || tok==T_SEMICOL || tok==T_COMMA || tok==T_EOS)) goto err;
		if (cc->result) step=getvalue(cc,cc->result);
		if (!(step && step->type==s_real))
			cc_error(cc,"<step> value must be real!");
		vstep=*realof(step);
	} else {
		vstep=1.0;
	}
	
	signum=(vstep>0);
	if (signum && rv.value>vend) { scan_end(cc); goto end; }
	else if (!signum && rv.value<vend) { scan_end(cc); goto end; }
	cc->newram=cc->endlocal;
	
	/* create the loop variable */
	kill_local(cc,name);
	cc->newram=cc->endlocal;
	hd=new_reference(cc,&rv.hd,name);
	cc->endlocal=cc->newram=(char *)hd+hd->size;
	
	jump=cc->next;
	vend=vend+cc->epsilon*vstep;
	while (1) {
		h=parse(cc);
		if (h==c_endfunction) cc_error(cc,"'end' missing!");
		if (h==c_return) break;
		if (h==c_break) { scan_end(cc); break; }
		if (h==c_end) {
			rv.value+=vstep;
			if (signum==1 && rv.value>vend) break;
			else if (!signum && rv.value<vend) break;
			else cc->next=jump;
			if (sys_test_key()==escape) cc_error(cc,"function interrupted!");
		}
	}
end :
	kill_local(cc,name);
	cc->result=NULL;
	return;
err:
	cc_error(cc,"'for' var=<start> 'to' <end> 'step' <step> do ...");
}

static void do_loop (Calc *cc)
/***** do_loop
	do a loop command in a UDF.
	loop value to value; .... ; end
*****/
{
	token_t tok;
	int h;
	char *jump;
	header *init=NULL,*end=NULL;
	long vend,oldindex;
	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'loop' only allowed in functions!");
	/* preserve englobing loop state */
	oldindex=cc->loopindex;
	/* parse the init value */
	tok=parse_expr(cc);
	if (tok!=T_TO) cc_error(cc,"'loop' <init> 'to' <end> 'do' <body> 'end'");
	if (cc->result) init=getvalue(cc,cc->result);
	if (!(init && init->type==s_real))
		cc_error(cc,"<init> value must be real!");
	cc->loopindex=(long)*realof(init);
	/* parse the end value */
	tok=parse_expr(cc);
	if (!(tok==T_DO || tok==T_SEMICOL || tok==T_COMMA))
		cc_error(cc,"'loop' <init> 'to' <end> 'do' <body> 'end'");
	if (cc->result) end=getvalue(cc,cc->result);
	if (!(end && end->type==s_real))
		cc_error(cc,"<end> value must be real!");
	vend=(long)*realof(end);
	if (cc->loopindex>vend) { scan_end(cc); goto end; }
	cc->newram=cc->endlocal;
	jump=cc->next;
	while (1) {
		h=parse(cc);
		if (h==c_endfunction) cc_error(cc,"'end' missing!");
		if (h==c_return) break;
		if (h==c_break) { scan_end(cc); break; }
		if (h==c_end) {
			cc->loopindex++;
			if (cc->loopindex>vend) break;
			else cc->next=jump;
			if (sys_test_key()==escape) {
				cc_error(cc,"function interrupted!");
			}
		}
	}
end:
	/* restore englobing loop state */
	cc->loopindex=oldindex;
}

static void do_repeat (Calc *cc)
/***** do_loop
	do a loop command in a UDF.
	for value to value; .... ; endfor
*****/
{	int h;
	char *jump;
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'end' only allowed in functions!");
	cc->newram=cc->endlocal;
	jump=cc->next;
	while (1) {
		h=parse(cc);
		if (h==c_endfunction) cc_error(cc,"'end' missing!");
		if (h==c_return) break;
		if (h==c_break) { scan_end(cc); break; }
		if (h==c_end) {
			cc->next=jump;
			if (sys_test_key()==escape) {
				cc_error(cc,"function interrupted!");
			}
		}
	}
}


static void do_end (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF)) cc_error(cc,"'end' only allowed in functions!");
}

static void do_break (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'break' only allowed in functions!");
}

static void do_to (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'to' only allowed in functions!");
}

static void do_do (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'do' only allowed in functions!");
}

static void do_step (Calc *cc)
{
	if (!CC_ISSET(cc,CC_EXEC_UDF))
		cc_error(cc,"'step' only allowed in functions!");
}

/***************************************************************************
 * parser
 ***************************************************************************/
int xor (char *n);

static token_t parse_matrix(Calc *cc)
{
	token_t tok;
	header* hd;
	int r=0, c=0, col=0;
	real *m, *ms;
	hd=new_matrix(cc,0,0,"");
	ms=m=matrixof(hd);
	while (1) {
		tok=parse_expr(cc);
		if (cc->result) {
			if (cc->result->type==s_reference) cc->result=getvalue(cc,cc->result);
			switch (hd->type) {
			case s_matrix:
				switch (cc->result->type) {
				case s_real:
					if (!r || (r && !c) || (r && (c+1>col))) {
						stack_realloc(cc, hd, (r*col+(c+1>col?c+1:col))*sizeof(real));
					}
					*ms++=*realof(cc->result);
					break;
				case s_complex:	/* upgrade the matrix to complex */
					*ms++=*realof(cc->result); *ms++=*imagof(cc->result);
					{
						real *p=m+2*(r*col+c+1);
						stack_realloc(cc,hd,(r*col+(c+1>col?c+1:col))*2*sizeof(real));
						hd->type=s_cmatrix;
						*--p=*--ms;*--p=*--ms;	/* copy the cplx values from end */
						while (p!=m) {
							*--p = 0.0;			/* imag value */
							*--p=*--ms;			/* real value */
						}
						ms=m+2*(r*col+c+1);		/* set the ptr at the end */
					}
					break;
				default:
					cc_error(cc, "Illegal type for building matrix");
				}
				break;
			case s_cmatrix:
				if (!r || (r && !c) || (r && (c+1>col))) {
					stack_realloc(cc, hd, (r*col+(c+1>col?c+1:col))*2*sizeof(real));
				}
				switch (cc->result->type) {
				case s_real:
					*ms++=*realof(cc->result);*ms++=0.0;
					break;
				case s_complex:
					*ms++=*realof(cc->result);*ms++=*imagof(cc->result);
					break;
				default:
					cc_error(cc, "Illegal type for building matrix");
				}
				break;
			}
		
			switch (tok) {
			case T_COMMA:				/* new value on the row */
				c++;
				break;
			case T_SEMICOL:				/* new row */
			case T_RBRACKET:			/* ']' so finished */
				c++;
				if (r) {
					if (c<col) {		/*    incomplete row : add extra 0 to not defined values */
						for (int i=0; i<col-c;i++) {
							switch (hd->type) {
							case s_matrix:
								*ms++=0.0;
								break;
							case s_cmatrix:
								*ms++=0.0; *ms++=0.0;
								break;
							}
						}
					} else if (c>col) { /*    expand number of columns for all rows */
						switch (hd->type) {
						case s_matrix:
							{
								stack_realloc(cc,hd,(r+1)*c*sizeof(real));
								real *d=m+r*c, *o=m+r*col;
								memmove(d,o,c*sizeof(real));	// move last row to make room for others
								d-=c;o-=col;
								for (int k=1;k<r;k++) {			// move and add extra 0
									memmove(d,o,col*sizeof(real));
									for (real*p=d+col;p<d+c;p++) {
										*p=0.0;
									}
									d-=c;o-=col;
								}
								for (real*p=d+col;p<d+c;p++) {// add extra 0 for first line
									*p=0.0;
								}
								ms=m+(r+1)*c;
							}
							break;
						case s_cmatrix:
							{
								stack_realloc(cc,hd,(r+1)*c*2*sizeof(real));
								real *d=m+2*r*c, *o=m+2*r*col;
								memmove(d,o,c*2*sizeof(real));	// move last row to make room for others
								d-=2*c;o-=2*col;
								for (int k=1;k<r;k++) {			// move and add extra 0
									memmove(d,o,col*2*sizeof(real));
									for (real*p=d+2*col;p<d+2*c;p++) {
										*p=0.0;
									}
									d-=2*c;o-=2*col;
								}
								for (real*p=d+2*col;p<d+2*c;p++) {// add extra 0 for first line
									*p=0.0;
								}
								ms=m+(r+1)*c;
							}
							break;
						}
						col=c;
					}
				}
				if (!col) col=c;
				r++; c=0;
				break;
			default:
				cc_error(cc, "Bad value or row delimiter or closing bracket");
			}
		}
		
		if (tok==T_RBRACKET) {
			dims *d=dimsof(hd);
			if (col) {
				d->r=r; d->c=col;	/* well formed matrix */
			} else if (r==0 && c==0) {
				d->r=1; d->c=0;		/* empty matrix */
			} else {
				d->r=1; d->c=c;		/* ill formed vector, but ok */
			}
//			d->r= (r>0) ? r : 1; d->c=col;
			cc->result = hd;
			cc->newram=(char*)hd+hd->size; /* set the matrix the last object 
											  on the stack */
			break;
		}
	}
	return tok;
}

header* get_mat_elt (Calc *cc, header *var)
/* get an element of a matrix in the form: var[i,j] */
{	header *st=(header*)cc->newram,*result,*hd,*hd1;
	token_t tok;

	while (var->type==s_reference) var=referenceof(var);
	CC_SET(cc,CC_PARSE_INDEX);
	tok=parse_expr(cc);
	hd=cc->result;
	if (tok==T_COMMA) {	/* two indexes: var[r,c] */
		tok=parse_expr(cc);
		hd1=cc->result;
	} else {			/* only one index: var[r] or var[c] for row vectors */
		if (dimsof(var)->r==1) {
			hd1=hd; hd=new_real(cc,1.0,"");
		} else {
			hd1=new_command(cc,c_allv);
		}
	}
	CC_UNSET(cc,CC_PARSE_INDEX);
	if (tok!=T_RBRACKET) {
		cc_error(cc,"Closing ']' missing");
	}
	
	/* else, get an element of a variable */
	hd=getvalue(cc,hd);
	hd1=getvalue(cc,hd1);
	if (var->type==s_matrix || var->type==s_real) {
		result=new_submatrix(cc,var,hd,hd1,"");
	} else if (var->type==s_cmatrix || var->type==s_complex) {
		result=new_csubmatrix(cc,var,hd,hd1,"");
	} else {
		cc_error(cc,"can't use [] to access this variable type!");
	}

	return moveresult(cc,st,result);
}


header* get_mat_elt1 (Calc *cc, header *var)
/* get an element of a matrix in the form: var{i} */
{	header *st=(header*)cc->newram,*result;
	token_t tok;
	int n,l;
	int r,c;
	real *m;
	while (var->type==s_reference) var=referenceof(var);
	tok=parse_expr(cc);
	cc->result=getvalue(cc,cc->result);
	if (!(cc->result && cc->result->type==s_real)) cc_error(cc,"Index must be a number!");
	if (tok!=T_RBRACE) cc_error(cc,"Closing '}' missing");
	
	if (var->type==s_real) {
		result=new_reference(cc,var,"");
	} else if (var->type==s_complex) {
		result=new_reference(cc,var,"");
	} else if (var->type==s_matrix) {
		getmatrix(var,&r,&c,&m);
		l=(int)(*realof(cc->result));
		n=r*c;
		if (n==0) cc_error(cc,"Matrix is empty!");
		if (l>n) l=n;
		if (l<1) l=1;
		if (CC_ISSET(cc,CC_NOSUBMREF)) result=new_real(cc,*(m+l-1),"");
		else result=new_subm(cc,var,l,"");
	} else if (var->type==s_cmatrix) {
		getmatrix(var,&r,&c,&m);
		l=(int)(*realof(cc->result));
		n=r*c;
		if (n==0) cc_error(cc,"Matrix is empty!");
		if (l>n) l=n;
		if (l<1) l=1;
		if (CC_ISSET(cc,CC_NOSUBMREF)) {
			m+=2*(l-1);
			result=new_complex(cc,*m,*(m+1),"");
		} else result=new_csubm(cc,var,l,"");
	} else {
		cc_error(cc,"can't use {} for this variable type!");
	}
	return moveresult(cc,st,result);
}

header* parse_list(Calc *cc)
{
	header *st=(header*)cc->newram;
	token_t tok;
	
	do {
		tok=parse_expr(cc);
		/* eliminate all references when returning multiple values
		   in user functions
		 */
		if (CC_ISSET(cc,CC_EXEC_RETURN)) {
			header *hd=getvalue(cc,cc->result);
			moveresult(cc,cc->result,hd);
		}

	} while (tok==T_COMMA);
	
	if (tok!=T_RBRACE) cc_error(cc,"Closing '}' missing!");
	
	return st;
}

header* parse_func_call(Calc *cc, char *name)
{
	header *st=(header*)cc->newram, *var=NULL, *res;
	token_t tok;
	unsigned int oldflags=cc->flags;
	int count=0;		/* classical value parameter counter*/
	int nocount=0;		/* boolean: 1 when a named parameter is given; do no count 
						   anymore, classical value parameter no more allowed */
	int epos=0;			/* index of the first extra parameter on the stack */
	int is_binfunc=1;	/* boolean: buitin/udf */
	int is_binfuncref=0;
	binfunc_t *func;
	char funcname[LABEL_LEN_MAX+1];
	
	/* make a copy of the name because it is, in fact cc->str, which may be
	   overwritten by subsequent calls to parse_expr */
	strcpy(funcname,name);
	
	if ((func=binfunc_find(funcname))!=NULL) {
		is_binfunc=1;
	}
	if (funcname[0]!='$' && (var=searchudf(cc,funcname))!=NULL) {
		is_binfuncref = (var->flags & FLAG_BINFUNC) ? 1 : 0;
		if (is_binfuncref) {
			strcpy(funcname,binfuncof(var)->name);
		}
#ifndef PRIO_TO_UDF
	}
#else
	/* better to search for strings in search udf */
	} else if ((var=searchvar(cc,cc->str))!=NULL) {
		if (var->type!=s_string) var=NULL;	/* not a string, so not a potential light user function */
	}
#endif
	if (!func && !var && !is_binfuncref) cc_error(cc, "Function '%s' not defined!", name);
	
	/* parse the parameter list, allow named parameters (which acts like local
	   variables for the function) */
	CC_SET(cc,CC_PARSE_PARAM_LIST|CC_NOSUBMREF);
//	CC_SET(cc,CC_PARSE_PARAM_LIST);
	do {
		tok=parse_expr(cc);
		if (tok==T_ASSIGN) {	/* a named parameter! */
			if (cc->result && cc->result->type==s_reference) {
				header*hd=cc->result;
				tok=parse_expr(cc);
				if (cc->result) {
					strcpy(cc->result->name,hd->name);
					cc->result->xor=hd->xor;
					moveresult(cc,hd,cc->result);
					nocount=1;
					epos++;
				}
			} else {
				cc_error(cc,"parameter name expected!");
			}
		} else if (nocount) {
			cc_error(cc,"unamed value parameter after named parameter not allowed!");
		} else if (cc->result) {
			/* count all the results provided by expression evaluation */
			header*hd=cc->result;
			while (hd<(header *)cc->newram) {
				if (!nocount) count++;
				epos++;
				hd=nextof(hd);
			}
			if (count>MAXARGS) cc_error(cc,"Too many arguments!");
		}
	} while (tok==T_COMMA);
	
	if (tok==T_SEMICOL) {
		do {
			tok=parse_expr(cc);
			if (tok==T_ASSIGN) {	/* a named parameter! */
				cc_error(cc,"only unamed parameters allowed here");
			}
			if (cc->result->type==s_reference && !referenceof(cc->result)) {
				cc_error(cc,"empty parameters not allowed here");
			}
		} while (tok==T_COMMA);
	}
	cc->flags=oldflags;

	if (tok!=T_RPAR) cc_error(cc,"Closing ')' missing!");
	
	cc->stack = st;		/* used by pushresults */
	
	if (var && var->type==s_string) {
		return interpret_luf(cc,var,st,count,epos);
	} else if (var && !is_binfuncref && funcname[0]!='$') {
		return interpret_udf(cc,var,st,count,epos);
	} else if ((is_binfunc || is_binfuncref) && (res=binfunc_exec(cc, funcname, count, st))!=NULL) {
		return res;
	} else cc_error(cc,"wrong number of parameters in binary function");
	return NULL;
}

#if ASSIGN_OP
#define shift_by_offset(hd,offset) ((header *)((char *)(hd)+offset))

static header *do_assignment(Calc *cc, header *hd)
{
	token_t tok;
	int varcount=0,rscount=0,i,j;
	header *variable[8],*rightside[8],*rs,*var=hd,*v,*mark, *st=hd;
	LONG offset,oldoffset,dif;
	char *oldendlocal;
	while (var<(header *)cc->newram) {
		if (var->type!=s_reference && var->type!=s_submatrixref
		&& var->type!=s_csubmatrixref) cc_error(cc,"Illegal assignment!");
		variable[varcount++]=var; var=nextof(var);
		if (varcount>=8) cc_error(cc,"Too many commas!");
	}
	
	CC_SET(cc,CC_NOSUBMREF); tok=parse_expr(cc); CC_UNSET(cc,CC_NOSUBMREF);
	if (tok==T_RBRACKET || tok==T_RBRACE || tok==T_RPAR) cc_error(cc,"Illegal separator: only ';', ',' or '\\n' allowed");
	/* count and note the values, that are assigned to the
	   variables */
	rscount=0; rs=cc->result;
	while (rs<(header *)cc->newram) {
		rightside[rscount++]=rs; rs=nextof(rs);
		if (rscount>=8) cc_error(cc,"Too many commas!");
	}
	
	/* cannot assign 2 values to 3 variables , e.g. */
	if (rscount>1 && rscount<varcount) cc_error(cc,"Illegal assignment: less results than variables to be assigned!");
	oldendlocal=cc->endlocal;
	offset=0;
	/* do all the assignments */
	if (varcount==1) {
		/* simple assignment */
		header *var=variable[0];
		if (var->type==s_reference && referenceof(var)) {
			header *h=referenceof(var);
			if (h->type==s_reference || h->type==s_submatrixref || h->type==s_csubmatrixref) {
				/* left value is a reference to a variable or a submatrix belonging to the
				   caller context. make it a fake type so that the variable will be replaced
				 */
				h->type=s_real;
			}
		}
		hd=assign(cc,variable[0],rightside[0]);
	} else {
		for (i=0; i<varcount; i++) {
			header *var=variable[i];
			oldoffset=offset;
			if (var->type==s_reference && referenceof(var)) {	
				header *h=referenceof(var);
				if (h->type==s_reference || h->type==s_submatrixref || h->type==s_csubmatrixref) {
					/* left value is a reference to a variable or a submatrix belonging to the
					   caller context. make it a fake type so that the variable will be replaced
					 */
					h->type=s_real;
				}
			}
			/* assign a variable */
			var=assign(cc,shift_by_offset(variable[i],offset),
				shift_by_offset(rightside[(rscount>1)?i:0],offset));
			if (i==0) hd=var;
			offset=cc->endlocal-oldendlocal;
			if (oldoffset!=offset) {
				/* size of var changed */
				v=shift_by_offset(variable[i],offset);
				if (v->type==s_reference) mark=referenceof(v);
				else mark=submrefof(v);
				/* now shift all references of the var.s */
				if (mark) {
					/* not a new variable */
					for (j=i+1; j<varcount; j++) {
						v=shift_by_offset(variable[j],offset);
						dif=offset-oldoffset;
						if (v->type==s_reference && referenceof(v)>mark)
							referenceof(v)=shift_by_offset(referenceof(v),dif);
						else if (submrefof(v)>mark)
							submrefof(v)=shift_by_offset(submrefof(v),dif);
					}
				}
			}
		}
	}
	return hd;
}
#endif

/* precedence (p) up to 15, is binary operator (b)
   is right associative (r), is ending symbol (e)
 */
#define P(e,r,b,p)	(0x80 | ((e)<<6) | ((r)<<5) | ((b)<<4) | (p))

typedef header* (*op_func_t)();
typedef header* (*op_func1_t)(Calc *cc, header *hd);
typedef header* (*op_func2_t)(Calc *cc, header *hd1, header *hd2);

static const struct Operator {
	unsigned int	flag;
	op_func_t		func;
} ops[] = {
	{0,	NULL},								/* T_NONE not used */
	{P(0,0,1,8), (op_func_t)add},			/* T_ADD, '+' */
	{P(0,0,1,8), (op_func_t)subtract},		/* T_SUB, '-' */
	{P(0,0,1,9), (op_func_t)dotmultiply},	/* T_MPY, '*' */
	{P(0,0,1,9), (op_func_t)dotdivide},		/* T_DIV, '/' */
	{P(0,1,1,10), (op_func_t)mpower},		/* T_POW, '^' (right-associative) */
	{P(0,0,1,9), (op_func_t)multiply},		/* T_DOT, '.' */
	{P(0,0,1,8), (op_func_t)msolve},		/* T_SOLVE, '\\' */
	{P(0,0,0,2), NULL},						/* T_LPAR, '(' */
	{P(0,0,0,2), NULL},						/* T_RPAR, ')' */
	{P(0,0,0,8), (op_func_t)opposite},		/* T_NEG, '-' (unary minus) */
	{P(0,0,0,2), (op_func_t)NULL},			/* T_LBRACKET, '[' */
	{P(1,0,0,2), NULL},						/* T_RBRACKET, ']' */
	{P(0,0,0,2), NULL},						/* T_LBRACE, '{' */
	{P(1,0,0,2), NULL},						/* T_RBRACE, '}' */
	{P(0,0,1,5), (op_func_t)mequal},		/* T_EQ, '==' */
	{P(0,0,1,5), (op_func_t)maboutequal},	/* T_ABOUTEQ, '~=' */
	{P(0,0,1,5), (op_func_t)munequal},		/* T_NE, '!=' */
	{P(0,0,1,5), (op_func_t)mgreater},		/* T_GT, '>' */
	{P(0,0,1,5), (op_func_t)mgreatereq},	/* T_GE, '>=' */
	{P(0,0,1,5), (op_func_t)mless},			/* T_LT, '<' */
	{P(0,0,1,5), (op_func_t)mlesseq},		/* T_LE, '<=' */
	{P(0,0,1,4), (op_func_t)mand},			/* T_AND, '&&' */
	{P(0,0,1,3), (op_func_t)mor},			/* T_OR, '||' */
	{P(0,1,0,10), (op_func_t)mnot},			/* T_NOT, '!' (right-associative) */
	{P(0,0,1,6), (op_func_t)mhconcat},		/* T_HCONCAT, '|' */
	{P(0,0,1,6), (op_func_t)mvconcat},		/* T_VCONCAT, '_' */
	{P(0,0,0,10), (op_func_t)transpose},		/* T_TRANSPOSE, '\'' */
	{P(0,1,0,7), NULL},						/* T_COL, ':' (binary/ternary operator: special handling) */
#ifndef ASSIGN_OP
	{P(1,0,0,1), NULL},						/* T_ASSIGN, '=' */
#else
	{P(0,1,1,1), (op_func_t)do_assignment},	/* T_ASSIGN, '=' */
#endif
	{P(1,0,0,1), NULL},						/* T_EOS - End Of Input String */
	{P(1,0,0,1), NULL},						/* T_SEMICOL, ';' */
	{P(1,0,0,1), NULL},						/* T_COMMA, ',' */
	{P(1,0,0,1), NULL},						/* T_THEN 'then' */
	{P(1,0,0,1), NULL},						/* T_ELSE, 'else' */
	{P(1,0,0,1), NULL},						/* T_ELSIF, 'elsif' */
	{P(1,0,0,1), NULL},						/* T_ENDIF, 'endif' */
	{P(1,0,0,1), NULL},						/* T_TO, 'to' */
	{P(1,0,0,1), NULL},						/* T_DO, 'do' */
	{P(1,0,0,1), NULL},						/* T_END, 'end' */
	{P(1,0,0,1), NULL},						/* T_STEP, 'step' */
};

#define PREC(op)	(ops[op].flag & 0xF)		/* precedence: priority grows higher */
#define IS_BIN(op)	(ops[op].flag & 0x10)		/* is binary operator */
#define IS_RASS(op)	(ops[op].flag & 0x20)		/* is right associative, else left associative */
#define IS_END(op)	(ops[op].flag & 0x40)		/* is ending symbol */
#define IS_OP(op)	(ops[op].flag & 0x80)		/* is an operator */

/* max operand stack to evaluate an expression */
#define OP_STACK_MAX		20
#define DATA_STACK_MAX		20

token_t parse_expr(Calc *cc)
{
	binfunc_t *fn;
	int usign=0;
	token_t tok;

	int      d_top=-1;				/* top of the data stack */
	int      o_top=0;				/* top of the operand stack */
	header*  data[DATA_STACK_MAX];	/* data stack */
	token_t  op[OP_STACK_MAX]={0};	/* operand stack */

	while (1) {
		/* get an operand */
		tok=scan(cc);
		
		switch (tok) {
		case T_REAL:
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=new_real(cc,cc->val,"");
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_IMAG:
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=new_complex(cc,0.0,cc->val,"");
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_LBRACKET:
			tok=parse_matrix(cc);
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=cc->result;
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_RBRACKET:	/* only when [] (empty matrix) */
			cc->result=NULL;
			return tok;
		case T_LBRACE:
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=parse_list(cc);
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_SUB:
			if (!usign) {
				usign=1;
				if (o_top<OP_STACK_MAX-1) {
					op[++o_top] = T_NEG;
				} else {
					cc_error(cc, "Operator stack overflow"); goto err;
				}
				continue;
			} else {
				goto bad_operand;
			}
			break;
		case T_ADD:
			if (!usign) {
				usign=1;
				continue;
			} else {
				goto bad_operand;
			}
			break;
		case T_NOT:
			if (o_top<OP_STACK_MAX-1) {
				op[++o_top] = tok;
			} else {
				cc_error(cc, "Operator stack overflow"); goto err;
			}
			usign=0;
			continue;
		case T_LPAR:
			if (o_top<OP_STACK_MAX-1) {
				op[++o_top] = tok;
			} else {
				cc_error(cc, "Operator stack overflow"); goto err;
			}
			usign=0;
			continue;
		case T_LABEL:
		{
			header* var;
			if ((var=searchvar(cc,cc->str))!=NULL) {
				/* variable exists, push a reference to it (will be
				   dereferenced by getvalue */
				data[++d_top] = new_reference(cc,var,cc->str);
			} else if ((var=searchudf(cc,cc->str))!=NULL) {
				data[++d_top] = new_funcref(cc,var,"");
			} else if ((fn=binfunc_find(cc->str))!=NULL) {
				data[++d_top] = new_binfuncref(cc,fn,"");
			} else {
				/* variable does not exist, push a NULL reference
				   (will be used in assignment when back to parse()) */
				data[++d_top]=new_reference(cc,NULL,cc->str);
			}
			usign=0;
			break;
		}
		case T_COL:
			if (CC_ISSET(cc,CC_PARSE_INDEX)) {
				if (d_top<DATA_STACK_MAX-1) {
					data[++d_top]=new_command(cc,c_allv);
					if (op[o_top]==T_NEG) cc_error(cc,"Bad index '-:'");
				} else {
					cc_error(cc, "reg file overflow"); goto err;
				}
			} else cc_error(cc,"':' only allowed as matrix index or to generate a vector");
			break;
		case T_FUNCREF:
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=parse_func_call(cc,cc->str);
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_RPAR:
			if (CC_ISSET(cc,CC_PARSE_PARAM_LIST)) {
				cc->result=NULL;
				return tok;
			} else goto bad_operand;
		case T_COMMA:
			if (CC_ISSET(cc,CC_PARSE_PARAM_LIST)) {
				cc->result=new_reference(cc,NULL,"");
				return tok;
			} else goto bad_operand;
		case T_MATREF:		/* var[i] */
			if (d_top<DATA_STACK_MAX-1) {
				header *var=searchvar(cc,cc->str);
				if (!var) cc_error(cc,"%s not a variable!",cc->str);
				data[++d_top]=get_mat_elt(cc,var);
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_MATREF1:		/* var{i} */
			if (d_top<DATA_STACK_MAX-1) {
				header *var=searchvar(cc,cc->str);
				if (!var) cc_error(cc,"%s not a variable!",cc->str);
				data[++d_top]=get_mat_elt1(cc,var);
			} else {
				cc_error(cc, "Data stack overflow"); goto err;
			}
			usign=0;
			break;
		case T_DQUOTE:
			if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=new_string(cc);
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		case T_HASH:	/* # (loopindex) */
				if (d_top<DATA_STACK_MAX-1) {
				data[++d_top]=new_real(cc,(real)cc->loopindex,"");
			} else {
				cc_error(cc, "Reg file overflow"); goto err;
			}
			usign=0;
			break;
		default:
		bad_operand:
			/* bad syntax */
			cc_error(cc, "Bad statement or lacks an operand!");
			goto err;
		}
		
		/* get an operator */
expect_operator:
		tok=scan(cc);
		if (!IS_OP(tok)) {
			/* bad syntax */
			cc_error(cc, "Bad statement or lacks an operator!");
			break;
		} else if (tok==T_LPAR) {
			cc_error(cc,"'*' operator missing?");
		}

		while (PREC(tok)<=PREC(op[o_top])) {
			/* reduce the stack by executing ops of greater 
			   precedence than tok */
			if (tok==T_ASSIGN) {
				cc_error(cc,"assignment: illegal op in expression");
			}
			if (op[o_top]==T_LPAR) {
				if (tok==T_RPAR) {
//					op[o_top]=T_NONE;
					o_top--;
					goto expect_operator;
				} else if (tok==T_EOS) {
					cc_error(cc, "Missing ')'"); goto err;
				}
			}
			/* deal with right-associative operators */
			if (IS_RASS(tok) && op[o_top]==tok) break;
			
			/* execute operator defined in table ops[] */
			header *b=data[d_top--], *a;
			
			if (op[o_top]==T_COL) {
				header* step;
				if (o_top && op[o_top-1]==T_COL) {
					step=data[d_top--];
					o_top-=2;
				} else {
					step=new_real(cc,1.0,"");
					o_top-=1;
				}
				a=data[d_top--];
				cc->stack=a;
				a=vectorize(cc,a,step,b);
			} else if (IS_BIN(op[o_top])) {
				a=data[d_top--];
				cc->stack=a;
				a=((op_func2_t)ops[op[o_top--]].func)(cc,a,b);
			} else {
				cc->stack=b;
				a=((op_func1_t)ops[op[o_top--]].func)(cc,b);
			}

			/* put the result on the data stack */
			data[++d_top] = a;
		}
		if (!CC_ISSET(cc,CC_PARSE_PARAM_LIST) && tok==T_RPAR) { cc_error(cc, "Missing '('"); goto err; }
//		if (!(tok==T_EOS || tok==T_SEMICOL || tok==T_COMMA || tok==T_RBRACKET  || tok==T_RBRACE || tok==T_ASSIGN || (tok==T_RPAR && CC_ISSET(cc,CC_PARSE_PARAM_LIST)))) {
		if (!(IS_END(tok) || (tok==T_RPAR && CC_ISSET(cc,CC_PARSE_PARAM_LIST)))) {
			/* push tok operator on stack */
			if (o_top<OP_STACK_MAX-1) {
				if (tok==T_COL && op[o_top]==T_COL && (o_top-1) && op[o_top-1]==T_COL) {
					cc_error(cc, "Too many ':' for vector generation"); goto err;
				}
				if (tok==T_LBRACKET) { /* index result as a matrix */
					header *var=data[d_top--];
					data[++d_top]=get_mat_elt(cc,var);
					goto expect_operator;
				} else if (tok==T_LBRACE) { /* index result as a matrix */
					header *var=data[d_top--];
					data[++d_top]=get_mat_elt1(cc,var);
					goto expect_operator;
				} else {
					op[++o_top]=tok;
				}
				if (tok==T_TRANSPOSE) goto expect_operator;
			} else {
				cc_error(cc, "Operator stack overflow"); goto err;
			}
		} else {
			/* finished, return the result */
			cc->result = data[0];
			break;
		}
	}
	
err:
	return tok;
}

#ifndef ASSIGN_OP 
#define shift_by_offset(hd,offset) ((header *)((char *)(hd)+offset))
#endif
int parse(Calc *cc)
{	
	token_t tok;
	int cmd=-1;
	
	while(1) {
		char c;
		/* skip space */
		while ((c=*cc->next)=='\t' || c==' ' || c==';') cc->next++;	
		/* and comments */
		while ((*cc->next=='#' && *(cc->next+1)=='#') || (*cc->next=='.' && *(cc->next+1)=='.')) next_line(cc);
		if (*cc->next && *cc->next!='\n') break;
		else next_line(cc);
	}
	
	if ((cmd=cmd_parse_and_exec(cc))!=c_none) {
		cc->newram=cc->endlocal;
		return cmd;
	}
	tok=parse_expr(cc);
	switch (tok) {
	case T_END: cmd=c_end; break;
	case T_ENDIF: cmd=c_endif; break;
	default: break;
	}
	if (!cc->result) return cmd;
	switch (cc->result->type) {
	case s_real:
	case s_complex:
	case s_matrix:
	case s_cmatrix:
	case s_string:
	case s_funcref:
		switch (tok) {
		case T_EOS:
		case T_COMMA:
		case T_ENDIF:
		case T_END:
			if (cc->result) give_out(cc, cc->result);
			break;
		case T_SEMICOL:
			break;
		case T_ASSIGN:
			if (cc->result->type==s_funcref)
				cc_error(cc,"function names have const status (can't write)!");
			else
				cc_error(cc,"illegal assignment");
			break;
		default:
			cc_error(cc,"bad statement!");
			break;
		}
		break;
		if (tok==T_RBRACKET || tok==T_RBRACE || tok==T_RPAR) cc_error(cc,"Illegal separator: only ';', ',' or '\\n' allowed");
		if (tok!=T_SEMICOL) give_out(cc, cc->result);
		break;
	case s_reference:
	case s_submatrixref:
	case s_csubmatrixref:
		if (tok==T_ASSIGN) {
			int varcount=0,rscount=0,i,j;
			header *variable[8],*rightside[8],*rs,*var=cc->result,*v,*mark;
			LONG offset,oldoffset,dif;
			char *oldendlocal;
			while (var<(header *)cc->newram) {
				if (var->type!=s_reference && var->type!=s_submatrixref
				&& var->type!=s_csubmatrixref) cc_error(cc,"Illegal assignment!");
				variable[varcount++]=var; var=nextof(var);
				if (varcount>=8) cc_error(cc,"Too many commas!");
			}
			
			CC_SET(cc,CC_NOSUBMREF); tok=parse_expr(cc); CC_UNSET(cc,CC_NOSUBMREF);
			if (tok==T_RBRACKET || tok==T_RBRACE || tok==T_RPAR) cc_error(cc,"Illegal separator: only ';', ',' or '\\n' allowed");
			/* count and note the values, that are assigned to the
			   variables */
			rscount=0; rs=cc->result;
			while (rs<(header *)cc->newram) {
				rightside[rscount++]=rs; rs=nextof(rs);
				if (rscount>=8) cc_error(cc,"Too many commas!");
			}
			
			/* cannot assign 2 values to 3 variables , e.g. */
			if (rscount>1 && rscount<varcount) cc_error(cc,"Illegal assignment: less results than variables to be assigned!");
			/* do all the assignments */
			if (varcount==1) {
				/* simple assignment */
				cc->result=assign(cc,variable[0],rightside[0]);
			} else {
#if 0
				for (i=0; i<varcount; i++) {
					char* oldnewram=cc->newram;
					oldendlocal=cc->endlocal;
					header *var=variable[i];
					if (var->type==s_reference && referenceof(var)) {	
						header *h=referenceof(var);
						if (h->type==s_reference || h->type==s_submatrixref || h->type==s_csubmatrixref) {
							/* left value is a reference to a variable or a submatrix belonging to the
							   caller context. make it a fake type so that the variable will be replaced
							 */
							h->type=s_real;
						}
					}
					/* assign a variable */
					var=assign(cc,variable[i],rightside[i]);
					if (i==0) cc->result=var;
					offset=cc->newram-oldnewram;
//					offset=cc->endlocal-oldendlocal;
					if (offset) {
						/* size of var changed, offset all references above var */
						for (j=i+1; j<varcount; j++) {
							/* shift references to var and submarices if necessary */
							v=variable[j]=shift_by_offset(variable[j],offset);
							if (v->type==s_reference && referenceof(v)>=var)
								referenceof(v)=shift_by_offset(referenceof(v),offset);
							else if (submrefof(v)>var)
								submrefof(v)=shift_by_offset(submrefof(v),offset);
							/* shift references in rightsides */
							v=rightside[j]=shift_by_offset(rightside[j],offset);
							if (v->type==s_reference && referenceof(v)>=var)
								referenceof(v)=shift_by_offset(referenceof(v),offset);
						}
					}
				}
#else
				offset=0;
				oldendlocal=cc->endlocal;
				for (i=0; i<varcount; i++) {
					header *var=variable[i];
					oldoffset=offset;
					/* assign a variable */
					var=assign(cc,shift_by_offset(variable[i],offset),
						shift_by_offset(rightside[(rscount>1)?i:0],offset));
					if (i==0) cc->result=var;
					offset=cc->endlocal-oldendlocal;
					if (oldoffset!=offset) {
						/* size of var changed */
						v=shift_by_offset(variable[i],offset);
						if (v->type==s_reference) mark=referenceof(v);
						else mark=submrefof(v);
						/* now shift all references of the var.s */
						if (mark) {
							/* not a new variable */
							for (j=i+1; j<varcount; j++) {
								v=shift_by_offset(variable[j],offset);
								dif=offset-oldoffset;
								if (v->type==s_reference && referenceof(v)>mark)
									referenceof(v)=shift_by_offset(referenceof(v),dif);
								else if (submrefof(v)>mark)
									submrefof(v)=shift_by_offset(submrefof(v),dif);
							}
						}
					}
				}

#endif
			}
		} else {
			/* expression is just a reference to a variable: expand reference */
			cc->result=getvalue(cc,cc->result);
			if (tok==T_RBRACKET || tok==T_RBRACE || tok==T_RPAR) cc_error(cc,"Illegal separator: only ';', ',' or '\\n' allowed");
		}
		switch (tok) {
		case T_EOS:
		case T_COMMA:
		case T_ENDIF:
		case T_END:
			if (cc->result) give_out(cc, cc->result);
			break;
		case T_SEMICOL:
			break;
		default:
			cc_error(cc,"bad statement");
			break;
		}
		break;
	default:
		break;
	}
	
	cc->newram=cc->endlocal;	/* remove the result from the calculus stack */
	
	return cmd;
}

/***************************************************************************
 * utils
 ***************************************************************************/

int extend(Calc *cc, char* start, char extend[16][LABEL_LEN_MAX+1])
/* extend
Extend a start string in up to 16 ways to a command or function.
This function is called from the line editor, whenever the HELP
key is pressed.
*/
{	int count=0,ln,l;
	header *hd=(header *)cc->udfstart;
#ifndef EMBED
	binfunc_t *b=binfunc_list;
	cmd_t *c=cmd_list;
#else
	binfunc_t const *b=binfunc_list;
	cmd_t const *c=cmd_list;
#endif
	ln=strlen(start);
	for (int i=0;i<BINFUNCS;i++) {
		if (!strncmp(start,b[i].name,ln)) {
			l=(int)strlen(b[i].name)-ln;
			if (l>0 && l<16) {
				strcpy(extend[count],b[i].name+ln);
				count++;
			}
			if (count>=16) return count;
		}
	}
	while (hd<(header *)cc->udfend) {
		if (!strncmp(start,hd->name,ln)) {
			l=(int)strlen(hd->name)-ln;
			if (l>0 && l<16) {
				strcpy(extend[count],hd->name+ln);
				count++;
			}
			if (count>=16) return count;
		}
		hd=nextof(hd);
	}
	for (int i=0;i<CMDS;i++) {
		if (!strncmp(start,c[i].name,ln)) {
			l=(int)strlen(c[i].name)-ln;
			if (l>0 && l<16) {
				strcpy(extend[count],c[i].name+ln);
				count++;
			}
			if (count>=16) return count;
		}
	}
	return count;
}

void cc_warn(Calc *cc, char *fmt, ...)
{
	char text [256];
	va_list v;
	int i=0;
	text_mode();
	sys_out_mode(CC_WARN);
	va_start(v,fmt);
	if (CC_ISSET(cc,CC_EXEC_UDF)) {
		output(cc,"  ");
		type_udfline(cc,cc->line);
	} else {
		if (cc->infile) outputf(cc,"  %s\n",cc->line);
		for (i=0;i<cc->next-cc->line+1;i++)
			text[i]=' ';
		text[i++]='^'; text[i++]='\n';
	}
	strncpy(text+i,"  warning: ",255-i);
	i+=11;
	i+=vsnprintf(text+i,255-i,fmt,v);
	text[i++]='\n';
	text[i]=0;
	va_end(v);
	output(cc,text);
}

void cc_error(Calc *cc, char *fmt, ...)
{
	char text [256];
	va_list v;
	int i=0;
	text_mode();
	sys_out_mode(CC_ERROR);
	va_start(v,fmt);
	if (CC_ISSET(cc,CC_EXEC_UDF)) {
		output(cc,"  ");
		type_udfline(cc,cc->line);
	} else {
		if (cc->infile) outputf(cc,"  %s\n",cc->line);
		for (i=0;i<cc->next-cc->line+1;i++)
			text[i]=' ';
		text[i++]='^'; text[i++]='\n';
	}
	i+=vsnprintf(text+i,255-i,fmt,v);
	text[i++]='\n';
	text[i]=0;
	va_end(v);
	output(cc,text);
	longjmp(*cc->env, 1);
}


void main_loop (Calc *cc, int argc, char *argv[])
{
	jmp_buf env;
	char input[LINEMAX]="";
	int i;
	header *hd;
	cc->globalstart=cc->globalend=cc->ramstart;
	cc->newram=cc->startlocal=cc->endlocal=cc->ramstart;
	cc->udfstart=cc->udfend=cc->ramend;
	cc->epsilon=EPSILON;
	cc->xstart=NULL;
	cc->xend=NULL;
	clear_fktext();
#ifndef EMBED
	binfunc_init(); cmd_init(); 
#endif
	cc->line = cc->next = input;	/* setup input line */
	cc->result = NULL;
	/* setup formats */
	cc->disp_mode=0;
	cc->disp_digits=6;
	cc->disp_fieldw=14;
	cc->disp_eng_sym=0;
	cc->maxexpo=1.0e6,
	cc->minexpo=1.0e-4;
	strcpy(cc->expoformat,"%0.4E");
	strcpy(cc->fixedformat,"%0.5G");
	/* setup display mode */
	cc->flags=CC_OUTPUTING;
	cc->infile = NULL;
	cc->outfile = NULL;
	cc->env=&env;
	
	hd=new_real(cc,M_PI,"pi");
	hd->flags=FLAG_CONST;
//	hd=new_real(cc,exp(1),"E");
//	hd->flags=FLAG_CONST;
	cc->globalend=cc->endlocal=cc->newram;
	
	strcpy(cc->line,"load first;");
	for (i=1; i<argc; i++)
	{	strncat(cc->line," load ",LINEMAX-1);
		strncat(cc->line,argv[i],LINEMAX-1);
		strncat(cc->line,";",LINEMAX-1);
	}

	switch (setjmp(env)) {
	case 0:
		/* register restore point */
		break;
	case 1:
		/* output(cc,"exception from global context\n"); */
		cc->line = cc->next = input;
		input[0]=0;
		break;
	case 2:
		/* output(cc,"back to global env from exception in function context\n"); */
		cc->line = cc->next = input;
		input[0]=0;
		break;
	}
	
	cc->result = NULL;
	cc->flags=CC_OUTPUTING;
	cc->loopindex=0;
	cc->level=0;
	
	/* interpret until "quit" */
	while (!cc->quit) {
		/* reset global context for commands evaluated in the 
	       lower level context (global scope) */
	    CC_UNSET(cc,CC_NOSUBMREF|CC_PARSE_INDEX|CC_PARSE_PARAM_LIST|CC_PARSE_UDF);
		cc->newram=cc->globalend=cc->endlocal;
		parse(cc);
		if (cc->trace<0) cc->trace=0;
	}
}
