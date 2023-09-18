#ifndef DSP_H
#define DSP_H

#include "calc.h"

/* FFT */
header* mfft (Calc *cc, header *hd);
header* mifft (Calc *cc, header *hd);

/* accelerometer */
header* maccel (Calc *cc, header *hd);

/* power quad */
header* mpqcos(Calc* cc, header* hd);
header* mpqfft(Calc* cc, header* hd);
header* mpqifft(Calc* cc, header* hd);

#endif
