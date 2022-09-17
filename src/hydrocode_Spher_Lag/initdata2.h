#ifndef _INITDTATA_H
#define _INITDTATA_H


#define Ncell 6000 // Number of computing cells in r direction
#define Tcell 400  // Number of computing cells in theta direction
#define Diaph1 (10.)
#define Diaph2 (10.2)
#define Diaph3 (10.25) // Domain length
#define Domlen (15) // Domain length

#define Tcell_plot 100 // Output zoom
#define Timeout (0.8) // Output time
#define D_PLOT_T (0.01) // Output time interval

#define Alpha (1.999) // GRP limiter parameter
#define LIMITER_CONF -2  /* LIMITER<0, add VIP limiter; LIMITER>0, only minmod limiter;
			   abs(LIMITER)=1, original minmod limiter; abs(LIMITER)=2, VIP-like minmod limiter */

#define GAMMAL (1.4)
#define GAMMAR (3.0) // Ratio of special heats Gamma=1.4 or 3.0
#define DL0 (0.129)
#define DR0 (19.237)
#define UL0 (0.)
#define UR0 (-10.)
#define PL0 (1.01325)
#define PR0 (1.01325)
#define DL1 (0.000129)
#define UL1 (0.)
#define PL1 (0.101325)

#define DATAOUT "../data_out2" //data out folder


#endif
