#include "fdtd-macro.h"
#include <stdio.h>
#include <complex.h>
#include <math.h>

/* update magnetic field */
void updateH(Grid *g) {
  int mm, nn, pp;
  complex double CurlE;

  if (Type == oneDGrid) {
    
    for (mm = 0; mm < SizeX; mm++)
      Hy1(mm) = Chyh1(mm) * Hy1(mm) 
	+ Chye1(mm) * (Ez1(mm + 1) - Ez1(mm));
    
  } else if (Type == tmZGrid) { 
    
    for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++)
	Hx2(mm, nn) = Chxh2(mm, nn) * Hx2(mm, nn) 
	  - Chxe2(mm, nn) * (Ez2(mm, nn + 1) - Ez2(mm, nn));
    
    for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++)
	Hy2(mm, nn) = Chyh2(mm, nn) * Hy2(mm, nn) 
	  + Chye2(mm, nn) * (Ez2(mm + 1, nn) - Ez2(mm, nn));
    
  } else if (Type == teZGrid) {
    
    for(mm = 0; mm < SizeX; mm++)
      for(nn = 0; nn < SizeY; nn++)
	Hz2(mm, nn) = Chzh2(mm, nn) * Hz2(mm, nn) -
	  Chze2(mm, nn) * ((Ey2(mm + 1, nn) - Ey2(mm, nn)) -
			   (Ex2(mm, nn + 1) - Ex2(mm, nn)));

  } else if (Type == threeDGrid) {
    /*Calculate the Hx field*/
    for (mm = 0; mm < SizeX; mm++){
      for (nn = 0; nn < SizeY - 1; nn++){
	for (pp = 0; pp < SizeZ - 1; pp++){
	  CurlE = Ey(mm, nn, pp + 1) - Ey(mm, nn, pp) -
				(Ez(mm, nn + 1, pp) - Ez(mm, nn, pp));
	  Hx(mm, nn, pp) = Chxh(mm, nn, pp) * Hx(mm, nn, pp) +
	    Chxe(mm, nn, pp) * CurlE;
    }
   }
  }
        
    for (mm = 0; mm < SizeX; mm++){
      for (pp = 0; pp < SizeZ -1; pp++){
	CurlE = Ey(mm, SizeY - 1, pp + 1) - Ey(mm, SizeY -1, pp) -
				(Phiy*Ez(mm, 0, pp) - Ez(mm, SizeY - 1, pp));
	  Hx(mm, SizeY - 1, pp) = Chxh(mm, SizeY - 1, pp) * Hx(mm, SizeY - 1, pp) +
	    Chxe(mm, SizeY - 1, pp) * CurlE;
      }
    }

    for (mm = 0; mm < SizeX; mm++){
      for (nn = 0; nn < SizeY -1; nn++){
	CurlE = (Phiz*Ey(mm, nn, 0) - Ey(mm, nn,SizeZ -1)) -
	        (Ez(mm, nn + 1, SizeZ -1) - Ez(mm, nn, SizeZ -1));
	  Hx(mm, nn, SizeZ -1) = Chxh(mm, nn, SizeZ -1) * Hx(mm, nn, SizeZ -1) +
	    Chxe(mm, nn, SizeZ -1) * CurlE;
      }
    }

    for (mm = 0; mm < SizeX; mm++){
      CurlE = (Phiz*Ey(mm, SizeY - 1, 0) - Ey(mm, SizeY - 1, SizeZ - 1)) -
				(Phiy*Ez(mm, 0, SizeZ - 1) - Ez(mm, SizeY - 1, SizeZ - 1));
	  Hx(mm, SizeY - 1, SizeZ - 1) = Chxh(mm,SizeY - 1, SizeZ - 1) * Hx(mm, SizeY - 1, SizeZ - 1) +
	    Chxe(mm, SizeY - 1, SizeZ - 1) * CurlE;
}
   
    /* Calculate the Hy field */	  

        
    for (mm = 0; mm < SizeX - 1; mm++){
      for (nn = 0; nn < SizeY; nn++){
        for (pp = 0; pp < SizeZ - 1; pp++){
	  CurlE = (Ez(mm + 1, nn, pp) - Ez(mm, nn, pp)) -
	    (Ex(mm, nn, pp + 1) - Ex(mm, nn, pp));
	  Hy(mm, nn, pp) = Chyh(mm, nn, pp) * Hy(mm, nn, pp) +
	    Chye(mm, nn, pp) *CurlE;
      }
    }
   }
    
    for (nn = 0; nn < SizeY; nn++){
       for (pp = 0; pp < SizeZ - 1; pp++){
	 CurlE = (Phix*Ez(0, nn, pp) - Ez(SizeX - 1, nn, pp)) -
	   (Ex(SizeX - 1, nn, pp + 1) - Ex(SizeX - 1, nn, pp));
	  Hy(SizeX - 1, nn, pp) = Chyh(SizeX - 1, nn, pp) * Hy(SizeX - 1, nn, pp) +
	    Chye(SizeX - 1, nn, pp) *CurlE;
	}
      }

    for (mm = 0; mm < SizeX - 1; mm++){
      for (nn = 0; nn < SizeY; nn++){
	  CurlE = (Ez(mm +1, nn, SizeZ - 1) - Ez(mm, nn, SizeZ - 1)) -
	    (Phiz*Ex(mm, nn, 0) - Ex(mm, nn, SizeZ - 1));
	  Hy(mm, nn, SizeZ - 1) = Chyh(mm, nn, SizeZ - 1) * Hy(mm, nn, SizeZ - 1) +
	    Chye(mm, nn, SizeZ - 1)*CurlE;	
      }
    }

      for (nn = 0; nn < SizeY; nn++){
	CurlE = (Phix*Ez(0, nn, SizeZ - 1) - Ez(SizeX - 1, nn, SizeZ - 1)) -
	  (Phiz*Ex(SizeX - 1, nn, 0) - Ex(SizeX - 1, nn, SizeZ - 1));
	  Hy(SizeX - 1, nn, SizeZ - 1) = Chyh(SizeX - 1, nn, SizeZ - 1) * Hy(SizeX - 1, nn, SizeZ - 1) +
	    Chye(SizeX - 1, nn, SizeZ - 1) *CurlE;
      }

   /* Calculate the Hz field */   
   for (mm = 0; mm < SizeX - 1; mm++){
     for (nn = 0; nn < SizeY - 1; nn++){
       for (pp = 0; pp < SizeZ; pp++){
	 CurlE = (Ex(mm, nn + 1, pp) - Ex(mm, nn, pp)) -
	   (Ey(mm + 1, nn, pp) - Ey(mm, nn, pp));
	  Hz(mm, nn, pp) = Chzh(mm, nn, pp) * Hz(mm, nn, pp) +
	    Chze(mm, nn, pp) * CurlE; 
       }
     }
   }


     for (nn = 0; nn < SizeY - 1; nn++){
       for (pp = 0; pp < SizeZ; pp++){
	 CurlE = (Ex(SizeX - 1, nn + 1, pp) - Ex(SizeX - 1, nn, pp)) -
	   (Phix*Ey(0, nn, pp) - Ey(SizeX - 1, nn, pp));
	  Hz(SizeX - 1, nn, pp) = Chzh(SizeX - 1, nn, pp) * Hz(SizeX - 1, nn, pp) +
	    Chze(SizeX - 1, nn, pp) * CurlE; 
       }
     }

   for (mm = 0; mm < SizeX - 1; mm++){
     for (pp = 0; pp < SizeZ; pp++){
       CurlE = (Phiy*Ex(mm, 0, pp) - Ex(mm, SizeY  - 1, pp)) -
	 (Ey(mm + 1, SizeY  - 1, pp) - Ey(mm, SizeY  - 1, pp));
	  Hz(mm, SizeY  - 1, pp) = Chzh(mm, SizeY  - 1, pp) * Hz(mm, SizeY  - 1, pp) +
	    Chze(mm, SizeY  - 1, pp) * CurlE; 
     }
   }

   for (pp = 0; pp < SizeZ; pp++){
     CurlE = (Phiy*Ex(SizeX - 1, 0, pp) - Ex(SizeX - 1, SizeY  - 1, pp)) -
       (Phix*Ey(0, SizeY - 1, pp) - Ey(SizeX - 1, SizeY - 1, pp));
	  Hz(SizeX - 1, SizeY  - 1, pp) = Chzh(SizeX - 1, SizeY  - 1, pp) * Hz(SizeX - 1, SizeY  - 1, pp) +
	    Chze(SizeX - 1, SizeY  - 1, pp) * CurlE; 
   }

  } else {
    fprintf(stderr, "updateH: Unknown grid type.  Terminating...\n");
  }
  
  return;
}  /* end updateH() */


/* update electric field */
void updateE(Grid *g) {
  int mm, nn, pp;
  complex double CurlH;
  complex double Psix = conj(Phix);
  complex double Psiy = conj(Phiy);
  complex double Psiz = conj(Phiz);
  
  if (Type == oneDGrid) {
    
    for (mm = 1; mm < SizeX - 1; mm++)
      Ez1(mm) = Ceze1(mm) * Ez1(mm) 
	+ Cezh1(mm) * (Hy1(mm) - Hy1(mm - 1));
    
  } else if (Type == tmZGrid) {
    
    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
	Ez2(mm, nn) = Ceze2(mm, nn) * Ez2(mm, nn) +
	  Cezh2(mm, nn) * ((Hy2(mm, nn) - Hy2(mm - 1, nn)) -
			   (Hx2(mm, nn) - Hx2(mm, nn - 1)));

  } else if (Type == teZGrid) {
    
    for(mm = 1; mm < SizeX - 1; mm++)
      for(nn = 1; nn < SizeY - 1; nn++)
	Ex2(mm, nn) = Cexe2(mm, nn) * Ex2(mm, nn) +
	  Cexh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm, nn - 1));
    
    for(mm = 1; mm < SizeX - 1; mm++)
      for(nn = 1; nn < SizeY - 1; nn++)
	Ey2(mm, nn) = Ceye2(mm, nn) * Ey2(mm, nn) -
	  Ceyh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm - 1, nn));
    
  } else if (Type == threeDGrid) {

    /* Calculate Ex field */
    for (mm = 0; mm < SizeX; mm++){
      for (nn = 1; nn < SizeY; nn++){
	for (pp = 1; pp < SizeZ; pp++){
	  CurlH = (Hz(mm, nn, pp) - Hz(mm, nn - 1, pp)) -
	    (Hy(mm, nn, pp) - Hy(mm, nn, pp - 1));
	  Ex(mm, nn, pp) = Cexe(mm, nn, pp) * Ex(mm, nn, pp) +
	    Cexh(mm, nn, pp) * CurlH;
	}
      }
    }

    for (mm = 0; mm < SizeX; mm++){
      for( pp = 1; pp < SizeZ; pp++){
	CurlH = (Hz(mm, 0, pp) - Psiy*Hz(mm, SizeY -1, pp)
		 - Hy(mm, 0, pp) + Hy(mm, 0, pp - 1));
	Ex(mm, 0, pp) = Cexe(mm, 0, pp)*Ex(mm, 0, pp) +
	  Cexh(mm, 0, pp)*CurlH;
      }
    }

    for(mm = 0; mm < SizeX; mm++){
      for(nn = 1; nn < SizeY; nn++){
	CurlH = (Hz(mm, nn, 0) - Hz(mm, nn - 1, 0)
		 - Hy(mm, nn, 0) + Psiz*Hy(mm, nn, SizeZ - 1));
	Ex(mm, nn, 0) = Cexe(mm, nn, 0)*Ex(mm, nn, 0) +
	  Cexh(mm, nn, 0)*CurlH;
      }
    }

    for(mm = 0; mm < SizeX; mm++){
      CurlH = (Hz(mm, 0, 0) - Psiy*Hz(mm, SizeY - 1, 0)
	       - Hy(mm, 0, 0) + Psiz*Hy(mm, 0, SizeZ - 1));
      Ex(mm, 0, 0) = Cexe(mm, 0, 0)*Ex(mm, 0, 0) +
	Cexh(mm, 0, 0)*CurlH;
    }

     /* Calculate EY field */
    
     for (mm = 1; mm < SizeX; mm++){
       for (nn = 0; nn < SizeY; nn++){
	 for (pp = 1; pp < SizeZ; pp++){
	  CurlH = (Hx(mm, nn, pp) - Hx(mm, nn, pp - 1)) -
	    (Hz(mm, nn, pp) - Hz(mm - 1, nn, pp));
	  Ey(mm, nn, pp) = Ceye(mm, nn, pp) * Ey(mm, nn, pp) + 
	    Ceyh(mm, nn, pp) * CurlH;
	 }
       }
     }

     for (mm = 1; mm < SizeX; mm++){
       for (nn = 0; nn < SizeY; nn++){
	 CurlH = (Hx(mm, nn, 0) - Psiz*Hx(mm, nn, SizeZ - 1)) -
	    (Hz(mm, nn, 0) - Hz(mm - 1, nn, 0));
	  Ey(mm, nn, 0) = Ceye(mm, nn, 0) * Ey(mm, nn, 0) + 
	    Ceyh(mm, nn, 0) * CurlH;
       }
     }

     for (nn = 0; nn < SizeY; nn++){
	 for (pp = 1; pp < SizeZ; pp++){
	   CurlH = (Hx(0, nn, pp) - Hx(0, nn, pp - 1)) -
	    (Hz(0, nn, pp) - Psix*Hz(SizeX - 1, nn, pp));
	  Ey(0, nn, pp) = Ceye(0, nn, pp) * Ey(0, nn, pp) + 
	    Ceyh(0, nn, pp) * CurlH;
	 }
     }

     for (nn = 0; nn < SizeY; nn++){
       CurlH = (Hx(0, nn, 0) - Psiz*Hx(0, nn, SizeZ - 1)) -
	    (Hz(0, nn, 0) - Psix*Hz(SizeX - 1, nn, 0));
	  Ey(0, nn, 0) = Ceye(0, nn, 0) * Ey(0, nn, 0) + 
	    Ceyh(0, nn, 0) * CurlH;
     }

     /* Calculate Ez field */
     for (mm = 1; mm < SizeX; mm++){
       for (nn = 1; nn < SizeY; nn++){
	 for (pp = 0; pp < SizeZ; pp++){
	   CurlH = (Hy(mm, nn, pp) - Hy(mm - 1, nn, pp)) -
	     (Hx(mm, nn, pp) - Hx(mm, nn - 1, pp));
	  Ez(mm, nn, pp) = Ceze(mm, nn, pp) * Ez(mm, nn, pp) +
	    Cezh(mm, nn, pp) * CurlH;
	 }
       }
     }

     for (nn = 1; nn < SizeY; nn++){
	 for (pp = 0; pp < SizeZ; pp++){
	   CurlH = (Hy(0, nn, pp) - Psix*Hy(SizeX - 1, nn, pp)) -
	     (Hx(0, nn, pp) - Hx(0, nn - 1, pp));
	  Ez(0, nn, pp) = Ceze(0, nn, pp) * Ez(0, nn, pp) +
	    Cezh(0, nn, pp) * CurlH;
	 }
     }

     for (mm = 1; mm < SizeX; mm++){
       for (pp = 0; pp < SizeZ; pp++){ 
	 CurlH = (Hy(mm, 0, pp) - Hy(mm - 1, 0, pp)) -
	     (Hx(mm, 0, pp) - Psiy*Hx(mm, SizeY - 1, pp));
	  Ez(mm, 0, pp) = Ceze(mm, 0, pp) * Ez(mm, 0, pp) +
	    Cezh(mm, 0, pp) * CurlH;
       }
     }

     for (pp = 0; pp < SizeZ; pp++){
       CurlH = (Hy(0, 0, pp) - Psix*Hy(SizeX - 1, 0, pp)) -
	 (Hx(0, 0, pp) - Psiy*Hx(0, SizeY - 1, pp));
	  Ez(0, 0, pp) = Ceze(0, 0, pp) * Ez(0, 0, pp) +
	    Cezh(0, 0, pp) * CurlH;
     }
     
	  
  } else {
    fprintf(stderr, "updateE: Unknown grid type.  Terminating...\n");
  }
  
  return;
}  /* end updateE() */
