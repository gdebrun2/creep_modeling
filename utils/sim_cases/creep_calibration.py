"""Stress-controlled creep simulation-case fft.in writer."""

from __future__ import annotations
from write_fft import FftConfig


def build_fft(cfg: FftConfig) -> str:
    """Return fft.in contents for the SRJ case."""
    if cfg.load is None:
        raise ValueError("cfg.load is required for creep")
    load = float(cfg.load)
    NPROC = 17

    # strain rate per second
    STRAIN_RATE_0 = 1e-4
    STRAIN_RATE_1 = 1e-5
    STRAIN_RATE_2 = 2e-6
    STRAIN_RATE_3 = 5e-7
    STRAIN_RATE_4 = 8.5e-8
    STRAIN_RATE_5 = 3.5e-9

    NSEC0 = 1
    NSEC1 = 720
    NSEC2 = 3600
    NSEC3 = 14400
    NSEC4 = 64080
    NSEC5 = 82000

    NSTEPS0 = NSEC0
    NSTEPS1 = 50 // 8
    NSTEPS2 = 200 // 64
    NSTEPS3 = 200 // 64
    NSTEPS4 = 200 // 64
    NSTEPS5 = 37 // 8

    TDOT0 = 1
    TDOT1 = NSEC1 / NSTEPS1  # 14s
    TDOT2 = NSEC2 / NSTEPS2  # 18s
    TDOT3 = NSEC3 / NSTEPS3  # 72s
    TDOT4 = NSEC4 / NSTEPS4  # 320s
    TDOT5 = NSEC5 / NSTEPS5  # 2200s
    NSTEPS5 += 15

    TDOT0_MIN = TDOT0
    TDOT1_MIN = TDOT1 // 2
    TDOT2_MIN = TDOT2 // 2
    TDOT3_MIN = TDOT3 // 2
    TDOT4_MIN = TDOT4 // 2
    TDOT5_MIN = TDOT5

    DEVMMX = 0.5

    ################### MANUAL TROUBLESHOOTING ##################
    # TDOT0 = 1
    # TDOT1 = 1
    # TDOT2 = 1
    # TDOT3 = 1  # 36
    # TDOT4 = 1  # 320
    # TDOT5 = 1  # 2200

    # NSTEPS0 = 1
    # NSTEPS1 = 2  # 20
    # NSTEPS2 = 2  # 20
    # NSTEPS3 = 2  # 20
    # NSTEPS4 = 2  # 20
    # NSTEPS5 = 2  # 20

    # STRAIN_RATE_0 = 1e-3
    # STRAIN_RATE_1 = 1e-3
    # STRAIN_RATE_2 = 2e-3
    # STRAIN_RATE_3 = 1e-3
    # STRAIN_RATE_4 = 1e-3
    # STRAIN_RATE_5 = 1e-3

    ################### END MANUAL TROUBLESHOOTING ##################

    NWRITE0 = 2
    NWRITE1 = NSTEPS1 // 5
    NWRITE2 = NSTEPS2 // 5
    NWRITE3 = NSTEPS3 // 5
    NWRITE4 = NSTEPS4 // 5
    NWRITE5 = NSTEPS5 // 5

    NWRITE0 = 0
    NWRITE1 = 0
    NWRITE2 = 0
    NWRITE3 = 3
    NWRITE4 = 3
    NWRITE5 = 4

    MAX_ITER0 = 200
    MAX_ITER1 = 50
    MAX_ITER2 = 200
    MAX_ITER3 = 200
    MAX_ITER4 = 200
    MAX_ITER5 = 200

    ERR0 = 1e-3
    ERR1 = 1e-3
    ERR2 = 1e-3
    ERR3 = 1e-3
    ERR4 = 1e-3
    ERR5 = 1e-3
    # devmx is maximum von mises stress update/increment
    # start at 375, work way upward
    inc1 = 375 + (load - 375) / 5
    inc2 = 375 + 2 * (load - 375) / 5
    inc3 = 375 + 3 * (load - 375) / 5
    inc4 = 375 + 4 * (load - 375) / 5

    if cfg.nphase == 1:
        stub = f"""{cfg.nphase}                      number of phases (nph)
{cfg.nx} {cfg.ny} {cfg.nz}               nx,ny,nz - number of points
1.  1.  1.             RVE dimensions (delt) - keep as 1.0 1.0 1.0
* name and path of microstructure file (filetext)
{cfg.microstructure}
*INFORMATION ABOUT PHASE #1
0                          igas(iph)
* name and path of single crystal files (dummy if igas(iph)=1)
{cfg.cupl2}
{cfg.cuel}"""
    elif cfg.nphase == 2:
        stub = f"""{cfg.nphase}                      number of phases (nph)
{cfg.nx} {cfg.ny} {cfg.nz}               nx,ny,nz - number of points
1.  1.  1.             RVE dimensions (delt) - keep as 1.0 1.0 1.0
* name and path of microstructure file (filetext)
{cfg.microstructure}
*INFORMATION ABOUT PHASE #1
0                          igas(iph)
* name and path of single crystal files (dummy if igas(iph)=1)
{cfg.cupl2}
{cfg.cuel}
*INFORMATION ABOUT PHASE #2
{cfg.igas}                          igas(iph)
* name and path of single crystal files (dummy if igas(iph)=1)
{cfg.cupl2_gas}
{cfg.cuel_gas}"""

    return (
        stub
        + f"""
* restart analysis, remeshing and gas phase compliance
0 200                  ires,imicrosave
{cfg.complgas}                   complgas (% of stiffness)
* initial displacement field
0                      IDISPINI (if idispini=1, next line is filedispl)
{NPROC}                      number of processes (nproc)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 56                   @                               
* other                                                                   
{TDOT0} 1 {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 112                  @                               
* other                                                                   
{TDOT0} 1 {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 168                  @                               
* other                                                                   
{TDOT0} 1 {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 224                  @                               
* other                                                                   
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 280                  @                               
* other                                                                   
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 350                  @                               
* other                                                                   
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 2            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       {STRAIN_RATE_0}              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 375                  @                               
* other                                                                  
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}              itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1  {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       {STRAIN_RATE_0}              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {inc1}                  @                               
* other                                                                  
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}              itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1  {NWRITE0}          IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       {STRAIN_RATE_0}              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {inc2}                  @                               
* other                                                                  
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-9  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}              itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1  {NWRITE0}          IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       {STRAIN_RATE_0}              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {inc3}                  @                               
* other                                                                  
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}   eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}              itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1  {NWRITE0}          IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       {STRAIN_RATE_0}              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {inc4}                  @                               
* other                                                                  
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}              itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1  {NWRITE0}          IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT0} {TDOT0} {TDOT0_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS0}               nsteps (total number of time steps)     
{ERR0} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER0}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_1}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT1} {TDOT1} {TDOT1_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS1}               nsteps (total number of time steps)     
{ERR1} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER1}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE1}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_2}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |     2                         
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT2} {TDOT2} {TDOT2_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS2}               nsteps (total number of time steps)     
{ERR2} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER2}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE2}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_3}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT3} {TDOT3} {TDOT3_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS3}               nsteps (total number of time steps)     
{ERR3} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER3}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE3}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_4}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT4} {TDOT4} {TDOT4_MIN} {DEVMMX}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS4}               nsteps (total number of time steps)     
{ERR4} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER4}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE4}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAIN_RATE_5}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {load}                   @                               
* other                                                                   
{TDOT5} {TDOT5} {TDOT5_MIN} {1}  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
{NSTEPS5}               nsteps (total number of time steps)     
{ERR5} 1.0e-8  err,erral (tolerances for fields and Augm.Lagr.)
{MAX_ITER5}             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 {NWRITE5}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
"""
    )
