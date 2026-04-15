"""Tension simulation-case fft.in writer"""

from __future__ import annotations
from utils.write_fft import FftConfig


def build_fft(cfg: FftConfig) -> str:
    """Return fft.in contents for the Tension case."""

    return f"""{cfg.nphase}                      number of phases (nph)
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
{cfg.cuel_gas}
* restart analysis, remeshing and gas phase compliance
0 200                  ires,imicrosave
{cfg.complgas}                   complgas (% of stiffness)
* initial displacement field
0                      IDISPINI (if idispini=1, next line is filedispl)
6                      number of processes (nproc)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.0008              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10            nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
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
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.005               |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10              nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 10            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.0052              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10              nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 25            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.005              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10              nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 25            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.0052              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10              nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 25            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       1                     |            
										  |                               
	0.0     0.        0.          udot    |    vel.grad                   
	0.       0.0      0.                  |                               
	0.       0.       0.0053              |                                       
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 0                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 0.                   @                               
* other                                                                   
1 1 .1 .1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS                                         
10              nsteps (total number of time steps)     
1.0e-3 1.0e-7  err,erral (tolerances for fields and Augm.Lagr.)
150             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?      
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?                    
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 25            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
"""
