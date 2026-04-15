"""Offset SRJ calibration simulation-case `fft.in` writer."""

from __future__ import annotations

from dataclasses import dataclass

from write_fft import FftConfig
from config import SRJ_OFFSET, SRJ_OFFSET_STRESS
from calc_utils import convert_load
from data_utils import MicrostructureInfo


@dataclass(frozen=True)
class _Process:
    """One loading segment ("process") in `fft.in`."""

    vgrad33: float
    tdot: float
    nsteps: int


_SOURCE_PROCESSES: tuple[_Process, ...] = (
    _Process(vgrad33=0.000314447, tdot=1.0239016666666667, nsteps=6),
    _Process(vgrad33=0.000890707, tdot=1.09073, nsteps=3),
    _Process(vgrad33=0.00201499, tdot=0.824895, nsteps=2),
    _Process(vgrad33=0.00285462, tdot=0.9869220000000001, nsteps=5),
    _Process(vgrad33=0.000972793, tdot=1.0153230769230768, nsteps=13),
    _Process(vgrad33=0.00958007, tdot=1.0654299999999999, nsteps=5),
    _Process(vgrad33=0.000261499, tdot=1.025576923076923, nsteps=26),
    _Process(vgrad33=0.026901, tdot=1.1111166666666665, nsteps=3),
    _Process(vgrad33=0.00248578, tdot=1.0644933333333333, nsteps=15),
)


def _process_block(*, proc: _Process) -> str:
    """Return the `fft.in` text for a single loading process."""

    tdot = float(proc.tdot)
    tdotref = tdot
    tdotmin = max(tdot / 1, 1.0e-6)
    v33 = float(proc.vgrad33)
    nsteps = int(proc.nsteps)

    return f"""*INFORMATION ABOUT TEST CONDITIONS
* boundary conditions
\t0       1       1           iudot     |    flag for vel.grad.
\t1       0       1                     |    (0:unknown-1:known)
\t1       1       1                     |
\t\t\t\t\t\t\t  |
\t0.0     0.0      0.0          udot    |    vel.grad
\t0.0     0.0      0.0                  |
\t0.0     0.0      {v33:.10g}               |
\t\t\t\t\t\t\t  |
\t1       0        0           iscau    |    flag for Cauchy
\t\t\t1        0                    |
\t\t\t\t\t 0                    |
\t\t\t\t\t\t\t  |
\t0.      0.       0.          scauchy  |    Cauchy stress
\t\t\t0.       0.                   |
\t\t\t\t\t 0                   @
* other
{tdot:.10g} {tdotref:.10g} {tdotmin:.10g} 1  eqincr or tdot, tdotref, tdotmin, devmmx
-1                     ictrl (1-6: strain comp, 0: VM eq, -1: tdot)
*INFORMATION ABOUT RUN CONDITIONS
{nsteps}               nsteps (total number of time steps)
1.0e-3 1.0e-10  err,erral (tolerances for fields and Augm.Lagr.)
190             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 8            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
"""


def _offset_processes(offset_s: float) -> list[_Process]:
    """Return the SRJ schedule starting at `offset_s` in the source case."""

    remaining = float(offset_s)
    shifted: list[_Process] = []
    tol = 1.0e-12

    for proc in _SOURCE_PROCESSES:
        dt = float(proc.tdot)
        nsteps = int(proc.nsteps)
        duration = dt * nsteps

        if remaining >= duration - tol:
            remaining = max(0.0, remaining - duration)
            continue

        if remaining > tol:
            skipped_full = int(remaining // dt)
            offset_in_step = remaining - skipped_full * dt
            steps_left = nsteps - skipped_full

            if offset_in_step > tol:
                shifted.append(
                    _Process(
                        vgrad33=proc.vgrad33,
                        tdot=dt - offset_in_step,
                        nsteps=1,
                    )
                )
                steps_left -= 1

            if steps_left > 0:
                shifted.append(
                    _Process(vgrad33=proc.vgrad33, tdot=dt, nsteps=steps_left)
                )
            remaining = 0.0
            continue

        shifted.append(proc)

    if not shifted:
        raise ValueError("SRJ offset removes the entire source schedule")

    return shifted


def build_fft(cfg: FftConfig) -> str:
    """Return `fft.in` contents for the offset SRJ calibration case."""

    processes = _offset_processes(SRJ_OFFSET)
    nproc = len(processes)
    body = "".join(_process_block(proc=p) for p in processes)
    convert = False
    if cfg.igas > 0:
        convert = True
    NSTEPS0 = 1
    TDOT0 = SRJ_OFFSET / (4 * NSTEPS0)
    ERR0 = 1e-3
    MAX_ITER0 = 200

    STRAINRATE0 = 0.00029838530595147696
    STRESS0 = SRJ_OFFSET_STRESS / 4
    STRESS1 = 2 * SRJ_OFFSET_STRESS / 4
    STRESS2 = 3 * SRJ_OFFSET_STRESS / 4
    STRESS3 = SRJ_OFFSET_STRESS
    if convert:
        info = MicrostructureInfo.load(cfg.microstructure)
        STRESS0 = convert_load(STRESS0, info)
        STRESS1 = convert_load(STRESS1, info)
        STRESS2 = convert_load(STRESS2, info)
        STRESS3 = convert_load(STRESS3, info)

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
{nproc+4}                      number of processes (nproc)
*INFORMATION ABOUT TEST CONDITIONS                                                            
* boundary conditions                                                     
	0       1       1           iudot     |    flag for vel.grad. 
	1       0       1                     |    (0:unknown-1:known)        
	1       1       0                     |            
										  |                               
	0.0     0.0      0.0          udot    |    vel.grad                   
	0.0     0.0      0.0                  |                               
	0.0     0.0      {STRAINRATE0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {STRESS0}                   @                               
* other                                                                   
{TDOT0} 1 {TDOT0/1} {1}  eqincr or tdot, tdotref, tdotmin, devmmx
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
1 {100}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
	0.0     0.0      {STRAINRATE0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {STRESS1}                   @                               
* other                                                                   
{TDOT0} 1 {TDOT0/1} {1}  eqincr or tdot, tdotref, tdotmin, devmmx
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
1 {100}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
	0.0     0.0      {STRAINRATE0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {STRESS2}                   @                               
* other                                                                   
{TDOT0} 1 {TDOT0/1} {1}  eqincr or tdot, tdotref, tdotmin, devmmx
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
1 {100}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
	0.0     0.0      {STRAINRATE0}               |                                      
										  |                               
	1       0        0           iscau    |    flag for Cauchy            
			1        0                    |                               
					 1                    |                               
										  |                               
	0.      0.       0.          scauchy  |    Cauchy stress              
			0.       0.                   |                               
					 {STRESS3}                   @                               
* other                                                                   
{TDOT0} 1 {TDOT0/1} {1}  eqincr or tdot, tdotref, tdotmin, devmmx
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
1 {100}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
{body}"""
