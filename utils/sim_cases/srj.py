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


# _SOURCE_PROCESSES: tuple[_Process, ...] = (
#     _Process(vgrad33=0.000314447, tdot=1.0239016666666667, nsteps=6),
#     _Process(vgrad33=0.000890707, tdot=1.09073, nsteps=3),
#     _Process(vgrad33=0.00201499, tdot=0.824895, nsteps=2),
#     _Process(vgrad33=0.00285462, tdot=0.9869220000000001, nsteps=5),
#     _Process(vgrad33=0.000972793, tdot=1.0153230769230768, nsteps=13),
#     _Process(vgrad33=0.00958007, tdot=1.0654299999999999, nsteps=5),
#     _Process(vgrad33=0.000261499, tdot=1.025576923076923, nsteps=26),
#     _Process(vgrad33=0.026901, tdot=1.1111166666666665, nsteps=3),
#     _Process(vgrad33=0.00248578, tdot=1.0644933333333333, nsteps=15),
# )

_SOURCE_PROCESSES: tuple[_Process, ...] = (
    _Process(vgrad33=0.000314447, tdot=1.0239016666666667 * 2, nsteps=6 // 2),
    _Process(vgrad33=0.000890707, tdot=1.09073, nsteps=3),
    _Process(vgrad33=0.00201499, tdot=0.824895, nsteps=2),
    _Process(vgrad33=0.00285462, tdot=0.9869220000000001, nsteps=5),
    _Process(vgrad33=0.000972793, tdot=2.1998666666666664, nsteps=13 // 2),
    _Process(vgrad33=0.00958007, tdot=1.0654299999999999, nsteps=5),
    _Process(vgrad33=0.000261499, tdot=2.051153846153846, nsteps=26 // 2),
    _Process(vgrad33=0.026901, tdot=1.1111166666666665, nsteps=3),
    _Process(vgrad33=0.00248578, tdot=2.2810571428571427, nsteps=15 // 2),
)

# _SOURCE_PROCESSES: tuple[_Process, ...] = (
#     _Process(vgrad33=0.0011560664997171146, tdot=5.81539, nsteps=1),
#     _Process(vgrad33=0.00285462, tdot=0.9869220000000001, nsteps=5),
#     _Process(vgrad33=0.000972793, tdot=2.1998666666666664, nsteps=13 // 2),
#     _Process(vgrad33=0.00958007, tdot=1.0654299999999999, nsteps=5),
#     _Process(vgrad33=0.000261499, tdot=2.051153846153846, nsteps=26 // 2),
#     _Process(vgrad33=0.026901, tdot=1.1111166666666665, nsteps=3),
#     _Process(vgrad33=0.00248578, tdot=2.2810571428571427, nsteps=15 // 2),
# )

# _SOURCE_PROCESSES: tuple[_Process, ...] = (
#     _Process(vgrad33=0.0006381057880091507, tdot=1.00280514835637, nsteps=3),
#     _Process(vgrad33=0.0015991680311581201, tdot=1.2236278803328355, nsteps=2),
#     _Process(vgrad33=0.0028561516598506393, tdot=1.0338657588530438, nsteps=5),
#     _Process(vgrad33=0.0009722323547671349, tdot=1.0152390212963953, nsteps=13),
#     _Process(vgrad33=0.009564105377266763, tdot=1.0678442338538239, nsteps=5),
#     _Process(vgrad33=0.0002614072133704713, tdot=1.025479119761422, nsteps=26),
#     _Process(vgrad33=0.026904749668203522, tdot=1.1109127225803448, nsteps=3),
#     _Process(vgrad33=0.0024858694834516146, tdot=1.0644984181559825, nsteps=15),
# )

# _SOURCE_PROCESSES: tuple[_Process, ...] = (
#     _Process(vgrad33=0.0006539430515336075, tdot=1.0179965188506863, nsteps=3),
#     _Process(vgrad33=0.0015825343746332043, tdot=1.1672055436239424, nsteps=2),
#     _Process(vgrad33=0.0028462130291877024, tdot=1.0473198712400111, nsteps=5),
#     _Process(vgrad33=0.0009790256301252297, tdot=2.204152180849833, nsteps=6),
#     _Process(vgrad33=0.009619546921950061, tdot=1.0597035791324458, nsteps=5),
#     _Process(vgrad33=0.0002614209382667542, tdot=2.051722802286515, nsteps=13),
#     _Process(vgrad33=0.026852484172231628, tdot=1.113347321476984, nsteps=3),
#     _Process(vgrad33=0.0024851869604735435, tdot=2.28059008929759, nsteps=7),
# )

_SOURCE_PROCESSES: tuple[_Process, ...] = (
    _Process(vgrad33=0.0007258945672295456, tdot=1.1562850978158548, nsteps=3),
    _Process(vgrad33=0.0017626484033519462, tdot=1.0293425955402755, nsteps=2),
    _Process(vgrad33=0.002851100952763204, tdot=1.019491903094377, nsteps=5),
    _Process(vgrad33=0.0009725913800538728, tdot=4.3996634381880035, nsteps=3),
    _Process(vgrad33=0.00956438776347836, tdot=1.067757728308191, nsteps=5),
    _Process(vgrad33=0.0002613403641101477, tdot=4.44365855705927, nsteps=6),
    _Process(vgrad33=0.026903177634321575, tdot=1.1110041781152968, nsteps=3),
    _Process(vgrad33=0.0024857678152979423, tdot=5.322419055731174, nsteps=3),
)

_SOURCE_PROCESSES: tuple[_Process, ...] = (
    _Process(vgrad33=0.0006636534464627386, tdot=1.0506103997986986, nsteps=3),
    _Process(vgrad33=0.0016579851804988492, tdot=1.1761925440302523, nsteps=2),
    _Process(vgrad33=0.0028574290630637893, tdot=1.02415674250868, nsteps=5),
    _Process(vgrad33=0.0009708232909342548, tdot=2.1989370879096284, nsteps=6),
    _Process(vgrad33=0.00955808811136056, tdot=1.0688129746967874, nsteps=5),
    _Process(vgrad33=0.0002616120354535676, tdot=2.050943738142053, nsteps=13),
    _Process(vgrad33=0.02690001973632781, tdot=1.111168217869621, nsteps=3),
    _Process(vgrad33=0.0024853740783381236, tdot=2.2809341928003914, nsteps=7),
)

_SOURCE_PROCESSES = (
    _Process(vgrad33=0.000314447, tdot=1.0239016666666667 * 2, nsteps=3),
    _Process(vgrad33=0.000890707, tdot=1.09073, nsteps=3),
    _Process(vgrad33=0.00201499, tdot=0.824895, nsteps=2),
    _Process(vgrad33=0.00285462, tdot=0.9869220000000001, nsteps=5),
    _Process(vgrad33=0.000972793, tdot=2.1998666666666664, nsteps=6),
    _Process(vgrad33=0.00958007, tdot=1.0654299999999999, nsteps=5),
    _Process(vgrad33=0.000261499, tdot=2.051153846153846, nsteps=13),
    _Process(vgrad33=0.026901, tdot=1.1111166666666665, nsteps=3),
    _Process(vgrad33=0.00248578, tdot=2.2810571428571427, nsteps=7),
)

_SOURCE_PROCESSES = (
    _Process(vgrad33=0.0009360056617230113, tdot=1.2345102497040341, nsteps=4),
    _Process(vgrad33=0.0028157872916966466, tdot=1.1373918002367727 / 2, nsteps=5 * 2),
    _Process(vgrad33=0.0009717718642897617, tdot=4.399024595708375, nsteps=3),
    _Process(vgrad33=0.00956422575394275, tdot=1.067869090102004, nsteps=5),
    _Process(vgrad33=0.0002616264347262647, tdot=4.444054803472794, nsteps=6),
    _Process(vgrad33=0.02691914669527462, tdot=1.1101935970086885, nsteps=3),
    _Process(vgrad33=0.002486246779341785, tdot=5.3228903835006784, nsteps=3),
)

# _SOURCE_PROCESSES = (
#     _Process(vgrad33=0.0009360056617230113, tdot=1.2345102497040341 * 4, nsteps=1),
#     _Process(vgrad33=0.0028157872916966466, tdot=1.1373918002367727 / 2, nsteps=5 * 2),
#     _Process(vgrad33=0.0009717718642897617, tdot=4.399024595708375, nsteps=3),
#     _Process(vgrad33=0.00956422575394275, tdot=1.067869090102004, nsteps=5),
#     _Process(vgrad33=0.0002616264347262647, tdot=4.444054803472794, nsteps=6),
#     _Process(vgrad33=0.02691914669527462, tdot=1.1101935970086885, nsteps=3),
#     _Process(vgrad33=0.002486246779341785, tdot=5.3228903835006784, nsteps=3),
# )


def _smooth(proc_idx: int = 1, substeps: int = 6) -> None:
    global _SOURCE_PROCESSES

    src = list(_SOURCE_PROCESSES)
    if len(src) < 2:
        raise RuntimeError("srj_offset source schedule is shorter than expected")
    if substeps < 2:
        raise ValueError("substeps must be at least 2")
    if not (1 <= proc_idx < len(src)):
        raise ValueError("proc_idx must satisfy 1 <= proc_idx < len(_SOURCE_PROCESSES)")

    prev = src[proc_idx - 1]
    target = src[proc_idx]
    total_duration = float(target.tdot) * int(target.nsteps)
    sub_dt = total_duration / float(substeps)

    v_start = float(prev.vgrad33)
    v_avg = float(target.vgrad33)
    v_end = 2.0 * v_avg - v_start

    ramp: list[_Process] = []
    for i in range(substeps):
        alpha = (i + 0.5) / substeps
        v = (1.0 - alpha) * v_start + alpha * v_end
        ramp.append(_Process(vgrad33=v, tdot=sub_dt, nsteps=1))

    _SOURCE_PROCESSES = tuple([*src[:proc_idx], *ramp, *src[proc_idx + 1 :]])


# for idx in sorted((1), reverse=True):
#     _smooth(idx)
# _smooth(1)


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
100             itmax (maximum allowed number of iterations)
0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?
0               ISAVE write grain states in STRESS.OUT (1) or not (0)?
1               IUPDATE update tex & RVE dim (1) or not (0)?
1               IUPHARD
1               IWTEX (write texture at the end)
1 4            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
"""


def build_fft(cfg: FftConfig) -> str:
    """Return `fft.in` contents for the offset SRJ calibration case."""

    # processes = _offset_processes(SRJ_OFFSET)
    processes = _SOURCE_PROCESSES
    nproc = len(processes)
    body = "".join(_process_block(proc=p) for p in processes)
    convert = False
    if cfg.igas > 0:
        convert = True
    NSTEPS0 = 1
    TDOT0 = SRJ_OFFSET / (4 * NSTEPS0)
    ERR0 = 1e-3
    MAX_ITER0 = 100

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
1 {0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
1 {0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
1 {0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
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
1 {0}            IWFIELDS,IWSTEP (write fields, and frequency of writing)
2. 0            XC0, IGAMMA (0-continious,1-discrete central diff,2-rotated)
0               ITHERMO (if ithermo=1, next line is filethermo)
{body}"""
