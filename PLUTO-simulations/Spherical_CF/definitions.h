#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  COMPONENTS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            5

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  T_CRIT                         0
#define  MACH_CRIT                      1
#define  N_CRIT                         2
#define  CHI                            3
#define  SINGULAR                       4

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL BOUNDARY              YES
#define  VX1RHO_OUTFLOW                 YES
#define  SUBSONIC_INI					YES
#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    1e3*CONST_pc
#define  UNIT_VELOCITY                  1.e5

/* [End] user-defined constants (do not change this line) */
