#define  TOWNSEND                       9
#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  COMPONENTS                     1
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     NO
#define  FORCED_TURB                    NO
#define  COOLING                        TOWNSEND
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        1
#define  USER_DEF_PARAMETERS            2

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  INTERNAL_BOUNDARY              YES
#define  PARTICLES                      NO
#define  SHOW_TIMING                    YES
#define  SHOW_TIME_STEPS                YES
#define  HIGHT_FLAT                     YES
 
/* -- user-defined parameters (labels) -- */

#define  N_HOT                          0
#define  T_HOT                          1

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   CONST_mp
#define  UNIT_LENGTH                    CONST_pc
#define  UNIT_VELOCITY                  1.e5

/* [End] user-defined constants (do not change this line) */
#define  MULTIPLE_LOG_FILES             YES
#define  VERBOSE                        YES
#define  DEBUG                          1
