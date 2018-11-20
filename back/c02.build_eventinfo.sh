#!/bin/tcsh



set PWD = `pwd`
set command = $PWD/.command


cat << EOF> ! $command
use hongyu_db1;
SELECT 
STA,
NET,
DIST,
AZ,
BAZ,
STA_LAT,
STA_LON,
EQ_LAT,
EQ_LON,
EQ_DEP,
EQ_MAG,
EQ_NAME,
POLAR_FLAG,
QUALITY_FLAG,
PREM,
AMP,
CCC_SEW,
SNR_AVE,
S40RTStime,
PHASE,
STRETCH_CCC,
STRETCH_FACTOR,
MISFIT_SIG,
COMP,
100,
100,
-100,
80,
100,
WEIGHT,
SNR2,
MISFIT2,
ONSET,
ENDSET,
TSTAR_FACTOR,
TSTAR_CCC,
CCC_GEW,
MISFIT_PRE,
MISFIT_POST,
RECORD_GAU,
EW_GAU,
GAU_MISFIT,
POLARITY,
POLARITY_COFLAG,
TRAFFIC_NEARBY,
PAIR1,
PAIR1DEP

FROM EQTIME

WHERE S40RTStime is not null
;
EOF

set big_event = $PWD/eventinfo.S40RTS

hongyusql_command $command |awk 'NR>1 {print $0}' >! $big_event

