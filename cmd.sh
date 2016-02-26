
export PATH="/opt/intel/bin:${PATH}";

source /opt/intel/bin/compilervars.sh intel64

make clean
export APP=decomp
export GRAPH=rmat_200m.graph
export LC_CTYPE=C
export LC_ALL=en_US
export LANG=en_US.UTF-8
make APP=$APP GRAPH="/home/15-418/asst3_graphs/${GRAPH}" jobs
qsub jobs/${USER}_${APP}_${GRAPH}.job


watch qstat -u haibinl


make grade
qsub jobs/${USER}_grade_performance.job
watch qstat -u haibinl


soc-pokec_30m.graph
rmat_200m.graph
soc-livejournal1_68m.graph
com-orkut_117m.graph
