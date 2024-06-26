#!/bin/bash 

#SBATCH --job-name=rstudio
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=4gb
#SBATCH --output=logs/rstudio.out
#SBATCH --partition=amilan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.wells-wrasman@cuanschutz.edu

# Set up environment
export ALPINE_SCRATCH=/gpfs/alpine1/scratch/$USER
export APPTAINER_TMPDIR=$ALPINE_SCRATCH/apptainer_tmp
export APPTAINER_CACHEDIR=$ALPINE_SCRATCH/apptainer_cache
mkdir -pv $APPTAINER_CACHEDIR $APPTAINER_TMPDIR

mkdir -p logs

# path to directory on HPC for persistant storage of R packages
USER_R_LIB=/projects/${USER}/R/RNA_seq/4.2
LOGIN_HOST=login-ci.rc.colorado.edu
RSA_KEY=/Users/wellskr/.ssh/id_rsa

# path to sif file on HPC
SINGULARITY_IMAGE="rnaseq_r_v1.sif"

# Change home path so that rstudio saves files to projects instead
export HOME=/projects/${USER}

# add options for singularity exec
# e.g. "--bind /path/to/some/other/user/directory"
SINGULARITY_EXEC_OPTS=$1

max_n_cores=$(grep processor /proc/cpuinfo | wc -l)
# customize --output path as appropriate (to a directory readable only by the user!)

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment
cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=${max_n_cores}
export R_LIBS_USER=$USER_R_LIB
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export APPTAINER_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server,/pl/active/Anschutz_BDC"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0

export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD=$(openssl rand -base64 15)
# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:
   ssh -i ${RSA_KEY} -N -L 8787:${HOSTNAME}:${PORT} ${APPTAINERENV_USER}@${LOGIN_HOST}
   and point your web browser to http://localhost:8787
2. log in to RStudio Server using the following credentials:
   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}
When done using RStudio Server, terminate the job by:
1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:
      scancel ${SLURM_JOB_ID}
END

singularity exec $SINGULARITY_EXEC_OPTS --cleanenv $SINGULARITY_IMAGE \
    rserver --www-port ${PORT} \
            --server-user ${USER} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh
printf 'rserver exited' 1>&2
