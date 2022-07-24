#!/bin/bash
#An example for MPI job.
#SBATCH -J job_name
#SBATCH -o test-%j.log
#SBATCH -e test-%j.err
#SBATCH -p test 
#SBATCH --qos=testqos
#SBATCH -N 1 -n 20 --cpus-per-task=2


echo Time is `date`
echo Directory is $PWD
echo This job runs on the following nodes:
echo $SLURM_JOB_NODELIST
echo This job has allocated $SLURM_JOB_CPUS_PER_NODE cpu cores.


module purge && module load gcc/9.2.0 elpa/2021.05.002/intelmpi2018 intelmpi/2018.update4 2>&1
module load anaconda3_nompi
source activate pytorch110
#module load python/3.9.1
module list 2>&1


MPIRUN=mpirun #Intel mpi and Open MPI
#MPIRUN=mpiexec #MPICH
MPIOPT="-env I_MPI_FABRICS shm:ofi" #Intel MPI
#MPIOPT="--mca mtl_ofi_provider_include psm2" #Open MPI
#MPIOPT="-iface ib0" #MPICH3
timeout 10 $MPIRUN $MPIOPT hostname 
echo '------'


echo "python3:"
which python3
python3 -u ../SIAB.py ORBITAL_INPUT_DZP

