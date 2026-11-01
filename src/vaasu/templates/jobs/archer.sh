#!/bin/bash
#SBATCH --job-name=@job_name@
#SBATCH --nodes=@nodes@
#SBATCH --ntasks=@ntasks@
#SBATCH --ntasks-per-node=@tasks_per_node@
#SBATCH --cpus-per-task=@cpus_per_task@
#SBATCH --output=out.%x.%A.%a.out
#SBATCH --time=@time_limit@
#SBATCH --account=e280-biggin
#SBATCH --partition=standard
#SBATCH --qos=@qos@

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
complex_dir=$PWD

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-@cpus_per_task@}"
export SRUN_CPUS_PER_TASK=$OMP_NUM_THREADS
export GMX_MAXBACKUP=-1
export BASE_DIR=/work/e280/e280-Biggin/reub0138/gromacs_plumed_installation
export GROMACS_DIR="$BASE_DIR/gromacs-2020.6/2020.6+plumed"
export PLUMED_DIR="$BASE_DIR/plumed2-2.8.2"
export PATH="$GROMACS_DIR/bin:$PLUMED_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$GROMACS_DIR/lib:$PLUMED_DIR/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

windows=(@windows@)
continuation=@continuation@
if (( continuation )); then
    stages=(prod)
else
    stages=(enmin nvt npt_b npt_pr prod)
fi

previous_stage=""
for stage in "${stages[@]}"; do
    directories=()
    for window in "${windows[@]}"; do
        directories+=("$window/$stage")
        if (( ! continuation )); then
            window_name=${window##*/}
            phase=${window_name%.*}
            number=${window_name##*.}
            if [[ "$stage" == enmin ]]; then
                coordinates="$complex_dir/complex_alchemical.gro"
            else
                coordinates="$complex_dir/$window/$previous_stage/confout.gro"
            fi
            # Preprocess one replica on one rank, then run the stage together.
            (
                cd "$window/$stage"
                OMP_NUM_THREADS=1 srun --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
                    gmx_mpi grompp -f "$phase.$stage.$number.mdp" \
                    -c "$coordinates" -r "$coordinates" \
                    -n "$complex_dir/index_alchemical.ndx" \
                    -p "$complex_dir/complex_alchemical.top" \
                    -o topol.tpr -maxwarn 2
            )
        fi
    done

    options=()
    if [[ "$stage" == prod ]]; then
        options=(-replex 200 -hrex -plumed "$complex_dir/plumed.dat" -dlb no)
        if (( continuation )); then
            options+=(-cpi state.cpt)
        fi
    fi
    srun --ntasks=@ntasks@ --cpu-bind=cores \
        gmx_mpi mdrun -multidir "${directories[@]}" "${options[@]}"
    previous_stage=$stage
done
