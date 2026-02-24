#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"

# LAMMPS executable
LMP_EXE="/users/qzhu8/GitHub/lammps-stable_23Jun2022_update3/src/lmp_mpi"

# Parameter grids
layers=(10 20 30)
k1s=(100)
k2s=(1 10 100)
k3s=(5)
v=(2e-7)

for layer in "${layers[@]}"; do
  NCPUS=$((layer * 12 / 10))
  for k1 in "${k1s[@]}"; do
    for k2 in "${k2s[@]}"; do
      for k3 in "${k3s[@]}"; do
        for v in "${v[@]}"; do
            job_name="L${layer}-s${k1}-b${k3}-l${k2}-v${v}"
            job_dir="${BASE_DIR}/Benchmark/${job_name}"
            dump_file="bend-L${layer}-s${k1}-b${k3}-l${k2}-v${v}.dump"

            mkdir -p "$job_dir"
            cp "${BASE_DIR}/${layer}layer.data" "$job_dir/"
            cp "${BASE_DIR}/in.bend" "$job_dir/in.bend"

            sed -i "s/variable nlayers equal .*/variable nlayers equal ${layer}/" "$job_dir/in.bend"
            sed -i "s/variable k1 equal .*/variable k1 equal ${k1}    # inplane strech (Morse)/" "$job_dir/in.bend"
            sed -i "s/variable k2 equal .*/variable k2 equal ${k2}    # out-of-plane stretch (Morse)/" "$job_dir/in.bend"
            sed -i "s/variable k3 equal .*/variable k3 equal ${k3}    # bending (Morse)/" "$job_dir/in.bend"
            sed -i "s/variable vdown equal .*/variable vdown equal ${v}   # slower for better quasi-static behavior/" "$job_dir/in.bend"
            sed -i "s/variable dumpfile equal .*/variable dumpfile equal ${dump_file}/" "$job_dir/in.bend"
            cat > "${job_dir}/run.sbatch" <<EOF
#!/bin/bash
#SBATCH -J ${job_name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${NCPUS}
#SBATCH --time=96:00:00
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err

cd "${job_dir}"

echo "Running on node: \$(hostname)"
echo "Job Dir : \$PWD"
echo "Job ID/Name : \$SLURM_JOBID / \$SLURM_JOB_NAME"

echo "Starting: ${job_name}"
mpirun -np ${NCPUS} ${LMP_EXE} < in.bend > out.bend
python "${BASE_DIR}/plot_combined.py" -d "${dump_file}"

cp combined_plot.png "../bend-${job_name}.png"
EOF

            sbatch "${job_dir}/run.sbatch"
        done
      done
    done
  done
done
