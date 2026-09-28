#!/bin/bash
# ground-rule checks after a final_pipeline.py change
cd /home/user/arrhythmia-detection-pipeline
PY=${PY:-/root/venv313/bin/python}
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "--- pytest default env"; $PY -m pytest -q -p no:cacheprovider tests 2>&1 | tail -2
echo "--- pytest AVX-512 disabled"; NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL AVX512_SPR" $PY -m pytest -q -p no:cacheprovider tests 2>&1 | tail -1
echo "--- same-env replicability"; $PY -m experiments.phase3_lle.replicate_check --compare ${REP_BEFORE:-/tmp/claude-0/rep_before.jsonl}
echo "--- run_phase2e --replicate"; $PY -m experiments.phase2e.run_phase2e --replicate --workers 4 | tail -1
cmp experiments/phase2e/results/replicability.json ${REPJSON_BEFORE:-/tmp/claude-0/replicability_env_before.json} && echo "replicability.json identical to pre-change (same machine)"
git checkout -q experiments/phase2e/results/replicability.json
