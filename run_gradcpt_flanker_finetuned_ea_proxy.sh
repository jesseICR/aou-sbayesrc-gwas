#!/bin/bash
# run_gradcpt_flanker_finetuned_ea_proxy.sh - Fine-tune SES-EA boosters toward GradCPT/Flanker mean.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export G4_FINETUNE_ETA="${G4_FINETUNE_ETA:-0.03}"
export G4_FINETUNE_MAX_DEPTH="${G4_FINETUNE_MAX_DEPTH:-4}"
export G4_FINETUNE_MIN_CHILD_WEIGHT="${G4_FINETUNE_MIN_CHILD_WEIGHT:-20}"
export G4_FINETUNE_LAMBDA="${G4_FINETUNE_LAMBDA:-2}"
export G4_FINETUNE_MAX_ROUNDS="${G4_FINETUNE_MAX_ROUNDS:-1000}"
export G4_FINETUNE_EARLY_STOPPING_ROUNDS="${G4_FINETUNE_EARLY_STOPPING_ROUNDS:-50}"

exec bash "${SCRIPT_DIR}/run_g4_finetuned_ea_proxy.sh" --target gradcpt-flanker-mean "$@"
