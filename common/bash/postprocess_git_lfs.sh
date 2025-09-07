#!/bin/bash
#SBATCH --account=p32397
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=timothyruel2024@u.northwestern.edu
#SBATCH --job-name=git_lfs_commit
#SBATCH --output=experiments/%x/postprocess.log
#SBATCH --nodes=1
#SBATCH --ntasks=1

# ===============================
# ✅ Resolve project root
# ===============================
PROJECT_ROOT="$SLURM_SUBMIT_DIR"
cd "$PROJECT_ROOT" || exit 1

# ===============================
# ✅ Load GitHub PAT for non-interactive push
# ===============================
PAT_FILE="$HOME/.github_pat"
if [[ -f "$PAT_FILE" ]]; then
    export GITHUB_PAT=$(cat "$PAT_FILE")
    git config credential.helper "store --file=$HOME/.git-credentials"
    echo "https://$GITHUB_PAT:@github.com" > $HOME/.git-credentials
    echo "✅ GitHub PAT loaded into environment"
else
    echo "❌ GitHub PAT file not found at $PAT_FILE"
    exit 1
fi

# ===============================
# ✅ Configure Git LFS
# ===============================
module load git-lfs/3.3.0-gcc-10.4.0
git lfs install --local

# ===============================
# ✅ Get experiment ID
# ===============================
EXP_ID="$1"
if [[ -z "$EXP_ID" ]]; then
    echo "❌ ERROR: Missing experiment ID argument"
    exit 1
fi

EXPERIMENT_DIR="experiments/${EXP_ID}"

# ===============================
# ✅ Track large experiment folder if not already tracked
# ===============================
if ! git check-attr --all "$EXPERIMENT_DIR" | grep -q "filter: lfs"; then
    echo "📌 Tracking $EXPERIMENT_DIR with Git LFS"
    git lfs track "${EXPERIMENT_DIR}/**"
    git add .gitattributes
    git commit -m "Track experiment ${EXP_ID} with Git LFS" || true
else
    echo "✅ $EXPERIMENT_DIR already tracked with LFS"
fi

# ===============================
# ✅ Add untracked or modified files
# ===============================
git add -u "$EXPERIMENT_DIR"          # stage modified/deleted files
git add "$EXPERIMENT_DIR"             # stage new files

# Commit only if there are changes
if git diff --cached --quiet; then
    echo "✅ No new changes to commit for ${EXP_ID}"
else
    echo "📌 Committing new simulation outputs for ${EXP_ID}"
    git commit -m "Add completed simulation results for ${EXP_ID}"
fi

# ===============================
# ✅ Push commits and LFS objects
# ===============================
echo "📌 Pushing commits and LFS objects to origin/main"
git push origin main
