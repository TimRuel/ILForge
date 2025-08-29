# ===== Makefile for ILForge Project =====

# ---------- Default Parameters ----------
APP         ?= poisson
ESTIMAND    ?= group_rates_weighted_sum
MODEL				?= naive_group_rates
EXP_ID      ?= exp_v1.0.0
SIM_ID      ?= sim_01
ITER_ID     ?= iter_0001
CORES       ?= 12
FLAGS       ?=

# ---------- Paths ----------
MAIN               := common/scripts/main.R
MODEL_DIR          := applications/$(APP)/$(ESTIMAND)/$(MODEL)
EXPERIMENT_DIR     := experiments/$(EXP_ID)
LAUNCH_SCRIPT      := common/bash/launch_experiment.sh
SIM_DIR            := $(EXPERIMENT_DIR)/simulations/$(SIM_ID)
SIM_REPORT_SCRIPT  := common/bash/render_sim_report.sh
ITER_DIR           := $(SIM_DIR)/$(ITER_ID)
ITER_REPORT_SCRIPT := common/bash/render_iter_report.sh

# ---------- Commands ----------

## Run a local experiment (not SLURM)
run:
	Rscript $(MAIN) $(APP) $(ESTIMAND) $(MODEL) $(EXP_ID) $(SIM_ID) $(ITER_ID) $(CORES) $(FLAGS)

## Submit to SLURM cluster
slurm:
	sbatch $(LAUNCH_SCRIPT) $(APP) $(ESTIMAND) $(MODEL) $(EXP_ID) $(SIM_ID) $(FLAGS)

## Render iteration report (must exist and be renderable)
iter_report:
	./$(ITER_REPORT_SCRIPT) $(APP) $(ESTIMAND) $(MODEL) $(EXP_ID) $(SIM_ID) $(ITER_ID)

## Render simulation report (must exist and be renderable)
sim_report:
	./$(SIM_REPORT_SCRIPT) $(APP) $(ESTIMAND) $(MODEL) $(EXP_ID) $(SIM_ID)

## Clean a run directory
clean:
	rm -rf $(ITER_DIR)

## List available experiments
list:
	@echo "Available experiments:"
	@ls -1 experiments | grep '^exp_'

## Help
help:
	@echo "ILForge Makefile — Usage:"
	@echo ""
	@echo "  make run FLAGS=\"--skip-integrated --skip-profile\""
	@echo "  make slurm FLAGS=\"--skip-profile\""
	@echo "  make clean       Remove one run directory"
	@echo "  make iter_report Render iteration report"
	@echo "  make sim_report  Render simulation-level report"
	@echo "  make list        List available experiments"
	@echo ""
	@echo "Parameters (override as needed):"
	@echo "  APP=$(APP)"
	@echo "  ESTIMAND=$(ESTIMAND)"
	@echo "  MODEL=$(MODEL)"
	@echo "  EXP_ID=$(EXP_ID)"
	@echo "  SIM_ID=$(SIM_ID)"
	@echo "  ITER_ID=$(ITER_ID)"
	@echo "  CORES=$(CORES)"
