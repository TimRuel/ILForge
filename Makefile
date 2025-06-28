# ===== Makefile for ILForge Project =====

# ---------- Default Parameters ----------
APP         ?= multinomial_logistic_regression
ESTIMAND    ?= entropy
EXP_ID      ?= exp_v1.0.0
SIM_ID      ?= sim_01
ITER_ID     ?= iter_0001
CORES       ?= 12

# ---------- Paths ----------
MAIN                := common/scripts/main.R
ESTIMAND_DIR        := applications/$(APP)/$(ESTIMAND)
EXPERIMENT_DIR      := experiments/$(EXP_ID)
LAUNCH_SCRIPT       := common/bash/launch_experiment.sh
SIM_DIR             := $(EXPERIMENT_DIR)/simulations/$(SIM_ID)
SIM_REPORT_SCRIPT   := common/bash/render_sim_report.sh
ITER_DIR            := $(SIM_DIR)/$(ITER_ID)
ITER_REPORT_SCRIPT  := common/bash/render_iter_report.sh

# ---------- Commands ----------

## Run a local experiment (not SLURM)
run:
	Rscript $(MAIN) $(APP) $(ESTIMAND) $(EXP_ID) $(SIM_ID) $(ITER_ID) $(CORES)

## Submit to SLURM cluster
slurm:
	sbatch $(LAUNCH_SCRIPT) $(APP) $(ESTIMAND) $(EXP_ID) $(SIM_ID)

## Render iteration report (must exist and be renderable)
iter_report:
	./$(ITER_REPORT_SCRIPT) $(APP) $(ESTIMAND) $(EXP_ID) $(SIM_ID) $(ITER_ID)

## Render simulation report (must exist and be renderable)
sim_report:
	./$(SIM_REPORT_SCRIPT) $(APP) $(ESTIMAND) $(EXP_ID) $(SIM_ID)

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
	@echo "  make run         Run a local test iteration"
	@echo "  make slurm       Submit a SLURM job"
	@echo "  make clean       Remove one run directory"
	@echo "  make run_report  Render run-level report"
	@echo "  make sim_report  Render simulation-level report"
	@echo "  make list        List available experiments"
	@echo ""
	@echo "Parameters (override as needed):"
	@echo "  APP=$(APP)"
	@echo "  ESTIMAND=$(ESTIMAND)"
	@echo "  EXP_ID=$(EXP_ID)"
	@echo "  SIM_ID=$(SIM_ID)"
	@echo "  ITER_ID=$(ITER_ID)"
	@echo "  CORES=$(CORES)"
	@echo ""
	@echo "Examples:"
	@echo "  make run APP=poisson_regression ESTIMAND=weighted_sum EXP_ID=exp_v1.0.6 SIM_ID=sim_02 ITER_ID=iter_0042"
	@echo "  make slurm APP=multinomial_logistic_regression ESTIMAND=shared_effect EXP_ID=exp_v2.1.0 SIM_ID=sim_01"
