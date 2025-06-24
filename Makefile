# ===== Makefile for ILForge Project =====

# ---------- Default Parameters ----------
APP        ?= multinomial_logistic_regression
ESTIMAND   ?= entropy
EXPID      ?= exp_v1.0.0
SIMID      ?= sim_01
RUNID      ?= iter_0001
CORES      ?= 12

# ---------- Paths ----------
COMMON_MAIN    := applications/$(APP)/$(ESTIMAND)/scripts/main.R
LAUNCH_SCRIPT  := launch_experiment.sh
EXPERIMENT_DIR := experiments/$(EXPID)
SIM_DIR        := $(EXPERIMENT_DIR)/simulations/$(SIMID)
RUN_DIR        := $(SIM_DIR)/$(RUNID)
REPORT_QMD     := $(RUN_DIR)/report.qmd

# ---------- Commands ----------

## Run a local experiment (not SLURM)
run:
	Rscript $(COMMON_MAIN) $(APP) $(ESTIMAND) $(EXPID) $(CORES) $(SIMID) $(RUNID) --force --non-interactive

## Submit to SLURM cluster
slurm:
	sbatch $(LAUNCH_SCRIPT) $(APP) $(ESTIMAND) $(EXPID) $(SIMID)

## Preview report (must exist and be renderable)
report:
	quarto preview $(REPORT_QMD)

## Clean a run directory
clean:
	rm -rf $(RUN_DIR)

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
	@echo "  make report      Preview Quarto report"
	@echo "  make list        List available experiments"
	@echo ""
	@echo "Parameters (override as needed):"
	@echo "  APP=$(APP)"
	@echo "  ESTIMAND=$(ESTIMAND)"
	@echo "  EXPID=$(EXPID)"
	@echo "  SIMID=$(SIMID)"
	@echo "  RUNID=$(RUNID)"
	@echo "  CORES=$(CORES)"
	@echo ""
	@echo "Examples:"
	@echo "  make run APP=poisson_regression ESTIMAND=weighted_sum EXPID=exp_v1.0.6 SIMID=sim_02 RUNID=iter_0042"
	@echo "  make slurm APP=multinomial_logistic_regression ESTIMAND=shared_effect EXPID=exp_v2.1.0 SIMID=sim_01"
