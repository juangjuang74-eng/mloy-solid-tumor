.PHONY: help install data baseline survival mediation test all quick clean

PY ?= python3
N  ?=            # empty = full cohort; set e.g. N=50000 for a quick run
NARG = $(if $(N),--n $(N),)

help:
	@echo "Available targets:"
	@echo "  make install    - install dependencies"
	@echo "  make data       - generate the synthetic cohort (N=... for a subsample)"
	@echo "  make baseline   - Table 1 (characteristics and mLOY prevalence)"
	@echo "  make survival   - Table 2 (Cox HRs). Add SITES=1 STRAT=1 for the full set"
	@echo "  make mediation  - mediation analysis, smoking -> mLOY -> tumour"
	@echo "  make all        - run the whole pipeline"
	@echo "  make quick      - quick pipeline on a subsample (N=50000)"
	@echo "  make test       - run pytest"
	@echo "  make clean      - remove generated data and results"

install:
	$(PY) -m pip install -r requirements.txt

data:
	$(PY) scripts/01_generate_data.py $(NARG)

baseline:
	$(PY) scripts/02_baseline_characteristics.py $(NARG)

survival:
	$(PY) scripts/03_survival_analysis.py $(NARG) \
		$(if $(SITES),--sites,) $(if $(STRAT),--stratified,)

mediation:
	$(PY) scripts/04_mediation_analysis.py $(NARG)

all: data baseline survival mediation

quick:
	$(MAKE) all N=50000

test:
	$(PY) -m pytest tests/ -q

clean:
	rm -f data/*.csv data/*.parquet results/*.csv results/*.txt
	rm -rf __pycache__ src/__pycache__ tests/__pycache__ .pytest_cache
