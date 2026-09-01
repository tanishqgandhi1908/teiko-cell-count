# make setup / pipeline / dashboard

PYTHON ?= $(shell command -v python3 2>/dev/null || command -v python)
PORT   ?= 8501

.PHONY: setup pipeline dashboard test clean all

setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) run_analysis.py

dashboard:
	@test -f cell_counts.db || $(MAKE) pipeline
	$(PYTHON) -m streamlit run dashboard/app.py --server.port $(PORT) --server.address 0.0.0.0

test:
	$(PYTHON) -m pytest -q tests

clean:
	rm -f cell_counts.db
	rm -rf outputs/*.csv outputs/*.md outputs/figures/*.png

all: setup pipeline
