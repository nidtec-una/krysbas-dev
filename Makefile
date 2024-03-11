.DEFAULT_GOAL := help

BROWSER := python -c "$$BROWSER_PYSCRIPT"

# TODO make more general to use the local matlab version
MATLAB = /usr/local/MATLAB/R2017a/bin/matlab
MATLAB_ARG    = -nodisplay -nosplash -nodesktop

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

# determines what "make help" will show
define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

################################################################################
# 	GENERIC
.PHONY: help clean clean-test lint install_dev

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)
clean: clean-test clean-doc ## remove all build, test, coverage artifacts

clean-test: ## remove test and coverage artifacts
	rm -rf coverage_html
	rm -f test_report.log

clean-doc: ## remove doc buils
	rm -rf docs/build
version.txt: ## update version.txt from CITATION.CFF info
	grep -w "^version" CITATION.cff | sed "s/version: /v/g" > version.txt

validate_cff: CITATION.cff ## validate citation file
	cffconvert --validate

manual: ## create pdf of the doc
	cd docs && sh create_manual.sh

################################################################################
# 	MATLAB

.PHONY: lint coverage

lint: ## lint and checks matlab code
	mh_style --fix && mh_metric --ci && mh_lint

coverage: run_tests.m ## runs tests and display coverage
	$(MATLAB) $(MATLAB_ARG) -r "run_tests; exit()"
	$(BROWSER) coverage_html/index.html
