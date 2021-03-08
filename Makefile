PYTHON ?= python
NOSETESTS ?= nosetests

all: clean develop

clean:
	rm -rf *.egg-info
	rm -f `find . -type f -name \*.py[co]`

develop:
	$(PYTHON) setup.py develop
