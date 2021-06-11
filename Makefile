PKG_DIR    = linkbudget
PKG_NAME   = link-budget
PY_FILES   = $(shell find . -type f -name '*.py')
VERSION    = $(shell grep "__version__ =" $(PKG_DIR)/cli.py | cut -d '"' -f2)
SDIST      = dist/$(PKG_NAME)-$(VERSION).tar.gz
WHEEL      = dist/$(PKG_NAME)-$(VERSION)-py3-none-any.whl

.PHONY: all clean clean-py sdist wheel install pypi testpypi

all: sdist

clean: clean-py

clean-py:
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

$(SDIST): $(PY_FILES)
	python3 setup.py sdist

sdist: $(SDIST)

$(WHEEL): $(PY_FILES)
	python3 setup.py bdist_wheel

wheel: $(WHEEL)

install: $(SDIST)
	pip3 install $(SDIST)

pypi: clean sdist wheel
	python3 -m twine upload --repository pypi dist/*

testpypi: clean sdist wheel
	python3 -m twine upload --repository testpypi dist/*
