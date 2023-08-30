params := $(wordlist 2,100,$(MAKECMDGOALS))
pwd := $(shell pwd)

# Styles
isort:
	isort .

black:
	black .

fix-style: isort black

isort-check:
	isort --check-only --diff .

flake8-check:
	flake8 .

check-style: flake8-check isort-check
