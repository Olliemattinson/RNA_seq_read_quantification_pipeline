#!/bin/bash
set -euo pipefail
DIR=$(dirname ${BASH_SOURCE[0]})
cd "$DIR"/..

echo 'linting python files with black'
black .

echo 'linting python files with ruff'
ruff . --fix

echo 'linting snakemake files'
snakefmt .

echo 'formatting R files with styler'
Rscript -e 'styler::style_dir()'