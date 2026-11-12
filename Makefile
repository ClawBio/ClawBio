.PHONY: demo test list demo-all

demo:
	python clawbio.py run pharmgx --demo

test:
	python -m pytest -v

list:
	python clawbio.py list

demo-all:
	python clawbio.py run pharmgx --demo
	python clawbio.py run equity --demo
	python clawbio.py run nutrigx --demo
	python clawbio.py run metagenomics --demo
	python clawbio.py run compare --demo
	python clawbio.py run rnaseq --demo
