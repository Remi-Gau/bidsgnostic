clean:
	rm -fr results

manual:
	python bidsgnostic/run.py tests/data/ds001 results participant -c24 --manual

demo:
	bidsgnostic tests/data/ds001 results participant -c24
