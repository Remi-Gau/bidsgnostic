clean:
	rm -fr results

manual:
	python bidsgnostic/run.py tests/data/bids-examples/ds001 results participant -c24 --manual

demo:
	mkdir -p results
	bidsgnostic tests/data/bids-examples/ds001 results participant -c24

get_ds114_test1:
	wget https://raw.githubusercontent.com/bids-apps/maintenance-tools/main/utils/get_data_from_osf.sh
	bash get_data_from_osf.sh ds114_test1

test_ds114_test1:
	mkdir -p tmp
	cd tmp
	bidsgnostic \
		~/data/ds114_test1 \
		outputs1 \
		participant \
		-c2 --verbose --log_level 2 \
		--participant_label 01 02
	bidsgnostic \
		~/data/ds114_test1 \
		outputs1 \
		group \
		-c2 --verbose --log_level 2
