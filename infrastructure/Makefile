STATIC_OUTPUTS = target/lines_of_code.csv target/nodes_in_scripts.csv target/scripts_to_benchmark.csv target/cyclomatic.csv target/shellmetrics.sh target/dynamic_analysis.jsonl

static: $(STATIC_OUTPUTS)

target/dynamic_analysis.jsonl: dynamic_analysis.py
	python3 $< | sort > $@

target/scripts_to_benchmark.csv: scripts_to_benchmark.py
	python3 $< | sort > $@

target/lines_of_code.csv: count_lines_of_code.py 
	python3 $< | sort > $@
	
target/nodes_in_scripts.csv: count_nodes_in_scripts.py syntax_analysis.py
	python3 $< | sort > $@

static-test: tests/test_syntax_analysis.py
	python3 -m unittest $<

clean-static: 
	rm -f $(STATIC_OUTPUTS)

target/shellmetrics.sh:
	wget --quiet -O $@ https://raw.githubusercontent.com/shellspec/shellmetrics/b3bfff2af6880443112cdbf2ea449440b30ab9b0/shellmetrics
	chmod +x $@

target/cyclomatic.csv: get_cyclomatic.py target/shellmetrics.sh
	python3 get_cyclomatic.py | sort > $@

.PHONY: static clean-static static-test
