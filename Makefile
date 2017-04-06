data/infered-types.txt : bin/infer.py data/database_link.txt data/family.txt data/so-simple.obo data/manual-assignments.json
	$^ > $@

data/database_link.txt :
	wget -O - 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/database_link.txt.gz' | gzip -d > $@

data/so-simple.obo :
	wget -O - 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/so-simple.obo' > $@

data/family.txt : 
	wget -O - 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz' | gzip -d > $@

.PHONY : clean

clean:
	[ ! -e data/infered-types.txt ] || rm data/infered-types.txt
	[ ! -e data/family.txt ] || rm data/family.txt
	[ ! -e data/database_link.txt ] || rm data/database_link.txt
	[ ! -e data/so-simple.obo ] || rm data/so-simple.obo
