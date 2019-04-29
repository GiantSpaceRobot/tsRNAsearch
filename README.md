# tsRNAsearch

A pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

INSTALL

chmod 755 setup.sh
sudo ./setup.sh

RUN

./tsRNAsearch_DE.sh -d ExampleData/ -e additional-files/Example_ExperimentLayout.csv -o MyResults -t 1 
