# tsRNAsearch

A pipeline for the identification, quantification and analysis of ncRNAs (especially tRNA fragments) in small/miRNA-seq datasets

INSTALL
```
#chmod 755 setup.sh
sudo ./setup.sh -g human # (human/mouse/both)
```
Add tsRNAsearch.sh and tsRNAsearch_DE.sh to your path
RUN

./tsRNAsearch_DE.sh -d ExampleData/ -e additional-files/Example_ExperimentLayout.csv -o MyResults -t 1 
