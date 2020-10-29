#!/bin/bash

# Create absolute path for bin files
echo "Creating absolute path for tsRNAsearch 'bin' and 'DBs'..."
myPath=$(pwd)
sed -i -e "s~ bin~ ${myPath}\/bin~g" bin/tsRNAsearch_single.sh # using tilde as delimiter here instead of slash as myPath variable contains slashes
sed -i -e "s~ bin~ ${myPath}\/bin~g" tsRNAsearch
sed -i -e "s~DBs~${myPath}\/DBs~g" bin/tsRNAsearch_single.sh
sed -i -e "s~DBs~${myPath}\/DBs~g" tsRNAsearch
sed -i -e "s~bin\/trim~${myPath}\/bin\/trim~g" bin/tsRNAsearch_single.sh
sed -i -e "s~bin\/feat~${myPath}\/bin\/feat~g" bin/tsRNAsearch_single.sh
sed -i -e "s~additional~${myPath}\/additional~g" bin/tsRNAsearch_single.sh
sed -i -e "s~additional~${myPath}\/additional~g" tsRNAsearch
sed -i -e "s~bin\/tsRNAsearch~${myPath}\/bin\/tsRNAsearch~g" tsRNAsearch
