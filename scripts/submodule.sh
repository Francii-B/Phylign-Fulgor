#! /usr/bin/env bash

#Install Rust
if ! command -v rustup  2>&1 >/dev/null
then
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
fi

#Compile Fulgor
if [ ! -f external/modified-Fulgor/build/fulgor ]; 
then
	rm -rf external/modified-Fulgor/build
	mkdir external/modified-Fulgor/build
	cd external/modified-Fulgor/build
	cmake ..
	make -j
	if [ ! -f external/modified-Fulgor/build/fulgor ]; 
	then
		echo "Fulgor compiled correctly"
	else
		echo "Fulgor DID NOT compile correctly"
	fi
else
	echo "Fulgor already exist"
    
fi
