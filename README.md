sudo apt install libbz2-dev libbz2-1.0

sudo apt install zlib1g-dev zlib1g

git clone --recurse-submodules https://github.com/YoraiLevi/DNA_Noam.git

mkdir build

cd build

cmake ..

make

./fm_indexer -f sequence.fna -f sequence2.fasta
