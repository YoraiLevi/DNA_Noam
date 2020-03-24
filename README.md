sudo apt install libbz2-dev libbz2-1.0

sudo apt install zlib1g-dev zlib1g

git clone --recurse-submodules https://github.com/YoraiLevi/DNA_Noam.git

cd DNA_Noam/

mkdir Release

cd Release

cmake -DCMAKE_BUILD_TYPE=ls ..

make

./fm_indexer -f sequence.fna -f sequence2.fasta

./onlyhuman -f ~/Source/datasets/DNA/fake_tests/seq0.fasta -f ~/Source/datasets/DNA/fake_tests/seq1.fasta &> out.log
