# haplovalidate
Install requirements: rust ver 1.3 or later, clang
```
curl https://sh.rustup.rs -sSf | sh
echo 'export PATH=~/.cargo/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
which cargo
```
If the build fails on the htslib dependency you might need xz.
```
export CFLAGS='-I/path/to/xz/<version>/include'
or add that to your .bashrc and source it
```
Then you should be able to clone and install the project.
```
git clone https://github.com/wheaton5/haplovalidate
cd haplovalidate
cargo build --release
```
And the usage is.
```
./target/release/haplo_validate -h
het_snp_kmers 1.0
Haynes Heaton <whheaton@gmail.com>
Finds kmer pairs that are different in the middle base and each have roughly haploid coverage. Meant for illumina data
as an initial step for de novo phasing.

USAGE:
    haplo_validate [OPTIONS] --primary_fasta <primary_fasta> --secondary_fasta <secondary_fasta>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -k, --het_kmers <het_kmers>                file containing het kmer pairs as output by het_kmer_pair tool <kmer1>
                                               <count> <kmer2> <count> per file
    -i, --primary_fasta <primary_fasta>        primary assembly fasta
    -s, --secondary_fasta <secondary_fasta>    alt fasta representing contigs from the second haplotype
```
which will tell you how many kmer pairs are in each assembly in a number of configurations such as seeing one is the primary and the other in the secondary (true phased)
