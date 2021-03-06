# Setup the running environment and dependencies <br/>
## Tested on MacOS 10.13.6

# Check to see if Anaconda is installed or not
which conda
if [ $? -eq 1 ];then
echo "**ERROR**
Dependency in need: Anaconda3.
To install, follow: https://conda.io/projects/conda/en/latest/user-guide/
or troubleshooting.txt setup session.
************************************************
After installation, run:
bash (if using Ubuntu)
bash setup.sh" >&2
exit 1
fi

# Setup running environment with Anaconda3
conda create -y -n py2-dms python=2.7.14
conda create -y -n py3-dms python=3.6.4

conda config --add channels r
conda config --add channels bioconda
conda config --add channels anaconda

(source activate py2-dms) && (source activate py3-dms)
if [ $? -eq 1 ]; then
echo "**ERROR**
Python environments not ready.
Check \"Setup running environments with Anaconda3\" session in troubleshooting.txt.
Set correctly and rerun." >&2
exit 1
fi

# Install third-party tools with conda install
source activate py3-dms
conda install -y bbmap=37.90.0
conda install -y bwa=0.7.15
conda install -y samtools=1.7
conda install -y seqtk=1.2
conda install -y sra-tools=2.8.2
conda install -y wget=1.19.1

conda install -y pysam=0.14.1
conda install -y biopython=1.72
conda install -y pytables=3.4.2

#(which bbduk.sh) && (which bwa) && (which samtools) && (which seqtk) && (which fastq-dump) && (which wget) && (conda install -y pysam=0.14.1)
if [ $? -eq 1 ]; then
echo "**ERROR**
Required tools not ready.
Check \"Install third-party tools with conda install\" session in troubleshooting.txt.
Install them correctly and rerun" >&2
exit 1
fi
## Note, pysam can be correctly installed with pip on Linux but not on MacOS.

# Install Python packages with pip install
python -m pip install -I pandas==0.22.0 numpy==1.14.1 matplotlib==2.2.3 seaborn==0.8.1
if [ $? -eq 1 ]; then
echo "**ERROR**
Required python packages not ready.
Check \"Install Python packages with pip install\" session in README_setup.md.
Install them correctly and rerun" >&2
exit 1
fi

# Install Enrich2 under py2-dms environment
wget https://github.com/FowlerLab/Enrich2/archive/master.zip
(unzip master.zip) && (rm master.zip)

cp setup.py Enrich2-master/setup.py

cd Enrich2-master
source activate py2-dms

python setup.py install

rm -r ../Enrich2-master
which enrich_cmd
if [ $? -eq 1 ]; then
echo "**ERROR**
Enrich2 not ready.
Check \"Install Enrich2 under py2-dms environment\" session in troubleshooting.txt.
Install it correctly and rerun" >&2
exit 1
fi

echo "
*** Setup accomplished! ***
"