
## Installation
### Singularity
First you will need to install some dependencies and install Go.

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config
    
export VERSION=1.18 OS=linux ARCH=amd64 &&     # change this as you need
wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz
```

Next, add /usr/local/go/bin to the PATH environment variable:

```bash
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc
```

Clone the singularity repository with git in a location of your choice:

```bash
git clone --recurse-submodules https://github.com/sylabs/singularity.git
cd singularity
```

You can configure, build, and install SingularityCE using the following commands:

```bash
./mconfig
make -C builddir
sudo make -C builddir install
```

### Snakemake
Next, you need to install a Conda-based Python3 distribution. The recommended choice is Mambaforge. Download the installer using curl and run the script.

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh
```

Given that Mamba is installed, run

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

to install Snakemake in an isolated environment, that has to be activated with

```bash
conda activate snakemake
```

### AWSCLI
Install AWS Command Line Interface

```bash
conda install -c conda-forge awscli
```

then use the command (use access key and secret as provided)

```bash
$ aws configure
AWS Access Key ID [None]:                                    # change this as you need
AWS Secret Access Key [None]:                                # change this as you need
Default region name [None]: eu-central-1
Default output format [None]: json
```

to set up your AWS CLI installation with the keys stored in "collaborator_accessKeys.csv".
