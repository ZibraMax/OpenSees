cd ~
sudo apt update
sudo apt upgrade

# Instalar Python 3.11 con pip y hacerlo el por defecto
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install -y python3.11 python3.11-venv python3.11-distutils python3.11-dev
sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1
sudo update-alternatives --config python3
sudo apt remove -y python3-apt
sudo apt install -y python3-apt
python3.11 -m ensurepip --upgrade
# En lugar de zibramax, poner el nombre de usuario
PATH=$PATH:/home/zibramax/.local/bin



# Instalar build tools para poder compilar software
sudo apt install -y cmake
sudo apt install -y gcc g++ gfortran
sudo apt install -y python3-pip
sudo apt install -y liblapack-dev
sudo apt install -y libopenmpi-dev
sudo apt install -y libmkl-rt
sudo apt install -y libmkl-blacs-openmpi-lp64
sudo apt install -y libscalapack-openmpi-dev

# Saldrá un anuncio purpura, decir Yes, luego enter

# Instalar CONAN
sudo python3 -m pip install --upgrade setuptools
sudo apt update
sudo apt install -y build-essential
# Paquetes que estaban dando problemas
sudo apt remove python3-urllib3
python3 -m pip install --upgrade urllib3
sudo apt remove python3-distro
python3 -m pip install --upgrade distro

sudo python3 -m pip install 'conan<2.0'

# No se si esto se necesite

# sudo apt update
# sudo apt install mpich libscalapack-mpi-dev libmumps-dev

wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | \
  gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg >/dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | \
  sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt update
sudo apt install intel-basekit
sudo apt install intel-oneapi-mkl-devel
source /opt/intel/oneapi/setvars.sh


# Clonar el repo OpenSees de ZibraMax
git clone https://github.com/ZibraMax/OpenSees
cd OpenSees



mkdir build
cd build
conan install .. --build missing
cmake ..
cmake --build . --target OpenSees -j12
cmake --build . --target OpenSeesPy -j12
mv ./lib/OpenSeesPy.so ./opensees.so
cd ..
export PYTHONPATH="./build/"
python3 -c "import sys; print(sys.path)"
python3 ./EXAMPLES/ExamplePython/example_variable_analysis.py