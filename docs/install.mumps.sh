cd /tmp
git clone https://github.com/giavancini/mumps.git
cd mumps/
sudo rm -rf /opt/mumps/*
# ls -la /opt/mumps/
rm -rf build
sudo update
sudo apt install -y gfortran cmake
cmake -Bbuild \
    -DBUILD_DOUBLE=on \
    -DBUILD_SHARED_LIBS=off \
    -DMUMPS_find_SCALAPACK=off \
    -DMUMPS_intsize64=on \
    -DMUMPS_metis=on \
    -DMUMPS_openmp=on \
    -DMUMPS_parallel=false \
    -DMUMPS_scalapack=off \
    --install-prefix /opt/mumps
sudo cmake --build build
sudo cmake --install build
# cmake -S example -B example/build -DMUMPS_ROOT=build/local
# cmake --build example/build