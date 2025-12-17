cd /tmp
git clone https://github.com/giavancini/mumps.git
cd mumps/
sudo rm -rf /opt/mumps/*
# ls -la /opt/mumps/
rm -rf build
cmake -Bbuild -DBUILD_DOUBLE=on -DMUMPS_parallel=false -DMUMPS_openmp=on -DBUILD_SHARED_LIBS=on --install-prefix /opt/mumps
sudo cmake --build build
sudo cmake --install build
# cmake -S example -B example/build -DMUMPS_ROOT=build/local
# cmake --build example/build