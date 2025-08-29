```bash
mkdir build
cd build
cmake ..
make -j${nproc}

./CABD -c grid_4corners
# ./CABD -c grid_1edge
```