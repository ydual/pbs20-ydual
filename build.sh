mkdir -p build && cd build
cmake ../project -DCMAKE_BUILD_TYPE=Release
make -j8
cd ..

