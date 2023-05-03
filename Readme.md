```
mkdir build_u
cd build_u
cmake ../
make

./geotest ../macros/check_gdml.mac
```


```
mkdir build_w
cd build_w
cmake ../
cmake --build . --config Release

./geotest ../macros/check_gdml.mac
```