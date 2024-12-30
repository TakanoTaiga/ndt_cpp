# NDT C++

Minimum NDT


https://github.com/TakanoTaiga/ndt_cpp/assets/53041471/510656ef-73a8-49dd-b51f-698165e1922a



## RUN

g++
```
g++ -O2 ./main.cpp -o main.out && ./main.out

g++ -O2 ./main_downsample.cpp -o main2.out && ./main2.out
```

NVIDIA HPC SDK
```
nvc++ -fast -O2 ./main.cpp -o main.out && ./main.out

nvc++ -fast -O2 ./main_downsample.cpp -o main2.out && ./main2.out
```


## LICENSE

NDT C++ specific code is distributed under Apache2.0 License.
The following extensions are included (*.cpp *.md)


NDT C++ specific data are distributed under MIT License.
The following extensions are included (*.txt)
