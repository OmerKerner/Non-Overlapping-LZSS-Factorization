# Non-Overlapping LZSS Factorization
Calculating non-overlapping LZSS factorization using sdsl-lite based on a paper by Dominik Köppl - [Non-Overlapping LZ77 Factorization and LZ78 Substring Compression Queries with Suffix Trees](https://doi.org/10.3390/a14020044)

## Compile sdsl-lite by adding to the install.sh arguments to cmake:
```
cmake -DCMAKE_C_FLAGS="-fPIC" -DCMAKE_CXX_FLAGS="-fPIC" ...
```

## Compile the executable:
```
g++ -std=c++11 -DNDEBUG -O3 -I/home/user/include -L/home/user/lib noLZSS.cpp -lsdsl -ldivsufsort -ldivsufsort64 -o noLZSS
```

## Compile the shared object with:
```
g++ -shared -fPIC -std=c++11 -DNDEBUG -O3 -I/home/user/include -L/home/user/lib noLZSS.cpp -lsdsl -ldivsufsort -ldivsufsort64 -o noLZSS.so
```

## Using the shared object in python:
```python
import ctypes
noLZSS = ctypes.CDLL('./noLZSS.so')
noLZSS.nolzss.argtypes = [ctypes.c_char_p]
noLZSS.nolzss.restype = ctypes.c_int
noLZSS.nolzss('ABRACADABRA'.encode())
```
