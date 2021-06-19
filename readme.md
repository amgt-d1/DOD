## Introduction
* This repository provides codes for bulding a MRPG and a distance-based outlier detection algorithm on a MRPG.
* Our algorithm supports metric space.
    * Our codes implement L2 (Euclidean distance), L1 (Manhattan distance), Jaccard distance, Edit distance, angular distance, and L4 distance by default.
    * The other distance functions are free to add.
* The details about our algorithm can be seen from our SIGMOD2021 paper, Fast and Exact Outlier Detection in Metric Spaces: A Proximity Graph-based Approach.

## Requirement
* Linux OS (Ubuntu).
   * The others have not been tested.
* `g++ 7.4.0` (or higher version) and `Openmp`.

## How to use
* Before running our DOD algorithm, build an MRPG.
* Parameter configuration can be done via txt files in `parameter` directory.
* Data files have to be at `dataset` directory.
   * You can implement data input in as you like manner at input_data() function in data.hpp.
   * Now dataset directory contains a dummy file only.


### MRPG
* Create `result/graph` directory.
* Compile: `g++ -O3 -o mrpg.out main.cpp -std=c++11 -fopenmp -Wall`.
* Run: `./mrpg.out`.

### Greedy-pivot (DOD algorithm)
* Create `result` directory.
* Compile: `g++ -O3 -o greedy-pivot.out main.cpp -std=c++11 -fopenmp -Wall`.
* Run: `./greedy-pivot.out`.
* If you test a low-dimensional dataset, you may enable VP-tree based verification (by setting `mode = 1` in main.cpp).
    * By default, verification is done by a linear scan.


## Citation
If you use our implementation, please cite the following paper.
``` 
@inproceedings{amagata2021dod,  
    title={Fast and Exact Outlier Detection in Metric Spaces: A Proximity Graph-based Approach},  
    author={Amagata, Daichi and Onizuka, Makoto and Hara, Takahiro},  
    booktitle={SIGMOD},  
    pages={36--48},  
    year={2021}  
}
``` 

## License
Copyright (c) 2020 Daichi Amagata  
This software is released under the [MIT license](https://github.com/amgt-d1/DOD/blob/main/license.txt).
