#include "file_output.hpp"
#include <algorithm>
#include <omp.h>


// random graph generation
void knn_list_initialize() {

    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(static)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            // random id generator
            std::default_random_engine engine(i);
            std::uniform_int_distribution<> rnd(0, (unsigned int)dataset.size() - 1);

            float distance = 0;

            while (dataset[i].neighbor_list.size() < degree) {

                unsigned int id = rnd(engine);
                if (id != i) {

                    // knn list update
                    distance = compute_distance(dataset[i], dataset[id]);
                    dataset[i].update_neighbor_list(distance, id);
                }
            }

            // init update flag
            dataset[i].update_flag = 0;
        }
    }
}

// reverse kNN update
void reverse(const std::unordered_set<unsigned int> &not_updated, std::vector<bool> &flag_new_reverse) {

    // init
    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(static)
        for (unsigned int i = 0; i < dataset.size(); ++i) {
            dataset[i].reverse_knn_list.clear();
        }
    }

    // reverse edge computation
    for (unsigned int i = 0; i < dataset.size(); ++i) {

        if (not_updated.find(i) == not_updated.end() || flag_new_reverse[i] == 1) {

            auto it = dataset[i].neighbor_list.begin();
            while (it != dataset[i].neighbor_list.end()) {
                dataset[it->second].reverse_knn_list.insert(i);
                ++it;
            }
        }
    }
}

// get neighbors
void get_neighbors(std::unordered_set<unsigned int> &not_updated, std::vector<bool> &flag_new_reverse) {
  
    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            // init
            dataset[i].neighbors.clear();
        
            // init by reverse kNN
            dataset[i].neighbors = dataset[i].reverse_knn_list;

            // + kNN list
            auto it = dataset[i].neighbor_list.begin();
            while (it != dataset[i].neighbor_list.end()) {
                if (not_updated.find(it->second) == not_updated.end() || flag_new_reverse[i] == 1) dataset[i].neighbors.insert(it->second);
                ++it;
            }

            // init flag
            flag_new_reverse[i] = 0;
        }
    }

    // clear not_updated
    not_updated.clear();
}

// NN descent
void nn_descent() {

    start = std::chrono::system_clock::now();

    // prepare not-updated data
    std::unordered_set<unsigned int> not_updated;

    // prepare reverse kNN update flags
    std::vector<bool> flag_new_reverse;
    flag_new_reverse.resize(dataset.size());

    // iterative improvement
    while (1) {

        // init update counter
        double count_update = 0;

        // reverse kNN computation
        reverse(not_updated, flag_new_reverse);

        // get neighbors
        get_neighbors(not_updated, flag_new_reverse);

        #pragma omp parallel num_threads(core_no)
	    {
		    #pragma omp for schedule(dynamic) reduction(+:count_update)
            for (unsigned int i = 0; i < dataset.size(); ++i) {

                // init local update counter
                double count = 0;

                float distance = 0;

                // get already accessed objects
                std::unordered_set<unsigned int> accessed;
                dataset[i].get_accessed(accessed);

                // iterative kNN update
                auto it = dataset[i].neighbors.begin();
                while (it != dataset[i].neighbors.end()) {

                    auto it_nei = dataset[*it].neighbors.begin();
                    while (it_nei != dataset[*it].neighbors.end()){

                        // not accessed object yet
                        if (accessed.find(*it_nei) == accessed.end()) {

                            // mark as accessed
                            accessed.insert(*it_nei);

                            // distance computation
                            distance = compute_distance(dataset[i], dataset[*it_nei]);

                            // kNN update
                            if (dataset[i].update_neighbor_list(distance, *it_nei)) {
                                ++count;

                                // update reverse flag
                                flag_new_reverse[*it_nei] = 1;
                            }
                        }
                        ++it_nei;
                    }
                    ++it;
                }
                count_update = count_update + count;
            }
        }

        // check update-flag
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            if (dataset[i].update_flag == 0) not_updated.insert(i);

            // init update flag
            dataset[i].update_flag = 0;
        }

        // if count_update is sufficiently small, iteration is over
        if (count_update < 1) break;
    }

    end = std::chrono::system_clock::now();
	cpu_update_edge = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << " NNDescent time: " << cpu_update_edge << "[msec]\n\n";
}
