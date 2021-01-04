#include "vptree.hpp"


// ANN search
void ann_greedy(const unsigned int destination, const unsigned int start_point, float &threshold, unsigned int &result_id) {

    // init queue
    std::deque<std::pair<unsigned int, unsigned int>> queue;
    queue.push_back({start_point, 1});

    // visited nodes
    std::unordered_set<unsigned int> visited;
    visited.insert(start_point);

    // max hop count
    const unsigned int hop_max = 10;

    while (queue.size() > 0) {

        // pop the front
        const unsigned int id = queue[0].first;
        unsigned int hop = queue[0].second;
        queue.pop_front();

        // prepare distance array
        std::vector<float> distance_array;
        distance_array.resize(dataset[id].edges.size());

        // distance computation
        #pragma omp parallel num_threads(core_no)
        {
            #pragma omp for schedule(dynamic)
            for (unsigned int i = 0; i < dataset[id].edges.size(); ++i) {

                // init
                distance_array[i] = FLT_MAX;

                // get id_
                const unsigned int id_ = dataset[id].edges[i];

                // compute only when non-visited
                if (visited.find(id_) == visited.end()) distance_array[i] = compute_distance(dataset[id_], dataset[destination]);
            }
        }

        float distance_min = FLT_MAX;
        unsigned int idx = 0;
        for (unsigned int i = 0; i < dataset[id].edges.size(); ++i) {

            // update visited
            visited.insert(dataset[id].edges[i]);

            // get min dist
            if (distance_array[i] < distance_min) {
                idx = i;
                distance_min = distance_array[i];
            }
        }

        // result update
        if (distance_min < threshold) {
            threshold = distance_min;
            result_id = dataset[id].edges[idx];
            if (hop < hop_max) queue.push_back({dataset[id].edges[idx], hop + 1});
        }
    }
}

// BFS
void breadth_first_search_test() {

    for (unsigned int i = 0; i < dataset.size(); ++i) dataset[i].update_flag = 0;

    // init
    std::deque<unsigned int> queue;
    queue.push_back(0);
    dataset[0].update_flag = 1;

    // #accessed nodes
    unsigned int count = 1;

    while (queue.size() > 0) {

        // get front
        unsigned int id = queue[0];
        queue.pop_front();

        for (unsigned int i = 0; i < dataset[id].edges.size(); ++i) {

            // get id_
            const unsigned int id_ = dataset[id].edges[i];

            if (dataset[id_].update_flag == 0) {

                dataset[id_].update_flag = 1;
                queue.push_back(id_);
                ++count;
            }
        }
    }

    std::cout << " count: " << count << "\n\n";
}