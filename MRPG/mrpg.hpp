#include "graph_search.hpp"


std::unordered_set<unsigned int> pivot_set;


// VP-tree like split
void split(const std::vector<unsigned int>& dataset_sub, const unsigned int seed, const bool direction) {

    const unsigned int coefficient = 2;

    if (dataset_sub.size() > degree * coefficient) {

        // random id generator
        std::default_random_engine engine(seed);
        std::uniform_int_distribution<> rnd(0, (unsigned int)dataset_sub.size() - 1);

        // get idx
        const unsigned int idx = rnd(engine);

        // get id
        const unsigned int id = dataset_sub[idx];

        // check pivot
        if (dataset_sub.size() / 2 <= degree * coefficient) dataset[id].pivot_flag = 1;

        // compute distance to the vp
		std::vector<float> distances;
        distances.resize(dataset_sub.size());

        #pragma omp parallel num_threads(core_no)
        {
            #pragma omp for schedule(static)
            for (unsigned int i = 0; i < dataset_sub.size(); ++i) {
                if (id != dataset_sub[i]) {
                    distances[i] = compute_distance(dataset[id], dataset[dataset_sub[i]]);
                    distances[i] += add_noise(dataset[id], dataset[dataset_sub[i]]);
                }
            }
        }

		std::vector<unsigned int> dataset_sub_l, dataset_sub_r;

        while (1) {

            // get mean distance
            std::vector<float> distances_cpy = distances;
            std::sort(distances_cpy.begin(), distances_cpy.end());
            float mean = distances_cpy[distances_cpy.size() / 2];

            // make new subsets
            for (unsigned int i = 0; i < dataset_sub.size(); ++i) {

                if (distances[i] <= mean) {
                    dataset_sub_l.push_back(dataset_sub[i]);
                }
                else {
                    dataset_sub_r.push_back(dataset_sub[i]);
                }
            }

            if (dataset_sub_r.size() > 0) {
                break;
            }
            else{

                // reset
                dataset_sub_l.clear();

                // add a random value
                std::default_random_engine engine_(1);
                std::uniform_real_distribution<> rnd_(0, 0.999);
                for (unsigned int i = 0; i < distances.size(); ++i) distances[i] += rnd_(engine_); 
            }
        }

        // recursive split
        split(dataset_sub_l, seed + 1, 1);
		split(dataset_sub_r, seed + 2, 0);
    }
    else {

        // join only left node
        if (direction) {

            // prepare distance
            float distance = 0;

            std::vector<std::vector<std::pair<float, std::pair<unsigned int, unsigned int>>>> join;
            join.resize(core_no);

            // parallel distance computation
            #pragma omp parallel num_threads(core_no)
            {
                #pragma omp for schedule(static)
                for (unsigned int i = 0; i < dataset_sub.size() - 1; ++i) {

                    for (unsigned int j = i + 1; j < dataset_sub.size(); ++j) {

                        distance = compute_distance(dataset[dataset_sub[i]], dataset[dataset_sub[j]]);
                        join[omp_get_thread_num()].push_back({distance, {dataset_sub[i], dataset_sub[j]}});
                    }
                }
            }

            // neighbor list update
            for (unsigned int i = 0; i < core_no; ++i) {

                for (unsigned int j = 0; j < join[i].size(); ++j) {

                    // get id
                    unsigned int id = join[i][j].second.first;
                    unsigned int id_ = join[i][j].second.second;

                    // get distance
                    float distance = join[i][j].first;

                    dataset[id].update_neighbor_list(distance, id_);
                    dataset[id_].update_neighbor_list(distance, id);
                }
            }
        }
    }
}

// initial edge determination by vp-tree like partitioning
void init_by_split() {

    unsigned int iteration_num = 15;

     // make a set of identifier
    std::vector<unsigned int> id_set;
    for (unsigned int i = 0; i < dataset.size(); ++i) id_set.push_back(i);

    // iterative spliting
    while (iteration_num > 0) {

        // split
        split(id_set, iteration_num, 1);

        --iteration_num;
    }
}

// init MRPG
void initialize() {

    start = std::chrono::system_clock::now();

    // vp-tree like partitioning init
    init_by_split();

    // random edge
    knn_list_initialize();

    end = std::chrono::system_clock::now();
    cpu_initial_edge = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " Initialization time: " << cpu_initial_edge << "[msec]\n\n";
}

// bool function for sorting pairs
bool compare_by_first(const std::pair<float, int> a, const std::pair<float, int> b) {
    return a.first > b.first;
}

// exact kNN evaluation in top-x% objects
void evaluate_knn_exact() {

    start = std::chrono::system_clock::now();

    // compute sum of kNN distance
    std::vector<std::pair<float, unsigned int>> aggregate_distance_set;
    aggregate_distance_set.resize(dataset.size());

    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(static)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            float aggregate_distance = 0;
            auto it = dataset[i].neighbor_list.begin();
            while (it != dataset[i].neighbor_list.end()) {
                aggregate_distance += it->first;
                ++it;
            }

            aggregate_distance_set[i] = {aggregate_distance, i};

            // copy neighbors
            dataset[i].copy_neighbors();
        }
    }

    // sort by descending order of aggregate distance
    std::sort(aggregate_distance_set.begin(), aggregate_distance_set.end(), compare_by_first);

    // determine the number of exact kNN evaluations
    const unsigned int exact_size = 50000;  // this incurs trade-off on pre-processing time and online time 

    // vp-tree build
    vp_tree vpt;
    unsigned int mode = 0;
    if (dataset_id == 0 || dataset_id == 1 || dataset_id == 3 || dataset_id == 6) mode = 1;
    if (mode == 1) vpt.build();

    // exact kNN evaluation
    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(static)
        for (unsigned int i = 0; i < exact_size; ++i) {

            // get id
            unsigned int id = aggregate_distance_set[i].second;

            // flag update
            dataset[id].exact_flag = 1;

            if (mode == 0) {    // scan-base

                // prepare distance
                float distance = 0;

                // init threshold
                dataset[id].threshold = FLT_MAX;

                // scan
                for (unsigned int j = 0; j < dataset.size(); ++j) {

                    if (j != id) {

                        distance = compute_distance(dataset[id], dataset[j], dataset[id].threshold);
                        dataset[id].update_neighbor_list_(distance, j);
                    }
                }
            }
            else {  // vp-tree base

                // init threshold
                dataset[id].threshold = FLT_MAX;

                // kNN search
                vpt.range_search(dataset[id]);
            }

            // copy neighbors
            dataset[id].copy_neighbors();
        }
    }

    end = std::chrono::system_clock::now();
    cpu_exact_knn = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " Exact kNN computation time: " << cpu_exact_knn << "[msec]\n\n";
}

// convert to undirected graph
void convert_to_undirected_graph() {

    start = std::chrono::system_clock::now();

    for (unsigned int i = 0; i < dataset.size(); ++i) {

        // init update flag
        dataset[i].update_flag = 0;

        // get pivot set
        if (dataset[i].pivot_flag) pivot_set.insert(i);

        // get reverse neighbors
        if (dataset[i].exact_flag == 0) {
            auto it = dataset[i].neighbor_list.begin();
            while (it != dataset[i].neighbor_list.end()) {

                dataset[it->second].neighbors.insert(i);
                ++it;
            }
        }
        else {
            unsigned int count = 0;
            auto it = dataset[i].neighbor_list.begin();
            while (it != dataset[i].neighbor_list.end()) {

                dataset[it->second].neighbors.insert(i);
                ++it;
                ++count;
                if (count == degree) break;
            }
        }
    }

    // convert to vector
    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(static)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            auto it = dataset[i].neighbors.begin();
            while (it != dataset[i].neighbors.end()) {
                dataset[i].edges.push_back(*it);
                ++it;
            }
        }
    }

    end = std::chrono::system_clock::now();
    cpu_reverse = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " Reverse edge addition time: " << cpu_reverse << "[msec]\n\n";
}

// connect graph
void connect_graph() {

    start = std::chrono::system_clock::now();

    // sample size
    unsigned int seed = 0;

    // prepare samples
    std::vector<unsigned int> samples;

    // copy pivot_set & dataset
    std::unordered_set<unsigned int> pivot_set_copy = pivot_set;
    std::unordered_set<unsigned int> dataset_copy;
    for (unsigned int i = 0; i < dataset.size(); ++i) dataset_copy.insert(i);

    while (dataset_copy.size() > 0) {

        // get start point
        auto it = pivot_set_copy.begin();
        unsigned int sp = *it;

        // remove start point
        pivot_set_copy.erase(sp);
        dataset_copy.erase(sp);

        // BFS
        std::deque<unsigned int> queue;
        queue.push_back(sp);
        dataset[sp].update_flag = 1;

        while (queue.size() > 0) {

            // pop the front
            unsigned int id = queue[0];
            queue.pop_front();

            for (unsigned int i = 0; i < dataset[id].edges.size(); ++i) {

                // get id
                unsigned int id_ = dataset[id].edges[i];

                if (dataset[id_].update_flag == 0) {

                    // update flag
                    dataset[id_].update_flag = 1;

                    // erase from dataset copy
                    dataset_copy.erase(id_);

                    // erase from pivot set copy & sample
                    if (dataset[id_].pivot_flag) {
                        pivot_set_copy.erase(id_);
                        samples.push_back(id_);
                    }

                    // update queue
                    queue.push_back(id_);
                }
            }
        }

        // connection
        if (dataset_copy.size() > 0) {

            // if no more pivot, add it
            if (pivot_set_copy.size() == 0) {
                auto it = dataset_copy.begin();
                dataset[*it].pivot_flag = 1;
                pivot_set.insert(*it);
                pivot_set_copy.insert(*it);
            }

            // get non-connected pivot
            it = pivot_set_copy.begin();
            sp = *it;

            // min dist destination
            float distance_to_destination_min = FLT_MAX;
            unsigned int dest = sp;
            unsigned int connect = sp;

            // determine iteration number
            unsigned int iteration = 250;

            // random generator
            std::uniform_int_distribution<> rnd_idx(0, samples.size() - 1);

            while (iteration > 0) {

                // determine destination
                std::mt19937 mt(seed);
                const unsigned int destination = samples[rnd_idx(mt)];

                // update seed
                ++seed;

                // threshold
                float ths = compute_distance(dataset[sp], dataset[destination]);

                // Approximate NN Search
                unsigned int result_id = sp;
                ann_greedy(destination, sp, ths, result_id);

                // update destination
                if (distance_to_destination_min > ths) {
                    distance_to_destination_min = ths;
                    dest = destination;
                    connect = result_id;
                }

                --iteration;
            }

            // connect result to destination (+ reverse)
            dataset[connect].edges.push_back(dest);
            dataset[dest].edges.push_back(connect);
        }
    }

    end = std::chrono::system_clock::now();
    cpu_connect = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " Graph connection time: " << cpu_connect<< "[msec]\n\n";
}

// get non-monotonic path
void get_non_monotonic(const unsigned int id, const unsigned int st, const unsigned int max_hop, std::multimap<float, unsigned int> &update_list, std::unordered_set<unsigned int> &unnecessary) {

    // determine inital distance
    float init_distance = FLT_MAX;
    if (id == st) init_distance = 0;

    // init queue
    std::deque<std::pair<unsigned int, std::pair<float, unsigned int>>> queue;  // next id, <prev-dist., hop-count>
    queue.push_back({st, {init_distance, 0}});

    // init visit
    std::unordered_set<unsigned int> visit;

    // map of prev-distance and distance
    std::unordered_map<unsigned int, std::pair<float, float>> prev_dist_map;

    // traversal
    while (queue.size() > 0) {

        // get id
        const unsigned int id_ = queue[0].first;

        // get hop-count
        const unsigned int hop = queue[0].second.second;

        // get prev. dist.
        float prev_dist = queue[0].second.first;

        // pop front
        queue.pop_front();

        // update prev-distance (to min.)
        if (prev_dist_map.find(id_) == prev_dist_map.end()) prev_dist_map[id_] = {prev_dist, 0};
        if (prev_dist_map[id_].first > prev_dist) prev_dist_map[id_].first = prev_dist;

        // visit check
        if (visit.find(id_) == visit.end()) {

            // mark as visit
            visit.insert(id_);

            // compute distance
            float distance = compute_distance(dataset[id], dataset[id_]);
            prev_dist_map[id_].second = distance;

            // update queue
            if (hop < max_hop) {
                for (unsigned int j = 0; j < dataset[id_].edges.size(); ++j) queue.push_back({dataset[id_].edges[j], {distance, hop + 1}});
            }
        }
    }

    // get non-monotonic path
    auto itr_map = prev_dist_map.begin();
    while (itr_map  != prev_dist_map.end()) {
        if (itr_map->second.first <= itr_map->second.second) {

            // mark as unnecessary
            unnecessary.insert(itr_map->first);

            // erase
            prev_dist_map.erase(itr_map++);
        }
        else {
            update_list.insert({itr_map->second.second, itr_map->first});
            ++itr_map;
        }
    }
}

// approx. remove detours
void remove_detour() {

    start = std::chrono::system_clock::now();

    std::mt19937 mt(0);
    std::uniform_real_distribution<> prob(0, 1.0);

    // prepare pivot-set
    std::vector<unsigned int> pivots;

    // get start objects
    std::vector<unsigned int> start_set;
    for (unsigned int i = 0; i < dataset.size(); ++i) {

        if (dataset[i].pivot_flag) {

            if (i % 2 == 0) {
                if (dataset[i].exact_flag == 0) start_set.push_back(i);
            }
            else {
                pivots.push_back(i);
            }
        }
        else {

            // check pivot connection
            bool flag = 0;
            if (dataset[i].exact_flag == 0) {

                for (unsigned int j = 0; j < dataset[i].edges.size(); ++j) {

                    // get id
                    const unsigned int id = dataset[i].edges[j];
                    if (dataset[id].pivot_flag) {
                        flag = 1;
                        break;
                    }
                }
            }
            else {
                flag = 1;
            }

            // sample if no connection to pivots
            if (flag == 0) {

                start_set.push_back(i);

                // make new pivot
                if (prob(mt) <= 0.2) dataset[i].pivot_flag = 1;
            }
        }
    }

    // determine sample size for detour removing
    const unsigned int sample_size_test = std::min(degree * 100, (unsigned int)(pivots.size() * 0.75));
    const unsigned int sample_size = degree * 3;

    // determine constant for 3-hop traversal
    const float constant = 2.5;

    // necessary edge lists
    std::vector<std::multimap<float, unsigned int>> add_edge_list;
    add_edge_list.resize(start_set.size());

    // remove detours
    #pragma omp parallel num_threads(core_no)
    {
        #pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < start_set.size(); ++i) {

            // get id (start object)
            const unsigned int id = start_set[i];

            // unnecessary nodes
            std::unordered_set<unsigned int> unnecessary;

            /*** 3-hop investigation ***/
            get_non_monotonic(id, id, 3, add_edge_list[i], unnecessary);

            // reduce non-monotonic path to constant
            while (add_edge_list[i].size() > degree * constant) {
                add_edge_list[i].erase(--add_edge_list[i].end());
            }

            /*** get samples ***/

            // init sample by neighbors
            std::unordered_set<unsigned int> sample_set, original_set;
            for (unsigned int j = 0; j < dataset[id].edges.size(); ++j) {
                sample_set.insert(dataset[id].edges[j]);
                original_set.insert(dataset[id].edges[j]);
                unnecessary.insert(dataset[id].edges[j]);
            }

            // random sampling with replacement
            std::mt19937 mt(i);
            std::uniform_int_distribution<> rnd_idx(0, pivots.size() - 1);
            while (sample_set.size() < sample_size_test) {

                // get idx
                const unsigned idx = rnd_idx(mt);
                if (pivots[idx] != id) sample_set.insert(pivots[idx]);
            }

            // sort distance between start and samples
            std::multimap<float, unsigned int> sample_distance;
            auto it_ = sample_set.begin();
            while (it_ != sample_set.end()) {
                sample_distance.insert({compute_distance(dataset[id], dataset[*it_]), *it_});
                ++it_;
            }

            // get first prev. distance
            auto it = sample_distance.begin();
            for (unsigned int j = 1; j < dataset[id].edges.size(); ++j) ++it;
            float distance_prev = it->first;
            ++it;

            // get distance difference max position
            float distance_difference = 0;
            float distance_difference_max = 0;
            float distance_difference_avg= 0;
            unsigned int idx = dataset[id].edges.size() - 1;
            unsigned int pos = dataset[id].edges.size() - 1;
            while (idx < sample_size) {

                // update max gap
                distance_difference = it->first - distance_prev;
                distance_difference_avg += (it->first - distance_prev);
                if (distance_difference > distance_difference_max) {
                    distance_difference_max = distance_difference;
                    pos = idx;
                }

                // update prev. distance
                distance_prev = it->first;
                ++idx;
                ++it;
            }
            distance_difference_avg /= sample_size - dataset[id].edges.size() + 1;
            if (distance_difference_max < distance_difference_avg * 2) pos = dataset[id].edges.size() - 1;

            // convert to vector (& remove original edges)
            unsigned int iteration = 0;
            std::vector<unsigned int> samples;
            it = sample_distance.begin();
            while (iteration < pos + 1) {
                if (original_set.find(it->second) == original_set.end()) samples.push_back(it->second);
                ++it;
                ++iteration;
            }

            // get non-monotonic path
            for (unsigned int j = 0; j < samples.size(); ++j) get_non_monotonic(id, samples[j], 2, add_edge_list[i], unnecessary);

            // clean add_edge_list
            auto itr_edge = add_edge_list[i].begin();
            while (itr_edge != add_edge_list[i].end()) {

                bool f = 0;

                if (unnecessary.find(itr_edge->second) != unnecessary.end()) {
                    f = 1;
                }
                else{

                    auto itr_next = itr_edge;
                    ++itr_next;
                    if (itr_next != add_edge_list[i].end()) {
                        if (itr_edge->second == itr_next->second) f = 1;
                    }
                }

                if (f) {
                    add_edge_list[i].erase(itr_edge++);
                }
                else {
                    ++itr_edge;
                }
            }
        }
    }

    // init neighbors
    for (unsigned int i = 0; i < dataset.size(); ++i) {

        dataset[i].neighbors.clear();
        for (unsigned int j = 0; j < dataset[i].edges.size(); ++j) dataset[i].neighbors.insert(dataset[i].edges[j]);
    }

    // add edges
    for (unsigned int i = 0; i < start_set.size(); ++i) {

        // get id
        const unsigned int id = start_set[i];

        // prepare prev. id
        unsigned int id_prev = id;

        auto it = add_edge_list[i].begin();
        while (it != add_edge_list[i].end()) {

            if (dataset[id_prev].neighbors.find(it->second) == dataset[id_prev].neighbors.end()) {
                dataset[id_prev].neighbors.insert(it->second);
                dataset[id_prev].edges.push_back(it->second);
            }

            if (dataset[it->second].neighbors.find(id_prev) == dataset[it->second].neighbors.end()) {
                dataset[it->second].neighbors.insert(id_prev);
                dataset[it->second].edges.push_back(id_prev);
            }

            id_prev = it->second;

            ++it;
        }
    }

    end = std::chrono::system_clock::now();
    cpu_detour = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " detours removing time: " << cpu_detour << "[msec]\n\n";
}

// remove unnecessary edges
void remove_edges() {

    start = std::chrono::system_clock::now();

    // convert to hash tables
    #pragma omp parallel num_threads(core_no)
    {
        #pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            // clear
            dataset[i].neighbors.clear();
            
            // insert
            for (unsigned int j = 0; j < dataset[i].edges.size(); ++j) dataset[i].neighbors.insert(dataset[i].edges[j]);
        }
    }

    // removing edges
    #pragma omp parallel num_threads(core_no)
    {
        #pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < dataset.size(); ++i) {

            if (dataset[i].pivot_flag == 0 && dataset[i].exact_flag == 0) {

                // prepare unnecessary
                std::unordered_set<unsigned int> unnecessary;

                // prepare necessary
                std::unordered_set<unsigned int> necessary;

                // clear edges
                dataset[i].edges.clear();

                auto it = dataset[i].neighbors.begin();
                while (it != dataset[i].neighbors.end()) {

                    if (dataset[*it].pivot_flag) {

                        // pivot must be inserted
                        necessary.insert(*it);

                        auto itr = dataset[i].neighbors.begin();
                        while (itr != dataset[i].neighbors.end()) {

                            if (dataset[*it].neighbors.find(*itr) != dataset[*it].neighbors.end()) unnecessary.insert(*itr);
                            ++itr;
                        }
                    }
                    ++it;
                }

                // remove unnecessary
                it = unnecessary.begin();
                while (it != unnecessary.end()) {
                    dataset[i].neighbors.erase(*it);
                    ++it;
                }

                // keep pivots
                it = necessary.begin();
                while (it != necessary.end()) {
                    dataset[i].neighbors.insert(*it);
                    ++it;
                }

                // edge insertion
                it = dataset[i].neighbors.begin();
                while (it != dataset[i].neighbors.end()) {
                    dataset[i].edges.push_back(*it);
                    ++it;
                }
            }
        }
    }

    end = std::chrono::system_clock::now();
    cpu_remove = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " edge removing time: " << cpu_remove << "[msec]\n\n";
}


// build an MRPG
void build_mrpg() {

    // init graph
    initialize();

    // build AkNN graph
    nn_descent();

    // exact kNN computation
    evaluate_knn_exact();

    // convert to undirected edge
    convert_to_undirected_graph();

    // connect graph
    connect_graph();

    // detour removing
    remove_detour();

    // remove unnecessary edges
    remove_edges();
}