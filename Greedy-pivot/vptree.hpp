#include "file_output.hpp"
#include <algorithm>
#include <deque>
#include <cfloat>


// capacity of one node
const unsigned int capacity = 16;

// distance computation
float compute_distance(const data& d1, const data& d2) {

	float distance = 0;

	if (type == 0) {

		// L2 norm
		for (unsigned int i = 0; i < dimensionality; ++i) {
			const float temp = d1.pt[i] - d2.pt[i];
			distance += temp * temp;
		}
		distance = sqrt(distance);
	}
	else if (type == 1) {

		// L1 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += fabsf(d1.pt[i] - d2.pt[i]);
	}
	else if (type == 2) {

		// jaccard distance
		double overlap = 0;

		if (d1.set.size() <= d2.set.size()) {

			auto itr = d1.set.begin();
			while (itr != d1.set.end()) {

				if (d2.set.find(*itr) != d2.set.end()) ++overlap;

				++itr;
			}
		}
		else {

			auto itr = d2.set.begin();
			while (itr != d2.set.end()) {

				if (d1.set.find(*itr) != d1.set.end()) ++overlap;

				++itr;
			}
		}

		float jaccard = (float)(overlap / (d1.set.size() + d2.set.size() - overlap));
		distance = 1.0 - jaccard;
	}
	else if (type == 3) {

		// edit distance
		const size_t len1 = d1.str.size();
		const size_t len2 = d2.str.size();
		std::vector<std::vector<size_t>> d(len1 + 1, std::vector<size_t>(len2 + 1));

		d[0][0] = 0;
		for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
		for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

		for (unsigned int i = 1; i <= len1; ++i) {
			for (unsigned int j = 1; j <= len2; ++j) d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (d1.str[i - 1] == d2.str[j - 1] ? 0 : 1) });
		}

		distance = (float)d[len1][len2];
	}
	else if (type == 4) {

		// angular distance

		// inner product
		float ip = 0;

		// L2 norm
		float l2norm_1 = 0;
		float l2norm_2 = 0;

		for (unsigned int i = 0; i < dimensionality; ++i) {

			ip += d1.pt[i] * d2.pt[i];

			l2norm_1 += powf(d1.pt[i], 2.0);
			l2norm_2 += powf(d2.pt[i], 2.0);
		}

		l2norm_1 = sqrt(l2norm_1);
		l2norm_2 = sqrt(l2norm_2);

		float cosine_sim = ip / (l2norm_1 * l2norm_2);

		distance = acos(cosine_sim) / 3.1416;
	}
	else if (type == 5) {

		// L4 norm
		for (unsigned int i = 0; i < dimensionality; ++i) {
			const float temp = d1.pt[i] - d2.pt[i];
			distance += temp * temp * temp * temp;
		}
		distance = powf(distance, 0.25);
	}

	return distance;
}

float compute_distance(const data& d1, const data& d2, const float &threshold) {

	float distance = 0;

	if (type == 0) {

		// L2 norm
		for (unsigned int i = 0; i < dimensionality; ++i) {
			const float temp = d1.pt[i] - d2.pt[i];
			distance += temp * temp;
		}
		distance = sqrt(distance);
	}
	else if (type == 1) {

		// L1 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += fabsf(d1.pt[i] - d2.pt[i]);
	}
	else if (type == 2) {

		// jaccard distance
		double overlap = 0;

		if (d1.set.size() <= d2.set.size()) {

			auto itr = d1.set.begin();
			while (itr != d1.set.end()) {

				if (d2.set.find(*itr) != d2.set.end()) ++overlap;

				++itr;
			}
		}
		else {

			auto itr = d2.set.begin();
			while (itr != d2.set.end()) {

				if (d1.set.find(*itr) != d1.set.end()) ++overlap;

				++itr;
			}
		}

		float jaccard = (float)(overlap / (d1.set.size() + d2.set.size() - overlap));
		distance = 1.0 - jaccard;
	}
	else if (type == 3) {

		// edit distance
		const size_t len1 = d1.str.size();
		const size_t len2 = d2.str.size();
		std::vector<std::vector<size_t>> d(len1 + 1, std::vector<size_t>(len2 + 1));

		d[0][0] = 0;
		for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
		for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

		for (unsigned int i = 1; i <= len1; ++i) {
			for (unsigned int j = 1; j <= len2; ++j) d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (d1.str[i - 1] == d2.str[j - 1] ? 0 : 1) });
		}

		distance = (float)d[len1][len2];
	}
	else if (type == 4) {

		// angular distance

		// inner product
		float ip = 0;

		// L2 norm
		float l2norm_1 = 0;
		float l2norm_2 = 0;

		for (unsigned int i = 0; i < dimensionality; ++i) {

			ip += d1.pt[i] * d2.pt[i];

			l2norm_1 += powf(d1.pt[i], 2.0);
			l2norm_2 += powf(d2.pt[i], 2.0);
		}

		l2norm_1 = sqrt(l2norm_1);
		l2norm_2 = sqrt(l2norm_2);

		float cosine_sim = ip / (l2norm_1 * l2norm_2);

		distance = acos(cosine_sim) / 3.1416;
	}
	else if (type == 5) {

		// L4 norm
		for (unsigned int i = 0; i < dimensionality; ++i) {
			const float temp = d1.pt[i] - d2.pt[i];
			distance += temp * temp * temp * temp;
			
			if (i % 30 == 0) {
				if (threshold * threshold * threshold * threshold < distance) break;
			}
		}
		distance = powf(distance, 0.25);
	}

	return distance;
}

// declaration of "node"
class node;

// implementation of "node"
class node {

public:

	// identifier of vantage point
	int identifier = -1;

	// mean distane
	float distance_mean = 0;

	// max distance
	float distance_max = 0;

	// pointer to left node
	node* left = NULL;

	// pointer to right node
	node* right = NULL;

	// set of data object identifier (only for leaf node)
	std::vector<unsigned int> leaf;

	// #objects in left sub-tree
	unsigned int left_count = 0;

	// #objects in left sub-tree
	unsigned int right_count = 0;


	/***************/
	/* constructor */
	/***************/

	// standard
	node() {}

	// with identifier
	node(unsigned int id) {
		identifier = id;
	}

	/*******************/
	/* member function */
	/*******************/

	void update_id(unsigned int id) { identifier = id; }
};

class vp_tree{

	std::vector<node> vptree;

public:
	/***************/
	/* constructor */
	/***************/
	vp_tree(){}

private:
	/*******************/
	/* member function */
	/*******************/

	// select vantage point randomly
	unsigned int select_vp(const std::vector<unsigned int>& dataset_sub) {

		// uniform distribution
		std::mt19937 mt(0);
		std::uniform_int_distribution<> rnd(0, (int)dataset_sub.size() - 1);

		unsigned int idx = (unsigned int)rnd(mt);

		return idx;
	}

	// split
	void split(node *n, const std::vector<unsigned int>& dataset_sub) {

		if (dataset_sub.size() > capacity) {

			// decide vp
			unsigned int idx = select_vp(dataset_sub);

			// update node
			n->update_id(dataset[dataset_sub[idx]].identifier);

			// compute distance mean
			std::vector<float> distances;
			distances.resize(dataset_sub.size());

			#pragma omp parallel num_threads(core_no)
			{
				#pragma omp for schedule(static)
				for (unsigned int i = 0; i < dataset_sub.size(); ++i) {
					distances[i] = compute_distance(dataset[n->identifier], dataset[dataset_sub[i]]);
				}
			}
			std::vector<float> distances_cpy = distances;
			std::sort(distances_cpy.begin(), distances_cpy.end());

			n->distance_mean = distances_cpy[distances_cpy.size() / 2];
			n->distance_max = distances_cpy[distances_cpy.size() - 1];
			
			// make new subsets
			std::vector<unsigned int> dataset_sub_l, dataset_sub_r;
			for (unsigned int i = 0; i < dataset_sub.size(); ++i) {

				if (i != idx) {
					if (distances[i] < n->distance_mean) {
						dataset_sub_l.push_back(dataset_sub[i]);
					}
					else {
						dataset_sub_r.push_back(dataset_sub[i]);
					}
				}
			}

			// update left & right count
			n->left_count = (unsigned int)dataset_sub_l.size();
			n->right_count = (unsigned int)dataset_sub_r.size();

			// make nodes
			node left, right;
			vptree.push_back(left);
			vptree.push_back(right);

			// update children
			node* l = &vptree[vptree.size() - 2];
			node* r = &vptree[vptree.size() - 1];

			n->left = l;
			n->right = r;

			split(l, dataset_sub_l);
			split(r, dataset_sub_r);
		}
		else {

			// leaf node
			n->leaf = dataset_sub;
		}
	}

public:

	// build a vp-tree
	void build() {

		// make a set of identifier
		std::vector<unsigned int> id_set;
		for (unsigned int i = 0; i < dataset.size(); ++i) id_set.push_back(i);

		// reserve memory for vp-tree
		vptree.reserve(dataset.size());

		// make the root node
		node root;
		vptree.push_back(root);

		// recursive split
		split(&vptree[0], id_set);
	}

	// range search
	unsigned int range_search(data &d) {

		unsigned int result_size = 0;

		// prepare a queue
		std::deque<node*> queue;

		// insert root into the queue
		queue.push_back(&vptree[0]);

		// recursive search
		while (queue.size() > 0) {

			// pop the front
			node* n = queue[0];
			queue.pop_front();

			if (n->identifier == -1) {

				// leaf node
				float distance = 0;

				for (unsigned int i = 0; i < n->leaf.size(); ++i) {
                    
                    if (d.identifier != n->leaf[i]) {
                        distance = compute_distance(d, dataset[n->leaf[i]]);
                        ++d.distance_comp_count;
                        if (distance <= radius) ++result_size;
                    }

					if (result_size == k) break;
				}
				
				if (result_size == k) break;
			}
			else {

				// intermediate node

				// distance computation between query and vp
				float distance = compute_distance(d, dataset[n->identifier]);
                ++d.distance_comp_count;

				if (distance <= radius && (int)d.identifier != n->identifier) ++result_size;

                if (result_size == k) break;

				// children check
				if (distance + radius >= n->distance_mean) queue.push_back(n->right);
				if (distance - radius <= n->distance_mean) queue.push_back(n->left);
			}
		}

		return result_size;
	}
};
