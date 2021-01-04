#include "nndescent.hpp"
#include <deque>


// declaration of "node"
class node;

// implementation of "node"
class node {

public:

	// identifier of vantage point
	int identifier = -1;

	// mean distane
	float distance_mean = 0;

	// pointer to left node
	node* left = NULL;

	// pointer to right node
	node* right = NULL;

	// set of data object identifier (only for leaf node)
	std::vector<unsigned int> leaf;


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

	// capacity of one node
	const unsigned int capacity = 16;

    // seed
    unsigned int seed = 0;

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
		std::mt19937 mt(seed);
		std::uniform_int_distribution<> rnd(0, (int)dataset_sub.size() - 1);
        ++seed;

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

			for (unsigned int i = 0; i < dataset_sub.size(); ++i) distances[i] = compute_distance(dataset[n->identifier], dataset[dataset_sub[i]]);
			std::vector<float> distances_cpy = distances;
			std::sort(distances_cpy.begin(), distances_cpy.end());

			n->distance_mean = distances_cpy[distances_cpy.size() / 2];

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
	void range_search(data &d) {

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
                        d.update_neighbor_list_(distance, n->leaf[i]);
                    }
				}
			}
			else {

				// intermediate node

				// distance computation between query and vp
				float distance = compute_distance(d, dataset[n->identifier]);

				if (distance <= d.threshold && (int)d.identifier != n->identifier) d.update_neighbor_list_(distance, n->identifier);

				// children check
				if (distance + d.threshold >= n->distance_mean) queue.push_back(n->right);
				if (distance - d.threshold <= n->distance_mean) queue.push_back(n->left);
			}
		}
	}
};