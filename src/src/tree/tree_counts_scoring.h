#ifndef COUNTS_SCORING_H
#define COUNTS_SCORING_H

#include "../input_data/input_data.h"
#include "../parameters/parameters.h"
#include "event_tree.h"

template <class Real_t> class CountsDispersionPenalty {
private:
  const size_t DEFAULT_CACHE_SIZE = 1000;

  std::vector<std::vector<Real_t>> sum_counts;
  std::vector<std::vector<Real_t>> squared_counts;
  std::vector<Real_t> counts_score_length_of_bin;
  Real_t all_bins_count;

  // Specifies which cell bin pairs have already bin assigned to a cluster
  std::vector<std::vector<bool>> bin_bitmap;
  // Will store current clustering induced by node
  std::vector<size_t> event_clusters;
  size_t max_cluster_id{0};
  // Used for persisting clusters induced by nodes
  std::vector<std::vector<size_t>> clusters_cache;
  size_t cache_id = 0;

  void save_clustering_in_cache(std::vector<size_t> &cluster, size_t cache_id) {
    if (cache_id >= clusters_cache.size()) {
      clusters_cache.resize(clusters_cache.size() + DEFAULT_CACHE_SIZE);
    }
    clusters_cache[cache_id] = cluster;
  }

  std::vector<size_t> get_clustering_from_cache(size_t cache_id) {
    return clusters_cache[cache_id];
  }

  /**
   * Calculates 1/Z * sum_{i= 0}^n (x_i - m)^2
   * Where @expected_mean = m, @sum_of_values = sum x_i,
   * @summed_squares_of_values = sum x_i^2, @number_of_values = n, @normalizer =
   * Z
   */
  inline Real_t calculate_l2_penalty(Real_t expected_mean, Real_t sum_of_values,
                                     Real_t summed_squares_of_values,
                                     Real_t number_of_values,
                                     Real_t normalizer) {
    return ((summed_squares_of_values - 2 * expected_mean * sum_of_values) /
                (normalizer) +
            expected_mean * expected_mean * number_of_values / normalizer);
  }

  /**
   * @brief Move all cells attached to @child to @parent
   */
  void move_cells_to_parent(EventTree::NodeHandle child,
                            EventTree::NodeHandle parent,
                            std::map<Event, std::set<size_t>> &attachment) {
    if (attachment.find(child->label) == attachment.end()) {
      return;
    }
    if (attachment.find(parent->label) == attachment.end()) {
      attachment[parent->label] = std::set<size_t>();
    }
    attachment[parent->label].insert(attachment[child->label].begin(),
                                     attachment[child->label].end());
    attachment.erase(child->label);
  }

  Real_t calculate_penalty_for_bins_at_node(
      EventTree::NodeHandle node,
      std::map<Event, std::set<size_t>> &attachment) {
    if (attachment.find(node->label) == attachment.end()) {
      return 0.0;
    }
    std::map<size_t, Real_t> cluster_to_counts_sum;
    std::map<size_t, Real_t> cluster_to_bin_count;
    std::map<size_t, Real_t> cluster_to_squared_counts_sum;

    Real_t result = 0.0;
    for (size_t i = node->label.first; i < node->label.second; i++) {
      cluster_to_counts_sum[event_clusters[i]] = 0.0;
      cluster_to_bin_count[event_clusters[i]] = 0.0;
      cluster_to_squared_counts_sum[event_clusters[i]] = 0.0;
    }

    for (auto cell : attachment[node->label]) {
      for (size_t bin = node->label.first; bin < node->label.second; bin++) {
        if (!bin_bitmap[cell][bin]) {
          cluster_to_counts_sum[event_clusters[bin]] += sum_counts[cell][bin];
          cluster_to_squared_counts_sum[event_clusters[bin]] +=
              squared_counts[cell][bin];
          cluster_to_bin_count[event_clusters[bin]] +=
              counts_score_length_of_bin[bin];
        }
      }
      std::fill(bin_bitmap[cell].begin() + node->label.first,
                bin_bitmap[cell].begin() + node->label.second, true);
    }

    for (const auto &cluster : cluster_to_counts_sum) {
      if (cluster_to_bin_count[cluster.first] > 0.0) {
        auto bin_count = cluster_to_bin_count[cluster.first];
        Real_t mean_count = cluster.second / bin_count;
        result +=
            COUNTS_SCORE_CONSTANT_0 *
            calculate_l2_penalty(mean_count, cluster.second,
                                 cluster_to_squared_counts_sum[cluster.first],
                                 bin_count, all_bins_count);
        if (mean_count >= NEUTRAL_CN - 0.5 && mean_count < NEUTRAL_CN + 0.5) {
          result += COUNTS_SCORE_CONSTANT_1 * bin_count / all_bins_count;
        }
      }
    }
    return result;
  }

  void update_clusters(std::vector<size_t> &clusters, Event event) {
    max_cluster_id++;
    std::map<size_t, size_t> cluster_to_new_id;
    for (size_t i = event.first; i < event.second; i++) {
      if (cluster_to_new_id.find(clusters[i]) == cluster_to_new_id.end()) {
        cluster_to_new_id[clusters[i]] = max_cluster_id;
        max_cluster_id++;
      }
      clusters[i] = cluster_to_new_id[clusters[i]];
    }
  }

  Real_t calculate_penalty_for_non_root_bins(
      EventTree::NodeHandle node,
      std::map<Event, std::set<size_t>> &attachment) {
    Real_t result = 0.0;
    auto node_cache_id = cache_id;
    cache_id++;

    save_clustering_in_cache(event_clusters, node_cache_id);
    update_clusters(event_clusters, node->label);

    for (auto child : node->children) {
      result += calculate_penalty_for_non_root_bins(child, attachment);
      move_cells_to_parent(child, node, attachment);
    }
    result += calculate_penalty_for_bins_at_node(node, attachment);

    /* Restore clustering of parent */
    event_clusters = get_clustering_from_cache(node_cache_id);
    return result;
  }

  Real_t calculate_penalty_for_bins_at_root() {
    Real_t counts_sum = 0.0;
    Real_t bin_count = 0.0;
    Real_t squared_counts_sum = 0.0;

    for (size_t cell = 0; cell < bin_bitmap.size(); cell++) {
      for (size_t i = 0; i < bin_bitmap[cell].size(); i++) {
        if (!bin_bitmap[cell][i]) {
          counts_sum += sum_counts[cell][i];
          squared_counts_sum += squared_counts[cell][i];
          bin_count += counts_score_length_of_bin[i];
        }
      }
    }

    if (bin_count == 0) {
      return 0.0;
    }
    return COUNTS_SCORE_CONSTANT_1 *
           calculate_l2_penalty(NEUTRAL_CN, counts_sum, squared_counts_sum,
                                bin_count, all_bins_count);
  }

  // Initialize state used for calculations
  void init_state() {
    if (clusters_cache.size() > DEFAULT_CACHE_SIZE) {
      clusters_cache.resize(DEFAULT_CACHE_SIZE);
    }
    for (auto &b : bin_bitmap) {
      std::fill(b.begin(), b.end(), false);
    }
    std::fill(event_clusters.begin(), event_clusters.end(), 0);
    max_cluster_id = 0;
    cache_id = 0;
  }

  Real_t calculate_log_score__(EventTree &tree, Attachment &at) {
    init_state();
    std::map<Event, std::set<size_t>> attachment =
        at.get_node_label_to_cells_map();
    Real_t result = 0.0;
    for (auto node : tree.get_children(tree.get_root())) {
      result += calculate_penalty_for_non_root_bins(node, attachment);
    }
    return -(result + calculate_penalty_for_bins_at_root());
  }

public:
  CountsDispersionPenalty<Real_t>(CONETInputData<Real_t> &cells)
      : sum_counts{cells.get_summed_counts()},
        squared_counts{cells.get_squared_counts()},
        counts_score_length_of_bin{cells.get_counts_scores_regions()} {

    event_clusters.resize(cells.get_loci_count());

    bin_bitmap.resize(cells.get_cells_count());
    for (auto &b : bin_bitmap) {
      b.resize(cells.get_loci_count());
    }
    for (size_t b = 0; b < cells.get_loci_count(); b++) {
      all_bins_count += counts_score_length_of_bin[b];
    }

    all_bins_count *= cells.get_cells_count();
    clusters_cache.resize(DEFAULT_CACHE_SIZE);
  }

  Real_t calculate_log_score(EventTree &tree, Attachment &at) {
    if (COUNTS_SCORE_CONSTANT_0 == 0.0 && COUNTS_SCORE_CONSTANT_1 == 0) {
      return 0.0;
    }

    return calculate_log_score__(tree, at);
  }
};
#endif // !COUNTS_SCORING_H
