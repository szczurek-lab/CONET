#ifndef ATTACHMENT_H
#define ATTACHMENT_H
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

#include "../types.h"

/**
 * @brief Attachment of cells to Event Tree nodes
 *
 */
class Attachment {
  std::vector<TreeLabel> cell_to_tree_label;

public:
  Attachment(TreeLabel default_label, size_t cells) {
    for (size_t i = 0; i < cells; i++) {
      cell_to_tree_label.push_back(default_label);
    }
  }

  bool has_attached_cells(TreeLabel label) const {
    return std::find(cell_to_tree_label.begin(), cell_to_tree_label.end(),
                     label) != cell_to_tree_label.end();
  }

  void set_attachment(size_t cell, TreeLabel node) {
    cell_to_tree_label[cell] = node;
  }

  std::map<TreeLabel, std::set<size_t>> get_node_label_to_cells_map() const {
    std::map<TreeLabel, std::set<size_t>> node_to_cells;
    for (size_t cell = 0; cell < cell_to_tree_label.size(); cell++) {
      if (node_to_cells.find(cell_to_tree_label[cell]) == node_to_cells.end()) {
        node_to_cells[cell_to_tree_label[cell]] = std::set<size_t>();
      }
      node_to_cells[cell_to_tree_label[cell]].insert(cell);
    }
    return node_to_cells;
  }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const Attachment &attachment) {
    for (size_t j = 0; j < attachment.cell_to_tree_label.size(); j++) {
      stream << j << ";" << attachment.cell_to_tree_label[j].first << ";"
             << attachment.cell_to_tree_label[j].second << "\n";
    }
    return stream;
  }

  friend void swap(Attachment &a, Attachment &b) {
    std::swap(a.cell_to_tree_label, b.cell_to_tree_label);
  }
};
#endif