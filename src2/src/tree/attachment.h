#ifndef ATTACHMENT_H
#define ATTACHMENT_H
#include <vector> 
#include <algorithm>
#include <map>

#include "../types.h"

class Attachment {
    public:
        std::vector<TreeLabel> cell_to_tree_label;

    public:
        Attachment(std::vector<TreeLabel> cell_to_tree_label): cell_to_tree_label{cell_to_tree_label}{ }
        bool has_attached_cells(TreeLabel label) {
            return std::find(cell_to_tree_label.begin(), cell_to_tree_label.end(), label) != cell_to_tree_label.end();
        } 


        void set_attachment(size_t cell, TreeLabel node) {
            cell_to_tree_label[cell] = node;
        }

        std::map<Event, std::set<size_t>> get_node_label_to_cells_map() {
            std::map<Event, std::set<size_t>> node_to_cells;
            for (size_t cell = 0; cell < cell_to_tree_label.size(); cell++) {
			    if (node_to_cells.find(cell_to_tree_label[cell]) == node_to_cells.end()) {
				    node_to_cells[cell_to_tree_label[cell]] = std::set<size_t>();
			    }
			    node_to_cells[cell_to_tree_label[cell]].insert(cell);
            }
            return node_to_cells;
        }
};

#endif 