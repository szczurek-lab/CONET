#ifndef ATTACHMENT_H
#define ATTACHMENT_H
#include <vector> 
#include <algorithm>

#include "../types.h"

class Attachment {
    private:
        std::vector<TreeLabel> cell_to_tree_label;

    public:
        Attachment(std::vector<TreeLabel> cell_to_tree_label): cell_to_tree_label{cell_to_tree_label}{ }
        bool has_attached_cells(TreeLabel label) {
            return std::find(cell_to_tree_label.begin(), cell_to_tree_label.end(), label) != cell_to_tree_label.end();
        } 

};

#endif 