
#include "../test_utils.h"
#include "../../src/tree/event_tree.h"
#include "../../src/tree/tree_sampler.h"
#include "../../src/tree/tree_formatter.h"

void sample_tree_test() {
    std::vector<size_t> chromosome_markers{10, 20, 45};
    size_t max_locus = 44;
    VertexLabelSampler<double> sampler{max_locus, chromosome_markers};
    Random<double> random(123124);

    BEGIN_TEST;
    auto tree = sample_tree(10, sampler, random);

    IS_EQUAL(tree.get_size(), 10)

    auto nodes = tree.get_descendants(tree.get_root());
    IS_EQUAL(nodes.size(), 10);

    auto events = tree.get_all_events();

    for (auto n : nodes) {
        if (n != tree.get_root()) {
            IS_TRUE(std::find(events.begin(), events.end(), tree.get_node_event(n)) != events.end())
        }

    }
    END_TEST;
}

void basic_tree_ops_test() {
    BEGIN_TEST;
    std::string tree_str = "(0,0)-(1,10)\n(1,10)-(2,10)";
    auto tree = TreeFormatter::from_string_representation(tree_str);

    IS_EQUAL(tree.get_size(), 3);
        std::cout << tree.get_children(tree.get_root()).size();
    IS_EQUAL(tree.get_children(tree.get_root()).size(), 1);

    auto root_child = tree.get_children(tree.get_root()).front();
    std::cout << tree.get_children(root_child).size();
    IS_EQUAL(tree.get_children(root_child).size(), 1);

    auto root_grandchild = tree.get_children(root_child).front();
    tree.add_leaf(root_child, label_from_str("(2,5)"));
    IS_EQUAL(tree.get_size(), 4);
    auto parent_of_deleted = tree.delete_leaf(root_grandchild);

    IS_EQUAL(parent_of_deleted, root_child);
    IS_EQUAL(tree.get_size(), 3);

    auto rep = TreeFormatter::to_string_representation(tree);
    IS_EQUAL(rep, "(0,0)-(1,10)\n(1,10)-(2,5)\n");

    END_TEST;
}

void prunning_test() {
    BEGIN_TEST;
    std::ifstream t("./test_structure");
    std::stringstream buffer;
    buffer << t.rdbuf();
    auto tree_str = buffer.str();

    auto tree = TreeFormatter::from_string_representation(tree_str);

    std::vector<std::pair<size_t, size_t>> attachment {
        std::make_pair(87, 119),
        std::make_pair(121, 137),
        std::make_pair(110, 119),
        std::make_pair(179, 217),
        std::make_pair(45, 65),
        std::make_pair(58, 63),
        std::make_pair(26, 28),
        std::make_pair(47, 52),
        std::make_pair(148, 164),
        std::make_pair(105, 106)
    };
    Attachment at{attachment};
    tree.prune_tree(tree.get_root(), at);

    END_TEST;
}

int main(void) {
    sample_tree_test();
    basic_tree_ops_test();
    prunning_test();
}