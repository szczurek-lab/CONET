#ifndef TREE_FORMATTER_H
#define TREE_FORMATTER_H

#include <sstream>
#include <string>
#include <map>

#include "event_tree.h"

class TreeFormatter {
private:
	using NodeHandle = EventTree::NodeHandle;
	// C++ 2017 >= 
	inline static const std::string ROOT_STRING_REP = "(0,0)";

    static std::string get_node_label(EventTree &tree, NodeHandle node)
	{
		return node == tree.get_root() ? ROOT_STRING_REP : label_to_str(node->label);
	}

	static void toString(EventTree &tree, NodeHandle node, std::stringstream &ss) {
		for (auto &child : node->children) {
			ss << get_node_label(tree, node) << "-" << get_node_label(tree, child) << "\n";
			toString(tree, child, ss);
		}
	}

	static std::vector<std::string> split_into_lines(std::string &s) {
		std::stringstream ss(s);
		std::string buffer; 
  		std::vector<std::string> lines; 

		while(std::getline(ss, buffer, '\n')){
			lines.push_back(buffer);
		}
		return lines;
	}

public: 
	static std::string to_string_representation(EventTree &tree) {
		std::stringstream ss;
		toString(tree, tree.get_root(), ss);
		return ss.str();
	}

	static EventTree from_string_representation(std::string &rep) {
		auto lines = split_into_lines(rep);
		EventTree tree;
		std::map<std::string, NodeHandle> label_to_node; 
		label_to_node[ROOT_STRING_REP] = tree.get_root();

		for (auto line : lines) {
			auto parent_str = line.substr(0, line.find('-'));
			auto child_str = line.substr(line.find('-') + 1, line.length());
			label_to_node[child_str] = tree.add_leaf(label_to_node[parent_str], label_from_str(child_str));

		}
		return tree;
	}
};

#endif 