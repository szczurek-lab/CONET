#include "move_type.h"

std::string move_type_to_string(MoveType type) {
	switch (type) {
	case DELETE_LEAF:
		return "Delete leaf";
	case ADD_LEAF:
		return "Add leaf";
	case PRUNE_REATTACH:
		return "Prune and reattach";
	case SWAP_LABELS:
		return "Swap labels";
	case CHANGE_LABEL:
		return "Change label";
	case SWAP_SUBTREES:
		return "Swap subtrees";
	case SWAP_ONE_BREAKPOINT:
  default:
		return "Swap one breakpoint";
	}
}

