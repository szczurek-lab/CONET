#ifndef MOVE_TYPE_H
#define MOVE_TYPE_H
#include <string>

enum MoveType {
	DELETE_LEAF,
	ADD_LEAF,
	PRUNE_REATTACH,
	SWAP_LABELS,
	CHANGE_LABEL,
	SWAP_SUBTREES,
	SWAP_ONE_BREAKPOINT
};

std::string moveTypeToString(MoveType type);

#endif // !MOVE_TYPE_H

