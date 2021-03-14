#ifndef MOVE_DATA_H
#define MOVE_DATA_H

#include "../tree/pointer_tree.h"
/**
* Keeps data about move.
* This is a simple struct which stores data allowing roll-back of any move. 
* For instance if delet leaf move has been made, <code>oldLeftBrkp, oldRightBrkp</node> will store
* breakpoints of the deleted leaf while <code>oldParent</code>  will store handle to its parent. 
* <code>simpleSwap</code> is only used for <code>SWAP_SUBTREES</code> move and is equal to <code>false</code> 
* if and only if roots of swapped trees have been in parent - descendant relation.
*/
struct MoveData {
	using NodeHandle_ = PointerTree::NodeHandle;
	NodeHandle_ oldParent;
	NodeHandle_ node1;
	NodeHandle_ node2;
	size_t oldLeftBrkp;
	size_t oldRightBrkp;
	bool simpleSwap;
};
 


#endif // !MOVE_DATA_H
