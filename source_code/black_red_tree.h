
// implementation follows Algorithmen - Eine Einfuehrung, Cormen, Leiserson, Rivest, Stein

// black red tree rules:
// every node is either red or black
// the root is black
// every Leaf (NULL) is black
// if a node is red, his children must be black
// for every node, all paths starting at the node and ending at a leaf of the subtree have the same number of black nodes

#ifndef BLACK_RED_TREE_H
#define BLACK_RED_TREE_H

#include <cstdlib>

class Sim;

class BR_Node
{
public:

    bool m_isBlack;
    double m_x;
    double m_y;
    BR_Node* m_left;
    BR_Node* m_right;
    BR_Node* m_parent;

    // constructors
    BR_Node();

    BR_Node(bool isBlack);

    BR_Node(double y);

    BR_Node(bool isBlack, double y);

    // destructor
    virtual ~BR_Node()
    {};

    static BR_Node nil; ///< guard
    static BR_Node* p_nil; ///< pointer to the guard

    /// compare two nodes (used in insertSite)
    inline static bool compare2(BR_Node* a, BR_Node* b)
    {
        // sort by y-coordinate
        if (a->m_y < b->m_y)
        {
            return true;
        }
        if (a->m_y > b->m_y)
        {
            return false;
        }

        // same y-coordinates, sort by x-coordinate
        if (a->m_x > b->m_x)
        {
            return true;
        }

        if (a->m_x < b->m_x)
        {
            return false;
        }

        // same x- and y-coordinate (same point)
        return false;
    }

    /// compare two nodes (used in insert)
    inline static bool compare(BR_Node* a, BR_Node* b)
    {
        return (a->m_y < b->m_y);
    }

    /// left rotation
    static void leftRotate(BR_Node* x, Sim* sim);

    /// right rotation
    static void rightRotate(BR_Node* x, Sim* sim);

    /// insert a node into the tree
    static void insert(BR_Node* z, Sim* sim);

    /// insert a site into the tree (other compare function)
    static void insertSite(BR_Node* z, Sim* sim);

    /// fix the red black tree after a insertion
    static void insertFixup(BR_Node* z, Sim* sim);

    /// delete a node
    static BR_Node* remove(BR_Node* z, Sim* sim);

    /// fix the red black tree after a removal
    static void removeFixup(BR_Node* x, Sim* sim);

    /// search downwards for the node with minmum key, starting at node
    static BR_Node* treeMinimum(BR_Node* node);

    /// search downwards for the node with maximum key, starting at node
    static BR_Node* treeMaximum(BR_Node* node);

    /// identify the successor of a node in respect to his key
    static BR_Node* treeSuccessor(BR_Node* node);

    /// identify the predecessor of a node in respect to his key
    static BR_Node* treePredecessor(BR_Node* node);

    /// print the Tree in gnuplot readable format, call with 0,0,1
    static void printTreeKeys(BR_Node* node, double x_level, double y_level, double factor);

    /// delete the hole tree
    static void freeTree(BR_Node* node);
};

#endif