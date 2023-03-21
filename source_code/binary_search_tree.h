#ifndef BINARY_SEARCH_TREE_H
#define BINARY_SEARCH_TREE_H

#include "config.h"
#include <cstdlib>
#include <cstdio>
#include "graph.h"
#include "r2.h"
#include "voronoi_particle.h"

// Algorithmen, eine Einfuehrung, p.245ff
// adapted to Voronoi (Computational Geometry, chapter 7)
// warning: binary search tree attribute holds for inner nodes (breakpoints) only

class Event;

class CircleEvent;

class SiteEvent;

class Node
{
public:

    Node* left;
    Node* right;
    Node* parent;

    /// is the node a left child?
    bool isLeftChild() const;

    virtual double xcoord() const = 0;

    virtual bool isLeaf() const = 0;

    virtual ~Node()
    {};
};

// beachline
class Arc : public Node
{
public:

    VoronoiParticle* m_vp;
    CircleEvent* m_circleEvent;

    // constructors
    Arc(VoronoiParticle* vp, CircleEvent* circleEvent);

    Arc();

    // destructor
    ~Arc()
    {}

    /// calculate the x-coordinate
    double xcoord() const;

    /// check whether the node is a leaf
    inline bool isLeaf() const
    {
        return true;
    }
};

class Breakpoint : public Node
{
public:

    VoronoiParticle* m_vpLeft; // represents left parabola
    VoronoiParticle* m_vpRight; // represents right parabola
    HalfEdge* m_halfEdge; // traced half edge
    static double sweepline_y; // y-coord of the sweepline

    // constructors
    Breakpoint(VoronoiParticle* vpLeft, VoronoiParticle* vpRight, HalfEdge* halfEdge);

    Breakpoint();

    // destructor
    ~Breakpoint()
    {}

    /// calculate the x-coordinate
    double xcoord() const;

    /// check whether the node is a leaf
    inline bool isLeaf() const
    {
        return false;
    }

    /// print infos
    void info(FILE* pFile) const;
};

namespace tree
{
    /**
    search for the arc in the tree above a given site
    \param root the root of the tree
    \param site the new site
    \return the Arc above the site. NULL if the function fails
    */
    Arc* searchArc(Node* root, const VoronoiParticle* vp);

    //tree operations
    //

    /// print all keys in increasing value
    void InorderTreeWalk(const Node* root); // warning: recursive

    /**
    print the tree in gnuplot readable format into a file
    standard: fPrintTreeKeys(pFile,root,0,0,1)
    */
    void fPrintTreeKeys(FILE* pFile, const Node* root, double x_level, double y_level, double factor);
    // warning: recursive

    /// insert a Node into the tree, which is represented by his root
    void TreeInsert(Node** root, Node* node);

    /// search for the node with minmum key
    Node* treeMinimum(Node* root);

    /// search for the node with maximum key
    Node* treeMaximum(Node* root);

    /// identify the successor of a node in respect to his key
    Node* treeSuccessor(Node* node);

    /// identify the predecessor of a node in respect to his key
    Node* treePredecessor(Node* node);

    /**
    delete a node from the tree, which is represented by its root.
    the node to delete must have one child or less
    */
    Node* TreeDeleteOneChild(Node** root, Node* node);

    /// find the twin breakpoint and store it in twinBP
    void findTwinBreakpoint(Node* root, const Breakpoint* BP, Breakpoint** twinBP);

    /// delete the hole tree
    void freeTree(Node* root); //warning recursive
};

#endif