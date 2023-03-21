#include "config.h"
#include "binary_search_tree.h"
#include "debug.h"

// check whether the node is a left child
bool Node::isLeftChild() const
{
    Node* parent = this->parent;

    if (parent == NULL)
    {
        return false;
    }

    if (parent->left == this)
    {
        return true;
    }

    return false;
}

// constructor
Arc::Arc(VoronoiParticle* vp, CircleEvent* circleEvent)
{
    m_vp = vp;
    m_circleEvent = circleEvent;
    left = NULL;
    right = NULL;
    parent = NULL;
}

Arc::Arc()
{
    m_vp = NULL;
    m_circleEvent = NULL;
    left = NULL;
    right = NULL;
    parent = NULL;
}

// calculate the x-coord
double Arc::xcoord() const
{
    return (m_vp->m_location.x);
}

// calculate the x-coord of the breakpoint
double Breakpoint::xcoord() const
{
    return r2_geometry::breakp_x(&(m_vpLeft->m_location), &(m_vpRight->m_location), sweepline_y);
}

// constructor
Breakpoint::Breakpoint(VoronoiParticle* vpLeft, VoronoiParticle* vpRight, HalfEdge* halfEdge)
{
    m_vpLeft = vpLeft;
    m_vpRight = vpRight;
    m_halfEdge = halfEdge;
    left = NULL;
    right = NULL;
    parent = NULL;
}

Breakpoint::Breakpoint()
{
    m_vpLeft = NULL;
    m_vpRight = NULL;
    m_halfEdge = NULL;
    left = NULL;
    right = NULL;
    parent = NULL;
}

// y-coord of the sweep line
double Breakpoint::sweepline_y = 0;

void Breakpoint::info(FILE* pFile) const
{
    fprintf(pFile, "half edge: %p\n", m_halfEdge);

    if (m_vpLeft != NULL)
    {
        fprintf(pFile, "left Voronoi particle: x=%f y=%f\n", m_vpLeft->m_location.x, m_vpLeft->m_location.y);
    }

    if (m_vpRight != NULL)
    {
        fprintf(pFile, "right Voronoi particle: x=%f y=%f\n", m_vpRight->m_location.x, m_vpRight->m_location.y);
    }
}

#ifdef SEARCH_ARC_2
// an alternative version of searchArc that may be faster
Arc* tree::searchArc(Node* root,const VoronoiParticle* vp)
{
    if(root ==NULL)
    {
        fprintf(stderr,"ERROR in searchArc: root==NULL");
        return NULL;
    }

    Breakpoint* breakpoint;
    Point top_site = Point(0,0);
    Point bottom_site = Point(0,0);

    while(!(root->isLeaf())) // as long as root is a breakpoint (inner node)
    {
        breakpoint = static_cast<Breakpoint*>(root);

        // determine top and bottom site
        if(breakpoint->m_vpLeft->m_location.y > breakpoint->m_vpRight->m_location.y)
        {
            top_site = breakpoint->m_vpLeft->m_location;
            bottom_site = breakpoint->m_vpRight->m_location;
        }
        else
        {
            top_site = breakpoint->m_vpRight->m_location;
            bottom_site = breakpoint->m_vpLeft->m_location;
        }

        // if dist to top site is smaller than to bottom site
        if(pow(vp->m_location.x - top_site.x, 2) + pow(vp->m_location.y - top_site.y,2) < pow(vp->m_location.x - bottom_site.x, 2) + pow(vp->m_location.y - bottom_site.y,2))
        {
            if(top_site.x < bottom_site.x)
            {
                root = root->left;
                continue;
            }
            else
            {
                root = root->right;
            }
        }
        else
        {
            if (root->xcoord() < vp->m_location.x) // xcoord() calls function breakp_x because root is breakpoint
            {
                root=root->right;
            }
            else
            {
                root=root->left;
            }
        }
    }

    return static_cast<Arc*>(root);
}
#else

// Search for the Arc in the tree above a given Site.
Arc* tree::searchArc(Node* root, const VoronoiParticle* vp)
{
    if (root == NULL)
    {
        fprintf(stderr, "ERROR in searchArc: root==NULL");
        return NULL;
    }

    while (!(root->isLeaf())) // as long as root is a breakpoint (inner node)
    {
        if (root->xcoord() < vp->m_location.x) // xcoord() calls function breakp_x because root is a breakpoint
        {
            root = root->right;
        }
        else
        {
            root = root->left;
        }
    }

    return static_cast<Arc*>(root);
}

#endif

// tree operations
//

// print all keys in increasing value
void tree::InorderTreeWalk(const Node* root)
{
    if (root == NULL)
    { return; }
    else
    {
        InorderTreeWalk(root->left);
        SD_PRINTF("%f\n", root->xcoord());
        InorderTreeWalk(root->right);
    }
} // warning: recursive

void tree::fPrintTreeKeys(FILE* pFile, const Node* root, double x_level, double y_level, double factor)
{
    if (root == NULL)
    {
        return;
    }
    else
    {
        if (root->isLeaf())
        {
            fprintf(pFile, "%f\t%f\t%.2f%s\n", x_level, y_level, root->xcoord(), "Arc");
        }
        else
        {
            fprintf(pFile, "%f\t%f\t%.2f%s\n", x_level, y_level, root->xcoord(), "BP");
        }

        fPrintTreeKeys(pFile, root->left, (x_level - 1 * factor), y_level - 1, factor * 0.8);
        if (root->left == NULL)
        { fprintf(pFile, "\n"); }

        if (root->right != NULL)
        {
            if (root->isLeaf())
            {
                fprintf(pFile, "%f\t%f\t%.2f%s\n", x_level, y_level, root->xcoord(), "Arc");
            }
            else
            {
                fprintf(pFile, "%f\t%f\t%.2f%s\n", x_level, y_level, root->xcoord(), "BP");
            }

            fPrintTreeKeys(pFile, root->right, (x_level + 1 * factor), y_level - 1, factor * 0.8);
        }

    }
} // warning: recursive

// insert a Node into the tree, which is represented by his root
void tree::TreeInsert(Node** root, Node* node)
{
    Node* y = NULL; // temporary parent
    Node* x = *root; // temporary node

    while (x != NULL)
    {
        y = x;

        if (node->xcoord() < x->xcoord())
        {
            x = x->left;
        }
        else
        {
            x = x->right;
        }
    }

    node->parent = y;

    if (y == NULL) // tree was empty
    {
        *root = node;
    }
    else
    {
        if (node->xcoord() < y->xcoord())
        {
            y->left = node;
        }
        else
        {
            y->right = node;
        }
    }
}

// search for the node with minmum key
Node* tree::treeMinimum(Node* root)
{
    while (root->left != NULL)
    {
        root = root->left;
    }

    return root;
}

// search for the node with maximum key
Node* tree::treeMaximum(Node* root)
{
    while (root->right != NULL)
    {
        root = root->right;
    }

    return root;
}

// identify the successor of a node in respect to his key
Node* tree::treeSuccessor(Node* node)
{
    if (node->right != NULL)
    {
        return treeMinimum(node->right);
    }

    Node* y = node->parent;

    while (y != NULL && node == y->right)
    {
        node = y;
        y = y->parent;
    }

    return y;
}

// identify the predecessor of a node in respect to his key
Node* tree::treePredecessor(Node* node)
{
    if (node->left != NULL)
    {
        return treeMaximum(node->left);
    }

    Node* y = node->parent;

    while (y != NULL && node == y->left)
    {
        node = y;
        y = y->parent;
    }

    return y;
}


// delete a node from the tree, which is represented by its root
// the node to delete must have one child or less
Node* tree::TreeDeleteOneChild(Node** root, Node* node)
{
    Node* y = node; // node to cut
    Node* x = NULL; // child of y

    if (y->left != NULL) // y has a left child
    {
        x = y->left;
    }
    else // y has a right child
    {
        x = y->right;
    }

    if (x != NULL) // y has a child
    {
        x->parent = y->parent;
    }

    if (y->parent == NULL)
    {
        *root = x;
    }
    else
    {
        if (y == y->parent->left)
        {
            y->parent->left = x;
        }
        else
        {
            y->parent->right = x;
        }
    }

    return y;
}

// find the twin breakpoint and store it in twinBP
void tree::findTwinBreakpoint(Node* root, const Breakpoint* BP, Breakpoint** twinBP)
{
    if (root == NULL)
    { return; }
    if (!root->isLeaf()) // node is a breakpoint
    {
        Breakpoint* temp = static_cast<Breakpoint*>(root);

        if (temp->m_halfEdge == BP->m_halfEdge->m_twin)
        {
            *twinBP = temp;
            return;
        }
    }

    findTwinBreakpoint(root->left, BP, twinBP);
    findTwinBreakpoint(root->right, BP, twinBP);
} // warning recursive

// delete the whole tree
void tree::freeTree(Node* root)
{
    if (root == NULL)
    { return; }
    if (root->isLeaf())
    {
        delete root;
        return;
    }
    if (root->left == NULL && root->right != NULL)
    {
        freeTree(root->right);
        delete root;
        return;
    }
    if (root->left != NULL && root->right == NULL)
    {
        freeTree(root->left);
        delete root;
        return;
    }

    freeTree(root->left);
    freeTree(root->right);
    delete root;
    return;
} // warning: recursive