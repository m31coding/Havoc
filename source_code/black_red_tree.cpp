#include "config.h"
#include "black_red_tree.h"
#include <cstdio>
#include <assert.h>
#include "debug.h"
#include "sim.h"

// define static member variables
BR_Node BR_Node::nil = BR_Node(true);
BR_Node* BR_Node::p_nil = &(BR_Node::nil);

// constructors
//

BR_Node::BR_Node()
        :
        m_isBlack(false),
        m_x(0),
        m_y(0),
        m_left(p_nil),
        m_right(p_nil),
        m_parent(p_nil)
{

}

BR_Node::BR_Node(bool isBlack)
        :
        m_isBlack(isBlack),
        m_x(0),
        m_y(0),
        m_left(p_nil),
        m_right(p_nil),
        m_parent(p_nil)
{

}

BR_Node::BR_Node(double y)
        :
        m_isBlack(false),
        m_x(0),
        m_y(y),
        m_left(p_nil),
        m_right(p_nil),
        m_parent(p_nil)
{

}

BR_Node::BR_Node(bool isBlack, double y)
        :
        m_isBlack(isBlack),
        m_x(0),
        m_y(y),
        m_left(p_nil),
        m_right(p_nil),
        m_parent(p_nil)
{

}

// left rotation
void BR_Node::leftRotate(BR_Node* x, Sim* sim)
{
    RB_PRINTF("\t\tleft rotation: %u\n", (int) x->m_key_1);

    assert(x->m_right != p_nil);

    BR_Node* y = x->m_right;
    x->m_right = y->m_left;

    if (!(y->m_left == p_nil))
    {
        y->m_left->m_parent = x;
    }

    y->m_parent = x->m_parent;

    if (x->m_parent == p_nil)
    {
        sim->m_BR_root = y;
    }
    else
    {
        if (x == x->m_parent->m_left)
        {
            x->m_parent->m_left = y;
        }
        else
        {
            x->m_parent->m_right = y;
        }
    }

    y->m_left = x;
    x->m_parent = y;
}

// right rotation
void BR_Node::rightRotate(BR_Node* x, Sim* sim)
{
    RB_PRINTF("\t\tright rotation: %u\n", (int) x->m_key_1);

    assert(x->m_left != p_nil);

    BR_Node* y = x->m_left;
    x->m_left = y->m_right;

    if (!(y->m_right == p_nil))
    {
        y->m_right->m_parent = x;
    }

    y->m_parent = x->m_parent;

    if (x->m_parent == p_nil)
    {
        sim->m_BR_root = y;
    }
    else
    {
        if (x == x->m_parent->m_right)
        {
            x->m_parent->m_right = y;
        }
        else
        {
            x->m_parent->m_left = y;
        }
    }

    y->m_right = x;
    x->m_parent = y;
}


// a node into the tree
void BR_Node::insert(BR_Node* z, Sim* sim)
{
    RB_PRINTF("insert node: %u\n", (int) z->m_key_1);

    BR_Node* y = p_nil;
    BR_Node* x = sim->m_BR_root;

    while (x != p_nil)
    {
        y = x;

        if (compare(z, x))
        {
            x = x->m_left;
        }
        else
        {
            x = x->m_right;
        }
    }

    z->m_parent = y;

    if (y == p_nil)
    {
        sim->m_BR_root = z;
    }

    else if (compare(z, y))
    {
        y->m_left = z;
    }

    else
    {
        y->m_right = z;
    }

    z->m_left = p_nil;
    z->m_right = p_nil;
    z->m_isBlack = false;

    insertFixup(z, sim);
}

// insert a site into the tree (other compare function)
void BR_Node::insertSite(BR_Node* z, Sim* sim)
{
    RB_PRINTF("insert node: %u\n", (int) z->m_key_1);

    BR_Node* y = p_nil;
    BR_Node* x = sim->m_BR_root;

    while (x != p_nil)
    {
        y = x;

        if (compare2(z, x))
        {
            x = x->m_left;
        }
        else
        {
            x = x->m_right;
        }
    }

    z->m_parent = y;

    if (y == p_nil)
    {
        sim->m_BR_root = z;
    }

    else if (compare2(z, y))
    {
        y->m_left = z;
    }

    else
    {
        y->m_right = z;
    }

    z->m_left = p_nil;
    z->m_right = p_nil;
    z->m_isBlack = false;

    insertFixup(z, sim);
}

// fix the red black tree after an insertion
void BR_Node::insertFixup(BR_Node* z, Sim* sim)
{
    RB_PRINTF("\tfixup node: %u\n", (int) z->m_key_1);

    while (!z->m_parent->m_isBlack)
    {
        if (z->m_parent == z->m_parent->m_parent->m_left)
        {
            BR_Node* y = z->m_parent->m_parent->m_right;

            if (!y->m_isBlack)
            {
                z->m_parent->m_isBlack = true;
                y->m_isBlack = true;
                z->m_parent->m_parent->m_isBlack = false;
                z = z->m_parent->m_parent;
            }
            else
            {
                if (z == z->m_parent->m_right)
                {
                    z = z->m_parent;
                    leftRotate(z, sim);
                }

                z->m_parent->m_isBlack = true;
                z->m_parent->m_parent->m_isBlack = false;
                rightRotate(z->m_parent->m_parent, sim);
            }

        }

        else
        {
            BR_Node* y = z->m_parent->m_parent->m_left;

            if (!y->m_isBlack)
            {
                z->m_parent->m_isBlack = true;
                y->m_isBlack = true;
                z->m_parent->m_parent->m_isBlack = false;
                z = z->m_parent->m_parent;
            }
            else
            {
                if (z == z->m_parent->m_left)
                {
                    z = z->m_parent;
                    rightRotate(z, sim);

                }

                z->m_parent->m_isBlack = true;
                z->m_parent->m_parent->m_isBlack = false;
                leftRotate(z->m_parent->m_parent, sim);
            }
        }
    }

    sim->m_BR_root->m_isBlack = true;
}

// delete a node
BR_Node* BR_Node::remove(BR_Node* z, Sim* sim)
{
    BR_Node* y;
    BR_Node* x;

    if (z->m_left == p_nil || z->m_right == p_nil)
    {
        y = z;
    }
    else
    {
        y = treeSuccessor(z);
    }

    if (y->m_left != p_nil)
    {
        x = y->m_left;
    }
    else
    {
        x = y->m_right;
    }

    x->m_parent = y->m_parent;

    if (y->m_parent == p_nil)
    {
        sim->m_BR_root = x;
    }
    else
    {
        if (y == y->m_parent->m_left)
        {
            y->m_parent->m_left = x;
        }
        else
        {
            y->m_parent->m_right = x;
        }
    }

    bool isYblack = y->m_isBlack;

    if (y != z)
    {
        // remove z, insert y
        y->m_left = z->m_left;
        y->m_right = z->m_right;
        y->m_parent = z->m_parent;
        y->m_isBlack = z->m_isBlack;

        z->m_left->m_parent = y;
        z->m_right->m_parent = y;

        if (z == z->m_parent->m_left) // z is left child
        {
            z->m_parent->m_left = y;
        }

        if (z == z->m_parent->m_right) // z is right child
        {
            z->m_parent->m_right = y;
        }

        if (sim->m_BR_root == z)
        {
            sim->m_BR_root = y;
        }
    }

    if (isYblack)
    {
        removeFixup(x, sim);
    }

    return z; // the node to delete (free memory)
}

// fix the red black tree after a removal
void BR_Node::removeFixup(BR_Node* x, Sim* sim)
{
    while (x != sim->m_BR_root && x->m_isBlack)
    {
        if (x == x->m_parent->m_left)
        {
            BR_Node* w = x->m_parent->m_right;

            if (!w->m_isBlack)
            {
                w->m_isBlack = true;
                x->m_parent->m_isBlack = false;
                leftRotate(x->m_parent, sim);
                w = x->m_parent->m_right;
            }

            if (w->m_left->m_isBlack && w->m_right->m_isBlack)
            {
                w->m_isBlack = false;
                x = x->m_parent;
            }

            else
            {
                if (w->m_right->m_isBlack)
                {
                    w->m_left->m_isBlack = true;
                    w->m_isBlack = false;
                    rightRotate(w, sim);
                    w = x->m_parent->m_right;
                }

                w->m_isBlack = x->m_parent->m_isBlack;
                x->m_parent->m_isBlack = true;
                w->m_right->m_isBlack = true;
                leftRotate(x->m_parent, sim);
                x = sim->m_BR_root;
            }
        }

        else
        {
            BR_Node* w = x->m_parent->m_left;

            if (!w->m_isBlack)
            {
                w->m_isBlack = true;
                x->m_parent->m_isBlack = false;
                rightRotate(x->m_parent, sim);
                w = x->m_parent->m_left;
            }

            if (w->m_right->m_isBlack && w->m_left->m_isBlack)
            {
                w->m_isBlack = false;
                x = x->m_parent;
            }

            else
            {
                if (w->m_left->m_isBlack)
                {
                    w->m_right->m_isBlack = true;
                    w->m_isBlack = false;
                    leftRotate(w, sim);
                    w = x->m_parent->m_left;
                }

                w->m_isBlack = x->m_parent->m_isBlack;
                x->m_parent->m_isBlack = true;
                w->m_left->m_isBlack = true;
                rightRotate(x->m_parent, sim);
                x = sim->m_BR_root;
            }
        }
    }

    x->m_isBlack = true;
}

// search for the node with minimum key
BR_Node* BR_Node::treeMinimum(BR_Node* node)
{
    while (node->m_left != p_nil)
    {
        node = node->m_left;
    }

    return node;
}

// search for the node with maximum key
BR_Node* BR_Node::treeMaximum(BR_Node* node)
{
    while (node->m_right != p_nil)
    {
        node = node->m_right;
    }

    return node;
}

// identify the successor of a node with respect to his key
BR_Node* BR_Node::treeSuccessor(BR_Node* node)
{
    if (node->m_right != p_nil)
    {
        return treeMinimum(node->m_right);
    }

    BR_Node* y = node->m_parent;

    while (y != p_nil && node == y->m_right)
    {
        node = y;
        y = y->m_parent;
    }

    return y;
}

// identify the predecessor of a node with respect to his key
BR_Node* BR_Node::treePredecessor(BR_Node* node)
{
    if (node->m_left != p_nil)
    {
        return treeMaximum(node->m_left);
    }

    BR_Node* y = node->m_parent;

    while (y != p_nil && node == y->m_left)
    {
        node = y;
        y = y->m_parent;
    }

    return y;
}

// print the char corresponding to the color of the node
static char charColor(BR_Node* node)
{
    if (node->m_isBlack)
    {
        return 'B';
    }

    else
    {
        return 'R';
    }
}

// print the tree in gnuplot readable format, call with 0,0,1
void BR_Node::printTreeKeys(BR_Node* node, double x_level, double y_level, double factor)
{
    if (node == p_nil)
    {
        return;
    }
    else
    {
        printf("%f\t%f\t%d%c\n", x_level, y_level, (int) node->m_y, charColor(node));
        printTreeKeys(node->m_left, (x_level - 1 * factor), y_level - 1, factor * 0.8);
        if (node->m_left == p_nil)
        { printf("\n"); }

        if (node->m_right != p_nil)
        {
            printf("%f\t%f\t%d%c\n", x_level, y_level, (int) node->m_y, charColor(node));
            printTreeKeys(node->m_right, (x_level + 1 * factor), y_level - 1, factor * 0.8);
        }
    }
} // warning: recursive

// delete the whole tree
void BR_Node::freeTree(BR_Node* node)
{
    if (node == p_nil)
    { return; }

    if (node->m_left == p_nil && node->m_right != p_nil)
    {
        BR_Node::freeTree(node->m_right);
        delete node;
        return;
    }

    if (node->m_left != p_nil && node->m_right == p_nil)
    {
        BR_Node::freeTree(node->m_left);
        delete node;
        return;
    }

    BR_Node::freeTree(node->m_left);
    BR_Node::freeTree(node->m_right);
    delete node;
    return;
}