#include "config.h"
#include "sim_graph.h"
#include "debug.h"
#include "my_random.h"

// delete a half edge
void sim_graph::delete_edge(HalfEdge* half_edge)
{
    HalfEdge* twin = half_edge->m_twin;
    Vertex* vertex_kept = twin->m_origin;

    // calculate new vertex position
    vertex_kept->m_position = (vertex_kept->m_position + half_edge->m_origin->m_position) * 0.5;

    // set m_outer components
    half_edge->m_vp->m_outerComponent = half_edge->m_next;
    twin->m_vp->m_outerComponent = twin->m_next;

    // set origins
    twin->m_next->m_origin = vertex_kept;
    half_edge->m_prev->m_twin->m_origin = vertex_kept;

    // set half edges
    half_edge->m_prev->m_next = half_edge->m_next;
    half_edge->m_next->m_prev = half_edge->m_prev;

    twin->m_prev->m_next = twin->m_next;
    twin->m_next->m_prev = twin->m_prev;

    // deactivate
    half_edge->m_inUse = false;
    twin->m_inUse = false;
    half_edge->m_origin->m_inUse = false;
}

// paint the cells with random colors
void sim_graph::randomColors(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_color = myRandom::randomDouble(0, 10);
    }
}

// check whether all real cells are closed
void sim_graph::assertRealCellsClosed(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        if (sim->m_particles[i].m_outerComponent == NULL)
        {
            fprintf(stderr, "ERROR in sim_graph::assertRealCellsClosed: not every particle has an outer component\n");
            exit(1);
        }

        if (!graph::isClosed(sim->m_particles[i].m_outerComponent))
        {
            fprintf(stderr,
                    "ERROR in sim_graph::assertRealCellsClosed: face not closed around the Voronoi particle: x=%f, y=%f\n",
                    sim->m_particles[i].m_location.x, sim->m_particles[i].m_location.y);
            exit(1);
        }
    }
}

// print the sites
void sim_graph::printSites(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        fprintf(pFile, "%.20f %.20f\n", sim->m_particles[i].m_location.x, sim->m_particles[i].m_location.y);
    }

    fclose(pFile);
}

// print the positions of the ghost particles
void sim_graph::printGhosts(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    FILE* pFile2;
    pFile2 = fopen("graph/ghosts_fixed", "w");
    printf("warning: fixed ghosts are written to file ghosts_fixed\n");

    for (unsigned int i = 0; i < sim->m_ghosts.counter(); i++)
    {
        if (sim->m_ghosts[i].m_flag_1 <= 5 && sim->m_ghosts[i].m_flag_1 >= 1)
        {

            fprintf(pFile2, "%.20f %.20f\n", sim->m_ghosts[i].m_location.x, sim->m_ghosts[i].m_location.y);
        }
        else
        {
            fprintf(pFile, "%.20f %.20f\n", sim->m_ghosts[i].m_location.x, sim->m_ghosts[i].m_location.y);
        }
    }

    fclose(pFile);
    fclose(pFile2);
}

// print the graph in gnuplot readable format
void sim_graph::printFaces(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {

        if (sim->m_particles[i].m_outerComponent == NULL)
        {
            fprintf(stderr, "ERROR in sim_graph::printFaces: not every particle has an outer component\n");
            continue;
        }

        graph::printFace(pFile, sim->m_particles[i].m_outerComponent);
    }

    fclose(pFile);
}

// print the ghost faces of the graph
void sim_graph::printGhostFaces(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int j = 0; j < sim->m_ghosts.counter(); j++)
    {
        if (sim->m_ghosts[j].m_outerComponent == NULL)
        {
            fprintf(stderr, "ERROR in sim_graph::printFaces: not every particle has an outer component\n");
            continue;
        }

        if (graph::isClosed(sim->m_ghosts[j].m_outerComponent))
        {
            graph::printFace(pFile, sim->m_ghosts[j].m_outerComponent);
        }
    }

    fclose(pFile);
}

// print the completed faces with half edges in gnuplot readable format
void sim_graph::printHalfFaces(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        if (sim->m_particles[i].m_outerComponent == NULL)
        {
            fprintf(stderr, "ERROR in sim_graph::printHalfFaces: not every particle has an outer component\n");
            continue;
        }

        graph::printHalfFace(pFile, sim->m_particles[i].m_outerComponent);
    }

    fclose(pFile);
}

// print data of the half edges
void sim_graph::printHalfEdgesData(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_halfEdges.counter(); i++)
    {
        if (!sim->m_halfEdges[i].m_inUse) continue; // skip half edges not in use

        sim->m_halfEdges[i].info(pFile);
        fprintf(pFile, "\n");
    }

    fclose(pFile);
}

// print the vertices
void sim_graph::printVertices(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_vertices.counter(); i++)
    {
        if (!sim->m_vertices[i].m_inUse) continue; // skip vertices not in use

        fprintf(pFile, "%f\t%f\n", sim->m_vertices[i].m_position.x, sim->m_vertices[i].m_position.y);
    }

    fclose(pFile);
}

// print the breakpoint vertices
void sim_graph::printBreakPointVertices(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_breakpointVertices.size(); i++)
    {
        fprintf(pFile, "%f\t%f\n", sim->m_breakpointVertices[i]->m_position.x,
                sim->m_breakpointVertices[i]->m_position.y);
    }

    fclose(pFile);
}

// print the edges adjacent to the breakpoint vertices
void sim_graph::printBreakPointEdges(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_breakpointVertices.size(); i++)
    {
        fprintf(pFile, "%f\t%f\n",
                sim->m_breakpointVertices[i]->m_position.x,
                sim->m_breakpointVertices[i]->m_position.y);
        fprintf(pFile, "%f\t%f\n",
                sim->m_breakpointVertices[i]->m_edge->m_next->m_origin->m_position.x,
                sim->m_breakpointVertices[i]->m_edge->m_next->m_origin->m_position.y);
        fprintf(pFile, "\n");
        fprintf(pFile, "\n");
    }

    fclose(pFile);
}

static void printFaceColor(FILE* pFile, HalfEdge* halfEdge)
{
    double color = myRandom::randomDouble(0, 10);

    HalfEdge* left = halfEdge;
    HalfEdge* right = halfEdge->m_prev;

    do
    {
        fprintf(pFile, "%f %f %f\n",
                left->m_origin->m_position.x,
                left->m_origin->m_position.y, color);
        fprintf(pFile, "%f %f %f\n",
                left->m_next->m_origin->m_position.x,
                left->m_next->m_origin->m_position.y, color);
        fprintf(pFile, "\n");
        fprintf(pFile, "%f %f %f\n",
                right->m_next->m_origin->m_position.x,
                right->m_next->m_origin->m_position.y, color);
        fprintf(pFile, "%f %f %f\n",
                right->m_origin->m_position.x,
                right->m_origin->m_position.y, color);
        fprintf(pFile, "\n");
        fprintf(pFile, "\n");

        left = left->m_next;
        right = right->m_prev;

    } while (left != right && left->m_prev != right);
}

// print the face centers
void sim_graph::printFaceCenters(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    HalfEdge* current = NULL;
    HalfEdge* start = NULL;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        start = sim->m_particles[i].m_outerComponent;
        current = start;

        do
        {
            fprintf(pFile, "%f %f\n", current->m_center.x, current->m_center.y);
            current = current->m_next;

        } while (current != start);
    }

    fclose(pFile);
}

// print the graph for gnuplot pm3d map splot
void sim_graph::printGraphColor(Sim* sim, const char* filename)
{
    FILE* pFile;
    pFile = fopen(filename, "w");

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        if (sim->m_particles[i].m_outerComponent == NULL)
        {
            fprintf(stderr, "ERROR in sim_graph::printGraphColor: not every particle has an outer component\n");
            continue;
        }

        printFaceColor(pFile, sim->m_particles[i].m_outerComponent);
    }

    fclose(pFile);
}

// print the whole graph
void sim_graph::print_graph_all(Sim* sim)
{
    // print the graph
    sim_graph::printSites(sim, "graph/sites");
    sim_graph::printGhosts(sim, "graph/ghosts");
    sim_graph::printVertices(sim, "graph/vertices");
    sim_graph::printFaces(sim, "graph/faces");
    sim_graph::printGhostFaces(sim, "graph/ghost_faces");
    sim_graph::printHalfFaces(sim, "graph/half_faces");
    sim_graph::printFaceCenters(sim, "graph/face_centers");
    sim_graph::printGraphColor(sim, "graph/graph_color");
    sim_graph::printBreakPointVertices(sim, "graph/breakpoint_vertices");
    sim_graph::printBreakPointEdges(sim, "graph/breakpoint_edges");

    // print the cell centers
    sim->printCellCentersOfVolume("graph/centers");
    sim->printGhostCellCentersOfVolume("graph/ghost_centers");
}