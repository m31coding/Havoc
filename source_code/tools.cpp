#include "config.h"
#include "tools.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include "options.h"
#include "r2.h"
#include "sim.h"
#include "binary_input_output.h"
#include "sim_graph.h"
#include "sweepline.h"
#include "point_location.h"
#include <limits>
#include "sim_core.h"

// static helpers
//

static double calculate_turbulence_kinetic_energy(Sim* sim);

static Vector row_mean_velocity(double height, double dx, Sim* sim);

static double tke_row(double height, double dx, Sim* sim);

static double mass_standard_deviation_square(Sim* sim);

// use a havoc tool
bool tools::use_a_tool()
{
    // use havoc to downsample binary output files for gnuplot matrix plot
    if (options::downsampling_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::downsampling_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);
        sim_core::calculateConservedVariables(sim);

        sim->m_box->determineGhostParticles();
        sim->freeGraph();
        sim->m_box->createGhostParticles();
        sim->m_box->reset();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_core::updateCellVolumesAndCenters(sim);
        sim_core::calculatePrimitiveVariables(sim);
        sim->m_box->copyPrimitiveVariablesCellCentersAndVolume();
        sim_core::calculate_gradients(sim);
        (Constants::get()->function_slopeLimitGradients())(sim);

        tools::downsampling(options::pixel_per_length, options::downsampling_filename, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use variable to calculate variables like energy_y or vorticity of the hole box
    if (options::add_up_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::add_up_filename);
        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);
        sim_core::calculateConservedVariables(sim);

        if (options::mass_std_dev_square)
        {
            printf("%f\t%e\n", sim_core::simulation_time, mass_standard_deviation_square(sim));
        }
        else if (options::abs_vorticity)
        {
            sim->m_box->determineGhostParticles();
            sim->freeGraph();
            sim->m_box->createGhostParticles();
            sim->m_box->reset();
            sim->eventQueueIni();
            sim->eventQueueIniGhosts();
            sweepline::makeTesselation(sim);
            (Constants::get()->function_delete_zero_length_edges())(sim);
            sim_core::updateHalfEdgeLengthsAndCenters(sim);
            sim_core::updateCellVolumesAndCenters(sim);
            sim_core::calculatePrimitiveVariables(sim);
            sim->m_box->copyPrimitiveVariablesCellCentersAndVolume();
            sim_core::calculate_gradients(sim);
            (Constants::get()->function_slopeLimitGradients())(sim);

            double vorticity_sum = 0;

            for (unsigned int l = 0; l < sim->m_particles.counter(); l++)
            {
                vorticity_sum += fabs(sim->m_particles[l].vorticity());

                if (options::color_on_the_fly)
                {
                    sim->m_particles[l].m_color = fabs(sim->m_particles[l].vorticity());
                }
            }

            printf("%f\t%f\n", sim_core::simulation_time, vorticity_sum);
        }
        else if (options::energy)
        {
            double energy_sum = 0;

            for (unsigned int j = 0; j < sim->m_particles.counter(); j++)
            {
                energy_sum += sim->m_particles[j].m_energy;

                if (options::color_on_the_fly)
                {
                    sim->m_particles[j].m_color = sim->m_particles[j].m_energy;
                }
            }

            printf("%f	%f\n", sim_core::simulation_time, energy_sum);
        }
        else if (options::kinetic_energy)
        {
            double kinetic_energy_sum = 0;

            for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
            {
                kinetic_energy_sum += 0.5 * sim->m_particles[i].m_mass * sim->m_particles[i].m_velocity.length_square();

                if (options::color_on_the_fly)
                {
                    sim->m_particles[i].m_color =
                            0.5 * sim->m_particles[i].m_mass * sim->m_particles[i].m_velocity.length_square();
                }
            }

            printf("%f	%e\n", sim_core::simulation_time, kinetic_energy_sum);
        }
        else if (options::kinetic_energy_y)
        {
            double kinetic_energy_y_sum = 0;

            for (unsigned int a = 0; a < sim->m_particles.counter(); a++)
            {
                kinetic_energy_y_sum += 0.5 * sim->m_particles[a].m_mass * pow(sim->m_particles[a].m_velocity.y, 2);

                if (options::color_on_the_fly)
                {
                    sim->m_particles[a].m_color =
                            0.5 * sim->m_particles[a].m_mass * pow(sim->m_particles[a].m_velocity.y, 2);
                }
            }

            printf("%f	%f\n", sim_core::simulation_time, kinetic_energy_y_sum);
        }
        else if (options::turbulence_kinetic_energy)
        {
            printf("%f	%f\n", sim_core::simulation_time, calculate_turbulence_kinetic_energy(sim));
        }

        if (options::color_on_the_fly)
        {
            printf("coloring on the fly and save to file 99999.hc\n");
            sim_core::counter--;
            bin_io::binaryOutput("99999.hc", 0, sim);
        }

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;

        return true;
    }

    if (options::whatever)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::whatever_filename);
        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        // implement whatever here
        //

        FILE* pFile;
        FILE* pFile2;

        pFile = fopen("particles", "w");
        pFile2 = fopen("obstacle_particles", "w");

        for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
        {
            if (sim->m_particles[i].m_obstacle_particle)
            {
                fprintf(pFile, "%f\t%f\n", sim->m_particles[i].m_location.x, sim->m_particles[i].m_location.y);
            }
            else
            {
                fprintf(pFile2, "%f\t%f\n", sim->m_particles[i].m_location.x, sim->m_particles[i].m_location.y);
            }
        }

        fclose(pFile);

        //
        // whatever end

        sim->freeGraph();
        delete sim;

        options::clear();
        Constants::del();

        return true;
    }

    // use havoc to find the minimum/maximum of a primitive variable in binary data
    if (options::minmax_filename != NULL)
    {
        Sim* sim = NULL;

        if (!options::subsequent)
        {
            sim = bin_io::binaryInput(options::minmax_filename);
            tools::find_minimum_maximum(sim);
            delete sim;
        }
        else
        {

            int filenumber = bin_io::filename_to_int(options::minmax_filename);
            char filename[9];
            bin_io::int_to_filename(filename, filenumber);

            while (bin_io::file_exists(filename))
            {
                printf("%s\n", filename);
                sim = bin_io::binaryInput(filename);
                tools::find_minimum_maximum(sim);
                delete sim;
                filenumber++;
                bin_io::int_to_filename(filename, filenumber);
            }
        }

        tools::print_minimum_maximum(stdout);
        options::clear();
        return true;
    }

    if (options::zero_grad)
    {
        if (options::convert_filename == NULL)
        {
            fprintf(stderr, "ERROR: zero_grad option must be used with convert\n");
            exit(1);
        }

        Constants::readConstants();
        Constants::check_and_adapt();

        double delta_x = Constants::getBOX_XMAX() - Constants::getBOX_XMIN();
        double delta_y = Constants::getBOX_YMAX() - Constants::getBOX_YMIN();

        options::cut = true;

        options::x_min = Constants::getBOX_XMIN() + options::cut_percent * delta_x;
        options::x_max = Constants::getBOX_XMAX() - options::cut_percent * delta_x;
        options::y_min = Constants::getBOX_YMIN() + options::cut_percent * delta_y;
        options::y_max = Constants::getBOX_YMAX() - options::cut_percent * delta_y;

        printf("\tx: %f - %f\n", options::x_min, options::x_max);
        printf("\ty: %f - %f\n", options::y_min, options::y_max);

        bin_io::convert(options::convert_filename);

        options::clear();
        Constants::del();

        return true;
    }

    // use havoc to convert a binary file
    if (options::convert_filename != NULL)
    {
        bin_io::convert(options::convert_filename);
        options::clear();
        return true;
    }

    // use havoc to do a cross section
    if (options::section_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();
        tools::cross_section(options::section_filename);
        options::clear();
        Constants::del();
        return true;
    }

    // use havoc to calculate gradients
    if (options::gradient_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();
        tools::calculate_gradients(options::gradient_filename);
        options::clear();
        Constants::del();
        return true;
    }

    // use havoc to print the graph
    if (options::graph_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::graph_filename);
        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        sim_graph::print_graph_all(sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to generate ghost data
    if (options::convert_ghost_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::convert_ghost_filename);
        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);

        bin_io::ghost_convert(options::convert_ghost_filename, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to color a cell
    if (options::color_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::color_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        tools::color(options::p.x, options::p.y, options::color, sim);

        bin_io::binaryOutput(options::color_filename, 0, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to color the particles according to their volumes
    if (options::color_volume_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::color_volume_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        tools::color_volume(sim);

        bin_io::binaryOutput(options::color_volume_filename, 0, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to color the particles according to their roundness
    if (options::color_roundness_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::color_roundness_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        tools::color_roundness(sim);

        bin_io::binaryOutput(options::color_roundness_filename, 0, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to color the particles according to their number of neighbours
    if (options::color_nof_neighbours_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::color_nof_neighbours_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);

        tools::color_nof_neighbours(sim);

        bin_io::binaryOutput(options::color_nof_neighbours_filename, 0, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to color the particles according to their vorticity
    if (options::color_vorticity_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();

        Sim* sim = bin_io::binaryInput(options::color_vorticity_filename);

        sim_core::counter--;

        sim->create_box();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        sim_graph::assertRealCellsClosed(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_graph::assertRealCellsClosed(sim);
        sim_core::updateCellVolumesAndCenters(sim);
        sim_core::calculateConservedVariables(sim);

        sim->m_box->determineGhostParticles();
        sim->freeGraph();
        sim->m_box->createGhostParticles();
        sim->m_box->reset();
        sim->eventQueueIni();
        sim->eventQueueIniGhosts();
        sweepline::makeTesselation(sim);
        (Constants::get()->function_delete_zero_length_edges())(sim);
        sim_core::updateHalfEdgeLengthsAndCenters(sim);
        sim_core::updateCellVolumesAndCenters(sim);
        sim_core::calculatePrimitiveVariables(sim);
        sim->m_box->copyPrimitiveVariablesCellCentersAndVolume();
        sim_core::calculate_gradients(sim);
        (Constants::get()->function_slopeLimitGradients())(sim);

        tools::color_vorticity(sim);

        bin_io::binaryOutput(options::color_vorticity_filename, 0, sim);

        options::clear();
        Constants::del();
        sim->freeGraph();
        delete sim;
        return true;
    }

    // use havoc to find an indice (debug)
    if (options::indice_filename != NULL)
    {
        Constants::readConstants();
        Constants::check_and_adapt();
        Sim* sim = bin_io::binaryInput(options::indice_filename);

        if (options::object == 1)
        {
            sim_core::findParticleIndice(sim, options::p.x, options::p.y);
        }
        else if (options::object == 2)
        {
            sim_core::findGhostIndice(sim, options::p.x, options::p.y);
        }
        else if (options::object == 3)
        {
            fprintf(stderr, "doesnt work yet\n");
            exit(1);
        }

        options::clear();
        Constants::del();
        delete sim;
        return true;
    }

    return false;
}

// calculate the turbulence kinetic energy of a KH simulation
static double calculate_turbulence_kinetic_energy(Sim* sim)
{
    // delete the file
    FILE* pFile = NULL;
    pFile = fopen("tke_matrix", "w");
    fclose(pFile);

    double tke = 0; // turbulence kinetic energy

    double dx = (Constants::getBOX_XMAX() - Constants::getBOX_XMIN()) / (2 * Constants::getNOF_PARTICLES_X());
    double dy = (Constants::getBOX_YMAX() - Constants::getBOX_YMIN()) / (2 * Constants::getNOF_PARTICLES_Y());

    Point first = Point(Constants::getBOX_XMIN() + 5 * dx, Constants::getBOX_YMIN() + 0.5 * dy);
    Point current = first;

    double counter_y = 0;

    // calculate mean velocity of the row
    while (true)
    {
        current.y = first.y + counter_y * dy;

        if (!(current.y < Constants::getBOX_YMAX()))
        {
            break;
        }

        tke += tke_row(current.y, dx, sim);
        counter_y++;
    }

    return tke;
}

// calculate the mass standard deviation square
static double mass_standard_deviation_square(Sim* sim)
{
    double mass_sum = 0;
    double mass_sum_of_squares = 0;
    unsigned int nof_particles = sim->m_particles.counter();

    for (unsigned int i = 0; i < nof_particles; i++)
    {
        mass_sum += sim->m_particles[i].m_mass;
        mass_sum_of_squares += pow(sim->m_particles[i].m_mass, 2);
    }

    return (mass_sum_of_squares / nof_particles - pow(mass_sum / nof_particles, 2));
}

// helper for calculate_turbulence_kinetic_energy
static double tke_row(double height, double dx, Sim* sim)
{
    double tke = 0;

    Point first = Point(Constants::getBOX_XMIN() + 5 * dx, height);
    Point current = first;
    Vector current_velocity = Vector(0, 0);
    double current_rho = 0;

    unsigned int counter_x = 0;

    // calculate mean velocity of the row
    Vector mean_v = row_mean_velocity(height, dx, sim);

    VoronoiParticle* start = NULL;

    FILE* pFile = NULL;
    pFile = fopen("tke_matrix", "a");

    while (true)
    {
        current.x = first.x + counter_x * dx;

        if (!(current.x < Constants::getBOX_XMAX() - 5 * dx))
        {
            break;
        }

        current_velocity.x = sim_core::point_v_x(current, sim, &start);
        current_velocity.y = sim_core::point_v_y(current, sim, &start);
        current_rho = sim_core::point_density(current, sim, &start);

        tke += current_rho * (mean_v - current_velocity).length_square();

        if (options::color_on_the_fly)
        {
            fprintf(pFile, "%f\t", current_rho * (mean_v - current_velocity).length_square());
            start->m_color = current_rho * (mean_v - current_velocity).length_square();
        }

        counter_x++;
    }

    fprintf(pFile, "\n");
    fclose(pFile);

    tke *= 0.5;
    return tke;
}

// helper for calculate_turbulence_kinetic_energy
static Vector row_mean_velocity(double height, double dx, Sim* sim)
{
    Point mean_v = Point(0, 0);

    Point first = Point(Constants::getBOX_XMIN() + 5 * dx, height);
    Point current = first;

    unsigned int counter_x = 0;

    VoronoiParticle* start = NULL;

    // calculate mean velocity of the row
    while (true)
    {
        current.x = first.x + counter_x * dx;

        if (!(current.x < Constants::getBOX_XMAX() - 5 * dx))
        {
            break;
        }

        mean_v.x += sim_core::point_v_x(current, sim, &start);
        mean_v.y += sim_core::point_v_y(current, sim, &start);

        counter_x++;
    }

    mean_v /= counter_x;
    return mean_v;
}

// downsampling
void tools::downsampling(unsigned int pixel_per_length, char* filename, Sim* sim)
{
    int name_length = strlen(filename) - 3;

    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char dataFilename[64];
    strcpy(dataFilename, name);

    // var:
    // 1 for density, 2 for pressure, 3 for velocity_x, 4 for velocity_y

    if (options::var == 1)
    {
        const char ending[] = ".density_sample";
        strcat(dataFilename, ending);
    }
    else if (options::var == 2)
    {
        const char ending[] = ".pressure_sample";
        strcat(dataFilename, ending);
    }
    else if (options::var == 3)
    {
        const char ending[] = ".velocity_x";
        strcat(dataFilename, ending);
    }
    else if (options::var == 3)
    {
        const char ending[] = ".velocity_y";
        strcat(dataFilename, ending);
    }

    FILE* pFile;

    pFile = fopen(dataFilename, "w");

    double dn = 1. / pixel_per_length;

    Point first = Point(options::p.x + 0.5 * dn, options::p.y + 0.5 * dn);
    Point current = first;

    int counter = 0;
    int counter_y = 0;

    VoronoiParticle* start = NULL;

    while (true)
    {
        current.y = first.y + dn * counter_y;

        if (!(current.y < options::q.y))
        {
            break;
        }

        while (true)
        {
            current.x = first.x + dn * counter;

            if (!(current.x < options::q.x))
            {
                counter = 0;
                break;
            }

            if (options::var == 1)
            {
                fprintf(pFile, "%f\t", sim_core::point_density(current, sim, &start));
            }
            else if (options::var == 2)
            {
                fprintf(pFile, "%f\t", sim_core::point_pressure(current, sim, &start));
            }
            else if (options::var == 3)
            {
                fprintf(pFile, "%f\t", sim_core::point_v_x(current, sim, &start));
            }
            else if (options::var == 3)
            {
                fprintf(pFile, "%f\t", sim_core::point_v_y(current, sim, &start));
            }

            counter++;
        }

        fprintf(pFile, "\n");
        counter_y++;
    }

    fclose(pFile);
}

// color a particle
void tools::color(double x, double y, double color, Sim* sim)
{
    Point p = Point(x, y);
    VoronoiParticle* particle = point_location::findCell(&p, sim, NULL);
    particle->m_color = color;
}

// color the particles according to their volume
void tools::color_volume(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_color = sim->m_particles[i].m_cellVolume;
    }
}

// color the particles according to their roundness
void tools::color_roundness(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_color = sim->m_particles[i].roundness();
    }
}

// color the particles according to their vorticity
void tools::color_vorticity(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_color = sim->m_particles[i].vorticity();
    }
}

// color the particles according to their number of neighbours
void tools::color_nof_neighbours(Sim* sim)
{
    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        sim->m_particles[i].m_color = sim->m_particles[i].number_of_neighbours();
    }
}

// do a cross section and write ASCII data
void tools::cross_section(const char* filename)
{
    unsigned int name_length = strlen(filename) - 3;

    // check whether the filename length is > 3
    if (name_length <= 0)
    {
        fprintf(stderr, "ERROR in tools::cross_section: invalid filename %s\n", filename);
        exit(1);
    }

    // check whether the file ending is .hc
    if (!(filename[name_length] == '.' && filename[name_length + 1] == 'h' && filename[name_length + 2] == 'c'))
    {
        fprintf(stderr, "ERROR in tools::cross_section: filename %s has not the ending .hc\n", filename);
        exit(1);
    }

    // binary input and calculation of the tesselation
    Sim* sim = bin_io::binaryInput(filename);
    sim->create_box();
    sim->eventQueueIni();
    sim->eventQueueIniGhosts();
    sweepline::makeTesselation(sim);
    sim_graph::assertRealCellsClosed(sim);

    (Constants::get()->function_delete_zero_length_edges())(sim);
    sim_core::updateHalfEdgeLengthsAndCenters(sim);
    sim_graph::assertRealCellsClosed(sim);
    sim_core::updateCellVolumesAndCenters(sim);

    // initialization of the output file
    //

    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char sectionFilename[64];
    strcpy(sectionFilename, name);

    char ending[64];

    if (options::var == 1)
    {
        char ending1[32] = "_density.section";
        strcpy(ending, ending1);
    }
    else if (options::var == 2)
    {
        char ending2[32] = "_pressure.section";
        strcpy(ending, ending2);
    }
    else if (options::var == 3)
    {
        char ending3[32] = "_velocity_x.section";
        strcpy(ending, ending3);
    }
    else if (options::var == 4)
    {
        char ending4[32] = "_velocity_y.section";
        strcpy(ending, ending4);
    }

    strcat(sectionFilename, ending);

    FILE* pSection;

    pSection = fopen(sectionFilename, "w");

    fprintf(pSection, "#section: px = %f, py = %f, qx = %f, qy = %f, nofp = %d\n", options::p.x, options::p.y,
            options::q.x, options::q.y, options::nofp);

    assert(options::var != 0);

    if (options::var == 1)
    {
        fprintf(pSection, "#x		y		l		rho		cell_center_x	cell_center_y\n");
    }
    else if (options::var == 2)
    {
        fprintf(pSection, "#x		y		l		p		cell_center_x	cell_center_y\n");
    }
    else if (options::var == 3)
    {
        fprintf(pSection, "#x		y		l		v_x		cell_center_x	cell_center_y\n");
    }
    else if (options::var == 4)
    {
        fprintf(pSection, "#x		y		l		v_y		cell_center_x	cell_center_y\n");
    }

    Vector relative = options::q - options::p;
    relative /= (options::nofp - 1); // nofp >= 2, see options::check_options

    double dl = sqrt(relative.length_square()); // distance between two points
    double l = 0;

    Point start = options::p;
    Point current = Point(0, 0);

    VoronoiParticle* cell = NULL;
    VoronoiParticle* prev_cell = NULL;

    for (unsigned int i = 0; i < options::nofp; i++)
    {
        current = start + relative * i;
        l = i * dl;
        cell = point_location::findCell(&current, sim, cell);

        if (cell == prev_cell)
        {
            continue;
        }
        else
        {
            prev_cell = cell;
        }

        if (options::var == 1)
        {
            fprintf(pSection, "%f	%f	%f	%f	%f	%f\n", current.x, current.y, l, cell->m_density,
                    cell->m_cellCenterOfVolume.x, cell->m_cellCenterOfVolume.y);
        }
        else if (options::var == 2)
        {
            fprintf(pSection, "%f	%f	%f	%f	%f	%f\n", current.x, current.y, l, cell->m_pressure,
                    cell->m_cellCenterOfVolume.x, cell->m_cellCenterOfVolume.y);
        }
        else if (options::var == 3)
        {
            fprintf(pSection, "%f	%f	%f	%f	%f	%f\n", current.x, current.y, l, cell->m_velocity.x,
                    cell->m_cellCenterOfVolume.x, cell->m_cellCenterOfVolume.y);
        }
        else if (options::var == 4)
        {
            fprintf(pSection, "%f	%f	%f	%f	%f	%f\n", current.x, current.y, l, cell->m_velocity.y,
                    cell->m_cellCenterOfVolume.x, cell->m_cellCenterOfVolume.y);
        }
    }

    sim->freeGraph();
    delete sim;
    fclose(pSection);
}

// calculate gradients and write ASCII data
void tools::calculate_gradients(const char* filename)
{
    unsigned int name_length = strlen(filename) - 3;

    // check whether the filename length is > 3
    if (name_length <= 0)
    {
        fprintf(stderr, "ERROR in tools::calculate_gradients: invalid filename %s\n", filename);
        exit(1);
    }

    // check whether the file ending is .hc
    if (!(filename[name_length] == '.' && filename[name_length + 1] == 'h' && filename[name_length + 2] == 'c'))
    {
        fprintf(stderr, "ERROR in tools::calculate_gradients: filename %s has not the ending .hc\n", filename);
        exit(1);
    }

    // binary input and calculation of the tesselation
    Sim* sim = bin_io::binaryInput(filename);
    sim->create_box();
    sim->eventQueueIni();
    sim->eventQueueIniGhosts();
    sweepline::makeTesselation(sim);
    sim_graph::assertRealCellsClosed(sim);
    (Constants::get()->function_delete_zero_length_edges())(sim);
    sim_core::updateHalfEdgeLengthsAndCenters(sim);
    sim_graph::assertRealCellsClosed(sim);
    sim_core::updateCellVolumesAndCenters(sim);
    sim_core::calculateConservedVariables(sim);

    sim->m_box->determineGhostParticles();
    sim->freeGraph();
    sim->m_box->createGhostParticles();
    sim->m_box->reset();
    sim->eventQueueIni();
    sim->eventQueueIniGhosts();
    sweepline::makeTesselation(sim);
    (Constants::get()->function_delete_zero_length_edges())(sim);
    sim_core::updateHalfEdgeLengthsAndCenters(sim);
    sim_core::updateCellVolumesAndCenters(sim);
    sim_core::calculatePrimitiveVariables(sim);
    sim->m_box->copyPrimitiveVariablesCellCentersAndVolume();
    sim_core::calculate_gradients(sim);
    (Constants::get()->function_slopeLimitGradients())(sim);

    // initialization of the output file
    //

    char name[name_length + 1];
    name[0] = '\0';
    strncat(name, filename, name_length);

    char gradientsFilename[64];
    strcpy(gradientsFilename, name);

    char ending[64];

    if (options::var == 1)
    {
        char ending1[32] = "_density.gradients";
        strcpy(ending, ending1);

    }
    else if (options::var == 2)
    {
        char ending2[32] = "_pressure.gradients";
        strcpy(ending, ending2);
    }
    else if (options::var == 3)
    {
        char ending3[32] = "_velocity_x.gradients";
        strcpy(ending, ending3);
    }
    else if (options::var == 4)
    {
        char ending4[32] = "_velocity_y.gradients";
        strcpy(ending, ending4);
    }

    strcat(gradientsFilename, ending);
    FILE* pGradients;
    pGradients = fopen(gradientsFilename, "w");

    assert(options::var != 0);

    if (options::var == 1)
    {
        fprintf(pGradients, "#density gradients\n");
    }
    else if (options::var == 2)
    {
        fprintf(pGradients, "#pressure gradients\n");
    }
    else if (options::var == 3)
    {
        fprintf(pGradients, "#velocity_x gradients\n");
    }
    else if (options::var == 4)
    {
        fprintf(pGradients, "#velocity_y gradients\n");
    }

    fprintf(pGradients, "#x		y		grad_x		grad_y\n");

    VoronoiParticle* particle = NULL;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        if (options::var == 1)
        {
            fprintf(pGradients, "%f	%f	%f	%f\n", particle->m_location.x, particle->m_location.y,
                    particle->m_gradientDensity.x, particle->m_gradientDensity.y);
        }
        else if (options::var == 2)
        {
            fprintf(pGradients, "%f	%f	%f	%f\n", particle->m_location.x, particle->m_location.y,
                    particle->m_gradientPressure.x, particle->m_gradientPressure.y);
        }
        else if (options::var == 3)
        {
            fprintf(pGradients, "%f	%f	%f	%f\n", particle->m_location.x, particle->m_location.y,
                    particle->m_gradientVelocityX.x, particle->m_gradientVelocityX.y);
        }
        else if (options::var == 4)
        {
            fprintf(pGradients, "%f	%f	%f	%f\n", particle->m_location.x, particle->m_location.y,
                    particle->m_gradientVelocityY.x, particle->m_gradientVelocityY.y);
        }
    }

    sim->freeGraph();
    delete sim;
    fclose(pGradients);
}

static double minimum_pressure = std::numeric_limits<double>::max();
static double maximum_pressure = 0;

static double minimum_density = std::numeric_limits<double>::max();
static double maximum_density = 0;

static double minimum_velocity_x = std::numeric_limits<double>::max();
static double maximum_velocity_x = -std::numeric_limits<double>::max();

static double minimum_velocity_y = std::numeric_limits<double>::max();
static double maximum_velocity_y = -std::numeric_limits<double>::max();

static double minimum_color = std::numeric_limits<double>::max();
static double maximum_color = 0;

// find minimum / maximum value of a primitive variable
void tools::find_minimum_maximum(Sim* sim)
{
    VoronoiParticle* particle;

    for (unsigned int i = 0; i < sim->m_particles.counter(); i++)
    {
        particle = &sim->m_particles[i];

        if (particle->m_obstacle_particle)
        {
            continue;
        }

        if (particle->m_pressure < minimum_pressure)
        {
            minimum_pressure = particle->m_pressure;
        }

        if (particle->m_pressure > maximum_pressure)
        {
            maximum_pressure = particle->m_pressure;
        }

        if (particle->m_density < minimum_density)
        {
            minimum_density = particle->m_density;
        }

        if (particle->m_density > maximum_density)
        {
            maximum_density = particle->m_density;
        }

        if (particle->m_velocity.x > maximum_velocity_x)
        {
            maximum_velocity_x = particle->m_velocity.x;
        }

        if (particle->m_velocity.x < minimum_velocity_x)
        {
            minimum_velocity_x = particle->m_velocity.x;
        }

        if (particle->m_velocity.y > maximum_velocity_y)
        {
            maximum_velocity_y = particle->m_velocity.y;
        }

        if (particle->m_velocity.y < minimum_velocity_y)
        {
            minimum_velocity_y = particle->m_velocity.y;
        }

        if (particle->m_color < minimum_color)
        {
            minimum_color = particle->m_color;
        }

        if (particle->m_color > maximum_color)
        {
            maximum_color = particle->m_color;
        }
    }
}

// print minimum and maximum pressure
void tools::print_minimum_maximum(FILE* pFile)
{
    fprintf(pFile, "%f <= density <= %f\n", minimum_density, maximum_density);
    fprintf(pFile, "%f <= pressure <= %f\n", minimum_pressure, maximum_pressure);
    fprintf(pFile, "%f <= velocity_x <= %f\n", minimum_velocity_x, maximum_velocity_x);
    fprintf(pFile, "%f <= velocity_y <= %f\n", minimum_velocity_y, maximum_velocity_y);
    fprintf(pFile, "%f <= color <= %f\n", minimum_color, maximum_color);
}