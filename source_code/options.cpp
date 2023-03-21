#include "config.h"
#include "options.h"
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include "r2.h"

// see also manual for getopt_long
// the global variables are set to default values

// read in a binary file
char* options::input_filename = NULL;

// convert a file
char* options::convert_filename = NULL;
bool options::hybrid = false;
bool options::zero_grad = false;
double options::cut_percent = 0;
bool options::cut = false;
bool options::color_only = false;
bool options::density_only = false;
double options::x_min = 0;
double options::y_min = 0;
double options::x_max = 0;
double options::y_max = 0;
static bool x_min_set = false;
static bool y_min_set = false;
static bool x_max_set = false;
static bool y_max_set = false;

char* options::add_up_filename = NULL;

bool options::abs_vorticity = false;
bool options::kinetic_energy = false;
bool options::turbulence_kinetic_energy = false;
bool options::energy = false;
bool options::kinetic_energy_y = false;
bool options::mass_std_dev_square = false;
bool options::color_on_the_fly = false;

bool options::whatever = false;
char* options::whatever_filename = NULL;

// convert a file into ghost data
char* options::convert_ghost_filename = NULL;

// print the graph
char* options::graph_filename = NULL;

// do a section
char* options::section_filename = NULL;
Point options::p = Point(0, 0);
Point options::q = Point(0, 0);
unsigned int options::nofp = 0;

static bool px_set = false;
static bool py_set = false;
static bool qx_set = false;
static bool qy_set = false;
static bool nofp_set = false;

// calculate the gradient
char* options::gradient_filename = NULL;

// primitive variable, for section and gradient
unsigned int options::var = 0;

// find minimum and maximum primitive variable in binary data
char* options::minmax_filename = NULL;
bool options::subsequent = false;

// color a Voronoi particle
double options::color = 0;
bool options::color_set = false;
char* options::color_filename = NULL;

// color the particles according to their volume
char* options::color_volume_filename = NULL;
char* options::color_roundness_filename = NULL;
char* options::color_nof_neighbours_filename = NULL;
char* options::color_vorticity_filename = NULL;

// find indice
char* options::indice_filename = NULL;
unsigned int options::object = 0;
bool options::object_set = false;

// obstacle
bool options::obstacle = false;

// downsampling
char* options::downsampling_filename = NULL;
unsigned int options::pixel_per_length = 0;

// short options, semicolon after option means required argument
static const char short_options[] = "-hi:c:s:m:g:";

static const struct option long_options[] = {
        {"help",                      no_argument,       NULL, 'h'},
        {"input",                     required_argument, NULL, 'i'},
        {"convert",                   required_argument, NULL, 'c'},
        {"section",                   required_argument, NULL, 's'},
        {"gradients",                 required_argument, NULL, 'g'},
        {"px",                        required_argument, NULL, 300},
        {"py",                        required_argument, NULL, 301},
        {"qx",                        required_argument, NULL, 302},
        {"qy",                        required_argument, NULL, 303},
        {"nofp",                      required_argument, NULL, 304},
        {"var",                       required_argument, NULL, 305},
        {"minmax",                    required_argument, NULL, 'm'},
        {"sub",                       no_argument,       NULL, 306},
        {"graph",                     required_argument, NULL, 307},
        {"hybrid",                    no_argument,       NULL, 308},
        {"cut",                       no_argument,       NULL, 309},
        {"xmin",                      required_argument, NULL, 310},
        {"ymin",                      required_argument, NULL, 311},
        {"xmax",                      required_argument, NULL, 312},
        {"ymax",                      required_argument, NULL, 313},
        {"color",                     required_argument, NULL, 314},
        {"value",                     required_argument, NULL, 315},
        {"color_volume",              required_argument, NULL, 316},
        {"color_roundness",           required_argument, NULL, 317},
        {"color_nof_neighbours",      required_argument, NULL, 318},
        {"color_only",                no_argument,       NULL, 319},
        {"color_vorticity",           required_argument, NULL, 320},
        {"indice",                    required_argument, NULL, 321},
        {"object",                    required_argument, NULL, 322},
        {"convert_ghost",             required_argument, NULL, 323},
        {"obstacle",                  no_argument,       NULL, 324},
        {"zero_grad",                 required_argument, NULL, 325},
        {"whatever",                  required_argument, NULL, 326},
        {"add_up",                    required_argument, NULL, 327},
        {"abs_vorticity",             no_argument,       NULL, 328},
        {"turbulence_kinetic_energy", no_argument,       NULL, 329},
        {"kinetic_energy",            no_argument,       NULL, 330},
        {"energy",                    no_argument,       NULL, 331},
        {"kinetic_energy_y",          no_argument,       NULL, 332},
        {"downsampling",              required_argument, NULL, 333},
        {"pixel_per_length",          required_argument, NULL, 334},
        {"density_only",              no_argument,       NULL, 335},
        {"color_on_the_fly",          no_argument,       NULL, 336},
        {"mass_std_dev_square",       no_argument,       NULL, 337},
        {0,                           0,                 0,    0}
};

char* input_filename;

// parse program options
void options::parse_options(int argc, char** argv)
{
    for (;;)
    {
        int index;
        int c = 0;

        c = getopt_long(argc, argv, short_options, long_options, &index);

        if (c == EOF)
            break;

        switch (c)
        {
            case 0:
                break;

            case 'h':
                usage(argc, argv);
                exit(1);
                break;

            case 'i':
                input_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(input_filename, optarg);
                break;

            case 'c':
                convert_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(convert_filename, optarg);
                break;

            case 's':
                section_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(section_filename, optarg);
                break;

            case 'g':
                gradient_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(gradient_filename, optarg);
                break;

            case 300:
                p.x = strtod(optarg, NULL);
                px_set = true;
                break;

            case 301:
                p.y = strtod(optarg, NULL);
                py_set = true;
                break;

            case 302:
                q.x = strtod(optarg, NULL);
                qx_set = true;
                break;

            case 303:
                q.y = strtod(optarg, NULL);
                qy_set = true;
                break;

            case 304:
                nofp = atoi(optarg);
                nofp_set = true;
                break;

            case 305:
                if (strncmp("density", optarg, 10) == 0)
                {
                    var = 1;
                }
                else if (strncmp("pressure", optarg, 10) == 0)
                {
                    var = 2;
                }
                else if (strncmp("velocity_x", optarg, 10) == 0)
                {
                    var = 3;
                }
                else if (strncmp("velocity_y", optarg, 10) == 0)
                {
                    var = 4;
                }
                else
                {
                    fprintf(stderr, "ERROR in options::parse_options: unknown option %s\n", optarg);
                    fprintf(stderr, "valid options: pressure / density / velocity_x / velocity_y\n");
                    exit(1);
                }
                break;

            case 'm':
            {
                minmax_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(minmax_filename, optarg);

            }
                break;

            case 306:
            {
                subsequent = true;
            }
                break;

            case 307:
            {
                graph_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(graph_filename, optarg);
            }
                break;

            case 308:
            {
                hybrid = true;
            }
                break;

            case 309:
            {
                cut = true;
            }
                break;

            case 310:
            {
                x_min = strtod(optarg, NULL);
                x_min_set = true;
            }
                break;

            case 311:
            {
                y_min = strtod(optarg, NULL);
                y_min_set = true;
            }
                break;

            case 312:
            {
                x_max = strtod(optarg, NULL);
                x_max_set = true;
            }
                break;

            case 313:
            {
                y_max = strtod(optarg, NULL);
                y_max_set = true;
            }
                break;

            case 314:
            {
                color_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(color_filename, optarg);
            }
                break;

            case 315:
            {
                color = strtod(optarg, NULL);
                color_set = true;
            }
                break;

            case 316:
            {
                color_volume_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(color_volume_filename, optarg);
            }
                break;

            case 317:
            {
                color_roundness_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(color_roundness_filename, optarg);
            }
                break;

            case 318:
            {
                color_nof_neighbours_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(color_nof_neighbours_filename, optarg);
            }
                break;

            case 319:
            {
                color_only = true;
            }
                break;

            case 320:
            {
                color_vorticity_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(color_vorticity_filename, optarg);
            }
                break;

            case 321:
            {
                indice_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(indice_filename, optarg);
            }
                break;

            case 322:
            {
                if (strncmp("cell", optarg, 10) == 0)
                {
                    object = 1;
                    object_set = true;
                }

                else if (strncmp("ghost_cell", optarg, 10) == 0)
                {
                    object = 2;
                    object_set = true;
                }

                else if (strncmp("half_edge", optarg, 10) == 0)
                {
                    object = 3;
                    object_set = true;
                }

                else
                {
                    fprintf(stderr, "ERROR in options::parse_options: unknown option %s\n", optarg);
                    fprintf(stderr, "valid options: cell / ghost_cell / half_edge\n");
                    exit(1);
                }
            }
                break;

            case 323:
            {
                convert_ghost_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(convert_ghost_filename, optarg);
            }
                break;

            case 324:
            {
                obstacle = true;
            }
                break;

            case 325:
            {
                zero_grad = true;
                cut_percent = strtod(optarg, NULL);
            }
                break;

            case 326:
            {
                whatever = true;
                whatever_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(whatever_filename, optarg);
            }
                break;

            case 327:
            {
                add_up_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(add_up_filename, optarg);
            }
                break;

            case 328:
            {
                abs_vorticity = true;
            }
                break;

            case 329:
            {
                turbulence_kinetic_energy = true;
            }
                break;

            case 330:
            {
                kinetic_energy = true;
            }
                break;

            case 331:
            {
                energy = true;
            }
                break;

            case 332:
            {
                kinetic_energy_y = true;
            }
                break;

            case 333:
            {
                downsampling_filename = (char*) malloc(sizeof(char) * (strlen(optarg) + 1));
                strcpy(downsampling_filename, optarg);
            }
                break;

            case 334:
            {
                pixel_per_length = strtol(optarg, NULL, 0);

            }
                break;

            case 335:
            {
                density_only = true;
            }
                break;

            case 336:
            {
                color_on_the_fly = true;
            }
                break;

            case 337:
            {
                mass_std_dev_square = true;
            }
                break;

            default:
            {
                fprintf(stderr, "ERROR in options::parse_options: unknown option\n");
                exit(1);
            }
        }
    }
}

// check the options
void options::check_options(int argc, char** argv)
{
    if (downsampling_filename != NULL &&
        (pixel_per_length == 0 || var == 0 || !px_set || !py_set || !qx_set || !qy_set))
    {
        fprintf(stderr,
                "ERROR in options::check_options: specify pixel_per_length, var, px, py, qx and qy for downsampling\n");
        exit(1);
    }

    unsigned int add_up_counter = 0;

    if (mass_std_dev_square)
    {
        add_up_counter++;
    }

    if (abs_vorticity)
    {
        add_up_counter++;
    }

    if (kinetic_energy)
    {
        add_up_counter++;
    }

    if (energy)
    {
        add_up_counter++;
    }

    if (turbulence_kinetic_energy)
    {
        add_up_counter++;
    }

    if (kinetic_energy_y)
    {
        add_up_counter++;
    }

    if (add_up_counter > 1)
    {
        fprintf(stderr, "ERROR in options::check_options: only one variable can be added up.\n");
        exit(1);
    }

    if (section_filename != NULL)
    {
        if (!(px_set && py_set && qx_set && qy_set && nofp_set && (var != 0)))
        {
            fprintf(stderr, "ERROR in options::check_options: not all options for section are set.\n");
            usage(argc, argv);
            exit(1);

        }

        if (!(nofp >= 2))
        {
            fprintf(stderr, "ERROR in options::check_options: number of points (nofp) has to be at least 2.\n");
            exit(1);
        }
    }

    if (gradient_filename != NULL)
    {
        if (var == 0)
        {
            fprintf(stderr, "ERROR in options::check_options: not all options for calculate gradients are set.\n");
            usage(argc, argv);
            exit(1);
        }
    }

    if (hybrid && cut)
    {
        fprintf(stderr, "ERROR in options::check_options: only one option is allowed, hybrid or cut.\n");
        usage(argc, argv);
        exit(1);
    }

    if (cut)
    {
        if (!(x_min_set && y_min_set && x_max_set && y_max_set))
        {
            fprintf(stderr, "ERROR in options::check_options: not all coordinates of the box are set.\n");
            usage(argc, argv);
            exit(1);
        }
    }

    if (color_filename != NULL)
    {
        if (!(px_set && py_set && color_set))
        {
            fprintf(stderr,
                    "ERROR in options::check_options: you have to specify the position which has to be colored and the color (value)\n");
            usage(argc, argv);
            exit(1);
        }
    }

    if (indice_filename != NULL)
    {
        if (!(px_set && py_set && object_set))
        {
            fprintf(stderr,
                    "ERROR in options::check_options: you have to specify the object and its x-/ and y-coordinate\n");
            usage(argc, argv);
            exit(1);
        }
    }
}

// print program usage
void options::usage(int argc, char** argv)
{
    fprintf(stderr, "\nusage: %s [options] \n\n"

                    " -i | --input       continue the simulation from a file, *.hc (required)\n"
                    "\n"
                    " -c | --convert     convert a file, *.hc (required)\n"
                    "\t --color_only          no argument\n"
                    "\t --density_only        no argument\n"
                    "\t --obstacle            no argument\n"
                    "\t --zero_grad           percentage to cut, e.g. 0.01\n"
                    "\t --hybrid              no argument\n"
                    "\t --cut                 no argument\n"
                    "\t --xmin                minimum x box value\n"
                    "\t --ymin                minimum y box value\n"
                    "\t --xmax                maximum x box value\n"
                    "\t --ymax                maximum y box value\n"
                    "\n"
                    " --convert_ghost    generate ghost data from a file, *.hc (required)\n"
                    "\n"
                    " -s | --section     cross sect data, *.hc (required)\n"
                    "\t --px                  x-coord of first point (required)\n"
                    "\t --py                  y-coord of first point (required)\n"
                    "\t --qx                  x-coord of second point (required)\n"
                    "\t --qy                  y-coord of second point (required)\n"
                    "\t --nofp                number of points (required)\n"
                    "\t --var                 primitive variable, \"density\" or \"pressure\" or \"velocity_x\" or \"velocity_y\" (required)\n"
                    "\n"
                    " -g | --gradient    calculate the gradients for a file, *.hc (required)\n"
                    "\t --var                 primitive variable, \"density\" or \"pressure\" or \"velocity_x\" or \"velocity_y\" (required)\n"
                    "\n"
                    " -m | --minmax      evaluate the minimum and maximum values of a file: *.hc (required)\n"
                    "\t --sub                 Also include the subsequent files, no argument\n"
                    "\n"
                    " --graph            print the graph for a file, *.hc (required)\n"
                    "\n"
                    " --color            color a file, *.hc (required)\n"
                    "\t --value               color of the cell, double (required)\n"
                    "\t --px                  x-coord in a cell (required)\n"
                    "\t --py                  y-coord in a cell (required)\n"
                    "\n"
                    "--color_volume      color cells according to their volume, *.hc (required)\n"
                    "--color_roundness   color cells according to their roundness, *.hc (required)\n"
                    "--color_nof_neighbours color cells according to their number of neighbours, *.hc (required)\n"
                    "--color_vorticity   color cells according to their vorticity, *.hc (required)\n"
                    "\n"
                    "--indice            find an index in a file, *.hc (required)\n"
                    "\t --object              object, \"cell\" or \"ghost_cell\" or \"half_edge\" (required)\n"
                    "\t --px                  x-coordinate (required)\n"
                    "\t --py                  y-coordinate (required)\n"
                    "--whatever          run the temporary code in tools.cpp for the file, *.hc (required)\n"
                    "--add_up            add up a quantity, *.hc (required)\n"
                    "\t --energy              no argument\n"
                    "\t --kinetic_energy      no argument\n"
                    "\t --turbulence_kinetic_energy no argument\n"
                    "\t --kinetic_energy_y    no argument\n"
                    "\t --mass_std_dev_square no argument\n"
                    "\t --color_on_the_fly    no argument\n"
                    "--downsampling      execute a downsampling, *.hc (required)\n"
                    "\t --pixel_per_length    integer (required)\n"
                    "\t --var                 primitive variable, \"density\" or \"pressure\" or \"velocity_x\" or \"velocity_y\" (required)\n"
                    "\t --px --py --qx --qy   doubles determining the box (required)\n", argv[0]);
}

// free memory
void options::clear()
{
    free(input_filename);
    free(convert_filename);
    free(section_filename);
    free(minmax_filename);
    free(gradient_filename);
    free(graph_filename);
    free(color_filename);
    free(color_volume_filename);
    free(color_roundness_filename);
    free(color_nof_neighbours_filename);
    free(color_vorticity_filename);
    free(indice_filename);
    free(convert_ghost_filename);
    free(whatever_filename);
    free(add_up_filename);
    free(downsampling_filename);
}