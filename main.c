// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include <unistd.h>
#include <groan.h>

const char VERSION[] = "v2022/08/30";

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **atoms1,
        char **atoms2,
        float *cutoff) 
{
    int gro_specified = 0, atoms1_specified = 0, atoms2_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:o:a:b:u:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file name
        case 'o':
            *output_file = optarg;
            break;
        // selection A
        case 'a':
            *atoms1 = optarg;
            atoms1_specified = 1;
            break;
        // selection B
        case 'b':
            *atoms2 = optarg;
            atoms2_specified = 1;
            break;
        // cut-off for contact
        case 'u':
            if (sscanf(optarg, "%f", cutoff) != 1) {
                fprintf(stderr, "Could not read cut-off value.\n");
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !atoms1_specified || !atoms2_specified) {
        fprintf(stderr, "Gro file and atoms specification must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -a SELECTION1 -b SELECTION2 [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file name (default: contacts.dat)\n");
    printf("-a STRING        selection of atoms\n");
    printf("-b STRING        selection of atoms\n");
    printf("-u FLOAT         cut-off for contact [nm] (default: 0.5)\n");
    printf("\n");
}

/*
 * Prints parameters that the program will use for the calculation.
 */
void print_arguments(
        FILE *stream,
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *atoms1,
        const char *atoms2,
        const float cutoff)
{
    fprintf(stream, "\nParameters for Contact Matrix calculation:\n");
    fprintf(stream, ">>> gro file:         %s\n", gro_file);
    if (xtc_file == NULL) fprintf(stream, ">>> xtc file:         ----\n");
    else fprintf(stream, ">>> xtc file:         %s\n", xtc_file);
    fprintf(stream, ">>> ndx file:         %s\n", ndx_file);
    fprintf(stream, ">>> output file:      %s\n", output_file);
    fprintf(stream, ">>> atoms1:           %s\n", atoms1);
    fprintf(stream, ">>> atoms2:           %s\n", atoms2);
    fprintf(stream, ">>> cut-off:          %f nm\n\n", cutoff);
}

void matrix_frame(
        const atom_selection_t *selection1, 
        const atom_selection_t *selection2, 
        size_t **matrix,
        box_t box,
        const float cutoff)
{
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        for (size_t j = 0; j < selection2->n_atoms; ++j) {
            if ( distance3D(selection1->atoms[i]->position, selection2->atoms[j]->position, box) < cutoff) ++(matrix[i][j]); 
        }
    }
}

void matrix_destroy(size_t **matrix, const size_t maxi) 
{
    for (size_t i = 0; i < maxi; ++i) {
        free(matrix[i]);
    }

    free(matrix);
}

int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = "contacts.dat";
    char *atoms1 = NULL;
    char *atoms2 = NULL;
    float cutoff = 0.5;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &atoms1, &atoms2, &cutoff) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    print_arguments(stdout, gro_file, xtc_file, ndx_file, output_file, atoms1, atoms2, cutoff);

    // try opening output file
    FILE *output = NULL;
    output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Could not open output file %s.\n", output_file);
        return 1;
    }

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

     // try reading ndx file (ignore if this fails)
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select specified atoms
    atom_selection_t *selection1 = smart_select(all, atoms1, ndx_groups);
    if (selection1 == NULL || selection1->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", atoms1);
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(selection1);
        fclose(output);
        return 1;
    }

    atom_selection_t *selection2 = smart_select(all, atoms2, ndx_groups);
    if (selection2 == NULL || selection2->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", atoms2);
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(selection1);
        free(selection2);
        fclose(output);
        return 1;
    }

    free(all);
    all = NULL;

    // prepare contact matrix
    size_t **contact_matrix = calloc(selection1->n_atoms, sizeof(size_t *));
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        contact_matrix[i] = calloc(selection2->n_atoms, sizeof(size_t));
    }

    // if there is no xtc file supplied, just use gro file
    size_t n_frames = 0;
    if (xtc_file == NULL) {
        matrix_frame(selection1, selection2, contact_matrix, system->box, cutoff);
        ++n_frames;
    } else {
        XDRFILE *xtc = xdrfile_open(xtc_file, "r");
        if (xtc == NULL) {
            fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);
            dict_destroy(ndx_groups);
            free(system);
            matrix_destroy(contact_matrix, selection1->n_atoms);
            free(selection1);
            free(selection2);
            fclose(output);
            return 1;
        }

        // check that the gro file and the xtc file match each other
        if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
            fprintf(stderr, "Number of atoms in %s does not match %s.\n", xtc_file, gro_file);
            xdrfile_close(xtc);
            dict_destroy(ndx_groups);
            free(system);
            matrix_destroy(contact_matrix, selection1->n_atoms);
            free(selection1);
            free(selection2);
            fclose(output);
            return 1;
        }

        while (read_xtc_step(xtc, system) == 0) {
            // print info about the progress of reading
            if ((int) system->time % PROGRESS_FREQ == 0) {
                printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
                fflush(stdout);
            }

            matrix_frame(selection1, selection2, contact_matrix, system->box, cutoff);
            ++n_frames;
        }

        xdrfile_close(xtc);
    }

    // write output
    fprintf(output, "# Generated with contact (C Contact Matrix Calculator) %s.\n", VERSION);
    fprintf(output, "# Command line: ");
    for (int i = 0; i < argc; ++i) {
        fprintf(output, "%s ", argv[i]);
    }
    fprintf(output, "\n       ");
    for (size_t i = 0; i < selection1->n_atoms; ++i) {
        fprintf(output, "%6d ", selection1->atoms[i]->atom_number);
    }
    fprintf(output, "\n");

    for (size_t j = 0; j < selection2->n_atoms; ++j) {
        fprintf(output, "%6d ", selection2->atoms[j]->atom_number);
        for (size_t i = 0; i < selection1->n_atoms; ++i) {
            fprintf(output, "%6.3f ", (float) contact_matrix[i][j] / n_frames);
        }
        fprintf(output, "\n");
    }

    printf("\n");
    fclose(output);
    matrix_destroy(contact_matrix, selection1->n_atoms);
    free(selection1);
    free(selection2);
    free(system);
    dict_destroy(ndx_groups);

    return 0;
}