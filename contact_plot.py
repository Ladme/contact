"""
This script plots the .dat file generated by 'contact'.
The path to the file that should be plotted must be specified on the command line.
Example: contact_plot.py examples/contact_matrix.dat
The script generates output png file 'contact.png'.
Requires python3 and matplotlib library.

Released under MIT License.
Copyright (c) 2022 Ladislav Bartos
"""

def main():
    try:
        import sys
        import matplotlib.pyplot as plt
    except:
        print("This script requires the following libraries to be installed: sys, matplotlib")
        return

    if len(sys.argv) != 2:
        print("This script expects exactly one command line argument.")
        print("Please, provide the path to the .dat file that should be plotted.")
        return

    INPUT_DAT  = sys.argv[1]
    OUTPUT_PNG = "contact.png"

    xticks, yticks, matrix = [], [], []
    for line in open(INPUT_DAT):
        if line.strip()[0] == "#":
            continue
        if line.strip() == "":
            continue
        
        splitted = line.split()
        if line[0:7] == "       ":
            xticks = splitted
        else:
            yticks.append(splitted[0])
            matrix.append([float(x)*100 for x in splitted[1:]])

    plt.pcolormesh([x for x in range(len(xticks))], [y for y in range(len(yticks))], matrix, cmap="inferno", shading = "nearest")
    ax = plt.colorbar()

    ax.set_label("contact percentage [%]", fontsize = 14)
    plt.xticks([x for x in range(len(yticks))], fontsize = 5, labels = xticks)
    plt.xlabel("selection A", fontsize = 14)
    plt.yticks([x for x in range(len(xticks))], fontsize = 5, labels = yticks)
    plt.ylabel("selection B", fontsize = 14)

    plt.tight_layout()

    plt.savefig(OUTPUT_PNG, format = "png", dpi = 500)

if __name__ == "__main__":
    main()