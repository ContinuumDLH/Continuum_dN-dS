"""
                         a88888b.                     dP     oo
                        d8'   `88                     88
                        88        .d8888b. 88d888b. d8888P   dP   88d888b. dP    dP dP    dP 88d8b.d8b.
                        88        88'  `88 88'  `88   88     88   88'  `88 88    88 88    88 88'`88'`88
                        Y8.   .88 88.  .88 88    88   88     88   88    88 88.  .88 88.  .88 88  88  88
                         Y88888P' `88888P' dP    dP   dP     dP   dP    dP `88888P' `88888P' dP  dP  dP
                                           dP 888888ba       d8'       dP .d88888b
                                           88 88    `8b      d8'       88 88.    "'
                                     .d888b88 88     88      d8' .d888b88 `Y88888b.
                                     88'  `88 88     88      d8 '88'  `88       `8b
                                     88.  .88 88     88      d8 '88.  .88 d8'   .8P
                                     `88888P8 dP     dP      88  `88888P8  Y88888P
"""


# Imports Required =====================================================================================================

import tkinter as tk                                        # To use file dialog                    (Files)
from tkinter import filedialog as fd                        # To pull files for user ease           (Files)
import time                                                 # For time.sleep                        (Routine)
import subprocess                                           # To process Codeml                     (Routine)
import matplotlib.pylab as plt                              # To plot                               (Graph)
from matplotlib.font_manager import FontProperties as fp    # To adjust legend fonts                (Graph)
import statistics                                           # To find the Moo (mean)                (Graph)
import itertools                                            # To iterate colors and names           (Graph)

# Static Variables =====================================================================================================

title_of_graph  = "Pten-C dN/dS in 12 Drosophila Species"       # Insert string you want for TITLE      (User Needed)
y_plane         = "dN/dS Ratio (ω)"                             # Insert string you want for Y          (User Needed)
x_plane         = "10 Amino Acid / 30n Sliding Window"          # Insert string you want for X          (User Needed)
legend          = "Species"                                     # Insert string you want for X          (User Needed)
fontp           = fp().set_size("large")                        # Insert desired font size for legend   (User Needed)

window_size     = 30                                            # Windows for AA, not codons            (User Needed)
slice_start     = 0                                             # Where slicing start                   (User Needed)
slide           = 1                                             # Each window step +n                   (User Needed)

dnds_file    = open("dnds_output.txt", "w+")                    # Write dN/dS coordinates for slides    (DO NOT TOUCH)
codeml_path  = r"C:\Users\conti\bin\paml4.9j\bin\codeml.exe"    # Codeml location in system files       (User Needed)

writing_switch = False                                          # To know when to start writing         (DO NOT TOUCH)

fasta_titles      = []                                          # Holds FASTA file names                {LISTS}
dnds_parse        = []                                          # Post Parse before dnds_writing        {LISTS}
solo_positions    = []                                          # Holds dN/dS based on slide            {LISTS}
avgs_positions    = []                                          # Holds AVG dN/dS based on slide        {LISTS}

debug = False                                                   # Code review and T-shooting            (Optional)

# Pre-Checks ===========================================================================================================

if window_size % 3 != 0:
    print("The window size does not follow biological convention. Please update your window to a multiple of 3")
    exit(SyntaxError)

if codeml_path == "":
    print("A Codeml path must be provided in order to run this program. Please ensure Codeml is installed and you are"
          "pointing to the correct path. May I suggest the readme?")
    exit(OSError)

if title_of_graph == "":
    title_answer = input("Your graph will not have a title, is this okay [Y]/[N]?"
                         ">>> ").upper()
    if title_answer == "Y":
        print("Graph will not have a title, continuing...")
    else:
        exit(SyntaxError)


# Definitions ==========================================================================================================

def select_files():
    global file_to_list
    try:
        root = tk.Tk()
        file_to_list    = []
        requested_files = fd.askopenfilenames()

        for requests in requested_files:
            title       = open(requests, "r").readlines()[0]
            title       = title.replace("\n", "").replace(">", "")
            fasta_titles.append(title)
            contents    = open(requests, "r").readlines()[1:]
            delimited   = ""

            for lines in contents:
                delimited += lines.replace("\n", "")

            file_to_list.append([title, [delimited[i:i + window_size] for i in range(0, len(delimited), int(slide))]])

        for sequences in file_to_list:
            for windows in sequences[1][0:]:
                if len(windows) < window_size:
                    del sequences[1][sequences[1].index(windows)]

        root.destroy()

        if debug:
            print(file_to_list)
    except ImportError:
        print("Could not continue with the select_files definition. Aborting...")
        exit(ImportError)


def transform():
    try:
        print("Parse request " + str(slice_start + 1) +
              "/" + str(len(file_to_list[1][1]) - (window_size-(window_size-1))))
        codeml_input = open("tarin.phylip", "w+")
        codeml_input.write(" " + str(len(file_to_list)) + " " + str(window_size) + "\n")

        for sections in file_to_list:
            codeml_input.write(sections[0] + "\n" + sections[1][slice_start] + "\n")
        codeml_input.close()
        dnds_file.write("Parse request " + str(slice_start + 1) +
                        "/" + str(len(file_to_list[1][1]) - (window_size-(window_size-1))) + "\n")
        time.sleep(2)
        codeml_input.close()
    except:
        print("Am Error has occurred where writing the Codeml input file has failed. May I suggest the following:\n"
            "   1. Check if your initial files are from an aliment FASTA (From Jalview or even Ugene)\n"
            "   2. Check if your files are the same length, if they aren't - this WILL cause problems\n"
            "   3. Debug the code to print our codeml input file. Is this showing something we cannot explain?\n")
        exit(SyntaxError)


def codeml():
    try:
        print("Running Codeml on: " + str(slice_start + 1) +
              "/" + str(len(file_to_list[1][1]) - (window_size-(window_size-1))) + "\n")
        sprocess = subprocess.Popen(codeml_path, stdin=subprocess.PIPE)
        sprocess.stdin.write(b"\r\n")
        time.sleep(5)
    except OSError or OverflowError:
        print("No results have been found, Fatal Error. It may have occurred from one of the following..."
              "1. Subprocess has stopped responding"
              "2. Sprocess.stdin.write may have failed to read in BYTE b'' format"
              "3. The subprocess could not complete in a time suitable for the human lifespan...")


def find_results():
    global writing_switch, dnds_parse, fasta_titles, solo_positions, slice_start
    solo_hold = []
    aver_hold = []

    time.sleep(1)
    faucet = open("results.txt", "r")
    for lines in faucet.readlines():
        if "Use runmode = -2 for ML pai" in lines:
            lines.replace("Use runmode = -2 for ML pairwise comparison.", "")
            writing_switch = True
        if "pairwise comparison, codon frequencies:" in lines:
            lines.replace("pairwise comparison, codon frequencies:", "")
            writing_switch = False
        try:
            if writing_switch:
                for name in fasta_titles[1:]:
                    if name in lines:
                        aver_hold.append(float(lines[0:27].split()[-1].replace("\n", "")))
                        solo_hold.append([slice_start, lines[0:27].split()[-1]])
        except:
            print("Failed to append the correct information - continuing")
            continue


    try:
        print(aver_hold)
        aver_hold = statistics.mean(aver_hold)
        print(aver_hold)
        avgs_positions.append([slice_start, aver_hold])
        solo_positions.append(solo_hold)
    except:
        print('No value calculated - to keep consistency, no calculations performed.')

    faucet.close()
    slice_start += 1

    del aver_hold, solo_hold


def graphics():

    iter_titles = itertools.cycle(fasta_titles[1:])                 # fasta_titles that is iterative        {ITERTOOL}
    colors = itertools.cycle(
                ['red', 'darkorange', 'gold',
                'darkgreen', 'blue', 'darkviolet',
                'darkcyan', 'darkmagenta', 'chartreuse',
                'palevioletred', 'indigo'])                         # Colors needed for each title type     (User Needed)

    fig   = plt.figure(figsize=(20 , 5))
    axes  = fig.add_axes([0.1, 0.1, 0.72, 0.8])

    x_avg = []
    y_avg = []

    for positions in avgs_positions:
        x_avg.append(float(positions[0]))
        y_avg.append(float(positions[1]))

    title_count = 0

    for sequences in range(len(fasta_titles)-1):
        x_solo = []
        y_solo = []
        title_count += 1
        for positions in solo_positions:
            x_solo.append(float(positions[sequences][0]))
            y_solo.append(float(positions[sequences][1]))
        if title_count <= (len(fasta_titles) - 1):
            axes.plot(x_solo, y_solo, label=next(iter_titles), c=next(colors), clip_on=False, linewidth=2, snap='plot')
        else:
            axes.plot(x_solo, y_solo, c=title_count)
        del x_solo, y_solo

    axes.plot(x_avg, y_avg, color='black', label="Average", clip_on=False, linewidth=2, snap='plot')
    plt.title(title_of_graph)
    plt.grid(color='k', linestyle='-', linewidth=0.1)
    plt.ylim((-1, 70))
    plt.ylabel(y_plane)
    plt.xlabel(x_plane)
    axes.legend(title=legend, prop=fontp, fancybox=True, loc=2,
                bbox_to_anchor=(1.02, 1), borderaxespad=0., edgecolor='black')
    fig.savefig('image_output.tiff',
                dpi=100,
                format='png')
    plt.show()
    print(avgs_positions)



# Main =================================================================================================================


select_files()

while slice_start != (len(file_to_list[1][1]) - (window_size-(window_size-1))):
    transform()
    codeml()
    find_results()

graphics()
