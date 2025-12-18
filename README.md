# sidisFESAPV
Sidis Italia shared code


# Basic workflow:


```
git clone https://github.com/nicolovalle/sidisFESAPV
```

Only once:
```
git config user.name "Your name"
git config user.email "your@email"
```

Edit code, then:
```
git add file.cpp
git commit -m "Your comment"
git push
```

## Content

### Macros
1. **epic_studies.cpp** runs over generated files (Pythia8), loops over all events, and fills a dedicated tree for the electron (MC and reco) and for the hadrons (MC and reco), storing all the relevant information (e.g. px, py, pz, PDG, Q2, xB, ... ). It takes as input a list file (.txt): the first line is the name of the output file, and each subsequent line is the path to an input file. To run the code from the terminal: 'root -l epic_studies.cpp\(\"lista_epic_studies.txt\"\)'.

2. **pion.plot2.cpp** runs over the epic_studies.cpp output and produce several plots

3. **relevant.plots.cpp** runs over the output produced by epic_studies.cpp and generates a set of relevant plots for a chosen particle species, saving them in a single .root output file. The macro takes as input: the pdg code of the particle of interest (211 for positve pions, -211 for negative pions, 321 for positive kaons and -321 for negative kaons) and the input directory (e.g. 25.10_10x100) containing the .root files produced by epic_studies.cpp (and where the macro will also write its output).

### Workflow

To run it locally:
1. Move inside the environment ./eic-shell
2. Create the list file (.txt)
3. Run the epic_studies.cpp macro in order to store the information in the trees
4. Create a dedicated folder for the campaign, Q2 and energy of interest.
5. Move the output file produced by epic_studies.cpp into this folder.
6. Run the relevant_plots.cpp macro for each particle of interest to produce the plots (saved in a .root output file).

# Useful information

### MC cross sections
Pythia8 DIS NC cross sections available [here](https://docs.google.com/spreadsheets/d/1yCUZFJjMbE-Ly73juGKYvsA-VufUplaDMWs4Uv7PqAk/edit?usp=sharing). (from Ralf who got them from Tyler)
