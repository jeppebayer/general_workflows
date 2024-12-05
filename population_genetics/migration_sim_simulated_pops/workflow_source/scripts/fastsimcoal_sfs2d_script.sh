#!/bin/bash
#SBATCH --account EcoGenetics
#SBATCH --partition normal
#SBATCH --mem-per-cpu 12G
#SBATCH --cpus-per-task 8
#SBATCH --time 1-24:00:00

/faststorage/project/EcoGenetics/people/Sarah/fsc28/fsc28_linux64/fsc28 -t *.tpl -e *.est -n 100000 -m -M -L 40
