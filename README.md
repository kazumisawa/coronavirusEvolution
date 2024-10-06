# coronavirusEvolution

Python codeds for the following paper:
Misawa and Ootsuki (2024)
A simple method for estimating time-irreversible nucleotide substitution rates in the SARS-CoV-2 genome
NAR Genomics and Bioinformatics 6(1):1-7
https://doi.org/10.1093/nargab/lqae009

pairRate.py Manual
Function: Estimates the evolutionary rate by comparing two FASTA files: one representing an ancestral sequence and the other representing descendant sequences.

Usage:
python pairRate.py ancestor.fasta descendants.fasta

Arguments:
ancestor.fasta: The FASTA file containing the ancestral sequence.
descendants.fasta: The FASTA file containing the descendant sequences.
Description: This script takes two input FASTA files and calculates the evolutionary rate by comparing the ancestral sequence to the descendant sequences. The evolutionary rate is output based on the differences between the two sequences over a given period.

trajectory.py Manual
Function: Predicts the expected nucleotide composition of a sequence after a specified time, using an evolutionary rate.

Usage:
python trajectory.py ancestor.fasta time_in_years oneYearRate.txt

Arguments:
ancestor.fasta: The FASTA file containing the initial (ancestral) sequence.
time_in_years: A floating-point value representing the number of years for which the evolution is predicted.
oneYearRate.txt: A text file containing the evolutionary rate (per year) of the sequences.
Description: This script calculates the expected nucleotide composition of the sequence over a specified period, using the given evolutionary rate. It outputs the predicted sequence's composition after the provided number of years.

Summary of Key Features:
pairRate.py: Estimates evolutionary rates by comparing ancestral and descendant sequences.
trajectory.py: Predicts future sequence evolution using a given evolutionary rate and time.
These tools are valuable for modeling viral evolution and understanding the dynamics of sequence changes over time.
