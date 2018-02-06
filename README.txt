The project aims at simulating the thermodynamics of polymer sequences binding. The focus is on DNA-DNA
interaction, intending to be applied to a CRISPR-CAS9 system, however it can be easily generalized to
other barcoding problems.

What follows is a list of the files in the project, ordered by their relation, with a short description
of each.

Python scripts
+ DNA binding simulations
    - model1.py
    simplest model of the interaction. Ignores genome probabilities, different energies
        distributions other than a matching and an unmatching.
    - model3.py
    most refined model. Includes probabilities, an energy distribution matrix.

+ Simplified models
    - distribution.py
    evaluates the distribution of binding selectivities based on a simplified model.

+ Analysis
    - correlations.py
    calculates probability matrix prob_matrix.txt based on real-world genome data
    - fill_matrix.py
    fills the missing gaps in the santa_lucia04 matrix.
    - length_analysis.py
    investigates the accuracy of the mathematical model vs the computational model for different values of the
        targeting sequence length l.

+ Miscellaneous/incomplete
    - defect_approx.py      approximates gaussian defect distribution.
    - defect_distribution.py generates the sequences randomly and plots defect distribution
    - kinetic_model1.py     kinetic model. Highly incomplete.
    - with_pam.py           includes PAM condition
    - his_en.py             creates histogram of Santa Lucia distribution of energies from matrix.
    - minimization.py
    - model3_obj.py         same as previous, rephrased using OOP. For efficiency.
