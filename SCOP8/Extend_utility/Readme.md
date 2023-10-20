## Readme

This stand-alone utility program is for users to transform the original lattice via a transformation matrix(written in 'mat.inp')



1. input file list (refer to folder <input_example>)
   - fc2.dat
   - fc3.dat
   - fc4.dat
   - lat_fc.dat
   - structure.params
   - mat.inp(this is the file you input the transformation matrix $P$)
2. output files list (refer to folder <output_example>)
   - fc2_new.dat
   - fc3_new.dat
   - fc4_new.dat
   - lat_fc_new.dat
   - structure_new.params

Regarding the transformation matrix $\textbf{P}$, given origina lattice transvectors are $\textbf{a},\textbf{b},\textbf{c}$, new lattice transvectors are $\textbf{a}^{new},\textbf{b}^{new},\textbf{c}^{new}$ , we have
$\begin{pmatrix}
\textbf{a}^{new}\\ \textbf{b}^{new} \\ \textbf{c}^{new}
\end{pmatrix} = \textbf{P}\cdot \begin{pmatrix}\textbf{a}\\ \textbf{b} \\ \textbf{c} \end{pmatrix}$


