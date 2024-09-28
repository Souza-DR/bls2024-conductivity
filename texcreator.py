import numpy as np
import sys
from pdb import set_trace

ndec = 2

loaded_data = np.load("./results/data.npz")
loaded_data1 = np.load("./results/data1.npz")
finit_array1 = loaded_data1["finit_array"]
ffinal_array1 = loaded_data1["ffinal_array"]

nsites_array = (loaded_data["nsites_array"]).astype(int)
nsources_array = (loaded_data["nsources_array"]).astype(int)
nmesh_array = loaded_data["nmesh_array"]
noise_coeff_array = loaded_data["noise_coeff_array"]

noise_level_array = [100]*loaded_data["noise_level_array"]
finit_array = loaded_data["finit_array"]
ffinal_array = loaded_data["ffinal_array"]
normgpinit_array = loaded_data["normgpinit_array"]
normgpfinal_array  = loaded_data["normgpfinal_array"]

erroropt_array  = [100]*loaded_data["erroropt_array"]
errorinit_array  = [100]*loaded_data["errorinit_array"]
flagsol_array = (loaded_data["flagsol_array"]).astype(int)
iter_array = (loaded_data["iter_array"]).astype(int)
numevalf_array = (loaded_data["numevalf_array"]).astype(int)
CPU_time_array = loaded_data["CPU_time_array"]

ffinal1 = []
ffinal2 = []
ffinal3 = []
vencedor = []

for l in range(len(nsites_array)):
    if nsources_array[l] == 1:
        ffinal1.append(ffinal_array1[l])
    # if nsources_array[l] == 2:
    #     ffinal2.append(ffinal_array1[l])
    if nsources_array[l] == 3:
        ffinal3.append(ffinal_array1[l])

best_ffinal1 = np.min(np.array(ffinal1))
# best_ffinal2 = np.min(np.array(ffinal2))
best_ffinal3 = np.min(np.array(ffinal3))

# print(best_ffinal1)


for l in range(len(nsites_array)):
    if nsources_array[l] == 1 and ffinal_array1[l] <= best_ffinal1:
        vencedor.append(l)
    # if nsources_array[l] == 2 and ffinal_array1[l] <= best_ffinal2:
    #     vencedor.append(l)
    if nsources_array[l] == 3 and ffinal_array1[l] <= best_ffinal3:
        vencedor.append(l)

# print(vencedor)
# sys.exit()

alphabet = []
for i in range(97, 123):
    alphabet.append(chr(i))

with open("./results/Results.tex", "w") as f:
    f.write(r"""
\documentclass[a4paper,10pt]{article}
\usepackage[utf8]{inputenc}

\usepackage{amsmath, amsthm, amssymb}
\usepackage{fullpage}
\usepackage{dsfont}
\usepackage{bm}
\usepackage{xcolor}
\usepackage{stmaryrd}
\usepackage{graphicx}
\usepackage{mathrsfs}
\newcommand{\bsi}{{\bf a}}
\newcommand{\ck}{{\kappa_0}}
\newcommand{\im}{\alpha}
\usepackage{multirow}

\graphicspath{{./figures/}}

\begin{document}      
    

\begin{table}[ht!]
\centering
\resizebox{0.7\textheight}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c||c|c|c}
\hline
$\ck$ & $\bar \im$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ & $G(\bsi^0)$ & $G(\hat\bsi)$ & $\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time & Flag & $G(\hat\bsi)$ \\ \hline \hline
""")
    for i in range(len(nsites_array)):
#         if i == 1:
#             f.write(
# str(nsites_array[i]) + r""" & """+str('%4.2f' % noise_level_array[i])+r"""& """+str('%5.2f' % errorinit_array[i])+r""" & """+str('%5.2f' % erroropt_array[i])+r""" & """+str('%4.2f' % finit_array[i])+r""" & """+str('%7.5f' % ffinal_array[i])+r""" & """+str('%7.5f' % normgpinit_array[i])+r""" & """+str('%7.5f' % normgpfinal_array[i])+r""" & """+str(iter_array[i])+r""" & """+str(numevalf_array[i])+r""" & """+str('%7.2f' % CPU_time_array[i])+r""" & """+str(flagsol_array[i]) + r""" & """+str('%10.3e' % ffinal_array1[i])+r"""\\ \hline
#                 """)

#         f.write(r"""
# \multirow{10}{*}{"""+str(nsites_array[i])+r"""} & \multirow{10}{*}{"""+str('%4.2f' % noise_level_array[i])+r"""}"""+str('%5.2f' % errorinit_array[i])+r""" & """+str('%5.2f' % erroropt_array[i])+r""" & """+str('%4.2f' % finit_array[i])+r""" & """+str('%7.5f' % ffinal_array[i])+r""" & """+str('%7.5f' % normgpinit_array[i])+r""" & """+str('%7.5f' % normgpfinal_array[i])+r""" & """+str(iter_array[i])+r""" & """+str(numevalf_array[i])+r""" & """+str('%7.2f' % CPU_time_array[i])+r""" & """+str(flagsol_array[i]) + r""" & """+str('%10.3e' % ffinal_array1[i])+r"""\\ \hline
#                 """)
        
        if i == vencedor[0] or i == vencedor[1]:
        # if i == vencedor[0]:
            f.write(str(nsites_array[i]) + r""" & """+str(nsources_array[i]) + r""" & """+str('%4.2f' % noise_level_array[i])+r"""& """+str('%5.2f' % errorinit_array[i])+r""" & """+str('%5.2f' % erroropt_array[i])+r""" & """+str('%4.2f' % finit_array[i])+r""" & """+str('%7.5f' % ffinal_array[i])+r""" & """+str('%7.5f' % normgpinit_array[i])+r""" & """+str('%7.5f' % normgpfinal_array[i])+r""" & """+str(iter_array[i])+r""" & """+str(numevalf_array[i])+r""" & """+str('%7.2f' % CPU_time_array[i])+r""" & """+str(flagsol_array[i]) + r""" & \textbf{"""+str('%10.3e' % ffinal_array1[i])+r"""}\\ \hline
                """)
        else:
            f.write(str(nsites_array[i]) + r""" & """+str(nsources_array[i]) + r""" & """+str('%4.2f' % noise_level_array[i])+r"""& """+str('%5.2f' % errorinit_array[i])+r""" & """+str('%5.2f' % erroropt_array[i])+r""" & """+str('%4.2f' % finit_array[i])+r""" & """+str('%7.5f' % ffinal_array[i])+r""" & """+str('%7.5f' % normgpinit_array[i])+r""" & """+str('%7.5f' % normgpfinal_array[i])+r""" & """+str(iter_array[i])+r""" & """+str(numevalf_array[i])+r""" & """+str('%7.2f' % CPU_time_array[i])+r""" & """+str(flagsol_array[i]) + r""" & """+str('%10.3e' % ffinal_array1[i])+r"""\\ \hline
                """)
        
        if (i+1)%10 == 0:
            f.write(r"""\hline
  """)
        
        if (i+1)%40 == 0:
            f.write(r"""
\end{tabular}}
\end{table}""")
            
            f.write(r"""
\begin{table}[ht!]
\centering
\resizebox{0.7\textheight}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c||c|c|c}
\hline
$\ck$ & $\bar \im$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ & $G(\bsi^0)$ & $G(\hat\bsi)$ & $\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time & Flag & $G(\hat\bsi)$ \\ \hline
                    """)
                
    
    f.write(r"""
\end{tabular}}
\end{table}""")
    

    f.write(r"""\begin{figure}[ht]
    \begin{center}
    \begin{tabular}{|c|c|c|c|}
    \hline
    & Ground truth & Initialization & Reconstruction\\
    \hline
    \hline""")

    l = 0
    name = 'blsfig6'
    for i in vencedor:

        # if i % 4 == 0:
        #     l = 0
        #     name = 'blsfig'+str((i//4))

        f.write(r"""\multirow{3}{*}{\rotatebox{90}{$\ck=""" + str(nsites_array[i]) + r"""\;$}} & & & Noise = $"""+str('%4.2f' % noise_level_array[i])+r"""\%$  \\
        \cline{4-4}
        & & $E(\bsi^0) = """+str('%5.2f' % errorinit_array[i])+r"""\% $ & $E(\hat\bsi) = """+str('%5.2f' % erroropt_array[i])+r"""\% $ \\""")
        for j in range(3):
            f.write(r"""
            & \includegraphics[scale=0.8]{./figures/"""+name+str(alphabet[l])+r""".eps}""")
            l += 1
        f.write(r"""\\ \hline""")
        
        # if (i+1)%4 == 0:
        #     f.write(r"""
        # \end{tabular}
        # \end{center}
        # \end{figure}
        # """)
        #     f.write(r"""\begin{figure}[ht]
        # \begin{center}
        # \begin{tabular}{|c|c|c|c|}
        # \hline
        # & Ground truth & Initialization & Reconstruction\\
        # \hline
        # \hline""")
    
    f.write(r"""
        \end{tabular}
        \end{center}
        \end{figure}
        """)

    
    f.write(r"""\end{document}""")




# import sys
# import numpy as np

# loaded_data = np.load("./results/data.npz")

# nsites_array = loaded_data["nsites_array"]
# nsources_array = loaded_data["nsources_array"]
# nmesh_array = loaded_data["nmesh_array"]
# noise_coeff_array = loaded_data["noise_coeff_array"]

# noise_level_array = [100]*loaded_data["noise_level_array"]
# finit_array = loaded_data["finit_array"]
# ffinal_array = loaded_data["ffinal_array"]
# normgpinit_array = loaded_data["normgpinit_array"]
# normgpfinal_array  = loaded_data["normgpfinal_array"]

# erroropt_array  = [100]*loaded_data["erroropt_array"]
# errorinit_array  = [100]*loaded_data["errorinit_array"]
# flagsol_array = loaded_data["flagsol_array"]
# iter_array = loaded_data["iter_array"]
# numevalf_array = loaded_data["numevalf_array"]
# CPU_time_array = loaded_data["CPU_time_array"]

# alphabet = []
# for i in range(97, 123):
#     alphabet.append(chr(i))

# with open("./results/Results.tex", "w") as f:
#     f.write(r"""\documentclass[a4paper,10pt]{article}
#     \usepackage[utf8]{inputenc}

#     \usepackage{amsmath, amsthm, amssymb}
#     \usepackage{fullpage}
#     \usepackage{dsfont}
#     \usepackage{bm}
#     \usepackage{xcolor}
#     \usepackage{stmaryrd}
#     \usepackage{graphicx}
#     \usepackage{mathrsfs}
#     \newcommand{\bsi}{{\bf a}}
#     \newcommand{\ck}{{\kappa_0}}
#     \usepackage{multirow}

#     \graphicspath{{./figures/}}

#     \begin{document}      

#     \begin{table}[ht]
#     \centering
#     \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}
#     \hline 
#     $\ck$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ &$G(\bsi^0)$ &$G(\hat\bsi)$ &$\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time \\ \hline""")
#     for i in range(len(nsites_array)):
#         f.write(str(nsites_array[i]) + r""" & """+str('%4.2f' % noise_level_array[i])+r"""& """+str('%5.2f' % errorinit_array[i])+r""" & """+str('%5.2f' % erroropt_array[i])+r""" & """+str('%4.2f' % finit_array[i])+r""" & """+str('%7.5f' % ffinal_array[i])+r""" & """+str('%7.5f' % normgpinit_array[i])+r""" & """+str('%7.5f' % normgpfinal_array[i])+r""" & """+str(iter_array[i])+r""" & """+str(numevalf_array[i])+r""" & """+str('%7.2f' % CPU_time_array[i])+r""" \\ \hline""")

#         if (i+1)%40 == 0:
#             f.write(r"""
#             \end{tabular}
#             \end{table}""")
        
#             f.write(r"""\begin{table}[ht]
#         \centering
#         \begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|}
#         \hline 
#         $\ck$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ &$G(\bsi^0)$ &$G(\hat\bsi)$ &$\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time \\ \hline""")
        
#     f.write(r"""
#     \end{tabular}
#     \end{table}""")

#     f.write(r"""\begin{figure}[!htbp]
#     \begin{center}
#     \begin{tabular}{|c|c|c|c|}
#     \hline
#     & Ground truth & Initialization & Reconstruction\\
#     \hline
#     \hline""")

#     for i in range(len(nsites_array)):
#         if i % 4 == 0:
#             l = 0
#             name = 'blsfig'+str((i//4))

#         f.write(r"""\multirow{3}{*}{\rotatebox{90}{$\ck=""" + str(nsites_array[i]) + r"""\;$}} & & & Noise = $"""+str('%4.2f' % noise_level_array[i])+r"""\%$  \\
#         \cline{4-4}
#         & & $E(\bsi^0) = """+str('%5.2f' % errorinit_array[i])+r"""\% $ & $E(\hat\bsi) = """+str('%5.2f' % erroropt_array[i])+r"""\% $ \\""")
#         for j in range(3):
#             f.write(r"""
#             & \includegraphics[scale=0.8]{./figures/"""+name+str(alphabet[l])+r""".eps}""")
#             l += 1
#         f.write(r"""\\ \hline""")
        
#         if (i+1)%4 == 0:
#             f.write(r"""
#         \end{tabular}
#         \end{center}
#         \end{figure}
#         """)
#             f.write(r"""\begin{figure}[!htbp]
#         \begin{center}
#         \begin{tabular}{|c|c|c|c|}
#         \hline
#         & Ground truth & Initialization & Reconstruction\\
#         \hline
#         \hline""")
    
#     f.write(r"""
#         \end{tabular}
#         \end{center}
#         \end{figure}
#         """)
#     f.write(r"""\end{document}""")