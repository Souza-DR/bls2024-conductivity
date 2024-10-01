import numpy as np
import sys
from pdb import set_trace

ndec = 2
color = "blue"

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

nsites_list = np.unique(nsites_array)
ninit_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
nsources_list = np.unique(nsources_array)
nmesh_list = np.unique(nmesh_array)
noise_coeff_list = np.unique(noise_coeff_array).astype(int)

ind = 0
winner = np.zeros(len(nsites_list)*len(nmesh_list)*len(noise_coeff_list)*len(nsources_list), dtype=int)
pst = 0
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:          
                best_ffinal = ffinal_array1[ind]
                for ninit in ninit_list:
                    if ffinal_array1[ind] <= best_ffinal:
                        best_ffinal = ffinal_array1[ind]
                        winner[pst] = ind
                    if nsites_array[ind] == nsites:
                        ind+=1
                pst+=1

# print(winner)
# sites = np.unique(nsites_array)
# if len(sites) == 1:
#     print("Usar outros critérios")
# else:
#     winner = []
#     for s in sites:    
#         ffinal = []
#         for l in range(len(nsites_array)):
#             if nsites_array[l] == s:
#                 ffinal.append(ffinal_array1[l])
        
#         best_ffinal = np.min(np.array(ffinal))

#         for l in range(len(nsites_array)):
#             if nsites_array[l] == s and ffinal_array1[l] == best_ffinal:
#                 winner.append(l)
# print(winner)
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
\usepackage{placeins}

\graphicspath{{./figures/}}

\begin{document}      
    
\begin{table}
\centering
\resizebox{0.6\textheight}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c||c|c|c}
\hline
$\ck$ & $\bar \im$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ & $G(\bsi^0)$ & $G(\hat\bsi)$ & $\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time & Flag & $G_{\zeta_1=1}(\hat \bsi)$ \\ \hline \hline
""")
    for i in range(len(nsites_array)):
                
        if i in winner:
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
\end{table}
\FloatBarrier""")
            
            f.write(r"""
\begin{table}
\centering
\resizebox{0.5\textheight}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c||c|c|c}
\hline
$\ck$ & $\bar \im$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ & $G(\bsi^0)$ & $G(\hat\bsi)$ & $\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time & Flag & $G(\hat\bsi)$ \\ \hline
                    """)
                
    
    f.write(r"""
\end{tabular}}
\end{table}""")
    

    f.write(r"""

\begin{figure}
\begin{center}
\resizebox{0.6\textheight}{!}{
\begin{tabular}{|c|c|c|c|}
\hline
& Ground truth & Initialization & Reconstruction\\
\hline
\hline""")

    for i in range(len(nsites_array)):

        if i % 4 == 0:
            l = 0
            name = 'blsfig'+str((i//4))
        
        if i in winner:
            f.write(r"""\multirow{3}{*}{\rotatebox{90}{\textcolor{"""+color+r"""}{$\ck=""" + str(nsites_array[i]) + r"""\;$}}} & & & \textcolor{"""+color+r"""}{Noise = $"""+str('%4.2f' % noise_level_array[i])+r"""\%$}  \\
        \cline{4-4}
        & & \textcolor{"""+color+r"""}{$E(\bsi^0) = """+str('%5.2f' % errorinit_array[i])+r"""\% $} & \textcolor{"""+color+r"""}{$E(\hat\bsi) = """+str('%5.2f' % erroropt_array[i])+r"""\% $} \\""")
        else:
            f.write(r"""\multirow{3}{*}{\rotatebox{90}{$\ck=""" + str(nsites_array[i]) + r"""\;$}} & & & Noise = $"""+str('%4.2f' % noise_level_array[i])+r"""\%$  \\
        \cline{4-4}
        & & $E(\bsi^0) = """+str('%5.2f' % errorinit_array[i])+r"""\% $ & $E(\hat\bsi) = """+str('%5.2f' % erroropt_array[i])+r"""\% $ \\""")

        
        for j in range(3):
            f.write(r"""
            & \includegraphics[scale=0.8]{./figures/"""+name+str(alphabet[l])+r""".eps}""")
            l += 1
        f.write(r"""\\ \hline""")
        
        if (i+1)%4 == 0:
            f.write(r"""
        \end{tabular}}
        \end{center}
        \end{figure}
        \FloatBarrier
        """)
            f.write(r"""\begin{figure}
        \begin{center}
        \resizebox{0.6\textheight}{!}{
        \begin{tabular}{|c|c|c|c|}
        \hline
        & Ground truth & Initialization & Reconstruction\\
        \hline
        \hline""")
    
    f.write(r"""
        \end{tabular}}
        \end{center}
        \end{figure}
        """)

    
    f.write(r"""\end{document}""")