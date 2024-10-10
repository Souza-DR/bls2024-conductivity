import numpy as np
import sys
from pdb import set_trace

nsites_list = [5]
ninit_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
nsources_list = [1, 3]
nmesh_list = [16, 32]
noise_coeff_list = [0.005, 0.01]

winner = []
for nsites in nsites_list:
    for nsources in nsources_list:
        for noise_coeff in noise_coeff_list:
            for nmesh in nmesh_list:
                ffinal1_list = []
                for ninit in ninit_list:
                    input_files = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
                    loaded_data = np.load("./results/data/"+input_files+".npz")
                    ffinal1_list.append(loaded_data["ffinal1"])
                winner.append(np.argmin(np.array(ffinal1_list))+1)

ndec = 2
color = "blue"

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
    pst = 0
    counter = 0
    for nsites in nsites_list:
        for nsources in nsources_list:
            for noise_coeff in noise_coeff_list:
                for nmesh in nmesh_list:
                    for ninit in ninit_list:
                        
                        input_files = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                        loaded_data = np.load("./results/data/"+input_files+".npz")

                        ffinal1 = loaded_data["ffinal1"]
                        nsites = (loaded_data["nsites"]).astype(int)
                        nsources = (loaded_data["nsources"]).astype(int)
                        nmesh = loaded_data["nmesh"]
                        noise_coeff = loaded_data["noise_coeff"]
                        noise_level = [100]*loaded_data["noise_level"]
                        finit = loaded_data["finit"]
                        ffinal = loaded_data["ffinal"]
                        normgpinit = loaded_data["normgpinit"]
                        normgpfinal  = loaded_data["normgpfinal"]
                        erroropt  = [100]*loaded_data["erroropt"]
                        errorinit  = [100]*loaded_data["errorinit"]
                        flagsol = (loaded_data["flagsol"]).astype(int)
                        iter = (loaded_data["iter"]).astype(int)
                        numevalf = (loaded_data["numevalf"]).astype(int)
                        CPU_time = loaded_data["CPU_time"]

                        if ninit == winner[pst]:
                            f.write(
str(nsites) + r""" & """+str(nsources) + r""" & """+str('%4.2f' % noise_level)+r"""& """+str('%5.2f' % errorinit)+r""" & """+str('%5.2f' % erroropt)+r""" & """+str('%4.2f' % finit)+r""" & """+str('%7.5f' % ffinal)+r""" & """+str('%7.5f' % normgpinit)+r""" & """+str('%7.5f' % normgpfinal)+r""" & """+str(iter)+r""" & """+str(numevalf)+r""" & """+str('%7.2f' % CPU_time)+r""" & """+str(flagsol) + r""" & \textbf{"""+str('%10.3e' % ffinal1)+r"""}\\ \hline
""")
                        else:
                            f.write(
str(nsites) + r""" & """+str(nsources) + r""" & """+str('%4.2f' % noise_level)+r"""& """+str('%5.2f' % errorinit)+r""" & """+str('%5.2f' % erroropt)+r""" & """+str('%4.2f' % finit)+r""" & """+str('%7.5f' % ffinal)+r""" & """+str('%7.5f' % normgpinit)+r""" & """+str('%7.5f' % normgpfinal)+r""" & """+str(iter)+r""" & """+str(numevalf)+r""" & """+str('%7.2f' % CPU_time)+r""" & """+str(flagsol) + r""" & """+str('%10.3e' % ffinal1)+r"""\\ \hline
""")
                        
                        if ninit == ninit_list[-1]:
                            f.write(r"""\hline
""")
                        
                        if (counter == (len(nsites_list)*len(nmesh_list)*len(noise_coeff_list)*len(nsources_list)*len(ninit_list)) - 1) or ((counter+1)%40) == 0:
                            f.write(r"""
\end{tabular}}
\end{table}
\FloatBarrier 
""")
                            if counter < (len(nsites_list)*len(nmesh_list)*len(noise_coeff_list)*len(nsources_list)*len(ninit_list) - 1):
                                f.write(r"""
\begin{table}
\centering
\resizebox{0.6\textheight}{!}{
\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c||c|c|c}
\hline
$\ck$ & $\bar \im$ & Noise & $E(\bsi^0)$ &$E(\hat\bsi)$ & $G(\bsi^0)$ & $G(\hat\bsi)$ & $\| \nabla G(\bsi^0)\|_{2}$ &$\|\nabla G(\hat\bsi)\|_{2}$ & $\#$iter & $\#$feval & Time & Flag & $G_{\zeta_1=1}(\hat \bsi)$ \\ \hline
""")
                        counter += 1
                    pst += 1                                   
    
    f.write(r"""

\begin{figure}
\begin{center}
\resizebox{0.6\textheight}{!}{
\begin{tabular}{|c|c|c|c|}
\hline
& Ground truth & Initialization & Reconstruction\\
\hline
\hline""")

    pst = 0
    counter = 0
    for nsites in nsites_list:
        for nsources in nsources_list:
            for noise_coeff in noise_coeff_list:
                for nmesh in nmesh_list:
                    for ninit in ninit_list:
                        
                        input_files = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                        loaded_data = np.load("./results/data/"+input_files+".npz")

                        ffinal1 = loaded_data["ffinal1"]
                        nsites = (loaded_data["nsites"]).astype(int)
                        nsources = (loaded_data["nsources"]).astype(int)
                        nmesh = loaded_data["nmesh"]
                        noise_coeff = loaded_data["noise_coeff"]
                        noise_level = 100*loaded_data["noise_level"]
                        finit = loaded_data["finit"]
                        ffinal = loaded_data["ffinal"]
                        normgpinit = loaded_data["normgpinit"]
                        normgpfinal  = loaded_data["normgpfinal"]
                        erroropt  = 100*loaded_data["erroropt"]
                        errorinit  = 100*loaded_data["errorinit"]
                        flagsol = (loaded_data["flagsol"]).astype(int)
                        iter = (loaded_data["iter"]).astype(int)
                        numevalf = (loaded_data["numevalf"]).astype(int)
                        CPU_time = loaded_data["CPU_time"]

                        if counter % 4 == 0:
                            l = 0
                            name = 'blsfig'+str((counter//4))
                        
                        if ninit == winner[pst]:

                            f.write(r"""
\multirow{3}{*}{\rotatebox{90}{\textcolor{"""+color+r"""}{$\ck=""" + str(nsites) + r"""\;$}}} & & & \textcolor{"""+color+r"""}{Noise = $"""+str('%4.2f' % noise_level)+r"""\%$}  \\
\cline{4-4} & & \textcolor{"""+color+r"""}{$E(\bsi^0) = """+str('%5.2f' % errorinit)+r"""\% $} & \textcolor{"""+color+r"""}{$E(\hat\bsi) = """+str('%5.2f' % erroropt)+r"""\% $} \\""")
                            
                        else:
                            f.write(r"""
\multirow{3}{*}{\rotatebox{90}{$\ck=""" + str(nsites) + r"""\;$}} & & & Noise = $"""+str('%4.2f' % noise_level)+r"""\%$  \\
\cline{4-4}
& & $E(\bsi^0) = """+str('%5.2f' % errorinit)+r"""\% $ & $E(\hat\bsi) = """+str('%5.2f' % erroropt)+r"""\% $ \\""")

                        
                        for j in range(3):
                            f.write(r"""
& \includegraphics[scale=0.8]{./figures/"""+name+str(alphabet[l])+r""".eps}""")
                            l += 1
                        f.write(r"""
\\ \hline""")
                        
                        if (counter == (len(nsites_list)*len(nmesh_list)*len(noise_coeff_list)*len(nsources_list)*len(ninit_list) - 1)) or ((counter+1)%4 == 0):
                            f.write(r"""
\end{tabular}}
\end{center}
\end{figure}
\FloatBarrier
                        """)
                            if counter < (len(nsites_list)*len(nmesh_list)*len(noise_coeff_list)*len(nsources_list)*len(ninit_list) - 1):
                                f.write(r"""
\begin{figure}
\begin{center}
\resizebox{0.6\textheight}{!}{
\begin{tabular}{|c|c|c|c|}
\hline
& Ground truth & Initialization & Reconstruction\\
\hline
\hline""")
                        counter += 1
                    pst += 1            

    f.write(r"""
\end{document}""")