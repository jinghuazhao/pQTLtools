\input{epsf}

\newcommand{\fig}[4]
{\begin{figure}\centerline{\epsfxsize=#1 \epsffile{#2}}
\caption{#3}\label{#4}
\end{figure}}

\newcommand{\ffig}[3]
{\begin{figure}\centerline{\epsfxsize=#1 \epsffile{#2}}\label{#3}
\end{figure}}

\newcommand{\fffig}[4]
{\begin{figure}[!bp]
\scalebox{#1}{\includegraphics{#2}}
\caption{#3}\label{#4}
\end{figure}}
