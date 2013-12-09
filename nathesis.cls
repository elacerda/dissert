% ** Natalia@woodstock - 27/Jan/2006
% ** Load hepthesis.cls with my hacks, plus some
%    other useful tex packages, such as natbib.
\ProvidesClass{nathesis}

%% Standard printing:
\LoadClass[a4paper,oneside]{hepthesis}

%%%% Tips from hepthesis class:
%% For normal draft builds
%\LoadClass[hyperpdf,nobind,draft]{nathesis}
%% For short draft builds (breaks citations by necessity)
%\LoadClass[hyperpdf,nobind,draft,hidefrontback]{hepthesis}
%% For Cambridge soft-bound version
%\LoadClass[hyperpdf,bindnopdf]{hepthesis}
%% For Cambridge hard-bound version (must be one-sided)
%\LoadClass[hyperpdf,oneside]{hepthesis}

% Figure packages
\usepackage{graphicx}
\usepackage{subfigure}
\usepackage[section]{placeins}

% Table package
\usepackage{longtable}

% The packages below seem to make a pretty good pdf using ps2pdf
\usepackage{ae,aecompl}

% Bibtex and Portuguese packages - cannot change this order!
\usepackage{natbib}

% Fix comma spacing in math mode.
\def\magyarOptions{mathhucomma=fix,figurecaptions=unchanged,tablecaptions=unchanged,sectiondot=none}
\usepackage[magyar,brazil]{babel}

\citestyle{aa}

% My definitions
\usepackage{nadefinitions}

%%%% André %%%%
\usepackage{hyperref}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}

% Used for code formatting.
\usepackage{fancyvrb}
\usepackage[scaled]{beramono}

\usepackage{txfonts}

