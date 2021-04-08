# heat
Author: Dominik A. Kardell

Publication: Kardell, D.A., Zhao, Z., Ramos, E.J., Estep, J., Christeson, G.L., Reece, R.S., Hesse, M.A., 2021. Hydrothermal models constrained by fine-scale seismic velocities confirm hydrothermal cooling of 7-63 Ma South Atlantic crst. Journal of Geophysical Research - Solid Earth.

Code used in Kardell et al. (JGR, 2021) to model hydrothermal fluid flow along the CREST transect in the western South Atlantic

The main directory contains all scripts and files called by the hydrothermal scripts in the line folders.

Folders line1AB, line1BC, line1CD, line1DE, and line1EF each contain:
1) A binary file containing the seismic velocity distribution from Full-Waveform Inversion
2) A conversion script that estimates physical property distributions from seismic velocities
3) A .mat file including the estimated physical properties (can be re-created by running the conversion script)
4) A convection script that runs the coupled fluid and heat flow simulation
5) A movie file showing the thermal distribution through time (can be re-created using the convection script)

If you have any questions, either contact me directly through GitHub or send me an email at dkardell@utexas.edu with the subject line "Kardell et al 2021 JGR model question".
