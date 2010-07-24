#!/bin/bash

python ../../examples/plot_2d_distances.py dist2d.png

convert -scale 244 dist2d_D_A_ony.png dist2d_D_A_ony_small.png 
convert -scale 244 dist2d_D_A.png dist2d_D_A_small.png 
convert -scale 244 dist2d_D_M_ony.png dist2d_D_M_ony_small.png 
convert -scale 244 dist2d_D_M.png dist2d_D_M_small.png 

