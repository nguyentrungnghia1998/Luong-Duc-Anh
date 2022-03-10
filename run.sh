#!/bin/bash
export OMP_NUM_THREADS=2
nohup ./LBM > out.LBM.txt &
