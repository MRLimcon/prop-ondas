import os

os.system(f"cd fortran_bins && f2py -c \
    --opt='-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops=6' -lgomp \
    calculations.f95 -m calculations")
os.system(f"cd fortran_bins && f2py -c \
    --opt='-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops=6' -lgomp \
    utils.f95 -m utils")