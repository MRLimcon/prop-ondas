import os

os.system(f"cd fortran_bins && f2py -c \
    --opt='-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -pipe' \
    calculations.f95 -m calculations")