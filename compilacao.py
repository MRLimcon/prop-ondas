import os

os.system(f'cd fortran_bins && f2py -c \
    --opt="-O3 -ftree-vectorize -march=native -fno-range-check -ffast-math \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops={os.cpu_count()}" -lgomp \
    calculations.f95 -m calculations')
os.system(f'cd fortran_bins && f2py -c \
    --opt="-O3 -ftree-vectorize -march=native -fno-range-check -ffast-math \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops={os.cpu_count()}" -lgomp \
    utils.f95 -m utils')