import os

os.system(
    f'cd fortran_bins && f2py -c \
    --opt="-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops={os.cpu_count()}" -lgomp \
    calculations.f95 -m calculations'
)
os.system(
    f'cd fortran_bins && f2py -c \
    --opt="-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -pipe -ftree-parallelize-loops={os.cpu_count()}" -lgomp \
    utils.f95 -m utils'
)

"""
os.system(f'cd fortran_bins && NPY_DISTUTILS_APPEND_FLAGS=1 LDFLAGS=-acc f2py -c --fcompiler=nv\
    --opt="-stdpar=gpu -fast -O3 -Mnouniform -Mconcur " -lgomp\
    calculations.f95 -m calculations')
os.system(f'cd fortran_bins && NPY_DISTUTILS_APPEND_FLAGS=1 LDFLAGS=-acc f2py -c  --fcompiler=nv\
    --opt="-stdpar=gpu -fast -O3 -Mnouniform -Mconcur " -lgomp\
    utils.f95 -m utils')

"""
