import argparse
import numpy as np
import pandas as pd
from skimage import io
import json
from plyfile import PlyData
import pyvista
import pyacvd
from pathlib import Path
import shutil
import os
import mmap
import re
import shell
from shell import ShellException

#-------------------------------------------------------------------------------
# Fixed parameters
#-------------------------------------------------------------------------------
SURFACE_AREA_CUTOFF = 70
N_ITER_FAIRING = 20
SPH_DEGREE = 12


#-------------------------------------------------------------------------------
# arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

# INPUT_DIR = Path("/path/to/gc-50_0.2_1_0.2/each/t0001")
INPUT_DIR = Path(args.input).resolve()
# OUTPUT_DIR = INPUT_DIR.parents[4] / "output"
OUTPUT_DIR = Path(args.output).resolve()

experiment = str(INPUT_DIR).split("/")[-1]
CSV_DIR = INPUT_DIR.parents[1] / "summary" / (experiment + ".csv")
CSV = pd.read_csv(CSV_DIR)
TMP_DIR = OUTPUT_DIR / "tmp" / experiment
os.makedirs(TMP_DIR, exist_ok=True)
LOG_DIR = OUTPUT_DIR / "log" / experiment
os.makedirs(LOG_DIR, exist_ok=True)

os.chdir(Path(__file__).parent)


#-------------------------------------------------------------------------------
# tif.gz -> bin
#-------------------------------------------------------------------------------
TIFF_DIR = TMP_DIR / "tiff"
os.makedirs(TIFF_DIR, exist_ok=True)
BIN_DIR = OUTPUT_DIR / "bin" / experiment
os.makedirs(BIN_DIR, exist_ok=True)

# decompress && filter small noise by surfacearea
is_garbage = (CSV.SurfaceArea < SURFACE_AREA_CUTOFF)
for i, f in enumerate(sorted(INPUT_DIR.glob("*.tif.gz"))):
    name = f.name.split(".")[0]
    if is_garbage[i]:
        with open(LOG_DIR / (name + ".log"), "w") as log:
            log.write(name + " below threshold of surface area " + str(SURFACE_AREA_CUTOFF))
        continue
    g = TIFF_DIR / (name + ".tif")
    shell.exec("gzip -dkc < {} > {}".format(f, g))

for tiff in TIFF_DIR.iterdir():
    im = io.imread(str(tiff))
    im = (im==False)

    # NOTE: original input uses left-handed coordinate system
    dim_z, dim_y, dim_x = [i+2 for i in im.shape]
    # we need to add zero-filled pixels to separate the object the border of the image
    out = np.zeros((dim_x + 2) * (dim_y + 2) * (dim_z + 2), dtype=np.uint8)
    # here, converting to right-handed system
    idz, idy, idx = im.nonzero()
    out[(idx + 1) + dim_x*(idy + 1) + dim_x*dim_y*(idz + 1)] = 1

    bname = "{}_x{}_y{}_z{}.bin".format(tiff.name.rsplit(".", 1)[0], dim_x, dim_y, dim_z)
    out.tofile(str(BIN_DIR / bname))

shutil.rmtree(TIFF_DIR)


#-------------------------------------------------------------------------------
# marching cube && bin -> ply
#-------------------------------------------------------------------------------
PLY_DIR0 = OUTPUT_DIR / "ply0" / experiment
os.makedirs(PLY_DIR0, exist_ok=True)
with open(INPUT_DIR.parents[2] / "info.json", "r") as f:
    data = json.load(f)
    z_interval = data["z_interval"]

for f in BIN_DIR.iterdir():
    name = re.compile("(.*)_x").search(f.name).group(1) + ".ply"
    x = re.compile("_x(\d*)_y").search(f.name).group(1)
    y = re.compile("_y(\d*)_z").search(f.name).group(1)
    z = re.compile("_z(\d*)").search(f.name).group(1)
    g = PLY_DIR0 / name

    # isovalue is to separete inside and outside of the object
    # any value which is 0 < x < 1 is fine, because input is binary file
    shell.exec("./bin/IsoSurfaceExtraction --flip --iso 0.95 --in {} --out {} --res {} {} {} --dim 1 1 {}"
        .format(f, g, x, y, z, z_interval))
    # below does not make any change
    # shell.exec("./bin/IsoSurfaceExtraction --flip --quadratic --iso 0.9 --in {} --out {} --res {} {} {} --dim 1 1 {}"
    #     .format(f, g, x, y, z, z_interval))


#-------------------------------------------------------------------------------
# extract mesh with genus = 0
#-------------------------------------------------------------------------------
PLY_DIR1 = OUTPUT_DIR / "ply1" / experiment
os.makedirs(PLY_DIR1, exist_ok=True)

for p in PLY_DIR0.iterdir():
    with open(p, "rb") as f:
        f.readline(); f.readline()
        line = f.readline()
        nV = int(line.split()[-1])
        f.readline(); f.readline(); f.readline()
        line = f.readline()
        nF = int(line.split()[-1])
        nE = 3 * nF // 2
    if ((2 - nV + nE - nF) // 2) == 0:
        shutil.move(str(p), str(PLY_DIR1 / p.name))


#-------------------------------------------------------------------------------
# normalizing by size
#-------------------------------------------------------------------------------
vol = CSV.Volume

PLY_DIR2 = OUTPUT_DIR / "ply2" / experiment
os.makedirs(PLY_DIR2, exist_ok=True)
for f in PLY_DIR1.iterdir():
    i = int(f.name[4:].split(".")[0]) - 1
    data = PlyData.read(f)
    data['vertex']['x'] -= np.mean(data['vertex']['x'])
    data['vertex']['y'] -= np.mean(data['vertex']['y'])
    data['vertex']['z'] -= np.mean(data['vertex']['z'])
    data['vertex']['x'] *= vol[i]**(-1/3) * 4 # 4 is for making the volume close to a unit sphere
    data['vertex']['y'] *= vol[i]**(-1/3) * 4
    data['vertex']['z'] *= vol[i]**(-1/3) * 4
    data.write(PLY_DIR2 / f.name)


# #-------------------------------------------------------------------------------
# # taubin smoothing
# #-------------------------------------------------------------------------------
# PLY_DIR3 = OUTPUT_DIR / "ply3" / experiment
# os.makedirs(PLY_DIR3, exist_ok=True)
#
# for f in PLY_DIR2.iterdir():
#     # name = f.name.rsplit(".", 1)[0]
#     # g = PLY_DIR3 / (name + ".ply")
#     g = PLY_DIR3 / f.name
#     try:
#         shell.exec(
#             "meshlabserver -i {} -o {} -s {}".format(f, g, "./mlx/taubin_lambda0.7_iter40.mlx"),
#             timeout_s = 60)
#     except ShellException as e:
#         with open(LOG_DIR / (name + ".log"), "w") as log:
#             log.write(str(e))


#-------------------------------------------------------------------------------
# uniform remeshing
#-------------------------------------------------------------------------------
PLY_DIR3 = OUTPUT_DIR / "ply3" / experiment
os.makedirs(PLY_DIR3, exist_ok=True)

for f in PLY_DIR2.iterdir():
    ply = pyvista.read(str(f))
    clus = pyacvd.Clustering(ply)
    clus.subdivide(4)
    clus.cluster(ply.number_of_points * 11 // 10)
    remesh = clus.create_mesh(flipnorm=False)
    remesh.save(PLY_DIR3 / f.name)

    # remove pyvtk header
    fd = os.open(PLY_DIR3 / f.name, os.O_RDWR)
    size = os.path.getsize(PLY_DIR3 / f.name)
    mm = mmap.mmap(fd, size)
    mm.readline(), mm.readline(), mm.readline()
    byte = ( "comment " + "s" * (mm.find(b'\n') - mm.tell() - 8) ).encode("latin-1")
    mm.write(byte)
    mm.flush()
    mm.close()


#-------------------------------------------------------------------------------
# post acvd fix
#-------------------------------------------------------------------------------
OBJ_DIR0 = OUTPUT_DIR / "obj0" / experiment
os.makedirs(OBJ_DIR0, exist_ok=True)

for f in PLY_DIR3.iterdir():
    name = f.name.rsplit(".", 1)[0]
    shell.exec("./bin/post_acvd --input {} --output {}".format(f, OBJ_DIR0 / (name + ".obj")))


#-------------------------------------------------------------------------------
# fairing
#-------------------------------------------------------------------------------
OBJ_DIR1 = OUTPUT_DIR / "obj1" / experiment
os.makedirs(OBJ_DIR1, exist_ok=True)

for f in OBJ_DIR0.iterdir():
    # HACK: modify file path because spinxFairingFast is actually inside the docker container
    input = "/hostroot" + str(f)
    output = "/hostroot" + str(OBJ_DIR1 / f.name)

    name = f.name.split(".")[0]
    shell.exec("./bin/spinxFairingFast {input} -0.95 {n_iter} {output}"
               " > /dev/null 2>&1" # silence the error log
               .format(input=input, output=output, n_iter=N_ITER_FAIRING))


#-------------------------------------------------------------------------------
# Sphere Harmonics Expansion
#-------------------------------------------------------------------------------
COEF_DIR = OUTPUT_DIR / "coef" / experiment
os.makedirs(COEF_DIR, exist_ok=True)

for f in OBJ_DIR1.iterdir():
    name = f.name.split("obj")[0][:-1]
    n = f.name.split("obj")[1][:-1] # times processed by fairing
    if n != "20":
        continue
    else:
        xyz = list(OBJ_DIR0.glob(name + ".obj"))[0]
        shell.exec("./bin/calc_coef --original {} --sphere {} --output {} --sph_degree {}"
            .format(xyz, f, COEF_DIR / (name + ".coef"), SPH_DEGREE),
            timeout_s = 60,
            log = LOG_DIR / (name + ".log"))


#-------------------------------------------------------------------------------
# Test: recovery of shape
#-------------------------------------------------------------------------------
OBJ_DIR2 = OUTPUT_DIR / "obj2" / experiment
os.makedirs(OBJ_DIR2, exist_ok=True)

for f in COEF_DIR.iterdir():
    name = f.name.split(".")[0]
    sphere = list(OBJ_DIR1.glob(name + ".obj20.obj"))[0]

    shell.exec("./bin/recv_shape --coef {} --sphere {} --output {}"
            .format(f, sphere, OBJ_DIR2 / (name + ".obj")))
