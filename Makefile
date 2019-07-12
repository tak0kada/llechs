make = make --no-print-directory

all: kazhdan fairing decomposition

kazhdan:
	$(make) -C thirdparty/IsoSurfaceExtraction
	mv thirdparty/IsoSurfaceExtraction/Bin/Linux/IsoSurfaceExtraction bin
	$(make) -C thirdparty/IsoSurfaceExtraction clean

fairing:
	$(make) -C thirdparty/spinxFairingFast
	cp thirdparty/spinxFairingFast/spinxFairingFast bin

decomposition:
	$(make) -C cpp
	mv cpp/calc_coef bin
	mv cpp/recv_shape bin
