make = make --no-print-directory

all: kazhdan fairing decomposition

kazhdan:
	$(make) -C ext/IsoSurfaceExtraction
	mv ext/IsoSurfaceExtraction/Bin/Linux/IsoSurfaceExtraction bin
	$(make) -C ext/IsoSurfaceExtraction clean

fairing:
	$(make) -C ext/spinxFairingFast
	cp ext/spinxFairingFast/spinxFairingFast bin

decomposition:
	$(make) -C cpp
	mv cpp/calc_coef bin
	mv cpp/recv_shape bin
