This repository was created as part of the KIT Modeling Week 2024, where we applied inverse uncertainty quantification to the simulation of a 1D shock tube.

## What's in the notebooks?

The notebooks containing `ex1` use the exact solver (Riemann.py from the AREPO public code release), while the notebooks containing `ex2` use the actual AREPO simulator. After downloading this repository and installing necessary Python packages, it should be possible to run the notebooks using the exact solver, while the simulator notebooks require an UM-Bride server, the code of which is not yet available in this repository, unfortunately.

The notebook names also indicate the value of `sigma`, i.e. the standard derivation of the measurement errors. Taking a look at the notebooks, one can compare how the inverse problem behaves for variations of `sigma`.

## Build docker container

```bash
docker build -t arepo-docker .
```

(TODO: add docker container to repository)

## Run docker container

```bash
$ docker run -p 4242:4242 --log-driver none -it arepo-docker
root@5434100aa413:/opt/arepo# python3 Shockwave1D-server.py
```

(We disable logging to avoid filling up the disk.)

## Run the notebooks

You can reproduce the results in the notebooks by running the cells from top to bottom.
Perhaps you will need to install some packages using pip.
If you just want to tweak the diagnostics, you can skip re-sampling the Markov chains. The chains are stored in the notebooks, so you can just run the cells starting from "Diagnostics".

Note: Some of the diagnostics are interactive and require `plotly`. It seems that they are not being rendered when previewing the notebooks on GitHub.
