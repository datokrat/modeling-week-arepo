TODO

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