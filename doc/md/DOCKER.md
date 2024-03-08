# tRIBS & Docker
[Docker](https://www.docker.com/) is a containerization platform that enables developers to package applications and their dependencies into portable containers. These containers effectively run like a virtual machine and can be utilized across different computing environments. Docker images of tRIBS (both parallel and serial versions) and MeshBuilder are maintained at [Docker Hub](https://hub.docker.com/).
## Getting Started with Docker
To get started with Docker, first, install Docker on your system by downloading the appropriate version for your operating system from the official Docker [website](https://www.docker.com/products/docker-desktop/). Docker can be ran both through the desktop client or via command line. In this tutorial we only focus on pulling and running tRIBS and or MeshBuilder via the command line. Additional information about running Docker can be found [here](https://docs.docker.com/).
## Pulling and Running tRIBS and MeshBuider Docker Images
Once docker has successfully been installed you can pull either image from Docker Hub. We walk through both case here highlighting minor differences in this processes.
### tRIBS
From the command line execute the following:
```
docker pull tribs/tribs:latest
```
You can check to see if the image is now available locally by:
```
docker image
```
From here the tRIBS image can be ran using `docker run -it tribs/tribs:latest`, where the `-it` flag creates an interactive session from the command line. However, in most cases, in order to successfully run tRIBS through the docker image you will need to be able to mount a local volume where data to run tRIBS is located. Note: it is also possible to download data into the image. This can be accomplished by using `-v` flag, where the local directory is mapped to a directory in the image following this structure __path/in/local:path/in/image__. For example, one could run: 
```
docker run -it -v /local/path/to/data:/tribs/data tribs/tribs:latest
```
In this case `/tribs/data` represents the directory in the image to access your data stored locally. From here one can execute tRIBS normally. There are two tRIBS binaries stored in `/tribs/bin`, a serial version, `tRIBS` and a version for parallel simulation `tRIBS_par`.

### MeshBuilder
Following the above steps one can obtain a Docker image of MeshBuilder as follows:
```
docker pull tribs/meshbuilder:latest
```
MeshBuilder also requires data that should be mounted to the image. The MeshBuilder image contains a shell script `meshbuild_workflow.sh`. This script runs through the various [steps](https://github.com/tribshms/MeshBuilder) required to successfully run MeshBuilder, but is predicated on the mounted volume containing:
1) A points file for the TIN mesh
2) A .in file, where the POINTFILENAME refers to the point file in the volume, and OUTFILENAME is the same as the base name for the model simulation.

To run the MeshBuilder image:
```
docker run -it -v /local/path/to/data:/meshbuild/data tribs/meshbuilder:latest
```
This will open up the command line for the image. To initiate the MeshBuilder workflow run `./meshbuild_workflow.sh`. This will prompt the user to enter the following information:
1) Name of the .in file stored in the volume
2) Number of compute nodes to partition the mesh on
3) Method for partitioning
4) And the model base name.

These steps can be recreated manually with a more detailed description of the processes [here](https://github.com/tribshms/MeshBuilder).

Return to [README](../../README.md)
