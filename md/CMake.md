# CMake 

## instructions for compiling tRIBS on your machine using CMake

Note: these instructions are for using CMake via terminal, there is additional documentation [here](https://cmake.org/cmake/help/latest/guide/user-interaction/index.html#guide:User%20Interaction%20Guide) for using the CMake gui. 

1) Use [Homebrew](https://formulae.brew.sh/formula/cmake)to install CMake, alternatively you can download [CMake](https://cmake.org/download/), but Homebrew is preferred as it will catch additional dependencies.

2) You can check to see if CMake is on your path, by typing cmake into the command line. If it says its not found then you will need to set cmake to your path. For example if you downloaded CMake and its now in your application folder you can use:

```bash
sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install
```
3) Next change directory to the tRIBS source code, should look something like this but is depenendent on where the tRIBS source code is located:

```bash
cd ~/Documents/tRIBS
```
4) Make a build directory and change into that directory.

```bash
mkdir build
cd build
```

5) The run CMake (which generates a Makefile) followed by make to create the tRIBS executable

```bash
cmake ../
make
```

6) After you can check to see that the executable was made by using.
```bash
 ls 
 ```
 The executable will have a name specifed in the CMakeList.txt file

 ## Content of CMakeFile.txt

 In some instance you may want to modify the CMakeFile.txt, for example if you want to change the name of the executable, or change compilation parallel mode to serial, or add additional compiler flags. This section will be updated with more detail to demonstrate how this is possible here, but the CMakeFile.txt is documented with where these changes can be made.


## Return to [README](../README.md)