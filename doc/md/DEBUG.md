# Debugging tRIBS

tRIBS is an extensive program and not all options and variations have been tested in conjunction with newer compilers.
Consequently, debugging the program maybe necessary as tRIBS continues to evolve and new and old functionality is tested
under different circumstances. The text below
describes some approaches that maybe utilized to help debug potential issues.

## Address Sanitizer
The preferred method for debugging memory related issues is to use Address Sanitizer 
(ASan). ASan is bundled with newer versions of GCC and Clang and a full description of ASan is provided [here](https://github.com/google/sanitizers/wiki/AddressSanitizer). Note the default version of Clang in macOS does not support
memory leak detection, to enable this feature on macOS you will need to install llvm (i.e. clang) using brew and specify the path 
separately while compiling tRIBS. The option to compile with ASan and a specific compiler can be specified in CMakeList.txt.

To use flags in [ASan](https://github.com/google/sanitizers/wiki/AddressSanitizerFlags)
```ASAN_OPTIONS``` must be specified in the IDE or shell. For example, to detect memory leaks on macOS the following must be executed before running 
tRIBS compiled with ASan: 
```export ASAN_OPTIONS=detect_leaks=1```(note on linux this is the default). Additional options can be added as follows,
```export ASAN_OPTIONS=detect_leaks=1:option_2=1:option_3:...```. Note you must use ```export``` if you are using ASan in a shell.
Alternatively, if you are using and IDE these options can be specified according to the relevant documentation.

## Semantic Errors and Model Behavior
For  more complex bugs related to model behavior we suggest user to use an IDE with debug functionality. Or alternatively one can use [GDB](https://www.sourceware.org/gdb/).

