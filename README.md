If using git, clone the git repo, making sure to get the dependencies with it
(they are git submodules in vendor/), and enter the directory:
```
git clone https://github.com/rumajo/fdck.git --recurse-submodules
cd fdck
```

If using an archive, just extract it using the following (or GUI etc):
```
unzip fdck.zip
cd fdck
```

Create a build directory and enter it:
```
mkdir build
cd build
```

Run `cmake` from within build directory
```
cmake <flags> ..
```

`<flags>` contains extra flags that we pass to `cmake`. Without any flags, `cmake` will do a debug build by default.

For an optimised build use:
```
<flags> = -DCMAKE_BUILD_TYPE=Release
```

To cross-compile for Windows using mingw use (this requires the appropriate cross compilers to be installed):
```
<flags> = -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++
```

Note that link-time optimisation is disabled when cross-compiling since it seems to create executables that immediately crash.

Then run
```
make <options> <targets>
```

where `<targets>` contains the list of drivers etc. to build (just the library
and tests by default), and `<options>` usually contains just `-jN`, to compile
using N cpus. Add `VERBOSE=1` to `<options>` to see exactly which commands are
run by `make` (this is a `cmake` specific option).

For example, to build just the library and the `minimal_1_dimensional` driver
using 4 CPU cores, run
```
make -j4 minimal_1_dimensional
```
and the executable can be found in `build/drivers/minimal_1_dimensional`. (The
executables might also have `.exe` extensions on Windows.)

The drivers are statically linked by default, which means they require no other
library files etc. at runtime. This makes the executables pretty much
self-contained, except for any config files they might need.
