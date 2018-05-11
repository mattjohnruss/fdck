Create a build directory and enter it:
```
mkdir build
cd build
```

Run `cmake` from within build directory
```
cmake <flags> ..
```

`<flags>` contains extra flags that we pass to `cmake`.

For an optimised build use:
```
<flags> = -DCMAKE_BUILD_TYPE=Release
```

To cross-compile for Windows using mingw use:
```
<flags> = -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_C_COMPILER=x86_64-w64-mingw32-gcc -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++
```

Note that link-time optimisation is disabled when cross-compiling since it seems to create executables that immediately crash.

Then run
```
make <options> <targets>
```

where `<targets>` contains the list of drivers etc. to build (all by default), and `<options>` usually contains just `-jN`, to compile using N cpus. Add `VERBOSE=1` to `<options>` to see exactly which commands are run by `make` (this is a `cmake` specific option).