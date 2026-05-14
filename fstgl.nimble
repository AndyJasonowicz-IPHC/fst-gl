# Package metadata
version = "0.1.0"
author = "Andy Jasonowicz"
description = "A program to estimate fst from genotype likelihoods"
license       = "TO DO"  # <-- ADD THIS LINE
srcDir        = "src"
bin           = @["fstgl"] # Compiles src/my_cli_tool.nim to a binary
namedBin = {"fstgl": "bin/fst-gl"}.toTable

# Dependencies
requires "nim >= 1.0.0"
requires "hts >= 0.3.22"
requires "zip"
requires "therapist == 0.3.0"


import std/os 

proc preBuild() =
  mkDir "bin"

## ---- Obtain and build dependencies HTSlib and libdeflate ----
proc buildDeps() =
  # --- Prerequisite Check ---
  for tool in ["autoreconf", "automake", "cmake", "make", "git"]:
    if findExe(tool) == "":
      # Using 'exec' forces the message out to the shell immediately
      exec "echo ''"
      exec "echo '[ERROR] Required build tool \"" & tool & "\" is missing.'"
      exec "echo 'Please install it to continue with the bundled HTSlib build.'"
      exec "echo ''"
      quit(1)

  let root = "static-build"
  mkDir root
  withDir root:
    let absRoot = getCurrentDir() 
    if not dirExists("libdeflate"):
      echo "--- Building libdeflate ---"
      exec "git clone https://github.com/ebiggers/libdeflate.git"
      withDir "libdeflate":
        exec "git checkout v1.25"
        exec "cmake -B build -DCMAKE_INSTALL_PREFIX=" & absRoot
        exec "cmake --build build"
        exec "cmake --install build"
    
    if not dirExists("htslib"):
      echo "--- Building htslib ---"
      exec "git clone https://github.com/samtools/htslib.git"
      withDir "htslib":
        exec "git checkout 1.23.1"
        exec "git submodule update --init --recursive"
        exec "autoreconf -i"
        exec "./configure --prefix=" & absRoot & " --disable-libcurl --with-libdeflate " &
             "CPPFLAGS=\"-I" & absRoot & "/include\" LDFLAGS=\"-L" & absRoot & "/lib64\""
        exec "make && make install"



task make, "Standard build":
  echo "--- Starting Standard Build ---"
  mkDir "bin"
  echo "Linking against system HTSlib..."
  exec "nim c -f -o:bin/fst-gl src/fstgl.nim"
  echo "Successfully built bin/fst-gl"

task bundle, "Developer build: Builds dependencies then fst-gl":
  buildDeps()
  mkDir "bin"
  exec "nim c -f -d:bundle -o:bin/fst-gl src/fstgl.nim"

task release, "Portable build":
  echo "--- Starting Portable Build ---"
  mkDir "bin"
  echo "Compiling with production optimizations and static flags..."
  exec "nim c -f -d:release -o:bin/fst-gl src/fstgl.nim"
  echo "Release build finished!"

task wipe, "Cleanup":
  echo "Cleaning up build artifacts..."
  rmDir "bin"
  rmDir "static-build"
  rmDir "nimcache"
  echo "Done."