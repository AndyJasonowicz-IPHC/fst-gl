#
# NimScript build file for Project X
#

# Switches
switch("verbosity", "0")
switch("hints", "off")

when defined(testing) :
  switch("verbosity", "1")
  switch("hints", "on")

#
# Tasks
#
task tests, "run the test":
  exec "testament pat \"test*.nim\""

task build, "build project":
  ## Build new parallelized module
  exec "mkdir -p bin"
  exec "nimble install --depsOnly --verbose"
  exec "nim c -f --cc:clang --exceptions:setjmp --debugger:native --debuginfo:on --linedir:on --threads:on --exceptions:setjmp -d:useMalloc -d:openmp -d:release -d:danger -d:speed -o:bin/fst-gl  src/fstgl.nim"

task clean, "cleanup":
  exec "rm -rf nim_modules"
  exec "rm -rf bin"