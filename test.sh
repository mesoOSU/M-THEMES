cp src/forge test_run/
cd test_run
../mpich/bin/mpiexec -np 4 ./forge input.in
