cd ../build
cmake ..
make
cd ./src
LOG_NUM_CLASS=7
NUM_CLASS=$(( 1 << LOG_NUM_CLASS ))

LOG_OMEGA_SIZE=10
NUM_OMEGA=$(( 1 << LOG_OMEGA_SIZE ))

B=5
A=$(( NUM_OMEGA - (NUM_CLASS - 1) * B ))
echo $A $B

./vrr $NUM_CLASS $NUM_OMEGA $A $B