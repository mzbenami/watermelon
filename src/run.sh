#! /bin/bash
# run from src

MAP=${1:-map0}
LENGTH=${2:-50}
WIDTH=${3:-50}
GUI=${4:-true}
GROUP=${5:-group4}
S=${6:-.01}

javac watermelon/sim/watermelon.java
javac watermelon/$GROUP/Player.java

java -Xmx1g watermelon.sim.watermelon $MAP $LENGTH $WIDTH $GUI $GROUP $S
