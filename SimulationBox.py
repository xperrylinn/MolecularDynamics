import numpy as np
import numpy.random as rand

class SimulationBox:

    # Initialize simulation box specifications
    def __init__(self):
        self.size = 1000 # Angstrom, Volume = size^3
        self.numAtoms = 500 # Number of Argon atoms
        self.timeStep = 1 # fs, Time step of simulation in femto second
        self.atomPositions = [] # List to hold atoms positions. Each element is a tuple of size 3. One value for each coordinate
        self.atomVelocities = [] # List to hold atom velocities. Each element is a tuple of size 3. One Value for each coordinate

    # Initialize atom positions
    def initializePositions(self):
        for i in range(0, self.numAtoms):
            # Choose random x, y, and z coordinates for each atom
            randomCoordinates = rand.random_integers(-500, 500, 3)
            self.atomPositions.append(randomCoordinates)

    # Dump coordinates to file
    def dumpPositions(self):
        # Iterate through atom positions at write to file
        file = open("positions.xyz", "w+")
        file.write(str(self.numAtoms) + "\n")
        file.write("\n")
        # file.write("atom x y z\n")
        for i in range(0, self.numAtoms):
            currentCoordinates = self.atomPositions[i]
            currX = currentCoordinates[0]
            currY = currentCoordinates[1]
            currZ = currentCoordinates[2]
            line = " " + str(i) + " " + str(currX) + " " + str(currY) + " " + str(currZ) + '\n'
            file.write(line)
        file.close()

    # Create frames for a pseudo simulation
    def pseudoSimulation(self):
        for i in range(0, 500):
            self.pseudoSimulate()

    # Randomly adds -1, 0, or 1 to each atom's coordinates
    def pseudoSimulate(self):
        # Open the positions.xyz
        file = open("positions.xyz", "a")
        file.write(str(self.numAtoms) + "\n")
        file.write("\n")
        i = 0
        for atomPosition in self.atomPositions:
            # Generate three radom numbers from -1, 0, or 1 into a tuple
            randoms = rand.random_integers(-1, 1, 3)
            # Add to each of the current atoms coordinates
            atomPosition[0] += randoms[0]
            atomPosition[1] += randoms[1]
            atomPosition[2] += randoms[2]
            currX = atomPosition[0]
            currY = atomPosition[1]
            currZ = atomPosition[2]
            line = " " + str(i) + " " + str(currX) + " " + str(currY) + " " + str(currZ) + '\n'
            file.write(line)
            i += 1





    def main(self):
        self.initializePositions()
        self.dumpPositions()
        self.pseudoSimulation()

print("Hello World from SimulationBox.py")
helloWorld = SimulationBox().main()