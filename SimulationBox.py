import numpy.random as rand
import numpy as np

class SimulationBox:

    # Initialize simulation box specifications
    def __init__(self):
        self.size = 10.0 # Angstrom, Volume = size^3
        self.numAtoms = 1000 # Number of Argon atoms
        self.temperature = 1000 # Tempeature, T in Kelvin K
        self.molecularWeight = 39.948 # Molar mass of Argon, g / mol
        self.massPerAtom = 1 #np.float64(6.634e-26) # Kg/atom of Argon 
        self.boltzman = 1.38064852 * 10**(-23) # Boltzman constant, (m^2 kg) / (s^2 K)
        self.avogadrosNumber = 6.022 * 10**(23) # Avogadros number, atoms per mole
        self.timeStep = 0.01 #1*10**(-15) # s, Time step of simulation in femto second
        self.atomPositions = [] # List to hold atoms positions. Each element is a tuple of size 3. One value for each coordinate
        self.atomPreviousPositions = [] # Used in integrating equations of motion. 
        self.atomVelocities = None # List to hold atom velocities. Each element is a tuple of size 3. One Value for each coordinate
        self.atomForces = [] # List to hold atom forces at current time step

    # Initialize routine - positions, previous positions, velocities, scale, etc..
    def initializationRoutine(self):
        print("Hello World")
        self.initializePositions()
        self.initializeVelocities()
        self.initializePreviousPositions()

    # Initialize atom positions
    def initializePositions(self):
        for i in range(0, self.numAtoms):
            # Choose random x, y, and z coordinates for each atom
            randomCoordinates = rand.random_integers(0.0, self.size, 3)
            self.atomPositions.append(randomCoordinates)
    
    # Initialize atom previous positions
    def initializePreviousPositions(self):
        # previous position = position - vel * dT
        for i in range(0, self.numAtoms):
            x, y, z = self.atomPositions[i]
            vx, vy, vz = self.atomVelocities[i]
            xm = x - vx * self.timeStep 
            ym = y - vy * self.timeStep
            zm = z - vz * self.timeStep
            self.atomPreviousPositions.append((xm, ym, zm)) 

    # Initialize velocities
    #   Random velocities based on normal distribution, mean = 0 and variance = 1
    #   Scales velocities to corrent temperature
    def initializeVelocities(self):
        self.atomVelocities = rand.randn(self.numAtoms, 3) # Generate random velocities from normal distrbution, mean 0 variance 1
        initialTemp = self.instantaneousTemperature()
        scailingFactor = (self.temperature / initialTemp)**(1/2) # Calculate scailing factor
        self.atomVelocities *= scailingFactor # Scale atom velocitie

    # Calculates instantaneous temperature
    def instantaneousTemperature(self):
        squaredVelocities = self.atomVelocities**2 # Square velocities
        sumSquaredVelocities = squaredVelocities.sum(axis=1) # Sum squared velocities
        instantaneousTemp = (self.molecularWeight * self.avogadrosNumber**(-1) * sumSquaredVelocities.sum()) / (self.boltzman * 3 * self.numAtoms) # Calculate initial instantaneous temp
        print("Instantaneous temperature = ", instantaneousTemp)
        return instantaneousTemp

    # Calculates the velocity center of mass
    def velocityCenterOfMass(self):
        # Calculate vcom
        massPerAtom = self.molecularWeight * self.avogadrosNumber**(-1) # Units of grams, convert to kg?
        totalMass = self.numAtoms * massPerAtom
        vCOMX = vCOMY = vCOMZ = 0
        print("Calculating VCOM")
        for velocity in self.atomVelocities:
            vCOMX += velocity[0] * massPerAtom / totalMass
            vCOMY += velocity[1] * massPerAtom / totalMass
            vCOMZ += velocity[2] * massPerAtom / totalMass
        print("vCOMX = ", vCOMX, " vCOMY = ", vCOMY, " vCOMZ = ", vCOMZ)

    # Dump coordinates to file
    def dumpPositions(self):
        # Iterate through atom positions and write to file
        file = open("positions.xyz", "a")
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

    # Dump forces to file
    def dumpForces(self):
        # Iterate through atoms positions and 
        file = open("forces.txt", "w+")
        file.write("atom Fx Fy Fz\n")
        for i in range(0, self.numAtoms):
            currentForces = self.atomForces[i]
            currFx = currentForces[0]
            currFy = currentForces[1]
            currFz = currentForces[2]
            line = " " + str(i) + " " + str(currFx) + " " + str(currFy) + " " + str(currFz) + '\n'
            file.write(line)
        file.close()

    # Force calculation - calculate the forces on each atom at the current time step
    # andadds forces to self.atomForces 
    def forceCalculation(self):
        for i in range(0, self.numAtoms): # For each atom..
            fx = fy = fz = 0.0
            for j in range(0, self.numAtoms): # Calculate forces from all other atoms within cutoff
                if i != j: # Calculate force in each dirrection between i and j
                    if True: 
                    #if self.inRange(self.atomPositions[i], self.atomPositions[j]):
                        forceCompoents = self.calculateForce(self.atomPositions[i], self.atomPositions[j])
                        fx += forceCompoents[0]
                        fy += forceCompoents[1]
                        fz += forceCompoents[2]
                else: # Do nothing..
                    pass
            self.atomForces.append((fx, fy, fz)) 

    # Lennard-Jones Potential between two atoms
    def lennardJones(self, atomPosition1, atomPosition2):
        return 0

    # Integrate equation of motion
    def integrateEquationsOfMotion(self):
        # Calculate new position, according to current forces and previous position
        for i in range(0, self.numAtoms):
           x, y, z = self.atomPositions[i]
           xm, ym, zm = self.atomPreviousPositions[i]
           fx, fy, fz = self.atomForces[i]
           nx = 2 * x - xm + fx * self.timeStep**2 / self.massPerAtom 
           ny = 2 * y - ym + fy * self.timeStep**2 / self.massPerAtom
           nz = 2 * z - zm + fz * self.timeStep**2 / self.massPerAtom
           self.atomPreviousPositions[i] = (x % (self.size/2), y % (self.size/2), z % (self.size/2))
           self.atomPositions[i] = (nx % (self.size/2), ny % (self.size/2), nz % (self.size/2))

    # Force between two atoms given by derivative of Lennard-Jones
    def calculateForce(self, atomPosition1, atomPosition2):
        x = atomPosition1[0] - atomPosition2[0]
        y = atomPosition1[1] - atomPosition2[1]
        z = atomPosition1[2] - atomPosition2[2]
        r = (x**2 + y**2 + z**2)**(0.5)
        fx = 48.0 * x / r * (1/r**(13) - 1/(2*(r**(7))))
        fy = 48.0 * y / r * (1/r**(13) - 1/(2*(r**(7))))
        fz = 48.0 * z / r * (1/r**(13) - 1/(2*(r**(7))))
        return (fx, fy, fz)

    # Returns True if atoms are within the potential cut off Euclidean distance of each other,
    #   False if not. Take cutt to be 10 Angstrom.
    def inRange(self, atomPosition1, atomPosition2):
        euclideanDistance = ((atomPosition1[0] - atomPosition2[0])**2 +
                             (atomPosition1[1] - atomPosition2[1])**2 +
                             (atomPosition1[2] - atomPosition2[2])**2)**(0.5)
        if euclideanDistance <= 10.0:
            return True
        else:
            return False

    # Create frames for a pseudo simulation
    def pseudoSimulation(self):
        for i in range(0, 500):
            self.pseudoSimulate()

    # Simulation 
    def simulation(self):
        for i in range(0, 10):
            self.simulate()

    # At each time step update forces, atom previous position and current position
    def simulate(self):
        for i in range(0, self.numAtoms):
            self.forceCalculation()
            self.integrateEquationsOfMotion()
            self.dumpPositions()

        # Call forces; calculate forces
        # Call integrate; integrate equations of motion
        # increment time


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

    # Used for testing
    def main(self):
        self.initializationRoutine()
        self.simulation()
        # self.initializePositions()
        # self.initializeVelocities()
        #self.initializePreviousPositions()
        # self.dumpPositions()
        # self.pseudoSimulation()
        # self.velocityCenterOfMass()
        # self.instantaneousTemperature()
        # self.initializeVelocities()
        # self.forceCalculation()
        self.dumpForces()

# Run
print("Hello World from SimulationBox.py")
helloWorld = SimulationBox().main()
