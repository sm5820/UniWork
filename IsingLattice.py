import numpy as np
import scipy.constants as c

class IsingLattice:

    E = 0.0
    E2 = 0.0
    M = 0.0
    M2 = 0.0
    n_cycles = 0

    def __init__(self, n_rows, n_cols):
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.lattice = np.random.choice([-1,1], size=(n_rows, n_cols))
        self.N = np.exp(3.9140120246589096*np.log(self.n_cols) -2.17478240465023 + 0.5)
        # self.N is the cut-off. 

    def energy(self):
        "Return the total energy of the current lattice configuration."
        J = 1.0
        energy = 0.0
        energy += np.roll(self.lattice, 0, axis =1)*np.roll(self.lattice, 1, axis = 1)
        energy += np.roll(self.lattice, 0, axis =0)*np.roll(self.lattice, 1, axis = 0)
        energy_sum = -float(np.sum(energy))
        # for (i,j),spin in np.ndenumerate(self.lattice):
        #     i_1 = spin*self.lattice[i-1][j]
        #     j_1 = spin*self.lattice[i][j-1]
        #     energy_step = i_1 + j_1
        #     energy += -energy_step
    
        return energy_sum

    def magnetisation(self):
        "Return the total magnetisation of the current lattice configuration."
        # magnetisation = 0.0
        # for (i,j),spin in np.ndenumerate(self.lattice):
        #   magnetisation += spin
          #magnetisation is the sum of all spins in the lattice
        magnetisation = np.sum(self.lattice)*1.0
        return magnetisation

    def montecarlostep(self, T):
        'Performs a Monte Carlo step'      
        energy_0 = self.energy()#energy of initial state
        random_i = np.random.choice(range(0, self.n_rows))
        random_j = np.random.choice(range(0, self.n_cols))
        self.lattice[random_i][random_j] = -self.lattice[random_i][random_j] #random spin is flipped
        energy_1 = self.energy()#energy of new state
        delta_energy = energy_1 - energy_0#difference in energy 
        if delta_energy <= 0:
            energy = energy_1# new spin accepted
        if delta_energy > 0:
            random_number = np.random.random()#random number in range (0,1]
            exp_val = np.exp(-delta_energy/(T))#Boltzmann distribution of reduced units
            if random_number <= exp_val:
                energy = energy_1#new spin accepted
            if random_number> exp_val:
                self.lattice[random_i][random_j] = -self.lattice[random_i][random_j]
                energy = energy_0#new spin rejected
        
        self.n_cycles += 1
        if self.n_cycles> self.N: 
            self.E += energy 
            self.E2 += energy**2
            self.M += self.magnetisation()
            self.M2 += self.magnetisation()**2
        return energy, self.magnetisation()    
                            
    def statistics(self):
        # Calculates values for the averages of E, E*E (E2), M, M*M (M2), and returns them with Nsteps
        stat_E = self.E/(self.n_cycles-self.N)
        stat_E2 = self.E2/(self.n_cycles-self.N)
        stat_M = self.M/(self.n_cycles-self.N)
        stat_M2 = self.M2/(self.n_cycles-self.N)
        return stat_E, stat_E2, stat_M, stat_M2, (self.n_cycles-self.N)