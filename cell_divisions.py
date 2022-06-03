import random
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from functools import reduce
import numpy as np
np.seterr(all="ignore")

# Global Time variable
T = 36


class Cell(ABC):
    id = 1
    
    @abstractmethod
    def __init__(self):
        self.id = Cell.id
        Cell.id+=1
        self.age = 0

    def getAge(self):
        return self.age
    
    def getDivisionTime(self):
        return self.division_time
    
    '''
    Increment age by 1
    If age is >= division time, cell divides (return False)
    and indicates 2 new cells should be created
    '''
    def grow(self):
        self.age+=1
        return self.age <= self.getDivisionTime() or self.getDivisionTime() <= 0
    
    def __str__(self):
        return f"Cell_{self.id}:\t{self.type}\tAge={self.age}\tDivision_Time={self.division_time}"


class UncorrelatedCell(Cell):
    def __init__(self):
        super().__init__()
        self.type = "Uncorrelated"
        self.division_time = random.randint(2, 8)
    

class CorrelatedCell(Cell):
    def __init__(self, parent_div_time = 5):
        super().__init__()
        self.type = "Correlated"
        self.division_time = random.randint(parent_div_time - 3, parent_div_time + 3)
    

def uncorrelated_sim():
    P = [0 for _ in range(100)]
    num_potential_mothers = [0 for _ in range(100)]
    num_mothers = [0 for _ in range(100)]
    uncorrelated_cells = [UncorrelatedCell() for _ in range(5)]
    population = []
    total_num = 5
    
    for i in range(T):
        population.append(len(uncorrelated_cells)) 
        new_cells = []       
        for cell in uncorrelated_cells:
            num_potential_mothers[cell.getAge()]+=1
            P[cell.getAge()]+=1
            if(not cell.grow()):
                # cell divides  
                num_mothers[cell.getAge()-1]+=1
                del cell
                new_cells.append(UncorrelatedCell())
                new_cells.append(UncorrelatedCell())
                total_num+=2
        uncorrelated_cells+=new_cells
        print(f"Time {i}: Population: {population[i]}", end='\r')
    
    # clean up
    for cell in uncorrelated_cells:
        del cell
        
    plt.plot(range(1,len(population)+1), population, '+', label="Uncorrelated")
    m = (np.divide(np.array(num_mothers), np.array(num_potential_mothers))).tolist()
    return [(p/total_num) for p in P], m, population
 

def correlated_sim():
    P = [0 for _ in range(100)]
    num_potential_mothers = [0 for _ in range(100)]
    num_mothers = [0 for _ in range(100)]
    correlated_cells = [CorrelatedCell() for _ in range(5)]
    population = []
    total_num = 5
    
    for i in range(T):
        population.append(len(correlated_cells)) 
        new_cells = []       
        for cell in correlated_cells:
            num_potential_mothers[cell.getAge()]+=1
            P[cell.getAge()]+=1
            if(not cell.grow()):
                # cell divides
                num_mothers[cell.getAge()-1]+=1
                parent_div_time = cell.getDivisionTime()
                del cell
                new_cells.append(CorrelatedCell(parent_div_time))
                new_cells.append(CorrelatedCell(parent_div_time))
                total_num+=2
        correlated_cells+=new_cells
        print(f"Time {i}: Population: {population[i]}", end='\r')

    # clean up
    for cell in correlated_cells:
        del cell

    plt.plot(range(1,len(population)+1), population, 'x', label='Correlated')
    m = (np.true_divide(np.array(num_mothers), np.array(num_potential_mothers))).tolist()
    return [(p/total_num) for p in P], m, population


def euler_lotka_estimation(P, m):
    
    temp = [y for y in m if not np.isnan(y)]
    default = sum(temp)/len(temp)
    
    result = 0
    for i in range(len(P)):
        if(not np.isnan(P[i])) and (not np.isnan(m[i])):    
            if(m[i] <= 0):
                result+=(P[i]*default)
            else:
                result+=(P[i]*m[i])
    return result


if __name__ == "__main__":
    # Part 1: Simulation of Growth
    print("\033[1m")
    print("." * 50)
    print("Starting Simulation\033[0m\n")
    
    print("Uncorrelated:")
    P_uncorr, m_uncorr, pop_uncorr = uncorrelated_sim()
    print("\nCorrelated:")
    P_corr, m_corr, pop_corr = correlated_sim()
    
    plt.title("Cell Population Growth")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.legend()
    plt.savefig("cell_divisions.png")
    
    print("\n\n\033[1mSimulation Complete\n")
    print("." * 50)
    # Part 2: Comparison to Euler-Lotka Model predictions
    print("Starting Model Checking\033[0m\n")
    
    R0_uncorr = euler_lotka_estimation(P_uncorr, m_uncorr)
    R0_corr = euler_lotka_estimation(P_corr, m_corr)
    
    ratios = list(filter((1).__ne__, [pop_uncorr[i]/pop_uncorr[i-1] for i in range(1,len(pop_uncorr))])) 
    actual_uncorr = reduce(lambda a, b: a + b, ratios) / len(ratios)
    ratios = list(filter((1).__ne__, [pop_corr[i]/pop_corr[i-1] for i in range(1,len(pop_corr))]))   
    actual_corr = reduce(lambda a, b: a + b, ratios) / len(ratios)

    print(f"Uncorrelated Growth Rate\n\tEstimate:\t{R0_uncorr}\n\tActual:\t\t{actual_uncorr}")
    print(f"Correlated Growth Rate\n\tEstimate:\t{R0_corr}\n\tActual:\t\t{actual_corr}")
    
    print("\n\033[1mModel Checking Complete\n")
    print("." * 50)
    print("\033[0m")
