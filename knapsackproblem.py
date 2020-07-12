import cirq
import numpy as np
import scipy
import operator

# Input parameters. Feel free to change them.
L = 100                                         # Limit value
objs = [60, 35, 23, 55, 70, 13]                 # weights of the items
capas = 5                                       # number of layers
initial = 0                                     # inizialization value
methods = ['COBYLA','SLSQP','BFGS']             # methods that you want to use
# See https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html#scipy.optimize.minimize for available methods

n = len(objs) 
gamma = 2 / (n*L)  

class ParametrizedGate(cirq.Gate):
    
    def _decompose_(self, qs):        
        for i in range(n):
            Rz = cirq.Z(qs[n]) ** (gamma * objs[i])
            yield Rz.controlled_by(qs[i])

    def _num_qubits_(self):
        return n + 1
    
    def _unitary_(self):
        return cirq.unitary(
            cirq.Circuit(self._decompose_(cirq.LineQubit.range(n + 1)))) 
    
def ansatz(params):
    
    qs = cirq.LineQubit.range(n + 2)
    ansatz = cirq.Circuit()
    
    for j in range(capas):
        
        ansatz.append([cirq.X(qs[i]) ** params[i + 2*j*n] for i in range(n)])
        ansatz.append([cirq.CX(qs[aux[0]], qs[aux[1]]) 
                       for aux in [[i,i+1] for i in range(n - 1)]])
        ansatz.append([cirq.Z(qs[i]) ** params[i + (2*j+1)*n] for i in range(n)])
        ansatz.append([cirq.CX(qs[aux[0]], qs[aux[1]]) 
                       for aux in [[i,i+1] for i in range(n - 1)]])

    return ansatz

def expected_values(params):
    
    sol = []
    
    for i in range(2):

        qs = cirq.LineQubit.range(n + 2)
        circuit = cirq.Circuit()
        circuit.append(cirq.X(qs[n]))
        cg = cirq.ControlledGate(ParametrizedGate())
        circuit.append(cirq.H(qs[n + 1]))
        circuit.append(cirq.S(qs[n + 1]) ** i)
        circuit.append(cg(qs[n + 1], *[qs[i] for i in range(n + 1)])) 
        circuit.append(cirq.H(qs[n + 1]))

        circuit.append(cirq.measure(qs[n + 1]))
        circuit = ansatz(params) + circuit
        s = cirq.Simulator()
        rep = 1000
        sol.append(s.run(circuit, repetitions = rep))

    real = 2 * int(str(sol[0]).count('0')) / rep - 1
    img = - 2 * int(str(sol[1]).count('0')) / rep + 1
    
    if np.sin(L * gamma * np.pi) - img < 0:
        return 3
    else:
        return np.abs(np.cos(L * gamma * np.pi) - real) 
    
    
def solution(method,n):
    
    initial_params = np.array([initial]*(2*capas*n))
    minimum = scipy.optimize.minimize(expected_values,
                                      initial_params, method=method)
    final = cirq.Circuit()
    qs = cirq.LineQubit.range(n + 2)
    final =ansatz(minimum.x)
    final.append(cirq.measure(*[qs[i] for i in range(n)] , key = 'm'))

    s = cirq.Simulator()
    rep = 10000
    sol = s.run(final, repetitions = rep)

    result = max(sol.histogram(key = 'm').items(),
                   key=operator.itemgetter(1))[0]
    times = max(sol.histogram(key = 'm').items(), 
                key=operator.itemgetter(1))[1]

    sum = 0
    print("to approach the value ", L, " without going over, we must take: ")
    for i, n in enumerate(np.binary_repr(result,n)):
        sum += int(n) * objs[i]
        if int(n) == 1:
            print("obj ", i + 1, " with value ", objs[i])

    print("METHOD: ",  method)
    print("sum: ",sum, "percentage: ", 100* times / rep , "%")

    

for method in methods:
    try:
        solution(method,n)
    except Exception as e:
        print(method, " does not work")
        print('Exception: ' + str(e)) # Printing exact error for debugging if necessary
