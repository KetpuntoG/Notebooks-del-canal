from qiskit import *

class Termino_hamiltoniano():
    
    def __init__(self, paulis, n_qubits, constante = 1):
        
        self.paulis = paulis
        self.n_qubits = n_qubits
        self.constante = constante
        
    def simplificar_termino(self):
        
        simplificacion = []
        qubits_recorridos = []
        
        for pauli, qubit in self.paulis:

            if qubit not in qubits_recorridos:
                
                qubits_recorridos.append(qubit)
                simplificacion.append([pauli, qubit])

            else:
                
                pauli_actual = simplificacion[qubits_recorridos.index(qubit)][0]
                nueva_pauli, constante = self.reducir_paulis(pauli_actual, pauli)
                
                if nueva_pauli == 'I':
                
                    simplificacion.remove([pauli_actual, qubit])
                    qubits_recorridos.remove(qubit)
                    
                else:
                    
                    simplificacion[qubits_recorridos.index(qubit)][0] = nueva_pauli
                    self.constante *= constante
                    
        return simplificacion
    
    def reducir_paulis(self, pauli1, pauli2):
        
        paulis = ['X', 'Y', 'Z']
        
        if paulis.index(pauli1) < paulis.index(pauli2):
            return set(paulis).difference(set([pauli1, pauli2])).pop(), complex('j')
        
        elif paulis.index(pauli1) > paulis.index(pauli2):
            return set(paulis).difference(set([pauli1, pauli2])).pop(), complex('-j')
        
        else:
            return 'I', 1
        
    def calcular_valor_esperado(self, circ_estado):
        
        puerta_estado = circ_estado.to_gate()
        
        termino_simp = self.simplificar_termino()
        qubits_interes = len(termino_simp)
        
        circuito = QuantumCircuit(self.n_qubits, qubits_interes)
        circuito.append(puerta_estado, range(self.n_qubits))
        
        for i, termino in enumerate(termino_simp):
            
            pauli, qubit = termino
            
            if pauli == 'X':
                circuito.h(qubit)
                
                
            elif pauli == 'Y':
                circuito.sdg(i)
                circuito.h(i)
            
            circuito.measure(qubit, i)
            
        conteo = self.ejecutar_circuito(circuito)
        solucion = 0
        for num in conteo:
            solucion += conteo[num] * (-1) ** (num.count('1'))
        return solucion * self.constante
                
    def ejecutar_circuito(self, circ):
        
        # devuelve la probabilidad de cada estado
        
        shots = 5000
        backend = Aer.get_backend("qasm_simulator")
        job = execute(circ, backend, shots = shots)
        counts = job.result().get_counts()
        
        for clave in counts:
            counts[clave] /= shots
            
        return counts


class Hamiltoniano():
    
    def __init__(self, n_qubits):
        
        self.terminos = []
        self.n_qubits = n_qubits
        
    def incluir_termino(self, termino):
        self.terminos.append(termino)
        
    def calcular_valor_esperado(self, circuito_estado):
        
        solucion = 0
        
        for termino in self.terminos:
            
            solucion += termino.calcular_valor_esperado(circuito_estado)
            
        return solucion
