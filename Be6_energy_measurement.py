
from VQE.Circuit import Circuits_Composser

thetas = [0.5912482290924049, -0.7853862664802604]      
op_used = [[0, 3, 4, 5], [0, 3, 1, 2]]

composer = Circuits_Composser(nucleus='Be6',
                             n_qubits = 6,
                             ref_state = 0,
                             parameters = thetas,
                             operators_used = op_used,
                             nshots = 1000) 


Et_shots = composer.Qibo_measure_Energy()

