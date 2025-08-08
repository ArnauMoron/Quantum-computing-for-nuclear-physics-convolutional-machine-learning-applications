from VQE.Circuit import Circuits_Composser


thetas = [-0.7853917646646815]      
op_used = [[0, 3, 1, 2]]

composer = Circuits_Composser(nucleus='Be6_red',
                             n_qubits = 4,
                             ref_state = 0,
                             parameters = thetas,
                             operators_used = op_used,
                             nshots = 1000) 


composer.Connect_platform()

Et_shots = composer.Qibo_measure_Energy()

composer.Disconnect_platform()
