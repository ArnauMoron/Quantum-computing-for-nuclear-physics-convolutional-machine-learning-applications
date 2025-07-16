# Nuclear Energy Calculation (Be6)

This project calculates the ground state energy of the Beryllium-6 (`Be6`) atomic nucleus. The program builds and simulates the necessary quantum circuits using the **Qibo**, **OpenFermion**, and **Qililab** libraries.

---



##  File Structure

The program is divided into three collaborating Python files.

### `Be6_energy_measurement.py`

This is the **main script** and entry point of the program.
* **Defines parameters**:
    * `thetas`: A fixed set of rotation angles for the ansatz gates.
    * `op_used`: The list of operators that define the structure of the ansatz circuit.
* **Initiates the calculation**: It creates an instance of `Circuits_Composser` and calls the `Qibo_measure_Energy()` method to measure the energy.

### `Circuit.py`

This is the **engine of the program**. It contains the logic for building and executing the quantum circuits.
* **`Circuits_Composser` Class**:
* **Energy Measurement (`Qibo_measure_Energy`)**: Creates all the neceessary circuits (9 in total) to measure the expectation value of diferent terma in the Hamiltonian and sums them to get the total energy. 

### `Nucleus.py`

This file acts as the **data layer** for the physics problem.
* **`Nucleus` Class**: Reads pre-calculated data for the `Be6` nucleus from the `nuclei/Be6_data/` directory.
* **Functionality**: Loads the Hamiltonian, the basis states, and the one- and two-body excitation operators, translating the nuclear physics problem into a format the rest of the program can use.

---

##  Execution Flow

1.  The `Be6_energy_measurement.py` script is executed.
2.  The script creates a `Circuits_Composser` object, passing it the `thetas` and `op_used` parameters.
3.  The `Circuits_Composser` initializes a `Nucleus` object to load all the physical data for `Be6`.
4.  The `Qibo_measure_Energy()` method is called.
5.  This method generates a list of all necessary circuits: the main ansatz circuit and several modified versions for measuring each term of the Hamiltonian.
6.  Each circuit is then executed in the simulator (`ql.execute`) to obtain the measurement probabilities.
7.  The energy contribution of each Hamiltonian term is calculated from these probabilities.
8.  Finally, all contributions are summed, and the total energy (`Et`) is displayed on the screen and saved to `data.dat`.

---

##  Libraries Used

* **NumPy**
* **Qibo**
* **OpenFermion**
* **Qililab**

---

