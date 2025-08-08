#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
"""Jordan-Wigner transform on fermionic operators."""
import itertools


from openfermion.ops.operators import FermionOperator, QubitOperator


def jordan_wigner(operator):
    r"""Apply the Jordan-Wigner transform to a FermionOperator,
    InteractionOperator, or DiagonalCoulombHamiltonian to convert
    to a QubitOperator.

    Operators are mapped as follows:

    $$
        a_j^\dagger -> Z_0 .. Z_{j-1} (X_j - iY_j) / 2
        a_j -> Z_0 .. Z_{j-1} (X_j + iY_j) / 2
    $$

    Returns:
        transformed_operator: An instance of the QubitOperator class.

    Warning:
        The runtime of this method is exponential in the maximum locality
        of the original FermionOperator.

    Raises:
        TypeError: Operator must be a FermionOperator,
            DiagonalCoulombHamiltonian, or InteractionOperator.
    """
    if isinstance(operator, FermionOperator):
        return _jordan_wigner_fermion_operator(operator)
    
    
    raise TypeError(
        "Operator must be a FermionOperator, "
        "MajoranaOperator, "
        "DiagonalCoulombHamiltonian, or "
        "InteractionOperator."
    )


def _jordan_wigner_fermion_operator(operator):
    transformed_operator = QubitOperator()
    # Purpose is storing ladder terms already transformed.
    lookup_ladder_terms = dict()
    for term in operator.terms:
        # Initialize identity matrix.
        transformed_term = QubitOperator((), operator.terms[term])
        # Loop through operators, transform and multiply.
        for ladder_operator in term:
            if ladder_operator not in lookup_ladder_terms:
                z_factors = tuple((index, 'Z') for index in range(ladder_operator[0]))
                pauli_x_component = QubitOperator(z_factors + ((ladder_operator[0], 'X'),), 0.5)
                if ladder_operator[1]:
                    pauli_y_component = QubitOperator(
                        z_factors + ((ladder_operator[0], 'Y'),), -0.5j
                    )
                else:
                    pauli_y_component = QubitOperator(
                        z_factors + ((ladder_operator[0], 'Y'),), 0.5j
                    )
                lookup_ladder_terms[ladder_operator] = pauli_x_component + pauli_y_component
            transformed_term *= lookup_ladder_terms[ladder_operator]
        transformed_operator += transformed_term
    return transformed_operator




def jordan_wigner_one_body(p, q, coefficient=1.0):
    r"""Map the term a^\dagger_p a_q + h.c. to QubitOperator.

    Note that the diagonal terms are divided by a factor of 2
    because they are equal to their own Hermitian conjugate.
    """
    # Handle off-diagonal terms.
    qubit_operator = QubitOperator()
    if p != q:
        if p > q:
            p, q = q, p
            coefficient = coefficient.conjugate()
        parity_string = tuple((z, 'Z') for z in range(p + 1, q))
        for c, (op_a, op_b) in [
            (coefficient.real, 'XX'),
            (coefficient.real, 'YY'),
            (coefficient.imag, 'YX'),
            (-coefficient.imag, 'XY'),
        ]:
            operators = ((p, op_a),) + parity_string + ((q, op_b),)
            qubit_operator += QubitOperator(operators, 0.5 * c)

    # Handle diagonal terms.
    else:
        qubit_operator += QubitOperator((), 0.5 * coefficient)
        qubit_operator += QubitOperator(((p, 'Z'),), -0.5 * coefficient)

    return qubit_operator


def jordan_wigner_two_body(p, q, r, s, coefficient=1.0):
    r"""Map the term a^\dagger_p a^\dagger_q a_r a_s + h.c. to QubitOperator.

    Note that the diagonal terms are divided by a factor of two
    because they are equal to their own Hermitian conjugate.
    """
    # Initialize qubit operator.
    qubit_operator = QubitOperator()

    # Return zero terms.
    if (p == q) or (r == s):
        return qubit_operator

    # Handle case of four unique indices.
    elif len(set([p, q, r, s])) == 4:
        if (p > q) ^ (r > s):
            coefficient *= -1
        # Loop through different operators which act on each tensor factor.
        for ops in itertools.product('XY', repeat=4):
            # Get coefficients.
            if ops.count('X') % 2:
                coeff = 0.125 * coefficient.imag
                if ''.join(ops) in ['XYXX', 'YXXX', 'YYXY', 'YYYX']:
                    coeff *= -1
            else:
                coeff = 0.125 * coefficient.real
                if ''.join(ops) not in ['XXYY', 'YYXX']:
                    coeff *= -1
            if not coeff:
                continue

            # Sort operators.
            [(a, operator_a), (b, operator_b), (c, operator_c), (d, operator_d)] = sorted(
                zip([p, q, r, s], ops)
            )

            # Compute operator strings.
            operators = ((a, operator_a),)
            operators += tuple((z, 'Z') for z in range(a + 1, b))
            operators += ((b, operator_b),)
            operators += ((c, operator_c),)
            operators += tuple((z, 'Z') for z in range(c + 1, d))
            operators += ((d, operator_d),)

            # Add term.
            qubit_operator += QubitOperator(operators, coeff)

    # Handle case of three unique indices.
    elif len(set([p, q, r, s])) == 3:
        # Identify equal tensor factors.
        if p == r:
            if q > s:
                a, b = s, q
                coefficient = -coefficient.conjugate()
            else:
                a, b = q, s
                coefficient = -coefficient
            c = p
        elif p == s:
            if q > r:
                a, b = r, q
                coefficient = coefficient.conjugate()
            else:
                a, b = q, r
            c = p
        elif q == r:
            if p > s:
                a, b = s, p
                coefficient = coefficient.conjugate()
            else:
                a, b = p, s
            c = q
        elif q == s:
            if p > r:
                a, b = r, p
                coefficient = -coefficient.conjugate()
            else:
                a, b = p, r
                coefficient = -coefficient
            c = q

        # Get operators.
        parity_string = tuple((z, 'Z') for z in range(a + 1, b))
        pauli_z = QubitOperator(((c, 'Z'),))
        for c, (op_a, op_b) in [
            (coefficient.real, 'XX'),
            (coefficient.real, 'YY'),
            (coefficient.imag, 'YX'),
            (-coefficient.imag, 'XY'),
        ]:
            operators = ((a, op_a),) + parity_string + ((b, op_b),)
            if not c:
                continue

            # Add term.
            hopping_term = QubitOperator(operators, c / 4)
            qubit_operator -= pauli_z * hopping_term
            qubit_operator += hopping_term

    # Handle case of two unique indices.
    elif len(set([p, q, r, s])) == 2:
        # Get coefficient.
        if p == s:
            coeff = -0.25 * coefficient
        else:
            coeff = 0.25 * coefficient

        # Add terms.
        qubit_operator -= QubitOperator((), coeff)
        qubit_operator += QubitOperator(((p, 'Z'),), coeff)
        qubit_operator += QubitOperator(((q, 'Z'),), coeff)
        qubit_operator -= QubitOperator(((min(q, p), 'Z'), (max(q, p), 'Z')), coeff)

    return qubit_operator
