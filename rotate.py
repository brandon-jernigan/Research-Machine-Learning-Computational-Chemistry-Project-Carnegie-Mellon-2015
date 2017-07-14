
import math
 
import numpy
 
def get_axis_rotation_matrix(axis, theta):
    # http://stackoverflow.com/questions/6721544/circular-rotation-around-an-arbitrary-axis
    ct = math.cos(theta)
    nct = 1 - ct
    st = math.sin(theta)
    r = numpy.linalg.norm(axis)
    if r == 0.0:
        return numpy.matrix(numpy.eye(3))
    ux = axis[0] / r
    uy = axis[1] / r
    uz = axis[2] / r
    rot = numpy.matrix([
        [ct + ux ** 2 * nct, ux * uy * nct - uz * st, ux * uz * nct + uy * st],
        [uy * ux * nct + uz * st, ct + uy ** 2 * nct, uy * uz * nct - ux * st],
        [uz * ux * nct - uy * st, uz * uy * nct + ux * st, ct + uz ** 2 * nct],
    ])
    return rot
 
 
def get_new_operator(operator, bond_vector):
    e_x = numpy.array([1., 0., 0.])
    e_bond = bond_vector / numpy.linalg.norm(bond_vector)
    e_1 = numpy.cross(e_bond, e_x)
    e_2 = numpy.cross(e_x, e_1)
    theta_1 = math.atan2(numpy.dot(e_bond, e_2), numpy.dot(e_bond, e_x))
 
    R_1 = get_axis_rotation_matrix(e_1, theta_1)
    operator_mid = R_1 * operator * R_1.T
 
    operator_lower = operator_mid[1:, 1:]
    eigen_vals, eigen_vectors = numpy.linalg.eig(operator_lower)
    # Sort eigenvalues [High, ..., Low]
    idx = eigen_vals.argsort()[::-1]
    eigen_sorted = eigen_vectors[:,idx]
 
    # Lazy way to add 1 to the top of the diagonal
    R_2 = numpy.matrix(numpy.eye(3))
    R_2[1:,1:] = eigen_sorted
 
    return R_2 * operator_mid * R_2.T
 
def get_rotation_matrix(bond_vector, operator):
    e_x = numpy.array([1., 0., 0.])
    e_bond = bond_vector / numpy.linalg.norm(bond_vector)
    e_1 = numpy.cross(e_bond, e_x)
    e_2 = numpy.cross(e_x, e_1)
    theta_1 = math.atan2(numpy.dot(e_bond, e_2), numpy.dot(e_bond, e_x))
 
    R_1 = get_axis_rotation_matrix(e_1, theta_1)
    operator_mid = R_1 * operator * R_1.T
 
    operator_lower = operator_mid[1:, 1:]
    eigen_vals, eigen_vectors = numpy.linalg.eig(operator_lower)
    # Sort eigenvalues [High, ..., Low]
    idx = eigen_vals.argsort()[::-1]
    eigen_sorted = eigen_vectors[:,idx]
 
    # Lazy way to add 1 to the top of the diagonal
    R_2 = numpy.matrix(numpy.eye(3))
    R_2[1:,1:] = eigen_sorted

    return R_2 * R_1    


 


    
    
    
if __name__ == "__main__":
    numpy.set_printoptions(precision=5, suppress=True)
    vectors = [
        numpy.array([1., 0., 0.]),
        numpy.array([0., 1., 0.]),
        numpy.array([0., 0., 1.]),
        numpy.random.rand(3),
    ]
    for bond_vector in vectors:
        half_op = numpy.random.rand(3, 3)
        operator = half_op + half_op.T
        print "Operator:"
        print operator
        print "Vector:"
        print bond_vector

        print "Result:"
        print get_new_operator(operator, bond_vector)