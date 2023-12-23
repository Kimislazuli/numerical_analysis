def find_matrix_norm(matrix):
    res = list()
    for line in matrix:
        temp_res = 0
        for element in line:
            temp_res += abs(element)
        res.append(temp_res)
    return max(res)


def find_vector_norm(vector):
    vector = [abs(x) for x in vector]
    return max(vector)


def yakobi_method(matrix_b, vector_c, approx, epsilon):
    norm_matrix = find_matrix_norm(matrix_b)
    approx_norm = find_vector_norm(approx)
    c_norm = find_vector_norm(vector_c)
    iteration = 0
    prev = [0] * 3
    result = approx
    while norm_matrix ** iteration * approx_norm + (norm_matrix ** iteration * c_norm) / (1 - norm_matrix) >= epsilon:
    # while find_vector_norm([abs(prev[0] - result[0]), abs(prev[1] - result[1]), abs(prev[2] - result[2])]) > epsilon:
        iteration += 1
        prev = result
        result = list()
        result.append(vector_c[0] + matrix_b[0][1] * prev[1] + matrix_b[0][2] * prev[2])
        result.append(vector_c[1] + matrix_b[1][0] * prev[0] + matrix_b[1][2] * prev[2])
        result.append(vector_c[2] + matrix_b[2][1] * prev[1] + matrix_b[2][0] * prev[0])
    return result, iteration


def gauss_seidel_method(matrix_b, vector_c, approx, epsilon):
    iteration = 0
    norm_matrix = find_matrix_norm(matrix_b)
    approx_norm = find_vector_norm(approx)
    c_norm = find_vector_norm(vector_c)
    prev = [0] * 3
    result = approx
    while norm_matrix ** iteration * approx_norm + (norm_matrix ** iteration * c_norm) / (1 - norm_matrix) >= epsilon:
    # while find_vector_norm([abs(prev[0] - result[0]), abs(prev[1] - result[1]), abs(prev[2] - result[2])]) > epsilon:
        iteration += 1
        prev = result
        result = list()
        result.append(vector_c[0] + matrix_b[0][1] * prev[1] + matrix_b[0][2] * prev[2])
        result.append(vector_c[1] + matrix_b[1][1] * prev[1] + matrix_b[1][2] * prev[2])
        result.append(vector_c[2] + matrix_b[2][1] * prev[1] + matrix_b[2][2] * prev[2])
    return result, iteration


print(yakobi_method(
    [
        [0, 0.005917, -0.11834],
        [0.18519, 0, -0.469148],
        [-0.22, -0.74, 0]
    ],
    [1.34319532, 3.2222277, 4.7],
    [1.34319532, 3.2222277, 4.7],
    0.5 * 10 ** -4))

print(gauss_seidel_method([
    [0, 0.0059172, -0.118344],
    [0, 0.0010958, -0.4910526],
    [0, -0.00211265, 0.3894134]
],
    [1.3432044, 3.4709743, 1.83598465],
    [1.34319532, 3.2222277, 4.7],
    0.5 * 10 ** -4
))