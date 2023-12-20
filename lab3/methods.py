def find_matrix_norm(matrix):
    res = list()
    for line in matrix:
        temp_res = 0
        for element in line:
            temp_res += abs(element)
        res.append(temp_res)
    print(res)
    return max(res)

def find_vector_norm(vector):
    vector = [abs(x) for x in vector]
    return max(vector)

def yakobi_method(matrix, vector_b, epsilon):
    iteration = 0
    norm_matrix = find_matrix_norm(matrix)
    norm_vector = find_vector_norm(vector_b)
    prev = [0] * 3
    result = vector_b
    # while norm_matrix ** 0 * norm_vector + :
    # print(f'{prev=}')
    # print(f'{approx=}')
    while find_vector_norm([abs(prev[0] - result[0]), abs(prev[1] - result[1]), abs(prev[2] - result[2])]) >= epsilon:
        prev = result
        result = list()
        for i in range(len(matrix)):
            local_result = vector_b[i]
            print(local_result)
            multiplier = 1 / matrix[i][i]
            for j in range(len(matrix)):
                if i == j:
                    continue
                # print(-1 * (matrix[i][j]))
                # print(f'{prev[j]=}')
                local_result += -1 * matrix[i][j] * prev[j]
            local_result *= multiplier
            result.append(local_result)
        print(result)
    return result

# print(find_vector_norm([1.3432, 3.2222, 4.7]))
# print(find_matrix_norm([
#         [1.69, -0.01, 0.2],
#         [-0.15, 0.81, 0.38],
#         [-0.11, -0.37, -0.5]
# ]))
print(yakobi_method([
        [1.69, -0.01, 0.2],
        [-0.15, 0.81, 0.38],
        [-0.11, -0.37, -0.5]
],
[1.3432, 3.2222, 4.7], 0.5 * 10 ** -4))