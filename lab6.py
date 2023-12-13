def divided_difference(x_values, y_values):
    n = len(x_values)
    if n == 1:
        return y_values[0]
    else:
        return (divided_difference(x_values[1:], y_values[1:]) - divided_difference(x_values[:-1], y_values[:-1])) / (x_values[-1] - x_values[0])

def newton_interpolation(x_values, y_values, x_interpolate):
    result = y_values[0]
    for i in range(1, len(x_values)):
        _ = 1
        for j in range(i):
            _ *= (x_interpolate - x_values[j])
        result += _ * divided_difference(x_values[:i+1], y_values[:i+1])
    return result

# Пример использования
x_values = [0, 1, 2, 3, 4]
y_values = [1, 2, 4, 3, 1]

print(divided_difference(x_values,y_values))
x_interpolate = 1.5  # Точка, в которой мы хотим интерполировать

result = newton_interpolation(x_values, y_values, x_interpolate)
print(f"N3({x_interpolate}) = {result}")
