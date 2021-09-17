def piecewise2(x, data):
    #请在下面编码实现分段常函数
    n = len(data)
    def I(x, L, R):
        if L <= x and x < R:
            return 1
        else :
            return 0
    fx = 0
    for i in range(n - 1):
        fx += data[i][1] * I(x, data[i][0], data[i + 1][0])
    #请不要修改下面的代码
    return fx
def piecewise(x, data):
    n = len(data)
    for i in range(n):
        if data[i][0] >= x and i == n - 1:
            return data[i - 1][1]
        if data[i][0] > x:
                return data[i - 1][1]
def my_piecewise(x, data):
    return piecewise(x, data)