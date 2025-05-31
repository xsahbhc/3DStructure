import numpy as np

def cal_dist(a, b, ord=2):
    """计算两点之间的距离"""
    if ord == 1:
        return sum(abs(a - b))
    elif ord == 2:
        return np.sqrt(sum((a - b) ** 2))
    else:
        return pow(sum((a - b) ** ord), 1 / ord)

def area_triangle(v1, v2):
    """计算由两个向量形成的三角形的面积"""
    return 0.5 * np.sqrt(sum(np.cross(v1, v2) ** 2))

def angle(v1, v2):
    """计算两个向量之间的夹角"""
    return np.arccos(np.clip(np.dot(v1, v2) / (np.sqrt(sum(v1**2)) * np.sqrt(sum(v2**2))), -1.0, 1.0))

def one_of_k_encoding_unk(x, allowable_set):
    """将标签编码为one-hot向量，如果不在允许集合中则标记为未知"""
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def one_of_k_encoding(x, allowable_set):
    """将标签编码为one-hot向量"""
    if x not in allowable_set:
        raise ValueError(f"输入 {x} 不在允许集合 {allowable_set} 中")
    return list(map(lambda s: x == s, allowable_set)) 