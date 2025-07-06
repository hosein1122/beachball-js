# compare_plotMT.py
from beachball import plot_mt, PrincipalAxis

# pure explosion example
T = PrincipalAxis(1.0, 0, 0)
N = PrincipalAxis(1.0, 0, 0)
P = PrincipalAxis(1.0, 0, 0)

colors, patches = plot_mt(T, N, P, size=100)
print("Python plot_mt →")
print(" colors =", colors)
print(" patch0:", patches[0])

print(patches[0])










# # compare_plotDC.py
# from beachball import mt2plane, plot_dc

# # اول یک NodalPlane می‌سازیم
# np1 = mt2plane(mt=type("MT", (), {"mt": [[4,1,2],[1,5,3],[2,3,6]]})())

# # سپس plot_dc را اجرا و نتایج را چاپ می‌کنیم
# colours, patches = plot_dc(np1, size=200, xy=[0,0], width=200)
# print("Python plot_dc →")
# print(" colours:", colours)
# print(" number of patches:", len(patches))
# for i, p in enumerate(patches):
#     print(f" patch {i}.vertices[:3] = {p.vertices[:3]}")  # فقط ۳ تا اول








# from beachball import mt2plane, MomentTensor

# Mrr=4.0
# Mtt=5.0
# Mpp=6.0
# Mrt=1.0
# Mrp=2.0
# Mtp=3.0
# # همان مثال بالا
# mt = MomentTensor(Mrr,Mtt,Mpp,Mrt,Mrp,Mtp, 0)

# # خروجی Python
# plane_py = mt2plane(mt)
# print("Python mt2plane →")
# print(f" strike = {plane_py.strike:.6f}")
# print(f" dip    = {plane_py.dip:.6f}")
# print(f" rake   = {plane_py.rake:.6f}")





# import numpy as np
# from beachball import mt2axes, MomentTensor, PrincipalAxis

# Mrr=1.0  # Up–Up
# Mtt=2.0  # North–South
# Mpp=3.0  # East–West
# Mrt=0.0
# Mrp=0.0
# Mtp=0.0

# # 1) Build a simple MomentTensor: λmin=1, λmid=2, λmax=3 on the diagonal
# mt = MomentTensor(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp,0)

# # 2) Compute the principal axes
# T, N, P = mt2axes(mt)

# # 3) Print out the results
# for name, axis in zip(("T", "N", "P"), (T, N, P)):
#     print(f"{name}-axis:")
#     print(f"  val = {axis.val:.6f}")
#     print(f"  strike = {axis.strike:.6f}°")
#     print(f"  dip    = {axis.dip:.6f}°")












# from beachball import tdl

# examples = [
#     ([1, 0, 0], [0, 0, 1]),   # flat‐horizontal branch
#     ([0, 0, 1], [1, 0, 0]),   # singular branch (sd == 0)
# ]

# for an, bn in examples:
#     result = tdl(an, bn)
#     print(f"an={an}, bn={bn} -> {result}")
# # Output should be:
# # an=[1,0,0], bn=[0,0,1] -> [270.0, 90.0, -90.0]
# # an=[0,0,1], bn=[1,0,0] -> None   (or whatever your Python returns for singular)






# from beachball import aux_plane

# cases = [
#     (0, 0, 0),
#     (0, 90, 0),
#     (45, 30, 45),
# ]

# for s1, d1, r1 in cases:
#     s2, d2, r2 = aux_plane(s1, d1, r1)
#     print(f"aux_plane({s1}, {d1}, {r1}) → strike2={s2:.6f}, dip2={d2:.6f}, rake2={r2:.6f}")






# import math

# def xy2patch_py(x, y, res, xy):
#     # normalize res
#     try:
#         length = len(res)
#     except TypeError:
#         res = [res, res]
#     else:
#         if length != 2:
#             raise ValueError("res must contain exactly two elements")

#     # build vertices
#     vertices = [
#         (xi * res[0] + xy[0], yi * res[1] + xy[1])
#         for xi, yi in zip(x, y)
#     ]
#     return {
#         "vertices": vertices,
#         "res": res,
#         "center": xy,
#     }

# if __name__ == "__main__":
#     # x = [0, 1, 2]
#     # y = [0, 10, 20]
#     # patch = xy2patch_py(x, y, 2, [5, 5])
#     # print("Python vertices:", patch["vertices"])
#     # # → Python vertices: [(5, 5), (7, 25), (9, 45)]


#     x = [0]
#     y = [0]
#     patch = xy2patch_py(x, y, [2], [0,0])
#     print("Python vertices:", patch["vertices"])






# import math
# from beachball import strike_dip

# # لیستِ بردارهای نرمال برای تست
# tests = [
#     (0, 0, 1),                           # horizontal plane
#     (1, 0, 0),                           # vertical facing north
#     (0, 1, 0),                           # vertical facing east
#     (0, 0, -1),                          # flipped down
#     (1/math.sqrt(3), 1/math.sqrt(3), 1/math.sqrt(3)),  # diagonal
# ]

# print(" n       e       u    →  strike (deg), dip (deg)")
# print("-" * 50)
# for n, e, u in tests:
#     strike, dip = strike_dip(n, e, u)
#     print(f"{n:7.4f}, {e:7.4f}, {u:7.4f} →  {strike:8.4f}, {dip:8.4f}")









# import numpy as np

# A = np.array([[4, 1, 2],
#               [1, 5, 3],
#               [2, 3, 6]], float)

# # A = np.array([[5, 0, 0],
# #       [0, 2, 0],
# #       [0, 0, -1]], float)
# w, v = np.linalg.eigh(A)

# # نزولی (descending)
# # w_desc = w[::-1]
# print("Eigenvalues:", w)
# # print("Eigenvalues (desc):", w_desc)
# # → [9.41883268, 3.38677016, 2.19439717]

# # بردارهای ویژه برای نزولی (ستون‌ها)
# # v_desc = v[:, ::-1]
# print("Eigenvectors :\n", v)
# # print("Eigenvectors for desc order:\n", v_desc)
# # → [[-0.3744, -0.8156, -0.4412],
# #    [-0.5774,  0.5774, -0.5774],
# #    [-0.7256, -0.0386,  0.6870]]


# print("===================")
