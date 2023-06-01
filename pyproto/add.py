def ADD(Px, Py, Qx, Qy, a, b):
    dY = Qy - Py        # N + 1
    dX = Qx - Px        # N + 1
    h1 = dY / dX        # 2N     (there is inversion here.)
    h2 = h1 * h1        # 4N
    h3 = h2 * h1        # 6N
    