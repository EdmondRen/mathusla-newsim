import numpy as np

# ----------------------------------------------
# https://stackoverflow.com/questions/3320311/fill-outside-of-polygon-mask-array-where-indicies-are-beyond-a-circular-bounda
def concat(*arrs) -> np.ndarray:
    return np.concatenate(tuple(map(np.asarray, arrs)))

def insert_at(outer_arr, arr, n) -> np.ndarray:
    outer_arr = np.asarray(outer_arr)
    prev, post = np.split(outer_arr, (n,))
    return concat(prev, arr, post)

def cross2d(x1, y1, x2, y2):
    return x1*y2-x2*y1

def is_clockwise(x1, y1, x2, y2):
    cp = cross2d(x1, y1, x2, y2)
    return cp < 0 if cp != 0 else None

def fill_outside(x, y, ll, ur, counter_clockwise=None):
    """
    Creates a polygon where x and y form a crevice of an outer
    rectangle with lower left and upper right corners `ll` and `ur`
    respectively. If `counter_clockwise` is `None` then the orientation
    of the outer polygon will be guessed to be the opposite of the
    inner connecting points.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    xmin, ymin = ll
    xmax, ymax = ur
    xmin, ymin = min(xmin, min(x)), min(ymin, min(y))
    xmax, ymax = max(xmax, max(x)), max(ymax, max(y))
    corners = np.array([
        [xmin, ymin],
        [xmin, ymax],
        [xmax, ymax],
        [xmax, ymin],
        [xmin, ymin],
    ])
    lower_left = corners[0]
    # Get closest point to splicing corner
    x_off, y_off = x-lower_left[0], y-lower_left[1]
    closest_n = (x_off**2+y_off**2).argmin()
    # Guess orientation
    p = [x_off[closest_n], y_off[closest_n]]
    try:
        pn = [x_off[closest_n+1], y_off[closest_n+1]]
    except IndexError:
        # wrap around if we're at the end of the array
        pn = [x_off[0], y_off[0]]
    if counter_clockwise is None:
        counter_clockwise = not is_clockwise(*p, *pn)
    corners = corners[::-1] if counter_clockwise else corners
    # Join the arrays
    corners = concat(np.array([[x[closest_n], y[closest_n]]]), corners)
    xs, ys = np.transpose(corners)
    return insert_at(x, xs, closest_n), insert_at(y, ys, closest_n)