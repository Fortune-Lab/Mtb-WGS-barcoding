import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt
from itertools import tee
import numpy as np

plt.rcParams["figure.figsize"] = (10,6)

df = pd.read_csv('LIB065162_primary_reads_sort_count_remove_singletons.tsv', sep='\t', header=None)

df['2'] = np.arange(1, len(df) + 1)
df['log_col1'] = np.log10(df[1])

x = df['2']

y = df['log_col1']


#coordinate_list = list(zip(x, y))

#def second_derivative(y, x):
#    """
#    Calculates the second derivative of y with respect to x using NumPy's gradient function.
#
#    Args:
#        y (numpy.ndarray): Array of y-values.
#        x (numpy.ndarray): Array of x-values.
#
#    Returns:
#        numpy.ndarray: Array of second derivative values.
#    """
#    dy_dx = np.gradient(y, x)
#    d2y_dx2 = np.gradient(dy_dx, x)
#    return d2y_dx2


#def window(iterable, size):
#    iters = tee(iterable, size)
#    for i in range(500, size):
#        for each in iters[i:(size-1)]:
#            next(each, None)
#    return zip(*iters)

def get_slopes(x, y, size, shoulder):
    coef = []
    index = 0
    while index <= len(y) :
        #print(size, shoulder)
        if index >= size and index <= shoulder :
            x1 = x[index-(size - 1):index]
            y1 = y[index-(size - 1):index]
            b,m = np.polyfit(x1, y1, 1)
            #print(index, b, m, sep="\t")
            coef.append([index, b, m])
            index += 1
        else :
            index += 1
    return coef

coef = get_slopes(x, y, 500, 50000)

counts = y[499:50000]

resulty = [t[1] for t in coef]
resultx = [t[0] for t in coef]       

def local_min(y, start, stop):
    min_value = min(y[start:stop])
    min_indices = [i for i, x in enumerate(resulty) if x == min_value]
    return min_indices

min_indices = local_min(resulty, 5000, 50000)
#min_indices = [i for i, x in enumerate(resulty) if x == min_value]
print(min_indices)

fig, ax1 = plt.subplots()
plt.axvline(x=min_indices[0] + 500, color='g', linestyle='-', label='inflection')
# Plot the first graph on the first axes
ax1.plot(resultx, counts, color='red', label='barcode counts')
ax1.set_xlabel('barcode index')
ax1.set_ylabel('barcode counts', color='red')
ax1.tick_params(axis='y', labelcolor='red')

ax2 = ax1.twinx()

# Plot the second graph on the second axes
ax2.plot(resultx, resulty, color='blue', label='linear regression slope')
ax2.set_ylabel('linear regression slope', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')

plt.savefig("library_graph.png")

