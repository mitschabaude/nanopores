import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np
x = np.random.randn(500)

data = [
    go.Histogram(
        x=x
    )
]
plot_url = py.plot(data, filename='basic-histogram')