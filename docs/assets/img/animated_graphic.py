
import plotly.graph_objects as go
import plotly
import numpy as np

# Parameters
num_points = 256
max_x = 10 * np.pi
t = np.linspace(0, max_x, num_points)
frequencies = np.linspace(0.02, 6.0, 64)
random_amplitudes = np.random.uniform(0.01, 0.1, len(frequencies))

# Complex wave creation with random amplitudes
def create_complex_wave(t, frequencies, amplitudes, phase):
    wave = sum(amplitude * np.sin(frequency * t + phase) for frequency, amplitude in zip(frequencies, amplitudes))
    return wave

# Damping function
def apply_fixed_damping(wave, phase_shift):
    damping_start = 0.25 * max_x
    damping_end = 0.50 * max_x
    damping_factor = np.ones_like(t)

    for i, x in enumerate(t):
        if x > damping_start and x < damping_end:
            damping_factor[i] = 1 - (x - damping_start) / (damping_end - damping_start)
        elif x >= damping_end:
            damping_factor[i] = 0.01
    
    return wave * damping_factor

# Function to assign color based on x-coordinate
def get_color_for_point(x, max_x):
    dark_slate_blue = (44, 56, 99)
    white = (255, 255, 255)
    white_transition_start = 0.75 * max_x

    if x < white_transition_start:
        color_ratio = x / white_transition_start
        r = int(dark_slate_blue[0] + (white[0] - dark_slate_blue[0]) * color_ratio)
        g = int(dark_slate_blue[1] + (white[1] - dark_slate_blue[1]) * color_ratio)
        b = int(dark_slate_blue[2] + (white[2] - dark_slate_blue[2]) * color_ratio)
    else:
        r, g, b = white
    
    return f'rgb({r}, {g}, {b})'

# Generate colors for each point
colors = [get_color_for_point(x, max_x) for x in t]

# Initial plot
initial_wave = create_complex_wave(t, frequencies, random_amplitudes, 0)
initial_damped_wave = apply_fixed_damping(initial_wave, 0)

# Plotly setup with color gradient
fig = go.Figure(
    data=[go.Scatter(
        x=t, 
        y=initial_damped_wave, 
        mode='markers', 
        marker=dict(color=colors, size=2.5, showscale=False)
    )],
    layout=go.Layout(
        xaxis=dict(range=[0.1*max_x, max_x], autorange=False, showgrid=False, visible=False),
        yaxis=dict(range=[-3.0, 3.0], autorange=False, showgrid=False, visible=False),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    ),
    frames=[go.Frame(data=[go.Scatter(
        x=t,
        y=apply_fixed_damping(create_complex_wave(t, frequencies, random_amplitudes, -phase), -phase),
        mode='markers', 
        marker=dict(color=colors, size=2.5)
    )]) for phase in np.linspace(0, 20 * np.pi, 800)]
)

# Convert frames to a list, append the first frame, and update the frames attribute
frame_list = list(fig.frames)
frame_list.append(frame_list[0])
fig.frames = tuple(frame_list)

# Animation configuration
fig.update_layout(
    updatemenus=[dict(
        type="buttons",
        showactive=False,
        buttons=[dict(
            label="Play",
            method="animate",
            args=[None, {"frame": {"duration": 100, "redraw": True}, "fromcurrent": True, "transition": {"duration": 100, "easing": "linear"}}]
        )]
    )]
)

plotly.offline.plot(fig, filename='hydrochrono_animation.html', include_plotlyjs=True)

fig.show()