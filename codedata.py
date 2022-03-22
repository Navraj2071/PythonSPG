# Hourly mean wind speed
super_height = [0, 10, 15, 20, 30, 50, 60, 70, 80, 90, 100, 1000]
mean_wind_speed_plain = [27.8, 27.8, 29.2, 30.3, 31.4, 33.1, 33.6, 34.0, 34.40, 34.9, 35.3, 35.3]
mean_wind_speed_obs = [17.8, 17.8, 19.6, 21.0, 22.8, 24.9, 25.6, 26.2, 26.9, 27.5, 28.2, 28.2]

# Live Load lanes
carriageway = [0, 5.3, 9.6, 13.1, 16.6, 20.1, 23.6, 100]
lanes = [1, 2, 3, 4, 5, 6, 7, 7]

# ClassA
classA = [2.7, 2.7, 11.4, 11.4, 6.8, 6.8, 6.8, 6.8, 2.7, 2.7, 11.4, 11.4, 6.8, 6.8, 6.8, 6.8]
classA_dist = [0, 1.1, 3.2, 1.2, 4.3, 3, 3, 3, 20, 1.1, 3.2, 1.2, 4.3, 3, 3, 3]

# 70R wheel
class_70r = [8, 12, 12, 17, 17, 17, 17, 8, 12, 12, 17, 17, 17, 17]
class_70r_dist = [0, 3.96, 1.52, 2.13, 1.37, 3.05, 1.37, 30, 3.96, 1.52, 2.13, 1.37, 3.05, 1.37]

# SV
class_sv = [6, 9.5, 9.5, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18]
class_sv_dist = [0, 3.2, 1.37, 5.389, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
                 1.5, 1.5, 1.5, ]

# Congestion Factor
span = [10, 30, 40, 50, 60, 70, 100]
congestion_factor = [1.15, 1.15, 1.3, 1.45, 1.6, 1.7, 1.7]

# Load Combinations
load = ['Selfweight of Girder', 'Weight of deck slab', 'Weight of crash barrier', 'Weight of green concrete',
        'Weight of shuttering', 'Construction Live Load', 'Weight of surfacing', 'Wind load construction stage',
        'Wind load service stage without Live Load', 'Wind Load service stage with Live Load', 'Seismic Load',
        'Live Load Normal', 'Live Load SV Load', 'Live Load Seismic']
load_combination = ['ULS_1', 'ULS_2', 'ULS_3', 'ULS_4', 'ULS_5', 'ULS_6', 'ULS_7', 'ULS_8']
factor = [
    [1.35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1.35, 0, 0, 1.35, 1.35, 1.5, 0, 1.5, 0, 0, 0, 0, 0, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 1.5, 0, 0, 0, 0, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 0, 0.9, 0, 1.5, 0, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 0, 1.5, 0, 1.15, 0, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 0, 0, 0, 0, 1.15, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 0, 0, 0.75, 0, 0, 0],
    [1.35, 1.35, 1.35, 0, 0, 0, 1.75, 0, 0, 0, 1.5, 0.2, 0, 1.5]
]
