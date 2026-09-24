from vpython import *

# ---------------------------------------
# 3D SCENE
# ---------------------------------------

scene = canvas(
    title="Virtual Medical Nanite",
    width=1000,
    height=700,
    background=vector(0.02, 0.02, 0.04)
)

scene.camera.pos = vector(0, 0, 18)
scene.camera.axis = vector(0, 0, -18)

# ---------------------------------------
# NANITE BODY
# ---------------------------------------

body = sphere(
    pos=vector(0, 0, 0),
    radius=2,
    color=vector(0.35, 0.40, 0.45),
    shininess=1
)

# ---------------------------------------
# CENTRAL GLOWING SENSOR
# ---------------------------------------

sensor = sphere(
    pos=vector(0, 0, 1.85),
    radius=0.45,
    color=vector(0.1, 0.8, 1),
    emissive=True
)

# Outer sensor ring
sensor_ring = ring(
    pos=vector(0, 0, 1.9),
    axis=vector(0, 0, 1),
    radius=0.75,
    thickness=0.12,
    color=vector(0.1, 0.6, 0.9)
)

# ---------------------------------------
# TOP DOME
# ---------------------------------------

dome = sphere(
    pos=vector(0, 0.8, 0.5),
    radius=0.8,
    color=vector(0.55, 0.58, 0.62),
    shininess=1
)

# Make dome appear partially embedded
dome.size = vector(1.2, 0.6, 1.2)

# ---------------------------------------
# ROBOTIC BODY RINGS
# ---------------------------------------

ring1 = ring(
    pos=vector(0, 0, 0),
    axis=vector(1, 0, 0),
    radius=2.05,
    thickness=0.12,
    color=vector(0.15, 0.6, 0.8)
)

ring2 = ring(
    pos=vector(0, 0, 0),
    axis=vector(0, 1, 0),
    radius=2.05,
    thickness=0.12,
    color=vector(0.15, 0.6, 0.8)
)

# ---------------------------------------
# MECHANICAL ARMS
# ---------------------------------------

arm_data = [
    (vector(1, 0, 0), 0),
    (vector(-1, 0, 0), 0),
    (vector(0, 1, 0), 0),
    (vector(0, -1, 0), 0),
    (vector(0, 0, 1), 0),
    (vector(0, 0, -1), 0)
]

arms = []

for direction, unused in arm_data:

    # First arm segment
    start = direction * 1.6

    segment1 = cylinder(
        pos=start,
        axis=direction * 1.0,
        radius=0.13,
        color=vector(0.65, 0.68, 0.72),
        shininess=1
    )

    # Mechanical joint
    joint = sphere(
        pos=start + direction * 1.0,
        radius=0.25,
        color=vector(0.25, 0.28, 0.32),
        shininess=1
    )

    # Second arm segment
    segment2 = cylinder(
        pos=joint.pos,
        axis=direction * 0.8,
        radius=0.10,
        color=vector(0.7, 0.72, 0.75),
        shininess=1
    )

    # Arm tip
    tip = sphere(
        pos=joint.pos + direction * 0.8,
        radius=0.18,
        color=vector(0.1, 0.7, 0.9),
        emissive=True
    )

    arms.append((segment1, joint, segment2, tip))

# ---------------------------------------
# SMALL BODY DETAILS
# ---------------------------------------

# Small circular sensors around body

sensor_positions = [
    vector(1.4, 1.0, 1.0),
    vector(-1.4, 1.0, 1.0),
    vector(1.4, -1.0, 1.0),
    vector(-1.4, -1.0, 1.0)
]

for position in sensor_positions:

    sphere(
        pos=position,
        radius=0.16,
        color=vector(0.9, 0.25, 0.15),
        emissive=True
    )

# ---------------------------------------
# LABEL
# ---------------------------------------

label(
    pos=vector(0, -3.5, 0),
    text="VIRTUAL NANITE",
    height=20,
    box=False,
    color=color.white
)

# ---------------------------------------
# ANIMATION
# ---------------------------------------

angle = 0

while True:

    rate(60)

    angle += 0.02

    # Gentle floating movement
    body.pos.y = 0.3 * sin(angle)
    body.pos.z = 0.3 * cos(angle)

    # Move sensor with body
    sensor.pos = body.pos + vector(0, 0, 1.85)

    sensor_ring.pos = body.pos + vector(0, 0, 1.9)

    # Move dome
    dome.pos = body.pos + vector(0, 0.8, 0.5)

    # Move body rings
    ring1.pos = body.pos
    ring2.pos = body.pos